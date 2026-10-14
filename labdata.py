"""
Reading and processing spectra exported by the Raman simulator.

Free of Streamlit imports, so the same functions work from a notebook:

    from labdata import parse_csv, despike, als_baseline
    s = parse_csv(open("Maya_scan03.csv").read())
    y = despike(s.y)
    y = y - als_baseline(y)

Processing follows the order a spectrometer package applies it in: cosmic rays
first (they wreck every later step), then reference subtraction, baseline,
smoothing and normalisation.
"""

from __future__ import annotations

import io
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, medfilt, peak_widths, savgol_filter
from scipy.sparse.linalg import spsolve

NUMERIC_KEYS = {
    "scan", "stage_x_mm", "stage_y_mm", "focus_um", "NA",
    "laser_power_pct", "integration_s", "accumulations",
    "total_exposure_s", "estimated_SNR",
}
BOOL_KEYS = {"illuminator_open"}


class ParseError(ValueError):
    """Raised when a file is not a spectrum export from this instrument."""


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
@dataclass
class Spectrum:
    meta: dict
    x: np.ndarray       # Raman shift, cm-1
    y: np.ndarray       # mean counts per accumulation
    source: str = ""    # file name, for uploads

    @property
    def label(self) -> str:
        m = self.meta
        if self.source:
            parts = [self.source.rsplit(".", 1)[0]]
        elif m.get("scan") is not None:
            parts = [f"Scan {int(m['scan'])}"]
        else:
            parts = ["spectrum"]
        if m.get("location"):
            parts.append(str(m["location"]))
        if m.get("objective"):
            parts.append(str(m["objective"]))
        return " · ".join(parts)


def _coerce(key: str, raw: str):
    raw = raw.strip()
    if key in NUMERIC_KEYS:
        try:
            return float(raw)
        except ValueError:
            return None
    if key in BOOL_KEYS:
        return raw.lower() == "true"
    return raw


def parse_csv(text: str | bytes, source: str = "") -> Spectrum:
    """Read an export produced by the instrument's 'Save spectrum' button."""
    if isinstance(text, bytes):
        text = text.decode("utf-8", errors="replace")

    meta, body = {}, []
    for line in text.splitlines():
        if line.startswith("#"):
            piece = line.lstrip("#").strip()
            if "," in piece:
                key, _, val = piece.partition(",")
                key = key.strip()
                meta[key] = _coerce(key, val)
        elif line.strip():
            body.append(line)

    if not body:
        raise ParseError("No data rows found. Is this a spectrum export?")

    df = pd.read_csv(io.StringIO("\n".join(body)))
    if not {"wavenumber_cm-1", "mean_counts"}.issubset(df.columns):
        raise ParseError(
            "Missing expected columns. Found: " + ", ".join(df.columns) +
            ". This does not look like a spectrum export."
        )
    df = df.dropna(subset=["wavenumber_cm-1", "mean_counts"])
    if df.empty:
        raise ParseError("The file parsed but contains no usable rows.")
    df = df.sort_values("wavenumber_cm-1")

    return Spectrum(meta=meta,
                    x=df["wavenumber_cm-1"].to_numpy(float),
                    y=df["mean_counts"].to_numpy(float),
                    source=source)


# --------------------------------------------------------------------------
# Processing
# --------------------------------------------------------------------------
def despike(y: np.ndarray, threshold: float = 6.0) -> np.ndarray:
    """Remove cosmic-ray spikes.

    A spike and a sharp Raman band can be equally tall, so height alone cannot
    separate them. Shape can: a spike is one to three pixels wide and has
    almost dropped to the floor two pixels out, while a band of this
    instrument's width is still above half height there.
    """
    y = np.asarray(y, float)
    n = y.size
    if n < 9:
        return y.copy()
    med = medfilt(y, 5)
    resid = y - med
    mad = np.median(np.abs(resid - np.median(resid))) * 1.4826 or 1.0

    pad = np.pad(y, 2, mode="edge")
    flank = (pad[:-4] + pad[4:]) / 2           # y[i-2] and y[i+2]
    floor = pd.Series(y).rolling(17, center=True, min_periods=1).min().to_numpy()
    above = np.maximum(y - floor, 1e-12)
    narrow = (y - flank) / above > 0.8

    spike = (resid > threshold * mad) & narrow
    if not spike.any():
        return y.copy()
    spike = np.convolve(spike, np.ones(3), mode="same") > 0
    good = ~spike
    out = y.copy()
    idx = np.arange(n)
    out[spike] = np.interp(idx[spike], idx[good], y[good])
    return out


def poly_baseline(x: np.ndarray, y: np.ndarray, order: int = 3,
                  iterations: int = 100) -> np.ndarray:
    """Iterative modified polynomial fit (Lieber & Mahadevan-Jansen, 2003).

    Each pass clips the data down to the current fit, so the polynomial sinks
    under the bands and settles on the background.
    """
    xs = (x - x.mean()) / (x.std() or 1.0)
    work = np.asarray(y, float).copy()
    base = work
    for _ in range(iterations):
        base = np.polyval(np.polyfit(xs, work, order), xs)
        nxt = np.minimum(work, base)
        if np.allclose(nxt, work):
            break
        work = nxt
    return base


def als_baseline(y: np.ndarray, lam: float = 1e6, p: float = 0.01,
                 iterations: int = 10) -> np.ndarray:
    """Asymmetric least squares baseline (Eilers & Boelens, 2005).

    lam sets stiffness; p is the weight given to points above the baseline, so
    small p keeps the baseline under the bands.
    """
    y = np.asarray(y, float)
    n = y.size
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(n, n - 2))
    P = lam * D.dot(D.T)
    w = np.ones(n)
    z = y
    for _ in range(iterations):
        W = sparse.spdiags(w, 0, n, n)
        z = spsolve(sparse.csc_matrix(W + P), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def smooth(y: np.ndarray, window: int = 9, order: int = 3) -> np.ndarray:
    """Savitzky-Golay smoothing. Window is in pixels and forced odd."""
    window = max(int(window) | 1, order + 2 | 1)
    if window > y.size:
        return y.copy()
    return savgol_filter(y, window, order)


def normalise(x: np.ndarray, y: np.ndarray, mode: str,
              band: tuple[float, float] | None = None) -> np.ndarray:
    if mode == "max":
        d = np.max(y)
    elif mode == "area":
        d = np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x)
    elif mode == "band" and band:
        sel = (x >= band[0]) & (x <= band[1])
        d = np.max(y[sel]) if sel.any() else 0
    else:
        return y
    return y / d if d else y


# --------------------------------------------------------------------------
# Peaks
# --------------------------------------------------------------------------
def noise_level(y: np.ndarray) -> float:
    """Robust noise estimate from point-to-point differences, which bands and
    sloping backgrounds barely affect."""
    d = np.diff(np.asarray(y, float))
    return float(np.median(np.abs(d - np.median(d))) * 1.4826 / np.sqrt(2)) or 1.0


def pick_peaks(x: np.ndarray, y: np.ndarray, snr: float = 8.0,
               min_sep: float = 10.0, min_fwhm: float = 5.0) -> pd.DataFrame:
    """Peak positions, heights and FWHM.

    A peak must stand `snr` noise-widths above its surroundings and be at
    least `min_fwhm` cm-1 wide, so single-pixel noise excursions do not
    qualify. min_sep is in cm-1.
    """
    step = float(np.median(np.diff(x))) or 1.0
    idx, props = find_peaks(y, prominence=snr * noise_level(y),
                            width=min_fwhm / step,
                            distance=max(1, int(min_sep / step)))
    if idx.size == 0:
        return pd.DataFrame(columns=["shift_cm-1", "height", "prominence",
                                     "FWHM_cm-1"])
    widths = peak_widths(y, idx, rel_height=0.5)[0] * step
    return pd.DataFrame({
        "shift_cm-1": x[idx].round(1),
        "height": y[idx],
        "prominence": props["prominences"],
        "FWHM_cm-1": widths.round(1),
    })


def integrate(x: np.ndarray, y: np.ndarray, lo: float, hi: float,
              local_baseline: bool = True) -> float:
    """Band area between two shifts, optionally above a straight line drawn
    between the end points."""
    sel = (x >= lo) & (x <= hi)
    if sel.sum() < 2:
        return float("nan")
    xs, ys = x[sel], y[sel]
    if local_baseline:
        ys = ys - np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])
    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(trap(ys, xs))


# --------------------------------------------------------------------------
# Band fitting
# --------------------------------------------------------------------------
SHAPES = ("Lorentzian", "Gaussian", "Pseudo-Voigt")
_G = 2 * np.sqrt(2 * np.log(2))     # FWHM / sigma


def _lor(x, c, f):
    return 1 / (1 + ((x - c) / (f / 2)) ** 2)


def _gau(x, c, f):
    return np.exp(-0.5 * ((x - c) / (f / _G)) ** 2)


def _model(shape: str, n: int):
    """Sum of n bands on a straight background. Per band the parameters are
    (height, centre, FWHM), plus eta for pseudo-Voigt."""
    k = 4 if shape == "Pseudo-Voigt" else 3

    def band(x, q):
        h, c, w = q[:3]
        if shape == "Lorentzian":
            return h * _lor(x, c, w)
        if shape == "Gaussian":
            return h * _gau(x, c, w)
        return h * (q[3] * _lor(x, c, w) + (1 - q[3]) * _gau(x, c, w))

    def f(x, *p):
        out = p[-2] + p[-1] * x
        for i in range(n):
            out = out + band(x, p[k * i: k * (i + 1)])
        return out
    return f, band, k


def _area(shape, h, w, eta=None):
    lor = np.pi * h * w / 2
    gau = h * (w / _G) * np.sqrt(2 * np.pi)
    if shape == "Lorentzian":
        return lor
    if shape == "Gaussian":
        return gau
    return eta * lor + (1 - eta) * gau


@dataclass
class Fit:
    x: np.ndarray
    y: np.ndarray
    total: np.ndarray
    background: np.ndarray
    components: list
    table: pd.DataFrame
    r2: float


def fit_bands(x: np.ndarray, y: np.ndarray, lo: float, hi: float,
              shape: str = "Lorentzian", n_bands: int = 1,
              guesses: list[float] | None = None) -> Fit:
    """Least-squares fit of n bands plus a linear background over [lo, hi].

    Starting centres come from `guesses` if given, otherwise from the most
    prominent maxima in the window.
    """
    sel = (x >= lo) & (x <= hi)
    xs, ys = x[sel], y[sel]
    if xs.size < 6:
        raise ValueError("Fit window holds too few points; widen it.")

    centres = list(guesses or [])[:n_bands]
    if len(centres) < n_bands:
        idx, pr = find_peaks(ys, prominence=0)
        order = np.argsort(pr["prominences"])[::-1]
        for i in idx[order]:
            if len(centres) >= n_bands:
                break
            if all(abs(xs[i] - c) > 4 for c in centres):
                centres.append(float(xs[i]))
    while len(centres) < n_bands:
        centres.append(float(lo + (hi - lo) * (len(centres) + 1) / (n_bands + 1)))
    centres.sort()

    f, band, k = _model(shape, n_bands)
    base0 = float(np.min(ys))
    span = float(np.ptp(ys)) or 1.0
    p0, lb, ub = [], [], []
    for c in centres:
        h0 = max(float(np.interp(c, xs, ys)) - base0, 0.05 * span)
        p0 += [h0, c, 15.0]
        lb += [0, lo, 1.0]
        ub += [np.inf, hi, hi - lo]
        if k == 4:
            p0.append(0.5); lb.append(0); ub.append(1)
    p0 += [base0, 0.0]
    lb += [-np.inf, -np.inf]
    ub += [np.inf, np.inf]

    popt, pcov = curve_fit(f, xs, ys, p0=p0, bounds=(lb, ub), maxfev=20000)
    perr = np.sqrt(np.clip(np.diag(pcov), 0, None))

    bg = popt[-2] + popt[-1] * xs
    rows, comps = [], []
    for i in range(n_bands):
        h, c, w = popt[k * i: k * i + 3]
        dh, dc, dw = perr[k * i: k * i + 3]
        eta = popt[k * i + 3] if k == 4 else None
        comps.append(band(xs, popt[k * i: k * (i + 1)]))
        row = {"band": i + 1,
               "centre_cm-1": c, "± centre": dc,
               "FWHM_cm-1": w, "± FWHM": dw,
               "height": h, "± height": dh,
               "area": _area(shape, h, w, eta)}
        if eta is not None:
            row["eta (Lorentzian fraction)"] = eta
        rows.append(row)

    total = f(xs, *popt)
    ss_res = float(np.sum((ys - total) ** 2))
    ss_tot = float(np.sum((ys - ys.mean()) ** 2)) or 1.0
    return Fit(x=xs, y=ys, total=total, background=bg, components=comps,
               table=pd.DataFrame(rows), r2=1 - ss_res / ss_tot)
