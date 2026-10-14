"""
Parsing and checking for spectra exported by the Raman simulator.

Free of Streamlit imports, so it can be unit tested and reused from a notebook
or a batch marking script.

The slide holds five capillaries: three glycine polymorphs and the
chloroform/deuterochloroform pair. Which tube holds what is shuffled from the
student ID, so the arrangement differs for every student and is recoverable here
for marking. `assign_tubes` mirrors the shuffle in raman_sim.html exactly; if you
change one, change both.
"""

from __future__ import annotations

import io
import math
from dataclasses import dataclass, field

import pandas as pd

# --------------------------------------------------------------------------
# Slide geometry, mirroring TUBE and SAMPLES in raman_sim.html
# --------------------------------------------------------------------------
BORE_W = 0.80   # mm, bore width
WALL_W = 0.26   # mm, glass wall thickness
LEN_Y = 9.0     # mm, tube length

TUBES = [
    # label,    stage x, stage y
    ("Tube 1", -6.40, 0.0),
    ("Tube 2", -3.20, 0.0),
    ("Tube 3",  0.00, 0.0),
    ("Tube 4",  3.20, 0.0),
    ("Tube 5",  6.40, 0.0),
]

# Order must match CONTENTS in raman_sim.html — the shuffle depends on it.
CONTENTS = [
    ("alpha_glycine", "\u03b1-glycine", "powder"),
    ("beta_glycine",  "\u03b2-glycine", "powder"),
    ("gamma_glycine", "\u03b3-glycine", "powder"),
    ("chcl3",         "CHCl\u2083",     "liquid"),
    ("cdcl3",         "CDCl\u2083",     "liquid"),
]

# Strongest band, plus the bands needed to separate near-identical materials.
# The polymorphs sit ~4 cm-1 apart at the CCN stretch, which no teaching
# spectrometer resolves, so the lattice modes below 300 cm-1 do the work.
DIAGNOSTIC = {
    "alpha_glycine": {"strongest": 893.0, "key_bands": [165, 198, 1324]},
    "beta_glycine":  {"strongest": 889.0, "key_bands": [172, 215, 1340]},
    "gamma_glycine": {"strongest": 877.0, "key_bands": [175, 251, 500]},
    "chcl3":         {"strongest": 667.0, "key_bands": [262, 761, 3019]},
    "cdcl3":         {"strongest": 651.0, "key_bands": [262, 737, 2256]},
}

NUMERIC_KEYS = {
    "scan", "stage_x_mm", "stage_y_mm", "focus_um", "NA",
    "laser_power_pct", "integration_s", "accumulations",
    "total_exposure_s", "estimated_SNR",
}
BOOL_KEYS = {"illuminator_open"}
FULL_WELL = 65535
_M32 = 0xFFFFFFFF


class ParseError(ValueError):
    """Raised when a file is not a spectrum export from this instrument."""


# --------------------------------------------------------------------------
# The seeded shuffle, ported from the instrument
# --------------------------------------------------------------------------
def _imul(a: int, b: int) -> int:
    return (a * b) & _M32


def hash_str(s: str) -> int:
    """FNV-1a, matching hashStr() in the instrument."""
    h = 2166136261
    for ch in s:
        h ^= ord(ch)
        h = _imul(h, 16777619)
    return h & _M32


class Mulberry32:
    """Matching mulberry32() in the instrument."""

    def __init__(self, seed: int):
        self.a = seed & _M32

    def __call__(self) -> float:
        self.a = (self.a + 0x6D2B79F5) & _M32
        a = self.a
        t = _imul(a ^ (a >> 15), 1 | a)
        t = ((t + _imul(t ^ (t >> 7), 61 | t)) ^ t) & _M32
        return ((t ^ (t >> 14)) & _M32) / 4294967296


def assign_tubes(student: str) -> dict:
    """Which material is in which tube, for this student.

    Fisher-Yates driven by the same seed string the instrument uses, so this
    reproduces the arrangement the student actually saw.
    """
    rng = Mulberry32(hash_str(f"{student}|tubes"))
    pool = list(CONTENTS)
    for i in range(len(pool) - 1, 0, -1):
        j = int(rng() * (i + 1))
        pool[i], pool[j] = pool[j], pool[i]
    return {label: pool[k] for k, (label, _, _) in enumerate(TUBES)}


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------
@dataclass
class Spectrum:
    meta: dict
    data: pd.DataFrame

    @property
    def student(self) -> str:
        return str(self.meta.get("student", "")).strip()

    @property
    def tube(self) -> str:
        """Which tube the beam was in, from the stage coordinates. Deliberately
        not written into the student's file."""
        x, y = self.meta.get("stage_x_mm"), self.meta.get("stage_y_mm")
        if x is None or y is None:
            return "unknown"
        for label, tx, ty in TUBES:
            if abs(y - ty) > LEN_Y / 2:
                continue
            if abs(x - tx) <= BORE_W / 2:
                return label
            if abs(x - tx) <= BORE_W / 2 + WALL_W:
                return f"{label} (glass wall)"
        return "off tube"

    @property
    def material(self) -> tuple | None:
        """(key, display name, phase) the student was actually sitting on."""
        tube = self.tube.replace(" (glass wall)", "")
        if not self.student or not tube.startswith("Tube"):
            return None
        return assign_tubes(self.student).get(tube)

    def observed_band(self, window: int = 60) -> float | None:
        """Strongest band above a local background, matching the estimator the
        instrument uses."""
        if self.data.empty:
            return None
        y = self.data["mean_counts"].rolling(5, center=True, min_periods=1).mean()
        floor = y.rolling(2 * window + 1, center=True, min_periods=1).min()
        height = y - floor
        if height.isna().all():
            return None
        return float(self.data["wavenumber_cm-1"].iloc[int(height.idxmax())])

    def intensity_at(self, wn: float, tol: float = 6.0) -> float | None:
        """Peak height above a local floor near a given shift, for checking
        whether a diagnostic band was actually captured."""
        d = self.data
        near = d[(d["wavenumber_cm-1"] - wn).abs() <= tol]
        if near.empty:
            return None
        wide = d[(d["wavenumber_cm-1"] - wn).abs() <= 180]
        return float(near["mean_counts"].max() - wide["mean_counts"].min())


@dataclass
class Check:
    ok: bool
    label: str
    detail: str
    severity: str = "warn"  # "warn" or "info"


@dataclass
class Report:
    spectrum: Spectrum
    checks: list = field(default_factory=list)

    @property
    def problems(self):
        return [c for c in self.checks if not c.ok and c.severity == "warn"]


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


def parse_csv(text: str | bytes) -> Spectrum:
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

    return Spectrum(meta=meta, data=df.reset_index(drop=True))


# --------------------------------------------------------------------------
# Checking
# --------------------------------------------------------------------------
def review(spec: Spectrum, *, min_accumulations: int = 4,
           min_snr: float = 10.0) -> Report:
    """Apply the acquisition-hygiene checks a demonstrator would apply."""
    m, checks = spec.meta, []

    beam = str(m.get("beam_mode", "")).lower()
    checks.append(Check(
        beam == "excite", "Excitation beam",
        "Acquired with the excitation beam." if beam == "excite"
        else f"Beam was '{beam or 'unrecorded'}'. The alignment attenuator passes "
             "~0.2%, so this spectrum carries almost no signal.",
    ))

    lamp = m.get("illuminator_open")
    checks.append(Check(
        lamp is False, "Illuminator closed",
        "White light was closed." if lamp is False
        else "Illuminator was left open; the broadband background will dominate.",
    ))

    obj = str(m.get("objective", ""))
    low = obj in ("5\u00d7", "5x")
    checks.append(Check(
        not low, "Objective",
        f"Collected with {obj or 'an unrecorded objective'}." if not low
        else "Collected at 5\u00d7. Its long depth response reaches straight "
             "through the bore into the capillary wall, so expect glass bands on "
             "top of the sample.",
    ))

    tube = spec.tube
    on_wall = "glass wall" in tube
    in_bore = tube.startswith("Tube") and not on_wall
    checks.append(Check(
        in_bore, "Beam in the bore",
        f"Beam was inside {tube}." if in_bore
        else (f"Beam was on the wall of {tube.replace(' (glass wall)', '')}, not "
              "in the bore — this is a spectrum of the capillary."
              if on_wall else "Beam was not on any tube."),
    ))

    n = m.get("accumulations")
    checks.append(Check(
        n is not None and n >= min_accumulations, "Accumulations",
        f"{int(n)} accumulations (\u221an gain {math.sqrt(n):.2f}\u00d7)." if n
        else "Accumulation count missing.",
    ))

    snr = m.get("estimated_SNR")
    checks.append(Check(
        snr is not None and snr >= min_snr, "Signal to noise",
        f"S/N {snr:.1f}." if snr is not None else "No S/N recorded.",
    ))

    peak = float(spec.data["mean_counts"].max())
    checks.append(Check(
        peak < FULL_WELL * 0.99, "Detector headroom",
        f"Peak {peak:,.0f} counts." if peak < FULL_WELL * 0.99
        else f"Peak reached {peak:,.0f} counts — the detector saturated and the "
             "band shape is clipped.",
    ))

    mat = spec.material
    if mat and in_bore:
        key, name, _ = mat
        exp = DIAGNOSTIC[key]["strongest"]
        obs = spec.observed_band()
        if obs is not None:
            good = abs(obs - exp) <= 10
            checks.append(Check(
                good, "Band position",
                f"Strongest band at {obs:.0f} cm\u207b\u00b9; {tube} holds {name}, "
                f"whose strongest band is {exp:.0f} cm\u207b\u00b9."
                + ("" if good else "  The spectrum is probably dominated by glass "
                   "or background rather than by the sample."),
                severity="info" if good else "warn",
            ))

        # The polymorphs differ by ~4 cm-1 at the CCN stretch, so record whether
        # the lattice region that actually separates them was captured.
        if key.endswith("glycine"):
            lattice = [b for b in DIAGNOSTIC[key]["key_bands"] if b < 300]
            got = [b for b in lattice if (spec.intensity_at(b) or 0) > 0.05 * peak]
            checks.append(Check(
                len(got) == len(lattice), "Lattice modes",
                f"Captured {len(got)} of {len(lattice)} lattice bands "
                f"({', '.join(str(b) for b in lattice)} cm\u207b\u00b9). These "
                "separate the polymorphs; the CCN stretch alone cannot.",
                severity="info",
            ))

    return Report(spectrum=spec, checks=checks)


def summarise(reports: list) -> pd.DataFrame:
    """One row per submission, for an at-a-glance class view."""
    rows = []
    for rep in reports:
        s, m = rep.spectrum, rep.spectrum.meta
        mat = s.material
        rows.append({
            "student": s.student,
            "scan": m.get("scan"),
            "tube": s.tube,
            "contents (key)": mat[1] if mat else "\u2014",
            "objective": m.get("objective"),
            "power %": m.get("laser_power_pct"),
            "accumulations": m.get("accumulations"),
            "S/N": m.get("estimated_SNR"),
            "band cm-1": s.observed_band(),
            "issues": len(rep.problems),
        })
    return pd.DataFrame(rows)


def answer_key(students: list) -> pd.DataFrame:
    """Which tube holds what, for each student. The instructor's crib sheet."""
    rows = []
    for st in students:
        a = assign_tubes(st)
        rows.append({"student": st,
                     **{label: a[label][1] for label, _, _ in TUBES}})
    return pd.DataFrame(rows)
