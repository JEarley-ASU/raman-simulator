"""
Raman microspectrometer teaching simulator — Streamlit wrapper.

The instrument itself stays client-side. Canvas rendering, keyboard jogging and
live focus feedback all need sub-100 ms response, and routing those through
Streamlit's rerun cycle would make the instrument unusable. The instrument
hands each finished scan back to the app, and the Analysis tab does the
processing a spectrometer's software package would.

instrument/index.html is the same file students can open on its own; there is
no Streamlit-specific copy to keep in step.

    streamlit run app.py
"""

from __future__ import annotations

import hashlib
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components
from plotly.subplots import make_subplots

from labdata import (SHAPES, ParseError, als_baseline, despike, fit_bands,
                     integrate, normalise, parse_csv, pick_peaks, poly_baseline,
                     smooth)

HERE = Path(__file__).parent
WN_MIN, WN_MAX = 150, 3200

st.set_page_config(page_title="Raman microspectrometer",
                   page_icon="🔬", layout="wide")

_instrument = components.declare_component(
    "raman_instrument", path=str(HERE / "instrument"))


@st.cache_data(show_spinner=False)
def _parse(text: str | bytes, source: str):
    return parse_csv(text, source)


# --------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------
with st.sidebar:
    st.subheader("Session")
    name = st.text_input(
        "First name", placeholder="e.g. Maya",
        help="Seeds the detector gain, noise and focus offsets. Use the same "
                 "name each time to get the same instrument back.").strip()
    st.divider()
    st.caption(
        "Scans you collect appear in the Analysis tab. They last until you "
        "close or reload this page, so save the CSVs you need."
    )

st.title("Raman microspectrometer")
tab_run, tab_an = st.tabs(["Instrument", "Analysis"])

# --------------------------------------------------------------------------
# Instrument
# --------------------------------------------------------------------------
scans = []
with tab_run:
    if not name:
        st.info("Enter your first name in the sidebar to load the instrument.",
                icon="👈")
    else:
        st.caption(
            "Five capillaries: three crystalline powders and two liquids. Find "
            "a tube at 5×, switch up to measure, frame and focus on the "
            "alignment beam, then acquire."
        )
        scans = _instrument(name=name, key="instrument", default=[]) or []


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------
def collect(uploads) -> dict:
    """Every spectrum available this run, keyed by a stable id."""
    out, bad = {}, []
    for s in scans:
        try:
            out[f"scan:{s['started']}"] = _parse(s["text"], "")
        except Exception as exc:  # a scan that failed to serialise
            bad.append((s.get("name", "scan"), str(exc)))
    for f in uploads or []:
        raw = f.getvalue()
        key = f"file:{f.name}:{hashlib.md5(raw).hexdigest()[:8]}"
        try:
            out[key] = _parse(raw, f.name)
        except (ParseError, Exception) as exc:
            bad.append((f.name, str(exc)))
    for fname, err in bad:
        st.error(f"Could not read {fname}: {err}")
    return out


def pick(spectra: dict) -> list:
    """Checkbox table. Selection lives in session state so it survives new
    scans arriving; the editor key changes with the row set so stale edits
    never land on the wrong row."""
    sel = st.session_state.setdefault("selected", {})
    for k in spectra:
        sel.setdefault(k, True)

    c1, c2, _ = st.columns([1, 1, 6])
    if c1.button("Select all"):
        sel.update({k: True for k in spectra})
        st.session_state.picker_rev = st.session_state.get("picker_rev", 0) + 1
    if c2.button("Select none"):
        sel.update({k: False for k in spectra})
        st.session_state.picker_rev = st.session_state.get("picker_rev", 0) + 1

    ids = list(spectra)
    m = [spectra[k].meta for k in ids]
    df = pd.DataFrame({
        "Plot": [sel[k] for k in ids],
        "Spectrum": [spectra[k].label for k in ids],
        "Source": ["this session" if k.startswith("scan:") else "uploaded"
                   for k in ids],
        "Power %": [x.get("laser_power_pct") for x in m],
        "Time (s)": [x.get("integration_s") for x in m],
        "Acc.": [x.get("accumulations") for x in m],
        "S/N": [x.get("estimated_SNR") for x in m],
        "X (mm)": [x.get("stage_x_mm") for x in m],
        "Y (mm)": [x.get("stage_y_mm") for x in m],
    }, index=ids)
    rev = st.session_state.get("picker_rev", 0)
    edited = st.data_editor(
        df, hide_index=True, width="stretch",
        disabled=[c for c in df.columns if c != "Plot"],
        column_config={"Plot": st.column_config.CheckboxColumn(width="small")},
        key=f"picker-{rev}-{hashlib.md5('|'.join(ids).encode()).hexdigest()}",
    )
    for k, v in zip(ids, edited["Plot"]):
        sel[k] = bool(v)
    return [k for k in ids if sel[k]]


def controls(spectra: dict, chosen: list) -> dict:
    o = {}
    with st.expander("Range", expanded=True):
        o["range"] = st.slider("Raman shift (cm⁻¹)", WN_MIN, WN_MAX,
                               (WN_MIN, WN_MAX), 10)
    with st.expander("Cosmic rays"):
        o["despike"] = st.checkbox("Remove cosmic-ray spikes", value=False)
        o["spike_thr"] = st.slider("Threshold (× noise)", 3.0, 20.0, 6.0, 0.5,
                                   disabled=not o["despike"],
                                   help="Lower catches smaller spikes but "
                                        "risks clipping sharp bands.")
    with st.expander("Subtract a reference"):
        labels = {k: spectra[k].label for k in spectra}
        ref = st.selectbox("Reference spectrum", ["None"] + list(spectra),
                           format_func=lambda k: labels.get(k, k),
                           help="For example an empty capillary or glass, "
                                "to remove its bands from the others.")
        o["ref"] = None if ref == "None" else ref
        o["ref_scale"] = st.number_input("Scale factor", 0.0, 10.0, 1.0, 0.05,
                                         disabled=o["ref"] is None)
    with st.expander("Baseline", expanded=True):
        o["base"] = st.radio("Method", ["None", "Polynomial", "Asymmetric least squares"],
                             index=0)
        if o["base"] == "Polynomial":
            o["poly_order"] = st.slider("Order", 1, 8, 3)
        elif o["base"] == "Asymmetric least squares":
            o["als_lam"] = 10 ** st.slider("Stiffness, log₁₀ λ", 2.0, 9.0, 6.0, 0.5)
            o["als_p"] = st.select_slider("Asymmetry p",
                                          [0.001, 0.002, 0.005, 0.01, 0.02,
                                           0.05, 0.1], value=0.01)
        o["show_base"] = st.checkbox("Preview baseline on the raw data",
                                     disabled=o["base"] == "None")
    with st.expander("Smoothing"):
        o["smooth"] = st.checkbox("Savitzky–Golay", value=False)
        o["sg_win"] = st.slider("Window (points)", 5, 41, 9, 2,
                                disabled=not o["smooth"])
        o["sg_ord"] = st.slider("Polynomial order", 2, 5, 3,
                                disabled=not o["smooth"])
    with st.expander("Normalise"):
        mode = st.radio("Scale each spectrum to", ["None", "Maximum = 1",
                                                   "Area = 1", "A chosen band = 1"])
        o["norm"] = {"None": None, "Maximum = 1": "max", "Area = 1": "area",
                     "A chosen band = 1": "band"}[mode]
        o["band"] = None
        if o["norm"] == "band":
            o["band"] = st.slider("Band window (cm⁻¹)", WN_MIN, WN_MAX, (860, 920), 2)
    with st.expander("Display", expanded=True):
        o["stack"] = st.checkbox("Stack with an offset", value=False)
        o["offset"] = st.slider("Offset (fraction of tallest)", 0.0, 1.5, 0.3, 0.05,
                                disabled=not o["stack"])
    return o


def process(spec, o: dict, spectra: dict):
    """Run one spectrum through the chain. Returns x, processed y, and the
    pre-baseline y with its baseline for the preview."""
    x, y = spec.x, spec.y.copy()
    if o["despike"]:
        y = despike(y, o["spike_thr"])
    if o["ref"] is not None and spectra[o["ref"]] is not spec:
        r = spectra[o["ref"]]
        ry = despike(r.y, o["spike_thr"]) if o["despike"] else r.y
        y = y - o["ref_scale"] * np.interp(x, r.x, ry)

    lo, hi = o["range"]
    keep = (x >= lo) & (x <= hi)
    x, y = x[keep], y[keep]
    if x.size < 10:
        return x, y, y, None

    base = None
    if o["base"] == "Polynomial":
        base = poly_baseline(x, y, o["poly_order"])
    elif o["base"] == "Asymmetric least squares":
        base = als_baseline(y, o["als_lam"], o["als_p"])
    pre = y
    if base is not None:
        y = y - base
    if o["smooth"]:
        y = smooth(y, o["sg_win"], o["sg_ord"])
    if o["norm"]:
        y = normalise(x, y, o["norm"], o["band"])
    return x, y, pre, base


def ylabel(o: dict) -> str:
    return "Normalised intensity" if o["norm"] else "Intensity (counts)"


with tab_an:
    st.write("Tick the spectra to overlay. Scans from this session appear "
             "automatically; upload saved CSVs to bring back earlier work.")
    uploads = st.file_uploader("Saved spectra (.csv)", type=["csv", "txt"],
                               accept_multiple_files=True)
    spectra = collect(uploads)

    if not spectra:
        st.caption("No spectra yet. Acquire one on the Instrument tab or upload "
                   "a CSV.")
    else:
        chosen = pick(spectra)
        left, right = st.columns([1, 3], gap="medium")
        with left:
            o = controls(spectra, chosen)

        done = {k: process(spectra[k], o, spectra) for k in chosen}

        with right:
            fig = go.Figure()
            preview = o["show_base"] and o["base"] != "None"
            tallest = max((float(np.max(v[1])) for v in done.values()
                           if v[1].size), default=1.0)
            offsets = {}
            for i, (k, (x, y, pre, base)) in enumerate(done.items()):
                lab = spectra[k].label
                if preview:
                    fig.add_scatter(x=x, y=pre, name=lab, mode="lines")
                    if base is not None:
                        fig.add_scatter(x=x, y=base, name=f"{lab} baseline",
                                        mode="lines", line=dict(dash="dash"),
                                        showlegend=False)
                    continue
                off = i * o["offset"] * tallest if o["stack"] else 0.0
                offsets[k] = off
                fig.add_scatter(x=x, y=y + off, name=lab, mode="lines")

            peak_tables = {}
            with st.expander("Peak picking", expanded=False):
                pk_on = st.checkbox("Mark peaks", value=False)
                c1, c2 = st.columns(2)
                prom = c1.slider("Minimum height (× noise)", 3, 50, 6,
                                 help="How far a peak must stand above its "
                                      "surroundings, in units of the "
                                      "spectrum's own noise.")
                sep = c2.slider("Minimum separation (cm⁻¹)", 3, 100, 10)
            if pk_on and not preview:
                for k, (x, y, _, _) in done.items():
                    t = pick_peaks(x, y, prom, sep)
                    peak_tables[k] = t
                    if len(t):
                        fig.add_scatter(
                            x=t["shift_cm-1"], y=t["height"] + offsets[k],
                            mode="markers+text", text=t["shift_cm-1"].round(0).astype(int),
                            textposition="top center", marker=dict(size=6),
                            showlegend=False, hoverinfo="x")

            fig.update_layout(
                height=560, margin=dict(l=10, r=10, t=30, b=10),
                xaxis_title="Raman shift (cm⁻¹)",
                yaxis_title=("Intensity (counts)" if preview else ylabel(o))
                            + (" + offset" if o["stack"] and not preview else ""),
                legend=dict(orientation="h", y=-0.15), hovermode="closest",
            )
            if not chosen:
                st.caption("Tick at least one spectrum above to plot it.")
            st.plotly_chart(fig)

            if peak_tables:
                rows = [t.assign(spectrum=spectra[k].label) for k, t in peak_tables.items()]
                allpk = pd.concat(rows, ignore_index=True)[
                    ["spectrum", "shift_cm-1", "height", "prominence", "FWHM_cm-1"]]
                st.dataframe(allpk, hide_index=True, width="stretch")
                st.download_button("Download peak table (.csv)",
                                   allpk.to_csv(index=False).encode(),
                                   "peaks.csv", "text/csv")

        if chosen:
            t_fit, t_int, t_exp = st.tabs(["Band fitting", "Integration", "Export"])

            with t_fit:
                c1, c2, c3, c4 = st.columns([2, 2, 1, 1])
                target = c1.selectbox("Spectrum", chosen,
                                      format_func=lambda k: spectra[k].label)
                lo, hi = c2.slider("Fit window (cm⁻¹)", WN_MIN, WN_MAX, (840, 940), 2)
                shape = c3.selectbox("Line shape", SHAPES)
                nb = c4.number_input("Bands", 1, 6, 1)
                guess_txt = st.text_input(
                    "Starting centres, cm⁻¹ (optional, comma-separated)",
                    placeholder="e.g. 889, 893",
                    help="Useful when two bands overlap. Left blank, the fit "
                         "starts from the strongest maxima in the window.")
                try:
                    guesses = [float(g) for g in guess_txt.replace(";", ",").split(",")
                               if g.strip()]
                except ValueError:
                    st.error("Starting centres must be numbers.")
                    guesses = []
                x, y, _, _ = done[target]
                try:
                    fit = fit_bands(x, y, lo, hi, shape, int(nb), guesses)
                except Exception as exc:
                    st.warning(f"The fit did not converge: {exc}")
                else:
                    ff = make_subplots(rows=2, cols=1, shared_xaxes=True,
                                       row_heights=[0.75, 0.25], vertical_spacing=0.04)
                    ff.add_scatter(x=fit.x, y=fit.y, name="data", mode="markers",
                                   marker=dict(size=4), row=1, col=1)
                    ff.add_scatter(x=fit.x, y=fit.total, name="fit", mode="lines",
                                   row=1, col=1)
                    for i, comp in enumerate(fit.components, 1):
                        ff.add_scatter(x=fit.x, y=comp + fit.background,
                                       name=f"band {i}", mode="lines",
                                       line=dict(dash="dot"), row=1, col=1)
                    ff.add_scatter(x=fit.x, y=fit.background, name="background",
                                   mode="lines", line=dict(dash="dash"), row=1, col=1)
                    ff.add_scatter(x=fit.x, y=fit.y - fit.total, name="residual",
                                   mode="lines", row=2, col=1)
                    ff.update_layout(height=520, margin=dict(l=10, r=10, t=30, b=10),
                                     legend=dict(orientation="h", y=-0.12))
                    ff.update_xaxes(title_text="Raman shift (cm⁻¹)", row=2, col=1)
                    ff.update_yaxes(title_text=ylabel(o), row=1, col=1)
                    ff.update_yaxes(title_text="residual", row=2, col=1)
                    st.plotly_chart(ff)
                    st.caption(f"R² = {fit.r2:.4f}. Uncertainties are one standard "
                               "error from the fit covariance.")
                    st.dataframe(fit.table.round(3), hide_index=True, width="stretch")
                    st.download_button("Download fit results (.csv)",
                                       fit.table.to_csv(index=False).encode(),
                                       "fit.csv", "text/csv")

            with t_int:
                c1, c2 = st.columns([3, 1])
                ilo, ihi = c1.slider("Band limits (cm⁻¹)", WN_MIN, WN_MAX, (860, 920), 2)
                local = c2.checkbox("Straight-line baseline between limits",
                                    value=True)
                rows = [{"spectrum": spectra[k].label,
                         "area": integrate(x, y, ilo, ihi, local)}
                        for k, (x, y, _, _) in done.items()]
                it = pd.DataFrame(rows)
                st.dataframe(it, hide_index=True, width="stretch")
                st.download_button("Download areas (.csv)",
                                   it.to_csv(index=False).encode(),
                                   "areas.csv", "text/csv")

            with t_exp:
                frames = [pd.Series(y, index=np.round(x, 2), name=spectra[k].label)
                          for k, (x, y, _, _) in done.items()]
                wide = pd.concat(frames, axis=1).sort_index()
                wide.index.name = "wavenumber_cm-1"
                st.download_button("Processed spectra, as plotted (.csv)",
                                   wide.to_csv().encode(), "processed_spectra.csv",
                                   "text/csv")
                st.caption("One column per ticked spectrum, after every "
                           "processing step above. Stacking offsets are not "
                           "included.")

        if scans:
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
                for s in scans:
                    z.writestr(s["name"], s["text"])
            st.download_button(
                f"Download all {len(scans)} raw scans from this session (.zip)",
                buf.getvalue(), f"{name or 'raman'}_scans.zip", "application/zip")
