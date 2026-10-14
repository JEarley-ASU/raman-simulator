"""
Raman microspectrometer teaching simulator — Streamlit wrapper.

The instrument itself stays client-side. Canvas rendering, keyboard jogging and
live focus feedback all need sub-100 ms response, and routing those through
Streamlit's rerun cycle would make the instrument unusable. Streamlit's job here
is the part a static HTML file cannot do: assign identities from a roster and
collect what students produce.

raman-sim.html is the same file students can open on its own; there is no
Streamlit-specific copy to keep in step.

    streamlit run app.py
"""

from __future__ import annotations

import inspect
import json
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from labdata import (TUBES, answer_key, assign_tubes, ParseError,
                     parse_csv, review, summarise)

HERE = Path(__file__).parent
SIM = HERE / "raman-sim.html"   # the one and only instrument

# Streamlit renamed both of these and the old spellings are past their stated
# removal dates, but plenty of installs are still on older versions. Detect once
# and use whichever the running version provides.
_HAS_IFRAME = hasattr(st, "iframe")
_HAS_WIDTH = "width" in inspect.signature(st.dataframe).parameters


def _frame(html: str, height: int) -> None:
    if _HAS_IFRAME:
        st.iframe(html, height=height)
    else:
        components.html(html, height=height, scrolling=True)


def _table(df: pd.DataFrame) -> None:
    if _HAS_WIDTH:
        st.dataframe(df, width="stretch", hide_index=True)
    else:
        st.dataframe(df, use_container_width=True, hide_index=True)

st.set_page_config(page_title="Raman microspectrometer",
                   page_icon="🔬", layout="wide")


@st.cache_data
def load_sim() -> str:
    if not SIM.exists():
        st.error(f"Cannot find {SIM.name} next to app.py.")
        st.stop()
    return SIM.read_text(encoding="utf-8")


def embed(student: str, lock: bool, height: int, reveal: bool = False) -> None:
    """Drop the instrument in, with the operator identity pinned.

    The identity matters: it seeds the noise, the detector gain, the parfocality
    offsets and which capillary holds what, so a given ID reproduces its session
    exactly while no two students get the same slide.
    """
    preset = json.dumps({"student": student, "lock": lock, "reveal": reveal})
    html = load_sim().replace(
        "</head>", f"<script>window.LAB_PRESET={preset};</script></head>", 1
    )
    _frame(html, height)


# --------------------------------------------------------------------------
# Sidebar
# --------------------------------------------------------------------------
with st.sidebar:
    st.subheader("Session")
    student = st.text_input("Student ID", value="",
                            placeholder="e.g. your NetID",
                            help="Seeds the noise, detector gain and parfocality. "
                                 "Use the same ID every time to get the same slide "
                                 "back.").strip()
    height = st.slider("Instrument panel height (px)", 700, 1800, 1150, 50,
                       help="The instrument fills its frame. Raise this on a "
                            "tall monitor.")
    lock = st.checkbox("Lock the ID in the instrument", value=True,
                       help="Stops students editing the seed to fish for a "
                            "cleaner dataset.")
    reveal = st.checkbox("Name the tube contents", value=False,
                         help="Shows what each capillary holds in the readout. "
                              "The arrangement is shuffled per student ID.")
    if reveal and student:
        a = assign_tubes(student)
        st.caption("  ·  ".join(f"**{l}** {a[l][1]}" for l, _, _ in TUBES))
    st.divider()
    st.caption(
        "Everything in the instrument runs in the browser, so a full lab "
        "section costs the server nothing."
    )


st.title("Raman microspectrometer")

tab_run, tab_submit, tab_class = st.tabs(
    ["Instrument", "Submit a spectrum", "Class view"]
)

# --------------------------------------------------------------------------
# Instrument
# --------------------------------------------------------------------------
with tab_run:
    if not student:
        st.info("Enter your student ID in the sidebar to load the instrument.",
                icon="👈")
    else:
        st.caption(
            f"Operating as **{student}**. Five capillaries: three crystalline "
            "powders and two liquids. Find a tube at 5×, switch up to measure, "
            "frame and focus on the alignment beam, then acquire."
        )
        embed(student, lock, height, reveal)
        st.info(
            "Export with **Save spectrum (.csv)**. If your browser blocks the "
            "download inside this frame, use **Show CSV text** and copy it "
            "instead — the Submit tab accepts pasted text.",
            icon="💾",
        )

# --------------------------------------------------------------------------
# Submission
# --------------------------------------------------------------------------
def show_report(rep) -> None:
    spec = rep.spectrum
    m = spec.meta

    cols = st.columns(4)
    cols[0].metric("Objective", str(m.get("objective", "—")))
    cols[1].metric("Accumulations", f"{int(m['accumulations'])}"
                   if m.get("accumulations") else "—")
    cols[2].metric("S/N", f"{m['estimated_SNR']:.1f}"
                   if m.get("estimated_SNR") is not None else "—")
    band = spec.observed_band()
    cols[3].metric("Strongest band", f"{band:.0f} cm⁻¹" if band else "—")
    st.caption(f"Beam was in **{spec.tube}**.")

    st.line_chart(spec.data.set_index("wavenumber_cm-1")[["mean_counts"]],
                  height=320)

    problems = rep.problems
    if problems:
        st.warning(f"{len(problems)} thing(s) worth a second look before you "
                   "write this up.")
    else:
        st.success("Acquisition hygiene looks sound.")

    for c in rep.checks:
        icon = "✅" if c.ok else ("⚠️" if c.severity == "warn" else "ℹ️")
        st.markdown(f"{icon} **{c.label}** — {c.detail}")

    with st.expander("Full acquisition record"):
        _table(pd.DataFrame({"field": list(m.keys()),
                             "value": [str(v) for v in m.values()]}))
        st.caption(
            "The contents are deliberately absent from the export. Which tube "
            "the beam was in is inferred from the stage coordinates, and what "
            "that tube holds is recomputed from the student ID — so the file "
            "itself never gives the identification away."
        )


with tab_submit:
    st.write("Upload the CSV you exported, or paste the text from "
             "**Show CSV text**.")

    up = st.file_uploader("Spectrum export", type=["csv", "txt"])
    pasted = st.text_area("…or paste here", height=120,
                          placeholder="# Raman microspectrometer - simulated acquisition")

    raw = None
    if up is not None:
        raw = up.getvalue()
    elif pasted.strip():
        raw = pasted

    if raw is not None:
        try:
            rep = review(parse_csv(raw))
        except ParseError as exc:
            st.error(f"Could not read that file: {exc}")
        except Exception as exc:  # malformed numbers, truncated paste, etc.
            st.error(f"Unexpected problem reading that file: {exc}")
        else:
            if student and rep.spectrum.student and rep.spectrum.student != student:
                st.warning(
                    f"This file was recorded by **{rep.spectrum.student}**, but "
                    f"the session is set to **{student}**.",
                    icon="👤",
                )
            show_report(rep)
    else:
        st.caption("Nothing loaded yet.")

# --------------------------------------------------------------------------
# Class view
# --------------------------------------------------------------------------
with tab_class:
    st.write("Drop in a whole section's exports to see them side by side.")
    ups = st.file_uploader("Spectrum exports", type=["csv", "txt"],
                           accept_multiple_files=True, key="batch")

    if ups:
        reports, failed = [], []
        for f in ups:
            try:
                reports.append(review(parse_csv(f.getvalue())))
            except Exception as exc:
                failed.append((f.name, str(exc)))

        if reports:
            table = summarise(reports)
            _table(table)
            st.download_button(
                "Download summary (.csv)",
                table.to_csv(index=False).encode("utf-8"),
                file_name="raman_class_summary.csv",
                mime="text/csv",
            )

            clean = sum(1 for r in reports if not r.problems)
            st.caption(f"{clean} of {len(reports)} submissions passed every check.")

            with st.expander("Per-submission detail"):
                for rep in reports:
                    st.markdown(
                        f"**{rep.spectrum.student}** — scan "
                        f"{rep.spectrum.meta.get('scan')} in {rep.spectrum.tube}"
                    )
                    for c in rep.problems:
                        st.markdown(f"&nbsp;&nbsp;⚠️ {c.label} — {c.detail}")
                    if not rep.problems:
                        st.markdown("&nbsp;&nbsp;✅ no issues")

        for name, err in failed:
            st.error(f"{name}: {err}")
    else:
        st.caption("No files loaded.")

    st.divider()
    with st.expander("Answer key for a roster"):
        st.caption(
            "Each student gets a different arrangement of the five capillaries. "
            "Paste a roster, one ID per line, to recover who had what."
        )
        roster = st.text_area("Student IDs", height=120,
                              placeholder="student-01\nstudent-02\nstudent-03")
        ids = [r.strip() for r in roster.splitlines() if r.strip()]
        if ids:
            key = answer_key(ids)
            _table(key)
            st.download_button("Download answer key (.csv)",
                               key.to_csv(index=False).encode("utf-8"),
                               file_name="raman_answer_key.csv", mime="text/csv")
