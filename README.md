# Raman microspectrometer simulator — Streamlit app

    pip install -r requirements.txt
    streamlit run app.py

## Deploying to Streamlit Community Cloud

1. Push this folder to a GitHub repository.
2. At [share.streamlit.io](https://share.streamlit.io), **Create app** → pick the
   repo, branch `main`, main file `app.py`.
3. Deploy, then give students the URL. Each enters their first name in the
   sidebar.

## What students see

- **Instrument.** The simulated microspectrometer. Every finished scan is sent
  to the Analysis tab automatically; **Save spectrum (.csv)** still saves a file.
- **Analysis.** A checkbox table of every scan from the session, plus any CSVs
  they upload (for work from an earlier session). Ticked spectra are overlaid,
  and the processing applies to all of them:
  cosmic-ray removal, reference subtraction, baseline correction (iterative
  polynomial or asymmetric least squares), Savitzky–Golay smoothing,
  normalisation, stacked display, peak picking, band fitting (Lorentzian,
  Gaussian, pseudo-Voigt, with uncertainties), band integration, and CSV export
  of the processed spectra and every results table. A zip of the session's raw
  scans is at the bottom.

Scans live in the browser page. Reloading or closing it loses them, so students
should save CSVs (or the zip) before leaving.

## Names

The first name seeds detector gain, noise and focus offsets, so the same name
reproduces the same session and different names get different data. It is
case-insensitive. Two students with the same first name get identical data;
have one add an initial. Nothing is stored on the server.

## The slide

Five capillaries, numbered 1–5 left to right. The arrangement is the same for
everyone. The instrument never says which tube is under the beam, in the
readout, the log or the exports; students work it out by scanning across the
slide, and the files record only stage coordinates. The positions and contents,
for your reference:

| tube | stage X (mm) | contents |
|---|---|---|
| Tube 1 | −6.4 | α-glycine (powder) |
| Tube 2 | −3.2 | β-glycine (powder) |
| Tube 3 | 0.0 | γ-glycine (powder) |
| Tube 4 | 3.2 | CHCl₃ (liquid) |
| Tube 5 | 6.4 | CDCl₃ (liquid) |

To change the order, reorder `CONTENTS` in `instrument/index.html`; entry *k*
goes in tube *k*.

Things students can get out of it:

- **Polymorphs.** α and β glycine differ by ~4 cm⁻¹ at the CCN stretch (893 vs
  889), which this spectrometer cannot resolve. They separate on the lattice
  modes instead: α at 198, β at 215, γ at 251 cm⁻¹. The detector reaches down to
  150 cm⁻¹ so those are in range.
- **Isotope shift.** The C–H stretch at 3019 cm⁻¹ moves to 2256 on deuteration,
  a ratio of 1.338 against the harmonic prediction √(μ_D/μ_H) = 1.363. The gap
  is anharmonicity. CDCl₃ also carries a small residual 3019 band at ~99.8 atom % D.

## Files

| file | role |
|---|---|
| `instrument/index.html` | The instrument. One copy: open it directly in a browser to use it standalone, or let `app.py` serve it |
| `app.py` | Streamlit app: name, instrument, Analysis tab |
| `labdata.py` | Parsing and spectral processing. No Streamlit import, so it works from a notebook |
| `make_sprites.py` | Cuts crystal sprites out of powder micrographs for the renderer |

The instrument runs entirely in the browser; canvas rendering and keyboard
jogging need faster response than Streamlit's rerun cycle gives. It talks to the
app through Streamlit's component protocol (name in, scans out), written
directly in the HTML so there is no build step.
