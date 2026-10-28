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
have one add an initial. Nothing is stored on the server. Nothing the students
see says the data are generated or that the name controls them; the help text
only says the name labels their scans and files.

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

- **Polymorphs.** Glycine bands follow Table II of Surovtsev et al., *J. Chem.
  Phys.* **137**, 065103 (2012), shifted from 78 K to room temperature using the
  anharmonic fits in Surovtsev et al., *J. Chem. Phys.* **134**, 045102 (2011).
  The C–C stretch sits at ~893 cm⁻¹ in all three forms, so it cannot tell them
  apart. In the standard window (150–3200 cm⁻¹, edge filter) the CO₂ bend
  (698 / 706 / 687), C–N stretch (1034 / 1042 / 1050) and CH₂ stretch
  (2972 / 2951 / 2962) do. The clearest fingerprint is the lattice modes, which
  need the **Spectral range** control switched to the notch filters
  (10–610 cm⁻¹). Strongest lattice bands: α 53, 72, 111 cm⁻¹; β 84, 108,
  144, 156; γ 91, 105, 143, 153. The low window costs light per pixel (about
  5× longer exposures) and shows a laser wing below ~30 cm⁻¹, as on a real
  instrument.
- **Isotope shift.** The C–H stretch at 3019 cm⁻¹ moves to 2256 on deuteration,
  a ratio of 1.338 against the harmonic prediction √(μ_D/μ_H) = 1.363. The gap
  is anharmonicity. CDCl₃ also carries a small residual 3019 band at ~99.8 atom % D.
- **Depolarization ratios.** The laser is vertically polarized, and the
  **Analyzer** control sets Out, Parallel (V) or Perpendicular (H). With the
  analyzer out, spectra are exactly what they were before the analyzer was
  added. In the liquids, E bands have ρ = I⊥/I∥ = 0.75. Polarized A1 bands use
  the values in the `A1_RHO` block of `instrument/index.html`, which you can
  edit: CHCl₃ 366 0.12, 667 0.08, 3019 0.20; CDCl₃ 366 0.12, 651 0.08,
  2256 0.18, 3019 0.20. Peak-height ratios come out a little off these values
  where bands overlap (an E band next to an A1 band reads about 0.74), so a
  band fit is the better measurement. Fluorescence and white light are
  unpolarized. In the powders each grain is a single crystal in a random
  orientation, so a band's ratio changes when the stage moves more than about
  4 µm (`GRAIN_UM`) and says nothing about symmetry. That is the lesson:
  depolarization ratios assign symmetry only in an isotropic sample.
- **Anti-Stokes.** The third **Spectral range** setting, −300 to +300 cm⁻¹, has
  the same dispersion and light per pixel as the 10–610 window. The notch blocks
  only within ~14 cm⁻¹ of the laser on both sides, and a laser wing leaks
  through on each side. Each anti-Stokes band is its Stokes partner times
  ((ν₀+ν)/(ν₀−ν))⁴·exp(−hcν/kT), with ν₀ = 1/532 nm and T = 296 K
  (`CONFIG.sampleK`). There is no fluorescence on the anti-Stokes side. From the
  262 cm⁻¹ CHCl₃ pair (ratio ≈ 0.31) students can recover the sample
  temperature. The edge-filter window blocks everything below 150 cm⁻¹,
  anti-Stokes included.

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
