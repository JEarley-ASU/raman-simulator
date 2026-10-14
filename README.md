# Raman microspectrometer simulator — Streamlit wrapper

    pip install -r requirements.txt
    streamlit run app.py

## Deploying to Streamlit Community Cloud

1. Push this folder to a GitHub repository (`.gitignore` keeps
   `.streamlit/secrets.toml` out of it).
2. At [share.streamlit.io](https://share.streamlit.io), **Create app** → pick the
   repo, branch `main`, main file `app.py`.
3. Under **Advanced settings → Secrets**, paste
   `INSTRUCTOR_PASSWORD = "your-password"`, then deploy.
4. Give students the URL. Each enters their own ID in the sidebar.

The public app shows students only the instrument and the Submit tab. The tube
reveal, the ID-lock toggle, the Class view and the answer key appear after
signing in under **Instructor** in the sidebar. For students, Submit feedback
leaves out the checks that name a tube's contents or list a polymorph's bands.
Without the secret set, instructor tools are disabled entirely. For local use,
put the same line in `.streamlit/secrets.toml`.

The shuffle runs in the page's JavaScript, so a determined student reading the
source could work out their arrangement. The gate stops casual lookup, not that.

## What runs where

The instrument runs entirely in the browser. Canvas rendering, keyboard jogging
and live focus feedback need sub-100 ms response, and routing those through
Streamlit's rerun cycle would make it unusable. Streamlit assigns identities and
collects submissions; it does no simulation work, so a full lab section costs the
server essentially nothing.

## Files

| file | role |
|---|---|
| `raman-sim.html` | The instrument. **One copy, used two ways** — open it directly in a browser, or let `app.py` embed it |
| `app.py` | Streamlit wrapper: roster, embedding, submission review |
| `labdata.py` | Parsing and acquisition checks. No Streamlit import, so it is usable from a notebook or a grading script |
| `make_sprites.py` | Cuts crystal sprites out of powder micrographs for the renderer |

There is deliberately no separate Streamlit build of the instrument. When the
wrapper embeds it, it sets `window.LAB_PRESET` to pin the student ID; run
standalone that object is absent and the instrument uses the ID in its own
field. Keeping two copies in step by hand is the kind of thing that works until
it quietly doesn't, and `labdata.py` would only agree with one of them.

## Student identity

The ID seeds the noise, detector gain and parfocality offsets, so a given ID
reproduces its session exactly while no two students get identical data. Set it
in the sidebar; leave **Lock the ID** on so students cannot edit the seed to fish
for a cleaner dataset.

## The slide

Five capillary tubes, numbered left to right, at stage X = −6.4, −3.2, 0, 3.2,
6.4 mm. Three hold glycine polymorphs (α, β, γ); two hold CHCl₃ and CDCl₃.

**Which tube holds what is shuffled per student**, seeded from the ID, so no two
students get the same slide. The tubes are visually identical apart from powder
versus liquid — only the spectrum distinguishes them, which is the exercise.

Two things students can get out of it:

- **Polymorphs.** α and β glycine differ by ~4 cm⁻¹ at the CCN stretch (893 vs
  889), which this spectrometer cannot resolve. They separate on the lattice
  modes instead: α at 198, β at 215, γ at 251 cm⁻¹. The detector reaches down to
  150 cm⁻¹ so those are in range.
- **Isotope shift.** The C–H stretch at 3019 cm⁻¹ moves to 2256 on deuteration,
  a ratio of 1.338 against the harmonic prediction √(μ_D/μ_H) = 1.363. The gap
  is anharmonicity. CDCl₃ also carries a small residual 3019 band at ~99.8 atom % D.

## Grading

The export omits both which tube was measured and what was in it. `labdata.py`
infers the tube from the stage coordinates and recomputes the contents from the
student ID, so the file never gives the identification away.

Checks applied to each submission:

- excitation beam used, not the alignment attenuator
- white-light illuminator closed
- objective adequate for measurement (5× is for finding, not measuring)
- beam inside the bore rather than on the capillary wall
- enough accumulations, and adequate S/N
- detector not saturated
- strongest band consistent with what that tube holds for that student
- for glycine, whether the lattice modes that separate the polymorphs were captured

The **Class view** tab has an answer key: paste a roster and get each student's
arrangement. Or from Python:

```python
from labdata import answer_key
answer_key(["student-01", "student-02"]).to_csv("key.csv", index=False)
```

With 24 students drawing from 120 possible arrangements, expect a couple of
pairs to share one. That is harmless for marking but worth knowing if you seat
students together.

Use `labdata.py` directly for offline marking:

```python
from labdata import parse_csv, review, summarise
reports = [review(parse_csv(open(p).read())) for p in paths]
summarise(reports).to_csv("marks.csv", index=False)
```

## Editing the slide

Tube positions and contents live in **both** `raman-sim.html` (`TUBE`,
`SAMPLES`, `CONTENTS`) and `labdata.py` (`BORE_W`, `TUBES`, `CONTENTS`). Keep
them in step, or grading will decode the wrong material. The order of `CONTENTS`
matters on both sides — the shuffle depends on it.

This is the one duplication left, and it is unavoidable in this shape: one side
is JavaScript in the browser, the other Python on the server. If you expect to
edit the slide often, have both read a shared JSON file instead.

`labdata.py` ports the instrument's RNG so the shuffle reproduces exactly; that
equivalence is worth re-checking if you touch either copy.

## Downloads inside the frame

On Streamlit versions with `st.iframe` the instrument is embedded with
JavaScript and same-origin access, so the CSV download works normally. Older
versions fall back to `components.v1.html`, whose sandbox can block downloads —
the instrument has a **Show CSV text** button for that case, and the Submit tab
accepts pasted text.

`app.py` detects which API the installed version provides, so it runs on both.
