# MQCT Input Generator (local GUI)

A small local web app for building MQCT `.inp` files for **molecule + atom**
systems (`SYS_TYPE` 1-4). Fill a form, get a valid input file.

## What it covers

| SYS_TYPE | System                    |
|----------|---------------------------|
| 1        | rigid diatom + atom       |
| 2        | vibrating diatom + atom   |
| 3        | symmetric top + atom      |
| 4        | asymmetric top + atom     |

The form is reactive: only the keywords relevant to the chosen `SYS_TYPE`
are shown. It builds the three blocks `$BASIS`, `$SYSTEM`, `$POTENTIAL`
and counts `NMB_CHNLS` / `NMB_ENERGS` automatically.

## Setup (macOS / Linux)

```bash
cd mqct_gui
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py
```

Then open <http://127.0.0.1:5000> in a browser.

To stop the server: Ctrl+C in the terminal.

## Usage

1. Pick the **system type** — the form reshapes itself.
2. Enter rotational/vibrational constants (cm-1).
3. Choose how to specify channels:
   - **Explicit list** — one tuple per line (`j` / `j,v` / `j,k,eps` / `j,ka,kc`)
   - **Energy cutoff** — `EMAX`
   - **j range** — `JMIN` / `JMAX` (+ `VMIN` / `VMAX` for type 2)
4. Fill `$SYSTEM` and `$POTENTIAL` fields.
5. Press **Generate**, then **Copy** or **Download**.

Place the downloaded `.inp` in your MQCT run directory and point
`INPUT_NAME.inp` at it.

## Notes / gotchas the form reminds you of

- **Prolate symmetric top**: enter A < C (opposite the usual convention).
- **Asymmetric top A,B,C**: must match the x,y,z axes of the PES reference
  orientation, not magnitude order.
- **U_ENERGY** is the effective (Billing) energy, not raw kinetic E.
- **B_IMPCT** overrides `JTOTL` / `JTOTU` if set.
- For diatom + atom (types 1,2) only the beta grid `GRD_ANG2` is used.

## Files

```
mqct_gui/
  app.py                 Flask backend + .inp generator
  templates/index.html   the form UI
  requirements.txt
  README.md
```

## Extending to molecule + molecule (types 5-0)

The keyword logic in `build_input()` is already structured by SYS_TYPE.
To add types 5-0 you would: add them to `SYS_TYPES` / `CHANNEL_LABELS`,
add the per-molecule constant fields (`BE1`, `BE2`, `A1`, `B1`, `C1`, ...),
and handle the two-valued `GRD_ANG*` keywords.
