"""
MQCT Input File Generator -- local Flask GUI
Covers all ten system types (SYS_TYPE = 1..9, 0):
  molecule + atom      (1-4)
  molecule + molecule  (5-8, 9, 0)

Run:
    pip install flask
    python app.py
Then open http://127.0.0.1:5000 in a browser.
"""

from flask import Flask, render_template, request, Response, jsonify

app = Flask(__name__)

# ---------------------------------------------------------------------------
# System types. Order in dict controls the dropdown order.
# ---------------------------------------------------------------------------
SYS_TYPES = {
    "1": "Rigid diatom + atom",
    "2": "Vibrating diatom + atom",
    "3": "Symmetric top + atom",
    "4": "Asymmetric top + atom",
    "5": "Rigid diatom + rigid diatom",
    "6": "Vibrating diatom + vibrating diatom",
    "7": "Symmetric top + rigid diatom",
    "8": "Asymmetric top + rigid diatom",
    "9": "Asymmetric top + symmetric top",
    "0": "Asymmetric top + asymmetric top",
}

# Quantum-number labels that make up one channel entry, per SYS_TYPE.
CHANNEL_LABELS = {
    "1": ["j"],
    "2": ["j", "v"],
    "3": ["j", "k", "eps"],
    "4": ["j", "ka", "kc"],
    "5": ["j1", "j2"],
    "6": ["j1", "v1", "j2", "v2"],
    "7": ["j1", "k1", "eps1", "j2"],
    "8": ["j1", "ka1", "kc1", "j2"],
    "9": ["j1", "ka1", "kc1", "j2", "k2", "eps2"],
    "0": ["j1", "ka1", "kc1", "j2", "ka2", "kc2"],
}

# Sets of SYS_TYPEs grouped by partner geometry, for readability below.
MOL_ATOM = {"1", "2", "3", "4"}
MOL_MOL = {"5", "6", "7", "8", "9", "0"}
HAS_VIB = {"2", "6"}                 # vibrational constants / grids needed
IDENTICAL_OK = {"5", "6", "0"}       # identical-partner option available


def num(value, default=None):
    """Return a stripped string if non-empty, else default."""
    if value is None:
        return default
    value = value.strip()
    return value if value else default


def add(lst, form, key, kw=None):
    """Append 'KW=value' to lst if the form field is non-empty."""
    v = num(form.get(key))
    if v is not None:
        lst.append(f"{kw or key.upper()}={v}")


def parse_channels(form):
    """Return list of channel tuples, each a list of token strings."""
    raw = num(form.get("chnls_list"), "")
    tuples = [t.strip() for t in raw.replace("\n", ";").split(";") if t.strip()]
    out = []
    for t in tuples:
        parts = [p for p in t.replace(",", " ").split() if p]
        if parts:
            out.append(parts)
    return out


def is_number(s):
    """True if s parses as a float."""
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def is_int(s):
    """True if s parses as an integer."""
    try:
        int(s)
        return True
    except (TypeError, ValueError):
        return False


def validate(form):
    """
    Check the submitted form for problems.
    Returns (errors, warnings) -- two lists of human-readable strings.
    Errors would make MQCT crash or run wrong physics.
    Warnings are probably-wrong but possibly intentional.
    """
    errors, warnings = [], []
    st = form.get("sys_type", "1")
    width = len(CHANNEL_LABELS.get(st, ["j"]))

    # ---- required scalar fields ----
    if not num(form.get("mass_red")):
        errors.append("MASS_RED is required (reduced mass of the pair, amu).")
    elif not is_number(form.get("mass_red")):
        errors.append("MASS_RED must be a number.")

    energies = num(form.get("u_energy"), "")
    elist = [e for e in energies.replace(",", " ").split() if e]
    if not elist:
        errors.append("At least one collision energy (U_ENERGY) is required.")
    else:
        for e in elist:
            if not is_number(e):
                errors.append(f"U_ENERGY value '{e}' is not a number.")
            elif float(e) <= 0:
                errors.append(f"U_ENERGY value '{e}' must be positive.")

    grd_r = num(form.get("grd_r"))
    if not grd_r:
        errors.append("GRD_R (number of R-grid points) is required.")
    elif not is_int(grd_r) or int(grd_r) < 2:
        errors.append("GRD_R must be an integer of at least 2.")

    # ---- rotational constants present, per system type ----
    def need(field, label):
        if not num(form.get(field)):
            errors.append(f"{label} is required for SYS_TYPE={st}.")
        elif not is_number(form.get(field)):
            errors.append(f"{label} must be a number.")

    if st in ("1", "2"):
        need("be", "BE")
    if st == "2":
        need("we", "WE")
    if st in ("3", "4"):
        need("a", "A")
        need("c", "C")
        if st == "4":
            need("b", "B")
    if st in ("5", "6"):
        need("be1", "BE1")
        need("be2", "BE2")
    if st == "6":
        need("we1", "WE1")
        need("we2", "WE2")
    if st in ("7", "8"):
        need("a_td", "A")
        need("c_td", "C")
        if st == "8":
            need("b_td", "B")
        need("be2", "BE2")
    if st in ("9", "0"):
        need("a1", "A1")
        need("b1", "B1")
        need("c1", "C1")
        need("a2", "A2")
        need("c2", "C2")
        if st == "0":
            need("b2", "B2")

    # ---- prolate / oblate check for symmetric tops ----
    # MQCT infers the shape from the A vs C ordering:
    #   oblate  -> A > C   (standard)
    #   prolate -> A < C   (opposite the usual spectroscopic convention)
    # The form supplies 'top_shape' so we can verify the ordering matches.
    def shape_check(a_field, c_field, who):
        a, c = num(form.get(a_field)), num(form.get(c_field))
        if not (a and c and is_number(a) and is_number(c)):
            return
        a, c = float(a), float(c)
        if a == c:
            warnings.append(f"{who}: A equals C -- this is a spherical top. "
                            "Allowed, but check this is intended.")
            return
        shape = form.get("top_shape", "oblate")
        if shape == "oblate" and a < c:
            errors.append(f"{who}: declared oblate but A < C. An oblate top "
                           "needs A > C. Either swap the values or set the "
                           "shape to prolate.")
        if shape == "prolate" and a > c:
            errors.append(f"{who}: declared prolate but A > C. In MQCT a "
                           "prolate top needs A < C (opposite the usual "
                           "convention). Swap the values or set shape to "
                           "oblate.")
    if st == "3":
        shape_check("a", "c", "Symmetric top")
    if st == "7":
        shape_check("a_td", "c_td", "Symmetric top")

    # ---- channels ----
    chan_mode = form.get("chan_mode", "list")
    init_raw = num(form.get("init_chnl"))
    init = None
    if init_raw:
        init = [p for p in init_raw.replace(",", " ").split() if p]

    if chan_mode == "list":
        chans = parse_channels(form)
        if not chans:
            errors.append("No channels listed. Add at least one channel, "
                          "or switch to the EMAX / j-range mode.")
        for i, c in enumerate(chans, 1):
            if len(c) != width:
                lbls = ", ".join(CHANNEL_LABELS[st])
                errors.append(f"Channel #{i} '{','.join(c)}' has {len(c)} "
                              f"value(s); SYS_TYPE={st} needs {width} "
                              f"({lbls}).")
            for tok in c:
                if not is_int(tok):
                    errors.append(f"Channel #{i}: '{tok}' is not an integer.")
        # INIT_CHNL must be one of the listed channels
        if init is not None and chans:
            if len(init) != width:
                errors.append(f"INIT_CHNL has {len(init)} value(s); "
                              f"SYS_TYPE={st} needs {width}.")
            elif init not in chans:
                errors.append("INIT_CHNL is not among the listed channels. "
                              "It must match one CHNLS_LIST entry exactly.")
    elif chan_mode == "emax":
        if st in MOL_MOL:
            if not (num(form.get("emax1")) and num(form.get("emax2"))):
                errors.append("EMAX1 and EMAX2 are both required in "
                              "energy-cutoff mode for molecule + molecule.")
        else:
            if not num(form.get("emax")):
                errors.append("EMAX is required in energy-cutoff mode.")
        if init is not None and len(init) != width:
            errors.append(f"INIT_CHNL has {len(init)} value(s); "
                          f"SYS_TYPE={st} needs {width}.")
    elif chan_mode == "jrange":
        if st in MOL_MOL:
            req = ["jmin1", "jmax1", "jmin2", "jmax2"]
        else:
            req = ["jmin", "jmax"]
        missing = [r.upper() for r in req if num(form.get(r)) is None]
        if missing:
            errors.append("j-range mode needs: " + ", ".join(missing) + ".")
        if init is not None and len(init) != width:
            errors.append(f"INIT_CHNL has {len(init)} value(s); "
                          f"SYS_TYPE={st} needs {width}.")

    if init is None:
        warnings.append("INIT_CHNL not set -- MQCT needs an initial state. "
                        "Add it unless you set it elsewhere.")

    # ---- R range ----
    rmin, rmax = num(form.get("rmin")), num(form.get("rmax"))
    if rmin and rmax and is_number(rmin) and is_number(rmax):
        if float(rmin) >= float(rmax):
            errors.append("RMIN must be smaller than RMAX.")
    if not rmin or not rmax:
        warnings.append("RMIN / RMAX not both set -- trajectories need a "
                        "defined R range (typical RMAX ~20 Bohr).")

    # ---- partial-wave sampling ----
    if not num(form.get("b_impct")):
        if num(form.get("jtotu")) is None:
            warnings.append("Neither B_IMPCT nor JTOTU set -- one is needed "
                            "to bound the partial-wave sampling.")

    # ---- MQCT hard rules from the manual ----
    if form.get("monte_carlo") and form.get("diff_cross"):
        errors.append("MONTE_CARLO cannot be combined with differential "
                      "cross-section calculations.")
    if (st in IDENTICAL_OK and form.get("identical")
            and num(form.get("dl"))):
        dl = form.get("dl")
        if is_int(dl) and int(dl) % 2 == 0:
            errors.append("With IDENTICAL=YES, DL must be odd so that both "
                          "even and odd l values are sampled.")

    # ---- time step ----
    if not num(form.get("time_step")):
        warnings.append("TIME_STEP not set -- propagation needs a step size; "
                        "check convergence of energy/norm.")

    return errors, warnings


def build_input(form):
    """Assemble the MQCT .inp text from submitted form fields."""
    st = form.get("sys_type", "1")
    lines = []

    # ======================= $BASIS =======================
    basis = [f"SYS_TYPE={st}"]

    # --- rotational / vibrational constants, by system type ---
    if st in ("1", "2"):
        add(basis, form, "be", "BE")
        add(basis, form, "de", "DE")
    if st == "2":
        add(basis, form, "we", "WE")
        add(basis, form, "xe", "XE")
    if st in ("3", "4"):
        # Symmetric top (3): B is implied (B = A) and must NOT be written.
        # Asymmetric top (4): all three of A, B, C are written.
        add(basis, form, "a", "A")
        if st == "4":
            add(basis, form, "b", "B")
        add(basis, form, "c", "C")
    if st in ("5", "6"):
        add(basis, form, "be1", "BE1")
        add(basis, form, "de1", "DE1")
        add(basis, form, "be2", "BE2")
        add(basis, form, "de2", "DE2")
    if st == "6":
        add(basis, form, "we1", "WE1")
        add(basis, form, "xe1", "XE1")
        add(basis, form, "we2", "WE2")
        add(basis, form, "xe2", "XE2")
    if st in ("7", "8"):
        # molecule 1 is a top (A,B,C); molecule 2 is a rigid diatom (BE2,DE2)
        add(basis, form, "a_td", "A")
        if st == "8":
            add(basis, form, "b_td", "B")
        add(basis, form, "c_td", "C")
        add(basis, form, "be2", "BE2")
        add(basis, form, "de2", "DE2")
    if st in ("9", "0"):
        add(basis, form, "a1", "A1")
        add(basis, form, "b1", "B1")
        add(basis, form, "c1", "C1")
        add(basis, form, "a2", "A2")
        if st == "0":
            add(basis, form, "b2", "B2")
        add(basis, form, "c2", "C2")

    # --- channels: explicit list / energy cutoff / j-range ---
    chan_mode = form.get("chan_mode", "list")
    if chan_mode == "list":
        raw = num(form.get("chnls_list"), "")
        tuples = [t.strip() for t in raw.replace("\n", ";").split(";") if t.strip()]
        flat = []
        for t in tuples:
            parts = [p for p in t.replace(",", " ").split() if p]
            if parts:
                flat.append(",".join(parts))
        if flat:
            basis.append(f"NMB_CHNLS={len(flat)}")
            basis.append("CHNLS_LIST=" + ", ".join(flat))
    elif chan_mode == "emax":
        if st in MOL_MOL:
            add(basis, form, "emax1", "EMAX1")
            add(basis, form, "emax2", "EMAX2")
        else:
            add(basis, form, "emax", "EMAX")
    elif chan_mode == "jrange":
        if st in MOL_MOL:
            add(basis, form, "jmin1", "JMIN1")
            add(basis, form, "jmax1", "JMAX1")
            add(basis, form, "jmin2", "JMIN2")
            add(basis, form, "jmax2", "JMAX2")
            if st == "6":
                add(basis, form, "vmin1", "VMIN1")
                add(basis, form, "vmax1", "VMAX1")
                add(basis, form, "vmin2", "VMIN2")
                add(basis, form, "vmax2", "VMAX2")
        else:
            add(basis, form, "jmin", "JMIN")
            add(basis, form, "jmax", "JMAX")
            if st == "2":
                add(basis, form, "vmin", "VMIN")
                add(basis, form, "vmax", "VMAX")

    init = num(form.get("init_chnl"))
    if init:
        init = ",".join(p for p in init.replace(",", " ").split() if p)
        basis.append(f"INIT_CHNL={init}")

    # --- optional flags ---
    if form.get("cs_approx"):
        basis.append("CS_APPROX=YES")
    if form.get("symmetry"):
        basis.append("SYMMETRY=YES")
    if st in IDENTICAL_OK and form.get("identical"):
        basis.append("IDENTICAL=YES")
    if form.get("print_states"):
        basis.append("PRINT_STATES=YES")

    lines.append("$BASIS")
    lines.append("  " + ",\n  ".join(basis))
    lines.append("$END")
    lines.append("")

    # ======================= $SYSTEM =======================
    system = []
    label = num(form.get("label"))
    if label:
        system.append(f'LABEL="{label}"')
    add(system, form, "mass_red", "MASS_RED")
    add(system, form, "rmin", "RMIN")
    add(system, form, "rmax", "RMAX")

    energies = num(form.get("u_energy"), "")
    elist = [e.strip() for e in energies.replace(",", " ").split() if e.strip()]
    if elist:
        system.append(f"NMB_ENERGS={len(elist)}")
        system.append("U_ENERGY=" + ", ".join(elist))

    b_impct = num(form.get("b_impct"))
    if b_impct:
        system.append(f"B_IMPCT={b_impct}")
    else:
        add(system, form, "jtotl", "JTOTL")
        add(system, form, "jtotu", "JTOTU")

    add(system, form, "time_step", "TIME_STEP")
    add(system, form, "time_lim", "TIME_LIM")

    if form.get("monte_carlo"):
        system.append("MONTE_CARLO=YES")
        add(system, form, "nmb_traj", "NMB_TRAJ")

    if st in IDENTICAL_OK and form.get("identical"):
        add(system, form, "wght_pospar", "WGHT_POSPAR")

    lines.append("$SYSTEM")
    lines.append("  " + ",\n  ".join(system))
    lines.append("$END")
    lines.append("")

    # ======================= $POTENTIAL =======================
    pot = []
    pot.append(f"E_UNITS={form.get('e_units', 'A.U.')}")
    pot.append(f"R_UNITS={form.get('r_units', 'A.U.')}")
    add(pot, form, "grd_r", "GRD_R")

    # Angular grids. For molecule+atom each keyword is a single integer.
    # For molecule+molecule each keyword takes two comma-separated integers;
    # some entries are dummy per SYS_TYPE (manual Sec. IV) but must be present.
    if st in MOL_ATOM:
        add(pot, form, "grd_ang2", "GRD_ANG2")          # beta
        if st in ("3", "4"):
            add(pot, form, "grd_ang3", "GRD_ANG3")      # gamma
        if st == "2":
            add(pot, form, "grd_vib", "GRD_VIB")
    else:
        # molecule + molecule: two values per angular keyword
        g1a = num(form.get("grd_ang1_a"), "1")
        g1b = num(form.get("grd_ang1_b"), "1")
        g2a = num(form.get("grd_ang2_a"), "1")
        g2b = num(form.get("grd_ang2_b"), "1")
        g3a = num(form.get("grd_ang3_a"), "1")
        g3b = num(form.get("grd_ang3_b"), "1")
        pot.append(f"GRD_ANG1={g1a},{g1b}")
        pot.append(f"GRD_ANG2={g2a},{g2b}")
        pot.append(f"GRD_ANG3={g3a},{g3b}")
        if st == "6":
            gv1 = num(form.get("grd_vib_a"), "1")
            gv2 = num(form.get("grd_vib_b"), "1")
            pot.append(f"GRD_VIB={gv1},{gv2}")

    if form.get("save_mtrx"):
        pot.append("SAVE_MTRX=YES")
    if form.get("read_mtrx"):
        pot.append("READ_MTRX=YES")
    if form.get("prog_run_no"):
        pot.append("PROG_RUN=NO")

    lines.append("$POTENTIAL")
    lines.append("  " + ",\n  ".join(pot))
    lines.append("$END")

    return "\n".join(lines) + "\n"


@app.route("/")
def index():
    return render_template(
        "index.html",
        sys_types=SYS_TYPES,
        channel_labels=CHANNEL_LABELS,
    )


@app.route("/generate", methods=["POST"])
def generate():
    text = build_input(request.form)
    errors, warnings = validate(request.form)

    # Download: still produce the file even with warnings, but block on errors.
    if request.form.get("download"):
        if errors:
            return jsonify(ok=False, errors=errors, warnings=warnings,
                           text=text), 400
        fname = (num(request.form.get("label")) or "mqct") + ".inp"
        return Response(
            text,
            mimetype="text/plain",
            headers={"Content-Disposition": f"attachment; filename={fname}"},
        )

    # Preview: return the text plus validation results as JSON.
    return jsonify(ok=not errors, errors=errors,
                   warnings=warnings, text=text)


if __name__ == "__main__":
    app.run(debug=True, port=5000)
