# -*- coding: utf-8 -*-
"""
Vergleichsabbildungen V3–V13 der EBSD-Auswertung (ab Version 1.4, ohne MTEX).

Zusätzlich zu den Kennwerten der Auswertung werden je Datensatz Versuchsdaten verwendet, die im Reiter
„Datensätze / Vergleich“ eingetragen werden (Legierung, Position, Prüfbedingungen, lokale Werte aus der
FE-Simulation, Fließspannung). Daraus wird der Zener–Hollomon-Parameter

    Z = φ̇ · exp(Q / (R·T))

mit einem gemeinsamen Q für alle Datensätze berechnet (Eingabefeld Q, Standard 142 kJ/mol).

Datensätze, die vor 1.4 ausgewertet wurden, enthalten nicht alle zusätzlichen Kennwerte im JSON; sie werden dann
aus Kennwerte_und_Parameter.txt gelesen (deutsche und englische Ausgabe).
"""
from __future__ import annotations

import csv
import io
import math
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, NullFormatter, ScalarFormatter

from ebsd_i18n import L, dez, num_g

R_GAS = 8.314
B_AL = 0.286e-3                      # Burgers-Vektor Al in µm
comma = FuncFormatter(lambda v, p: num_g(v))

# ======================================================================== Optionen (von der Oberfläche gesetzt)
OPTIONS = {
    "Q_kJ": 142.0,            # gemeinsame Aktivierungsenergie für Z
    "z_basis": "lokal",       # "lokal": lokale T und φ̇ (FE), sonst Prüfbedingung; "nominell": immer Prüfbedingung
    "v3_key": "overview",     # Kennwert in V3 ("overview" = sechs Kennwerte)
    "same_scale": True,       # V6: Karten im gleichen µm-Maßstab
}

# ======================================================================== Versuchsdaten je Datensatz
# (Schlüssel, Beschriftung, Beschriftung engl., Einheit, Typ, Nachkommastellen)
TEST_FIELDS = [
    ("legierung", "Legierung", "Alloy", "", "str", None),
    ("fe_pct", "Fe-Gehalt", "Fe content", "Ma.-%", "float", 2),
    ("si_pct", "Si-Gehalt", "Si content", "Ma.-%", "float", 2),
    ("position", "Position", "Position", "", "str", None),
    ("T_C", "Prüftemperatur", "Test temperature", "°C", "float", 0),
    ("rate_s", "Umformgeschwindigkeit (Prüfung)", "Strain rate (test)", "s⁻¹", "float", 3),
    ("T_lok_C", "lokale Temperatur (Ende)", "Local temperature (end)", "°C", "float", 1),
    ("phi_lok", "lokaler Umformgrad φ", "Local strain φ", "", "float", 2),
    ("rate_lok_s", "lokale Umformgeschwindigkeit φ̇", "Local strain rate φ̇", "s⁻¹", "float", 3),
    ("kf_MPa", "Fließspannung k_f", "Flow stress k_f", "MPa", "float", 1),
]
CENTRE_WORDS = ("mitte", "centre", "center", "m")


def test(ds):
    """Versuchsdaten des Datensatzes (leeres dict, wenn keine eingetragen sind)"""
    return dict(ds.meta.get("versuch") or {})


def tval(ds, key):
    v = test(ds).get(key)
    if v in (None, ""):
        return None
    if key in ("legierung", "position"):
        return str(v)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def guess_test(ds):
    """Vorschlag aus Bezeichnung und Dateiname: Legierung, Prüftemperatur, Umformgeschwindigkeit, Position"""
    out = {}
    txt = " ".join([ds.label, ds.metrics.get("datei", ""), str(ds.path)])
    m = re.search(r"(\d{3,4}-\d+)[_ ,]+(\d{3})(?:\s*°C)?[_ /]+(\d+(?:[.,]\d+)?)", txt)
    if m:
        out["legierung"] = m.group(1)
        out["T_C"] = float(m.group(2))
        out["rate_s"] = float(m.group(3).replace(",", "."))
    low = txt.lower()
    for key, short in (("mitte", "M"), ("centre", "M"), ("schmiedekreuz", "SK"), ("forging cross", "SK"), ("rand", "Rand")):
        if key in low:
            out["position"] = short
            break
    return out


def is_centre(ds):
    p = (tval(ds, "position") or "").strip().lower()
    return p == "" or p in CENTRE_WORDS or p.startswith("mitte") or p.startswith("cent")


def alloy(ds):
    a = tval(ds, "legierung")
    if a:
        return a
    return ds.label.split(",")[0].strip()


def lnZ(ds, basis=None):
    """(ln Z, Quelle) mit dem gemeinsamen Q; Quelle "lokal", "Prüfung" oder None"""
    basis = basis or OPTIONS["z_basis"]
    Q = float(OPTIONS["Q_kJ"]) * 1e3
    cand = []
    if basis == "lokal":
        cand.append(("lokal", tval(ds, "T_lok_C"), tval(ds, "rate_lok_s")))
    cand.append(("Prüfung", tval(ds, "T_C"), tval(ds, "rate_s")))
    for src, T, r in cand:
        if T is not None and r is not None and r > 0:
            return math.log(r) + Q / (R_GAS * (T + 273.15)), src
    return None, None


def temperature(ds):
    """Temperatur für G(T): lokal, sonst Prüftemperatur (°C)"""
    return tval(ds, "T_lok_C") if tval(ds, "T_lok_C") is not None else tval(ds, "T_C")


def shear_modulus_MPa(T_C):
    """Schubmodul Al, linear in T (25,4 GPa bei 300 K, −50 % bei der Schmelztemperatur 933 K)"""
    return 25.4e3 * (1 - 0.5 * (T_C + 273.15 - 300) / 933)


def z_label():
    return (L("ln Z, lokal", "ln Z, local") if OPTIONS["z_basis"] == "lokal" else L("ln Z, Prüfbedingung", "ln Z, test condition")) \
        + f" (Q = {dez(OPTIONS['Q_kJ'], 0)} kJ mol⁻¹)"


# ======================================================================== zusätzliche Kennwerte (JSON oder Text)
_NUM = r"(-?\d+(?:[.,]\d+)?)"
_cache = {}


def _f(s):
    try:
        return float(str(s).replace(",", "."))
    except ValueError:
        return None


def _parse_text(ds):
    """liest Kennwerte aus Kennwerte_und_Parameter.txt (Datensätze vor 1.4); Ergebnis wird je Datensatz zwischengespeichert"""
    key = (str(ds.folder), ds.meta.get("geaendert"))
    if key in _cache:
        return _cache[key]
    t = ds.kennwerte()
    out = {}
    m = re.search(r"(?:rekristallisierter Anteil vs\. Abstand zu Partikeln|recrystallized fraction vs\. distance to particles)[^:]*\):(.*)", t)
    if m:
        items = re.findall(r"([^;:]+?) µm: " + _NUM + r" % \((?:Flächenanteil|area fraction) " + _NUM + r" %\)", m.group(1))
        if items:
            out["rx_abstand_klassen"] = [i[0].strip() for i in items]
            out["rx_abstand_pct"] = [_f(i[1]) for i in items]
            out["rx_abstand_flaeche_pct"] = [_f(i[2]) for i in items]
    m = re.search(r"(?:Sensitivität rekristallisierter Anteil|Sensitivity of the recrystallized fraction):(.*)", t)
    if m:
        items = re.findall(r"GOS <= " + _NUM + "°: " + _NUM + " %", m.group(1))
        if items:
            out["gos_sens"] = {"grenzen_deg": [_f(a) for a, _ in items], "rx_pct": [_f(b) for _, b in items]}
    m = re.search(r"\((?:3 größte Körner|3 largest grains): " + _NUM + ", " + _NUM + ", " + _NUM + " %", t)
    if m:
        out["groesste_koerner_pct"] = [_f(m.group(i)) for i in (1, 2, 3)]
    m = re.search(r"(?:Körner|Grains) \(MTEX calcGrains\): (\d+) \((?:klassifizierbar|classifiable) (\d+)\)", t)
    if m:
        out["n_koerner"], out["n_koerner_klassifizierbar"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"(?:Subkörner|subgrains) \(>= [^)]*\): (\d+)", t)
    if m:
        out["n_subkoerner"] = int(m.group(1))
    m = re.search(r"(?:Mittelwert|mean) " + _NUM + "°;(?:Median|median) " + _NUM + "°;90% ?-?(?:Quantil|quantile) " + _NUM + "°", t)
    if m:
        out["kam_median_deg"], out["kam_q90_deg"] = _f(m.group(2)), _f(m.group(3))
    m = re.search(r"(?:mittlere KAM|mean KAM): (?:rekristallisiert|recrystallized) " + _NUM + "°, (?:substrukturiert|substructured) "
                  + _NUM + "°, (?:verformt|deformed) " + _NUM + "°", t)
    if m:
        out["kam_klassen_deg"] = {"rx": _f(m.group(1)), "sub": _f(m.group(2)), "def": _f(m.group(3))}
    m = re.search(r"<110> \|\| CD \(15°\): (?:rekristallisiert|recrystallized) " + _NUM + " %, (?:substrukturiert|substructured) "
                  + _NUM + " %, (?:verformt|deformed) " + _NUM + " %", t)
    if m:
        out["faser110_klassen_pct"] = {"rx": _f(m.group(1)), "sub": _f(m.group(2)), "def": _f(m.group(3))}
    part = 0.0
    found = False
    for line in t.splitlines():
        f_ = line.split(";")
        if len(f_) >= 4 and f_[1].strip().lower() in ("partikel", "particle", "particles"):
            v = _f(f_[3])
            if v is not None:
                part += v
                found = True
    if found:
        out["partikel_pct"] = part
    _cache[key] = out
    return out


def xm(ds, key):
    """Kennwert aus dem JSON des Datensatzes, sonst aus dem Kennwerte-Text"""
    v = ds.metrics.get(key)
    if v is not None:
        return v
    return _parse_text(ds).get(key)


def subgrain_ecd(ds):
    """ECD der Subkörner in der Statistik (Kornliste_Subkoerner.csv bzw. Grain_list_subgrains.csv) oder None"""
    for name, text in ds.tables().items():
        low = name.lower()
        if "subk" not in low and "subgrain" not in low:
            continue
        sep = ";" if text.count(";") > text.count(",") else ","
        rows = list(csv.reader(io.StringIO(text), delimiter=sep))
        if len(rows) < 2:
            return None
        head = [h.strip().lower() for h in rows[0]]
        try:
            i_ecd = next(i for i, h in enumerate(head) if h.startswith("ecd"))
            i_in = next(i for i, h in enumerate(head) if h.startswith("in stat"))
        except StopIteration:
            return None
        vals = []
        for r in rows[1:]:
            if len(r) > max(i_ecd, i_in) and r[i_in].strip() == "1":
                v = _f(r[i_ecd]) if sep == ";" else (float(r[i_ecd]) if r[i_ecd] else None)
                if v is not None:
                    vals.append(v)
        return np.array(vals) if vals else None
    return None


# ======================================================================== Gestaltung
ALLOY_COLORS = ["#14466B", "#E3A35B", "#6A9F58", "#A23B2A", "#7FA7C9", "#8E6BB0", "#C9A227", "#4D4D4D"]


def style(datasets):
    """Farbe je Legierung, gefülltes Symbol = Probenmitte, offenes Symbol = andere Position"""
    alloys = []
    for d in datasets:
        if alloy(d) not in alloys:
            alloys.append(alloy(d))
    col = {a: ALLOY_COLORS[i % len(ALLOY_COLORS)] for i, a in enumerate(alloys)}
    out = {}
    for d in datasets:
        c = col[alloy(d)]
        out[d.id] = dict(color=c, mfc=c if is_centre(d) else "white", mec=c, marker="o" if is_centre(d) else "s")
    return out, col


def short_label(d):
    return d.label


def _finish(ax, xlabel=None, ylabel=None, grid=True):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if grid:
        ax.grid(alpha=0.3)
    ax.set_axisbelow(True)
    for a in (ax.xaxis, ax.yaxis):
        if a.get_scale() == "linear" and isinstance(a.get_major_formatter(), ScalarFormatter):   # nur Zahlenachsen
            a.set_major_formatter(comma)


def _note(fig, lines, y=0.005):
    lines = [l for l in lines if l]
    if lines:
        fig.text(0.01, y, "\n".join(lines), fontsize=8.5, color="0.35", va="bottom")


def _empty(fig, text):
    fig.text(0.5, 0.5, text, ha="center", va="center", color="0.4", fontsize=12)
    return fig


def _legend_alloys(ax, col, extra_pos=True, loc="best"):
    h = [Line2D([], [], color=c, marker="o", ls="", ms=7, label=a) for a, c in col.items()]
    if extra_pos:
        h += [Line2D([], [], color="0.3", marker="o", ls="", ms=7, label=L("Probenmitte (gefüllt)", "specimen centre (filled)")),
              Line2D([], [], color="0.3", marker="s", mfc="white", ls="", ms=7, label=L("andere Position (offen)", "other position (open)"))]
    ax.legend(handles=h, frameon=False, fontsize=8.5, loc=loc)


# ======================================================================== V3 Kennwerte über ln Z
V3_KEYS = [  # (Schlüssel, Beschriftung, Beschriftung engl., Einheit, log-Achse, Exponentialfit)
    ("subkorn_ecd_median_um", "Subkorngröße (Median-ECD)", "Subgrain size (median ECD)", "µm", True, True),
    ("hagb_anteil_pct", "HAGB-Anteil der Grenzlänge", "HAGB fraction of boundary length", "%", False, False),
    ("rx_pct", "rekristallisierter Anteil (GOS ≤ θc)", "recrystallized fraction (GOS ≤ θc)", "%", False, False),
    ("faser110_CD_pct", "<110> ∥ Stauchrichtung (±15°)", "<110> ∥ compression direction (±15°)", "%", False, False),
    ("kam_mittel_deg", "KAM Mittelwert", "KAM mean", "°", False, False),
    ("korn_ecd_median_um", "Korngröße (Median-ECD)", "Grain size (median ECD)", "µm", True, True),
    ("sehne_hagb_CD_um", "HAGB-Sehnenlänge in CD", "HAGB intercept along CD", "µm", True, True),
    ("_rx_nah", "rekristallisiert 0–1 µm von Partikeln", "recrystallized 0–1 µm from particles", "%", False, False),
    ("def_pct", "verformt (Subkorn-GOS > θc)", "deformed (subgrain GOS > θc)", "%", False, False),
    ("partikel_pct", "Partikelanteil (Fläche)", "Particle fraction (area)", "%", False, False),
]
V3_OVERVIEW = ["subkorn_ecd_median_um", "hagb_anteil_pct", "rx_pct", "faser110_CD_pct", "kam_mittel_deg", "_rx_nah"]


def v3_choices():
    return [("overview", L("Übersicht (6 Kennwerte)", "Overview (6 key figures)"))] + [(k, L(a, b)) for k, a, b, *_ in V3_KEYS]


def v3_value(ds, key):
    if key == "_rx_nah":
        v = xm(ds, "rx_abstand_pct")
        return v[0] if v else None
    v = xm(ds, key)
    return float(v) if isinstance(v, (int, float)) else None


def point_label(d, datasets):
    """Beschriftung eines Punktes: ohne Legierung, wenn alle Datensätze dieselbe Legierung haben"""
    if len({alloy(x) for x in datasets}) > 1:
        return d.label
    return d.label.split(", ", 1)[-1] if ", " in d.label else d.label


def _v3_panel(ax, datasets, key, st):
    k = next(e for e in V3_KEYS if e[0] == key)
    _, lab, lab_en, unit, logy, fit = k
    pts = {}
    for d in datasets:
        z, src = lnZ(d)
        v = v3_value(d, key)
        if z is None or v is None:
            continue
        s = st[d.id]
        ax.plot(z, v, marker=s["marker"], ms=8, mfc=s["mfc"], mec=s["mec"], mew=1.6, ls="", zorder=3)
        ax.annotate(point_label(d, datasets), (z, v), xytext=(5, 4),
                    textcoords="offset points", fontsize=7, color="0.35")
        pts.setdefault((alloy(d), is_centre(d)), []).append((z, v))
    if fit:
        for (a, centre), p in pts.items():
            if len(p) < 3 or not centre:
                continue
            z = np.array([q[0] for q in p])
            y = np.array([q[1] for q in p])
            if np.any(y <= 0) or np.ptp(z) < 0.5:
                continue
            c = np.polyfit(z, np.log(y), 1)
            r = np.corrcoef(z, np.log(y))[0, 1]
            zz = np.linspace(z.min() - 0.5, z.max() + 0.5, 50)
            col = st[next(d.id for d in datasets if alloy(d) == a)]["color"]
            ax.plot(zz, np.exp(np.polyval(c, zz)), "--", color=col, lw=1.2)
            ax.text(0.98, 0.97 - 0.07 * list(pts).index((a, centre)), f"{a}: ∝ Z^{dez(c[0], 3)}, r = {dez(r, 3)}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8, color=col)
    if logy:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(comma)
        ax.yaxis.set_minor_formatter(comma)
    ax.set_title(L(lab, lab_en), fontsize=11)
    _finish(ax, z_label(), unit)
    return sum(len(p) for p in pts.values())


def fig_v3(datasets):
    key = OPTIONS["v3_key"]
    st, col = style(datasets)
    missing = [d.label for d in datasets if lnZ(d)[0] is None]
    if key == "overview":
        fig, axs = plt.subplots(2, 3, figsize=(16, 9.4), dpi=100)
        n = 0
        for ax, k in zip(axs.ravel(), V3_OVERVIEW):
            n = max(n, _v3_panel(ax, datasets, k, st))
        _legend_alloys(axs[0, 0], col)
    else:
        fig, ax = plt.subplots(figsize=(10, 6.8), dpi=100)
        n = _v3_panel(ax, datasets, key, st)
        _legend_alloys(ax, col)
    fig.suptitle(L("Kennwerte über dem Zener–Hollomon-Parameter", "Key figures versus the Zener–Hollomon parameter"), fontsize=14)
    fits = {k for k, *_r, fit in V3_KEYS if fit}
    notes = [L("Fit: exponentiell, je Legierung über die Datensätze der Probenmitte (mind. 3)",
               "Fit: exponential, per alloy over the specimen-centre datasets (at least 3)")] \
        if (key == "overview" and fits & set(V3_OVERVIEW)) or key in fits else []
    if any(lnZ(d)[1] == "Prüfung" for d in datasets) and OPTIONS["z_basis"] == "lokal":
        notes.append(L("Ohne lokale Werte wird ln Z aus der Prüfbedingung berechnet.",
                       "Without local values, ln Z is calculated from the test condition."))
    if missing:
        notes.append(L("ohne Versuchsdaten (nicht dargestellt): ", "without test data (not shown): ") + "; ".join(missing))
    fig.tight_layout(rect=(0, 0.06 if missing else 0.04, 1, 0.96))
    _note(fig, notes)
    if n == 0:
        _empty(fig, L("Keine Datensätze mit Versuchsdaten.\n„Versuchsdaten …“ eintragen (T und φ̇).",
                      "No datasets with test data.\nEnter “Test data …” (T and φ̇)."))
    return fig


# ======================================================================== V4 Rekristallisation über Partikelabstand
def fig_v4(datasets):
    from ebsd_datasets import colors_for
    cols = colors_for(datasets)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), dpi=100, gridspec_kw=dict(width_ratios=[1.6, 1]))
    names, ratio, used = [], [], []
    for d in datasets:
        v = xm(d, "rx_abstand_pct")
        cls = xm(d, "rx_abstand_klassen")
        if not v or not cls:
            continue
        x = np.arange(len(v))
        ax.plot(x, v, "-o", color=cols[d.id], lw=2, ms=6, label=d.label)
        used.append(cls)
        names.append(d.label)
        ratio.append(v[0] / v[-1] if v[-1] else np.nan)
    if used:
        ax.set_xticks(np.arange(len(used[0])))
        ax.set_xticklabels([c.replace("-", "–") for c in used[0]])
        ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(0, None)
    _finish(ax, L("Abstand zum nächsten Partikel (µm)", "Distance to the nearest particle (µm)"),
            L("rekristallisierter Anteil (%)", "recrystallized fraction (%)"))
    ax.set_title(L("Rekristallisierter Anteil über dem Partikelabstand", "Recrystallized fraction versus particle distance"), fontsize=12)
    if ratio:
        y = np.arange(len(ratio))
        ax2.barh(y, np.nan_to_num(ratio), color=[cols[d.id] for d in datasets if xm(d, "rx_abstand_pct")])
        for yi, r in zip(y, ratio):
            ax2.text(np.nan_to_num(r) + 0.05, yi, dez(r, 1), va="center", fontsize=9)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=9)
        ax2.invert_yaxis()
        ax2.axvline(1, color="k", ls=":", lw=1)
    _finish(ax2, L("Verhältnis rx(0–1 µm) / rx(> 10 µm)", "ratio rx(0–1 µm) / rx(> 10 µm)"))
    ax2.set_title(L("Bevorzugung partikelnaher Bereiche", "Preference for regions near particles"), fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    _note(fig, [L("Partikel: Partikelphasen und Mg2Si-Cluster ab dem PSN-Mindestdurchmesser (Parameter „PSN: Partikel ab“)",
                  "Particles: particle phases and Mg2Si clusters above the PSN minimum diameter (parameter “PSN: particles from”)")])
    if not names:
        _empty(fig, L("Keine Rekristallisationsauswertung in den gewählten Datensätzen.", "No recrystallization analysis in the selected datasets."))
    return fig


# ======================================================================== V5 Zerlegung (Referenz = erster Datensatz)
def decompose(ref, d):
    """Symmetrische Zerlegung der Differenz des rekristallisierten Anteils (d − Referenz) in
    Partikelnähe (andere Abstandsverteilung) und Rekristallisation je Abstandsklasse"""
    r0, a0 = np.array(xm(ref, "rx_abstand_pct"), float), np.array(xm(ref, "rx_abstand_flaeche_pct"), float)
    r1, a1 = np.array(xm(d, "rx_abstand_pct"), float), np.array(xm(d, "rx_abstand_flaeche_pct"), float)
    a0, a1 = a0 / a0.sum(), a1 / a1.sum()
    tot0, tot1 = float((r0 * a0).sum()), float((r1 * a1).sum())
    near = float(((a1 - a0) * (r1 + r0) / 2).sum())
    per = float(((r1 - r0) * (a1 + a0) / 2).sum())
    return dict(tot_ref=tot0, tot=tot1, near=near, per=per, cf_d_with_ref_dist=float((r1 * a0).sum()),
                cf_ref_with_d_dist=float((r0 * a1).sum()))


def _has_dist(d):
    return bool(xm(d, "rx_abstand_pct")) and bool(xm(d, "rx_abstand_flaeche_pct"))


def fig_v5(datasets):
    from ebsd_datasets import colors_for
    cols = colors_for(datasets)
    ok = [d for d in datasets if _has_dist(d)]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, max(5.4, 2.0 + 0.75 * len(ok))), dpi=100)
    if len(ok) < 2:
        plt.close(fig)
        fig = plt.figure(figsize=(12, 5), dpi=100)
        return _empty(fig, L("Mindestens zwei Datensätze mit Rekristallisation über dem Partikelabstand nötig.",
                             "At least two datasets with recrystallization versus particle distance are needed."))
    ref = ok[0]
    others = ok[1:]
    y = np.arange(len(others))
    res = [decompose(ref, d) for d in others]
    near = np.array([r["near"] for r in res])
    per = np.array([r["per"] for r in res])
    ax.barh(y, near, color="#E3A35B", label=L("Partikelnähe (Abstandsverteilung)", "proximity to particles (distance distribution)"))
    ax.barh(y, per, left=np.where(np.sign(per) == np.sign(near), near, 0), color="#14466B",
            label=L("Rekristallisation je Abstandsklasse", "recrystallization per distance class"))
    for yi, r in zip(y, res):
        ax.text(max(r["near"], 0) + max(r["per"], 0) + 0.2, yi, L(f"Δ {dez(r['tot'] - r['tot_ref'], 1)} Pp.", f"Δ {dez(r['tot'] - r['tot_ref'], 1)} pp"),
                va="center", fontsize=9)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([d.label for d in others], fontsize=9)
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    _finish(ax, L(f"Differenz zur Referenz „{ref.label}“ (Prozentpunkte)", f"Difference to the reference “{ref.label}” (percentage points)"))
    ax.set_title(L("Zerlegung der Differenz des rekristallisierten Anteils", "Decomposition of the difference in recrystallized fraction"),
                 fontsize=12)
    w = 0.38
    ax2.barh(y - w / 2, [r["tot"] for r in res], w, color=[cols[d.id] for d in others], label=L("gemessen", "measured"))
    ax2.barh(y + w / 2, [r["cf_d_with_ref_dist"] for r in res], w, color=[cols[d.id] for d in others], alpha=0.4, hatch="//",
             label=L("mit der Abstandsverteilung der Referenz", "with the distance distribution of the reference"))
    ax2.axvline(res[0]["tot_ref"], color=cols[ref.id], ls="--", lw=1.5, label=L(f"Referenz {dez(res[0]['tot_ref'], 1)} %",
                                                                                 f"reference {dez(res[0]['tot_ref'], 1)} %"))
    for yi, r in zip(y, res):
        ax2.text(r["tot"] + 0.2, yi - w / 2, dez(r["tot"], 1), va="center", fontsize=8.5)
        ax2.text(r["cf_d_with_ref_dist"] + 0.2, yi + w / 2, dez(r["cf_d_with_ref_dist"], 1), va="center", fontsize=8.5)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.invert_yaxis()
    ax2.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    _finish(ax2, L("rekristallisierter Anteil (%)", "recrystallized fraction (%)"))
    ax2.set_title(L("Was-wäre-wenn: gleiche Partikelnähe wie die Referenz", "What-if: same particle proximity as the reference"), fontsize=12)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    _note(fig, [L("Referenz = erster Datensatz im Vergleich (Reihenfolge mit ▲/▼ ändern). Zerlegung symmetrisch: "
                  "Δ = Σ (a − a_ref)·(r + r_ref)/2 + Σ (r − r_ref)·(a + a_ref)/2; r = rekristallisierter Anteil je Abstandsklasse, "
                  "a = Flächenanteil der Klasse.",
                  "Reference = first dataset in the comparison (change the order with ▲/▼). Symmetric decomposition: "
                  "Δ = Σ (a − a_ref)·(r + r_ref)/2 + Σ (r − r_ref)·(a + a_ref)/2; r = recrystallized fraction per distance class, "
                  "a = area fraction of the class.")])
    return fig


# ======================================================================== V6 Kartenvergleich
# Maße der Kartenabbildungen (wie ebsd_core.MapLayout, Vorschau mit 110 dpi)
MAP_W_MAX, MAP_H_MAX, LEFT_MIN, M_L, TITLE_H, PREVIEW_DPI = 9.5, 7.0, 5.6, 0.3, 1.0, 110


def map_crop(ds, fid):
    """Kartenbereich aus der Vorschau der Abbildung fid (PIL-Bild) oder None"""
    b = ds.preview(fid)
    if not b:
        return None
    from PIL import Image
    im = Image.open(io.BytesIO(b)).convert("RGB")
    W, H = ds.metrics.get("flaeche_um") or [None, None]
    if not W or not H:
        return None
    s = min(MAP_W_MAX / W, MAP_H_MAX / H)
    mw, mh = W * s, H * s
    left_w = max(mw, LEFT_MIN)
    x0 = (M_L + (left_w - mw) / 2) * PREVIEW_DPI
    y0 = TITLE_H * PREVIEW_DPI
    box = (int(round(x0)) + 1, int(round(y0)) + 1, int(round(x0 + mw * PREVIEW_DPI)) - 1, int(round(y0 + mh * PREVIEW_DPI)) - 1)
    if box[2] > im.size[0] or box[3] > im.size[1]:
        return None
    return im.crop(box)


def fig_v6(datasets):
    rows = [("01", L("IPF ∥ CD", "IPF ∥ CD")), ("09", L("GOS-Klassen", "GOS classes"))]
    n = len(datasets)
    if n == 0:
        return _empty(plt.figure(figsize=(10, 4)), "")
    areas = [(d.metrics.get("flaeche_um") or [1, 1]) for d in datasets]
    wmax = max(a[0] for a in areas)
    hmax = max(a[1] for a in areas)
    cell_w = min(4.2, 17.0 / n)
    cell_h = cell_w * hmax / wmax
    fig = plt.figure(figsize=(0.7 + n * (cell_w + 0.15), 0.9 + len(rows) * (cell_h + 0.55) + 0.4), dpi=100)
    FW, FH = fig.get_size_inches()
    for j, d in enumerate(datasets):
        W, H = areas[j]
        if OPTIONS["same_scale"]:
            w_in, h_in = cell_w * W / wmax, cell_w * H / wmax
        else:
            k = min(cell_w / W, cell_h / H)
            w_in, h_in = W * k, H * k
        for i, (fid, rlab) in enumerate(rows):
            x = 0.55 + j * (cell_w + 0.15)
            ytop = 0.8 + i * (cell_h + 0.55)
            ax = fig.add_axes([x / FW, 1 - (ytop + h_in) / FH, w_in / FW, h_in / FH])
            im = map_crop(d, fid)
            if im is not None:
                ax.imshow(np.asarray(im), extent=(0, W, H, 0), interpolation="antialiased")
            else:
                ax.text(0.5, 0.5, L("keine Vorschau", "no preview"), ha="center", va="center", transform=ax.transAxes, color="0.5")
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_title(f"{d.label}\n{dez(W, 0)} × {dez(H, 0)} µm, {dez(d.metrics.get('schrittweite_um'), 2)} µm",
                             fontsize=8.5, loc="left")
            if j == 0:
                fig.text(0.12 / FW, 1 - (ytop + h_in / 2) / FH, rlab, rotation=90, va="center", ha="left", fontsize=11)
    fig.suptitle(L("Kartenvergleich", "Map comparison") + (L(" – gleicher Maßstab (µm)", " – same scale (µm)") if OPTIONS["same_scale"]
                                                         else L(" – jede Karte auf die Zellgröße skaliert", " – each map scaled to the cell")),
                 fontsize=13)
    _note(fig, [L("Ausschnitte aus den gespeicherten Vorschauen (Abbildungen 01 und 09); Farbschlüssel siehe dort.",
                  "Cut-outs from the stored previews (figures 01 and 09); see there for the colour keys.")])
    return fig


# ======================================================================== V7 GOS-Klassen
def _order_by_z(datasets):
    z = [lnZ(d)[0] for d in datasets]
    if all(v is not None for v in z):
        return [d for _, d in sorted(zip(z, datasets), key=lambda t: t[0])], True
    return list(datasets), False


def fig_v7(datasets):
    dss, by_z = _order_by_z(datasets)
    n = len(dss)
    fig, ax = plt.subplots(figsize=(max(9, 1.6 + 1.3 * n), 6.4), dpi=100)
    x = np.arange(n)
    parts = [("rx_pct", L("rekristallisiert", "recrystallized"), "#2B6FD6"),
             ("sub_pct", L("substrukturiert", "substructured"), "#F2C744"),
             ("def_pct", L("verformt", "deformed"), "#D63B2B")]
    bottom = np.zeros(n)
    for key, lab, c in parts:
        v = np.array([np.nan_to_num(_nanf(xm(d, key))) for d in dss])
        if key == "def_pct":
            v = np.array([np.nan_to_num(_nanf(xm(d, "def_pct") if xm(d, "def_pct") is not None else xm(d, "nicht_rx_pct"))) for d in dss])
        ax.bar(x, v, 0.7, bottom=bottom, color=c, label=lab, edgecolor="white", linewidth=1.5)
        for xi, b, vi in zip(x, bottom, v):
            if vi >= 4:
                ax.text(xi, b + vi / 2, dez(vi, 1), ha="center", va="center", fontsize=8.5, color="white" if c != "#F2C744" else "k")
        bottom += v
    labs = []
    for d in dss:
        z, _ = lnZ(d)
        labs.append(d.label.replace(", ", "\n") + (f"\nln Z = {dez(z, 1)}" if z is not None else ""))
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylim(0, 105)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.08))
    _finish(ax, None, L("Flächenanteil der Matrix (%)", "Area fraction of the matrix (%)"))
    ax.set_title(L("GOS-Klassen", "GOS classes") + (L(" (sortiert nach ln Z)", " (sorted by ln Z)") if by_z else ""), fontsize=12, pad=26)
    fig.tight_layout()
    return fig


def _nanf(v):
    return np.nan if v is None else float(v)


# ======================================================================== V8 Subkorngrößenverteilung
def fig_v8(datasets):
    from ebsd_datasets import colors_for
    cols = colors_for(datasets)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), dpi=100)
    miss = []
    shown = 0
    for d in datasets:
        e = subgrain_ecd(d)
        if e is None or len(e) < 10:
            miss.append(d.label)
            continue
        e = np.sort(e)
        y = np.arange(1, len(e) + 1) / len(e) * 100
        ax.plot(e, y, color=cols[d.id], lw=2, label=f"{d.label} (n = {len(e)})")
        bins = np.logspace(np.log10(max(e.min(), 0.1)), np.log10(e.max()), 40)
        h, b = np.histogram(e, bins=bins)
        ax2.plot(np.sqrt(b[:-1] * b[1:]), h / h.sum() * 100, "-", color=cols[d.id], lw=1.8)
        shown += 1
    for a in (ax, ax2):
        a.set_xscale("log")
        a.xaxis.set_major_formatter(comma)
        a.xaxis.set_minor_formatter(NullFormatter())
    ax.axhline(50, color="0.5", ls=":", lw=1)
    _finish(ax, L("Subkorn-ECD (µm)", "Subgrain ECD (µm)"), L("Summenhäufigkeit (Anzahl, %)", "Cumulative frequency (number, %)"))
    _finish(ax2, L("Subkorn-ECD (µm)", "Subgrain ECD (µm)"), L("Häufigkeit je logarithmischer Klasse (%)", "Frequency per logarithmic class (%)"))
    ax.set_title(L("Summenkurve", "Cumulative distribution"), fontsize=12)
    ax2.set_title(L("Verteilung", "Distribution"), fontsize=12)
    if shown:
        ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.suptitle(L("Subkorngrößenverteilung (Subkörner in der Statistik)", "Subgrain size distribution (subgrains in the statistics)"),
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    _note(fig, [L("ohne Subkorn-Kornliste: ", "without subgrain list: ") + "; ".join(miss) if miss else ""])
    if not shown:
        _empty(fig, L("Keine Subkorn-Kornlisten in den gewählten Datensätzen.", "No subgrain lists in the selected datasets."))
    return fig


# ======================================================================== V9 Subkorngröße und Fließspannung
def holt_K(d):
    dsub, kf, T = xm(d, "subkorn_ecd_median_um"), tval(d, "kf_MPa"), temperature(d)
    if dsub is None or kf is None or T is None:
        return None
    return dsub * kf / (shear_modulus_MPa(T) * B_AL)


def fig_v9(datasets):
    st, col = style(datasets)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), dpi=100)
    names, Ks, cs = [], [], []
    for d in datasets:
        dsub, kf, T = xm(d, "subkorn_ecd_median_um"), tval(d, "kf_MPa"), temperature(d)
        if dsub is None or kf is None or T is None:
            continue
        G = shear_modulus_MPa(T)
        s = st[d.id]
        ax.plot(kf / G, dsub / B_AL, marker=s["marker"], ms=8, mfc=s["mfc"], mec=s["mec"], mew=1.6, ls="")
        ax.annotate(point_label(d, datasets), (kf / G, dsub / B_AL), xytext=(5, 4), textcoords="offset points", fontsize=7, color="0.35")
        names.append(d.label)
        Ks.append(holt_K(d))
        cs.append(s["color"])
    if names:
        xx = np.logspace(np.log10(min(ax.get_xlim()[0], 1e-3)), np.log10(ax.get_xlim()[1] * 1.2), 20)
        for K, ls in ((10, ":"), (20, "--"), (30, ":")):
            ax.plot(xx, K / xx, ls, color="0.5", lw=1)
            ax.text(xx[-1], K / xx[-1], f"K = {K}", fontsize=8, color="0.4", va="bottom", ha="right")
        ax.set_xscale("log")
        ax.set_yscale("log")
        for a_ in (ax.xaxis, ax.yaxis):
            a_.set_major_formatter(comma)
            a_.set_minor_formatter(comma)
        y = np.arange(len(names))
        ax2.barh(y, Ks, color=cs)
        for yi, K in zip(y, Ks):
            ax2.text(K + 0.3, yi, dez(K, 1), va="center", fontsize=9)
        ax2.axvspan(10, 30, color="0.9", zorder=0)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=9)
        ax2.invert_yaxis()
        _legend_alloys(ax, col)
    _finish(ax, "k_f / G", "d / b")
    _finish(ax2, "K = d · k_f / (G · b)")
    ax.set_title(L("Subkorngröße über der normierten Fließspannung", "Subgrain size versus normalised flow stress"), fontsize=12)
    ax2.set_title(L("Holt-Konstante je Datensatz (grau: üblicher Bereich 10–30)", "Holt constant per dataset (grey: usual range 10–30)"),
                  fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _note(fig, [L("d = Median-ECD der Subkörner; k_f aus den Versuchsdaten; G(T) = 25,4 GPa·(1 − 0,5·(T − 300 K)/933 K) mit der lokalen "
                  "Temperatur (sonst Prüftemperatur); b = 0,286 nm",
                  "d = median ECD of the subgrains; k_f from the test data; G(T) = 25.4 GPa·(1 − 0.5·(T − 300 K)/933 K) with the local "
                  "temperature (else test temperature); b = 0.286 nm")])
    if not names:
        _empty(fig, L("Für V9 Fließspannung k_f und Temperatur unter „Versuchsdaten …“ eintragen.",
                      "For V9, enter the flow stress k_f and the temperature under “Test data …”."))
    return fig


# ======================================================================== V10 geometrische DRX
def fig_v10(datasets):
    st, col = style(datasets)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), dpi=100)
    n = 0
    names, vals, cs = [], [], []
    for d in datasets:
        h, s_ = xm(d, "sehne_hagb_CD_um"), xm(d, "subkorn_ecd_median_um")
        if h is None or not s_:
            continue
        r = h / s_
        names.append(d.label)
        vals.append(r)
        cs.append(st[d.id]["color"])
        z, _ = lnZ(d)
        if z is not None:
            s = st[d.id]
            ax.plot(z, r, marker=s["marker"], ms=8, mfc=s["mfc"], mec=s["mec"], mew=1.6, ls="")
            ax.annotate(point_label(d, datasets), (z, r), xytext=(5, 4), textcoords="offset points", fontsize=7, color="0.35")
            n += 1
    ax.axhline(1, color="k", ls=":", lw=1)
    ax.set_ylim(0, None)
    _finish(ax, z_label(), L("HAGB-Sehnenlänge CD / Subkorngröße", "HAGB intercept CD / subgrain size"))
    ax.set_title(L("Geometrische DRX: Abstand der HAGB in Stauchrichtung", "Geometric DRX: HAGB spacing along the compression direction"),
                 fontsize=12)
    if n:
        _legend_alloys(ax, col)
    else:
        ax.text(0.5, 0.5, L("ln Z: Versuchsdaten eintragen", "ln Z: enter test data"), transform=ax.transAxes, ha="center", color="0.5")
    if names:
        y = np.arange(len(names))
        ax2.barh(y, vals, color=cs)
        for yi, v in zip(y, vals):
            ax2.text(v + 0.05, yi, dez(v, 2), va="center", fontsize=9)
        ax2.axvline(1, color="k", ls=":", lw=1)
        ax2.set_yticks(y)
        ax2.set_yticklabels(names, fontsize=9)
        ax2.invert_yaxis()
    _finish(ax2, L("HAGB-Sehnenlänge CD / Subkorngröße", "HAGB intercept CD / subgrain size"))
    ax2.set_title(L("je Datensatz", "per dataset"), fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _note(fig, [L("Werte nahe 1: die gestauchten Ausgangskörner sind so dünn wie die Subkörner – Hinweis auf geometrische DRX "
                  "(HAGB-Anteil steigt ohne Keimbildung).",
                  "Values close to 1: the flattened initial grains are as thin as the subgrains – indication of geometric DRX "
                  "(HAGB fraction rises without nucleation).")])
    return fig


# ======================================================================== V11 KAM
def fig_v11(datasets):
    from ebsd_datasets import colors_for
    cols = colors_for(datasets)
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), dpi=100)
    any_hist = False
    stat = []
    for d in datasets:
        h = xm(d, "kam_hist")
        if h:
            e = np.asarray(h["kanten"], float)
            c = np.asarray(h["anzahl"], float)
            ax.plot(0.5 * (e[:-1] + e[1:]), c / c.sum() / np.diff(e), color=cols[d.id], lw=1.8, label=d.label)
            any_hist = True
        stat.append((d, xm(d, "kam_mittel_deg"), xm(d, "kam_median_deg"), xm(d, "kam_q90_deg"), xm(d, "kam_klassen_deg")))
    if any_hist:
        ax.legend(frameon=False, fontsize=9)
        ax.set_xlim(0, 2.5)
        _finish(ax, "KAM (°)", L("Häufigkeitsdichte (1/°)", "Frequency density (1/°)"))
        ax.set_title(L("KAM-Verteilung (Auswertungen ab 1.4)", "KAM distribution (analyses from 1.4)"), fontsize=12)
    else:
        n = len(stat)
        x = np.arange(n)
        for k, (lab, j, mk) in enumerate([(L("Median", "median"), 2, "o"), (L("Mittelwert", "mean"), 1, "D"), (L("90%-Quantil", "90% quantile"), 3, "^")]):
            ax.plot(x, [_nanf(s[j]) for s in stat], mk, ms=8, color=["#7FA7C9", "#14466B", "#A23B2A"][k], label=lab)
        ax.set_xticks(x)
        ax.set_xticklabels([s[0].label.replace(", ", "\n") for s in stat], fontsize=8)
        ax.legend(frameon=False, fontsize=9)
        _finish(ax, None, "KAM (°)")
        ax.set_title(L("KAM-Kennwerte (Verteilung erst ab Auswertungen mit 1.4)", "KAM statistics (distribution only for analyses from 1.4)"),
                     fontsize=12)
    n = len(stat)
    x = np.arange(n)
    w = 0.26
    for k, (cls, lab, c) in enumerate([("rx", L("rekristallisiert", "recrystallized"), "#2B6FD6"),
                                       ("sub", L("substrukturiert", "substructured"), "#F2C744"),
                                       ("def", L("verformt", "deformed"), "#D63B2B")]):
        v = [_nanf((s[4] or {}).get(cls)) for s in stat]
        ax2.bar(x + (k - 1) * w, np.nan_to_num(v), w, color=c, label=lab)
    ax2.set_xticks(x)
    ax2.set_xticklabels([s[0].label.replace(", ", "\n") for s in stat], fontsize=8)
    ax2.legend(frameon=False, fontsize=9)
    _finish(ax2, None, L("mittlere KAM (°)", "mean KAM (°)"))
    ax2.set_title(L("mittlere KAM je GOS-Klasse", "Mean KAM per GOS class"), fontsize=12)
    steps = {d.metrics.get("schrittweite_um") for d in datasets}
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    _note(fig, [L("Achtung: KAM hängt von der Schrittweite ab – nur Karten gleicher Schrittweite vergleichen (siehe V12).",
                  "Note: KAM depends on the step size – compare only maps with the same step size (see V12).") if len(steps) > 1 else ""])
    return fig


# ======================================================================== V12 Vergleichbarkeit
def comparability(datasets):
    """Zeilen (Beschriftung, Werte, Warnung je Datensatz) und Warnungen"""
    rows = []
    def add(lab, vals, nd, warn=None):
        rows.append((lab, vals, nd, warn or [False] * len(vals)))
    area = [((d.metrics.get("flaeche_um") or [np.nan, np.nan])[0] * (d.metrics.get("flaeche_um") or [np.nan, np.nan])[1] / 1e3)
            for d in datasets]
    step = [_nanf(d.metrics.get("schrittweite_um")) for d in datasets]
    null = [_nanf(d.metrics.get("nullloesung_roh_pct")) for d in datasets]
    big = [sum(xm(d, "groesste_koerner_pct") or [np.nan]) for d in datasets]
    ncls = [_nanf(xm(d, "n_koerner_klassifizierbar")) for d in datasets]
    part = [_nanf(xm(d, "partikel_pct")) for d in datasets]
    sub_res = []
    for d in datasets:
        v = d.metrics.get("subkorn_aufloesung_schritte")
        if v is None and d.metrics.get("subkorn_ecd_median_um") and d.metrics.get("schrittweite_um"):
            v = d.metrics["subkorn_ecd_median_um"] / d.metrics["schrittweite_um"]
        sub_res.append(_nanf(v))
    warns = []
    smin, smax = np.nanmin(step), np.nanmax(step)
    add(L("Schrittweite (µm)", "Step size (µm)"), step, 3, [s < smax / 1.05 for s in step] if smax > smin * 1.05 else None)
    if smax > smin * 1.05:
        warns.append(L(f"Schrittweiten unterschiedlich ({dez(smin, 2)}–{dez(smax, 2)} µm): „Schrittweite angleichen …“ verwenden "
                       f"(Ziel {dez(smax, 2)} µm).",
                       f"Different step sizes ({dez(smin, 2)}–{dez(smax, 2)} µm): use “Match step size …” (target {dez(smax, 2)} µm)."))
    amed = np.nanmedian(area)
    small = [a < 0.5 * amed for a in area]
    add(L("Kartenfläche (10³ µm²)", "Map area (10³ µm²)"), area, 1, small)
    if any(small):
        warns.append(L("Mindestens eine Karte ist kleiner als die Hälfte der mittleren Kartenfläche (Median).",
                       "At least one map is smaller than half of the median map area."))
    add(L("Nulllösungen Rohdaten (%)", "Zero solutions raw data (%)"), null, 1, [v > 20 for v in null])
    add(L("3 größte Körner (% der Fläche)", "3 largest grains (% of the area)"), big, 1, [v > 50 for v in big])
    if any(v > 50 for v in big):
        warns.append(L("In mindestens einer Karte decken 3 Körner mehr als die Hälfte der Fläche ab – Korngrenzanteile und "
                       "Rekristallisation sind dort wenig repräsentativ.",
                       "In at least one map, 3 grains cover more than half of the area – boundary fractions and recrystallization "
                       "are hardly representative there."))
    add(L("klassifizierbare Körner", "classifiable grains"), ncls, 0, [v < 500 for v in ncls])
    add(L("Partikelanteil (%)", "Particle fraction (%)"), part, 2)
    add(L("Subkorn-Median in Schrittweiten", "Subgrain median in step sizes"), sub_res, 1, [v < 5 for v in sub_res])
    if any(v < 5 for v in sub_res):
        warns.append(L("Subkörner in mindestens einer Karte kleiner als 5 Schrittweiten – kaum aufgelöst.",
                       "Subgrains smaller than 5 step sizes in at least one map – hardly resolved."))
    rs = [d.metrics.get("resample_um") for d in datasets]
    add(L("Schrittweite angeglichen von (µm)", "Step size resampled from (µm)"),
        [_nanf(d.metrics.get("schrittweite_original_um")) if r else np.nan for d, r in zip(datasets, rs)], 3)
    return rows, warns


def fig_v12(datasets):
    rows, warns = comparability(datasets)
    n = len(datasets)
    fw = max(10, 4.2 + 2.0 * n)
    import textwrap
    warns = [textwrap.fill("⚠ " + w, width=int(fw * 11.5), subsequent_indent="   ") for w in warns]
    nl = sum(w.count("\n") + 1 for w in warns)
    fig = plt.figure(figsize=(fw, 1.6 + 0.42 * len(rows) + 0.32 * nl + 0.6), dpi=100)
    ax = fig.add_axes([0.01, 0.0, 0.98, 1.0])
    ax.axis("off")
    FW, FH = fig.get_size_inches()
    y = FH - 0.6
    ax.text(0.01, y / FH, L("Vergleichbarkeit der Karten", "Comparability of the maps"), fontsize=14, weight="bold", transform=ax.transAxes)
    y -= 0.55
    xs = [0.01] + [(3.9 + 2.0 * j) / FW for j in range(n)]
    for j, d in enumerate(datasets):
        ax.text(xs[j + 1], y / FH, "\n".join(d.label.split(", ")[:3]), fontsize=8.5, weight="bold", transform=ax.transAxes, va="top")
    y -= 0.62
    for lab, vals, nd, wr in rows:
        ax.text(xs[0], y / FH, lab, fontsize=9.5, transform=ax.transAxes, va="center")
        for j, (v, w) in enumerate(zip(vals, wr)):
            ax.text(xs[j + 1], y / FH, ("⚠ " if w else "") + dez(v, nd), fontsize=9.5, transform=ax.transAxes, va="center",
                    color="#A23B2A" if w else "black")
        y -= 0.42
    y -= 0.15
    for w in warns:
        k = w.count("\n")
        ax.text(xs[0], (y + 0.16) / FH, w, fontsize=9.5, color="#A23B2A", transform=ax.transAxes, va="top", linespacing=1.35)
        y -= 0.32 * (k + 1)
    if not warns:
        ax.text(xs[0], y / FH, L("✓ keine Auffälligkeiten", "✓ no issues"), fontsize=10, color="#3C7A3C", transform=ax.transAxes, va="center")
    return fig


# ======================================================================== V13 Sensitivität gegenüber θc
def fig_v13(datasets):
    from ebsd_datasets import colors_for
    cols = colors_for(datasets)
    fig, ax = plt.subplots(figsize=(10, 6.2), dpi=100)
    n = 0
    for d in datasets:
        g = xm(d, "gos_sens")
        if not g:
            continue
        ax.plot(g["grenzen_deg"], g["rx_pct"], "-o", ms=4, lw=1.8, color=cols[d.id], label=d.label)
        n += 1
    thc = {d.metrics.get("theta_c") for d in datasets}
    for t in thc:
        if t:
            ax.axvline(t, color="0.5", ls="--", lw=1)
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    _finish(ax, L("GOS-Grenzwert θc (°)", "GOS threshold θc (°)"), L("rekristallisierter Anteil (%)", "recrystallized fraction (%)"))
    ax.set_title(L("Sensitivität des rekristallisierten Anteils gegenüber θc", "Sensitivity of the recrystallized fraction to θc"),
                 fontsize=12)
    if n:
        ax.legend(frameon=False, fontsize=9)
        fig.tight_layout(rect=(0, 0.05, 1, 1))
        _note(fig, [L("Kreuzen sich die Kurven nicht, ist die Reihenfolge der Zustände unabhängig vom gewählten θc.",
                      "If the curves do not cross, the order of the states does not depend on the chosen θc.")])
    else:
        _empty(fig, L("Keine Sensitivitätsdaten in den gewählten Datensätzen.", "No sensitivity data in the selected datasets."))
    return fig


# ======================================================================== V2 erweitert: Summenkurve und regellose Verteilung
_MACKENZIE = None


def mackenzie(width=2.0, n=400_000, seed=1):
    """Misorientierungswinkel regellos orientierter kubischer Kristalle (Monte-Carlo, 2°-Klassen, % je Klasse)"""
    global _MACKENZIE
    if _MACKENZIE is not None:
        return _MACKENZIE
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1)[:, None]
    # 24 Drehungen der kubischen Punktgruppe als Quaternionen
    s2 = np.sqrt(0.5)
    ops = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
    ops += [(s2, s2, 0, 0), (s2, -s2, 0, 0), (s2, 0, s2, 0), (s2, 0, -s2, 0), (s2, 0, 0, s2), (s2, 0, 0, -s2)]
    ops += [(0, s2, s2, 0), (0, s2, -s2, 0), (0, s2, 0, s2), (0, s2, 0, -s2), (0, 0, s2, s2), (0, 0, s2, -s2)]
    ops += [(0.5, a, b, c) for a in (0.5, -0.5) for b in (0.5, -0.5) for c in (0.5, -0.5)]
    S = np.array(ops, float)
    w = np.abs(q @ S.T)                       # |Skalarteil| von s·q (Betrag des Skalarprodukts)
    ang = 2 * np.degrees(np.arccos(np.clip(w.max(axis=1), -1, 1)))
    edges = np.arange(0, 64 + width, width)
    h, _ = np.histogram(ang, bins=edges)
    _MACKENZIE = (0.5 * (edges[:-1] + edges[1:]), h / n * 100)
    return _MACKENZIE


BUILDERS = {"V3": fig_v3, "V4": fig_v4, "V5": fig_v5, "V6": fig_v6, "V7": fig_v7, "V8": fig_v8, "V9": fig_v9,
            "V10": fig_v10, "V11": fig_v11, "V12": fig_v12, "V13": fig_v13}
