# -*- coding: utf-8 -*-
"""
Gespeicherte Datensätze und Vergleichsgrafiken der EBSD-Auswertung (ohne MTEX, damit die Oberfläche
die Datensätze sofort beim Start anzeigen kann).

Jeder ausgewertete Datensatz (ctf-Datei + Ausschnitt) wird in einem eigenen Ordner abgelegt:
    datensatz.json        Bezeichnung, Datei, Rollen, Einstellungen, Kennwerte für den Vergleich
    Kennwerte_und_Parameter.txt, Kornliste_*.csv
    vorschau/<Kennung>.png  Vorschauen der Abbildungen (wie im Reiter „Ergebnisse“)
"""
from __future__ import annotations

import hashlib
import io
import json
import re
import shutil
import textwrap
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

from ebsd_settings import APP_NAME, APP_VERSION, compare_figures, fig_fname, kennwerte_name
from ebsd_i18n import L, csv_join, dez, get_language, num_g

comma = FuncFormatter(lambda v, p: num_g(v))
KENNWERTE_FILE = "Kennwerte_und_Parameter.txt"          # interner Dateiname in der Ablage (sprachunabhängig)
TABLE_GLOBS = ("Kornliste_*.csv", "Grain_list_*.csv")
PALETTE = ["#7FA7C9", "#14466B", "#E3A35B", "#A23B2A", "#6A9F58", "#8E6BB0", "#C9A227", "#4D4D4D",
           "#D17AA5", "#3FA7A3", "#9C755F", "#B7B7B7"]


def _nan(x):
    return np.nan if x is None else float(x)


# ======================================================================== Bezeichnung
def default_label(path, crop_um=None):
    """Kurzbezeichnung aus Datei- und Ordnernamen, z. B. „750-3, 350 °C / 0,1 s⁻¹, Mitte“"""
    p = Path(path)
    parts = [p.stem] + [q.name for q in list(p.parents)[:3]]
    lab = None
    for t in parts:
        m = re.search(r"(\d{3,4}-\d+)_(\d{3})_(\d+(?:[.,]\d+)?)(?:_|$)", t)
        if m:
            lab = f"{m.group(1)}, {m.group(2)} °C / {m.group(3).replace('.', ',')} s⁻¹"
            break
    if lab is None:
        stem = re.sub(r"[_ -]*EBSD-Daten$", "", p.stem)
        return stem + (L(" (Ausschnitt)", " (crop)") if crop_um else "")
    low = p.stem.lower()
    for key, short in (("probenmitte", L("Mitte", "centre")), ("mitte", L("Mitte", "centre")),
                       ("schmiedekreuz", L("Schmiedekreuz", "forging cross")), ("rand", L("Rand", "edge")),
                       ("kern", L("Kern", "core"))):
        if key in low:
            lab += f", {short}"
            break
    else:
        if not re.search(r"\d{3,4}-\d+_\d{3}_", p.parent.name) and len(p.parent.name) <= 12:
            lab += f", {p.parent.name}"              # Unterordner wie „A4“ als Position
    if crop_um:
        lab += L(" (Ausschnitt)", " (crop)")
    return lab


def dataset_key(path, crop_um):
    crop = "" if crop_um is None else ",".join(f"{float(v):.1f}" for v in crop_um)
    return hashlib.sha1(f"{str(Path(path)).lower()}|{crop}".encode("utf-8")).hexdigest()[:12]


# ======================================================================== Datensätze
class Dataset:
    def __init__(self, folder: Path, meta: dict):
        self.folder = Path(folder)
        self.meta = meta

    id = property(lambda self: self.meta["id"])
    label = property(lambda self: self.meta.get("bezeichnung") or self.meta["kennwerte"].get("datei", self.id))
    metrics = property(lambda self: self.meta.get("kennwerte", {}))
    path = property(lambda self: Path(self.meta["kennwerte"].get("pfad", "")))
    crop_um = property(lambda self: self.meta["kennwerte"].get("ausschnitt_um"))
    in_compare = property(lambda self: bool(self.meta.get("im_vergleich", True)))

    @property
    def figures(self):
        return [tuple(f) for f in self.meta.get("abbildungen", [])]

    def preview(self, fid):
        f = self.folder / "vorschau" / f"{fid}.png"
        return f.read_bytes() if f.exists() else None

    def kennwerte(self):
        f = self.folder / KENNWERTE_FILE
        return f.read_text(encoding="utf-8") if f.exists() else ""

    def tables(self):
        return {f.name: f.read_text(encoding="utf-8-sig") for g in TABLE_GLOBS for f in sorted(self.folder.glob(g))}

    @property
    def strokes(self):
        """manuelle Kratzermarkierungen der Auswertung (ab Version 1.3)"""
        return list(self.metrics.get("kratzer_manuell") or [])

    def info(self):
        k = self.metrics
        t = time.strftime(L("%d.%m.%Y %H:%M", "%Y-%m-%d %H:%M"), time.localtime(self.meta.get("geaendert", 0)))
        crop = k.get("ausschnitt_um")
        pts = k.get("punkte", ["?", "?"])
        sub = subgrains_on(self)
        n_str = len(self.strokes)
        return (f"{k.get('datei', '')}"
                + (L("  |  Ausschnitt ", "  |  crop ") + f"{'; '.join(dez(v, 1) for v in crop)} µm" if crop else "")
                + "\n" + L(f"{pts[0]} × {pts[1]} Punkte, Schrittweite {dez(k.get('schrittweite_um'), 3)} µm  |  ausgewertet {t}",
                           f"{pts[0]} × {pts[1]} points, step size {dez(k.get('schrittweite_um'), 3)} µm  |  analysed {t}")
                + ("" if sub is not False else L("  |  ohne Subkörner", "  |  without subgrains"))
                + (L(f"  |  {n_str} manuelle Kratzermarkierung(en)", f"  |  {n_str} manual scratch mark(s)") if n_str else ""))

    def save_meta(self):
        tmp = self.folder / "datensatz.json.tmp"
        tmp.write_text(json.dumps(self.meta, indent=1, ensure_ascii=False), encoding="utf-8")
        tmp.replace(self.folder / "datensatz.json")


class DatasetStore:
    """Ordner mit gespeicherten Datensätzen; Reihenfolge in reihenfolge.json"""

    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _order(self):
        f = self.root / "reihenfolge.json"
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            return []

    def set_order(self, ids):
        (self.root / "reihenfolge.json").write_text(json.dumps(list(ids)), encoding="utf-8")

    def list(self):
        out = {}
        for d in self.root.iterdir():
            f = d / "datensatz.json"
            if d.is_dir() and f.exists():
                try:
                    out[d.name] = Dataset(d, json.loads(f.read_text(encoding="utf-8")))
                except Exception:  # noqa: BLE001
                    continue
        order = [i for i in self._order() if i in out]
        rest = sorted((i for i in out if i not in order), key=lambda i: out[i].meta.get("erstellt", 0))
        return [out[i] for i in order + rest]

    def get(self, ds_id):
        f = self.root / ds_id / "datensatz.json"
        return Dataset(f.parent, json.loads(f.read_text(encoding="utf-8"))) if f.exists() else None

    def save_result(self, res, previews, label=None):
        """legt das Ergebnis einer Auswertung ab (gleiche Datei und gleicher Ausschnitt: wird ersetzt)"""
        k = dict(res.metrics)
        ds_id = dataset_key(k.get("pfad", res.project.path), k.get("ausschnitt_um"))
        folder = self.root / ds_id
        old = self.get(ds_id)
        if folder.exists():
            shutil.rmtree(folder / "vorschau", ignore_errors=True)
            for g in TABLE_GLOBS:
                for f in folder.glob(g):
                    f.unlink()
        (folder / "vorschau").mkdir(parents=True, exist_ok=True)
        for fid, b in previews.items():
            (folder / "vorschau" / f"{fid}.png").write_bytes(b)
        (folder / KENNWERTE_FILE).write_text(res.kennwerte, encoding="utf-8")
        for name, text in res.tables.items():
            (folder / name).write_text(text, encoding="utf-8-sig")
        now = time.time()
        meta = {"id": ds_id, "programm": f"{APP_NAME} {APP_VERSION}", "sprache": get_language(),
                "bezeichnung": label or (old.meta.get("bezeichnung") if old else None)
                               or default_label(k.get("pfad", ""), k.get("ausschnitt_um")),
                "im_vergleich": old.in_compare if old else True,
                "erstellt": old.meta.get("erstellt", now) if old else now, "geaendert": now,
                "abbildungen": [list(a) for a in res.available], "kennwerte": k}
        ds = Dataset(folder, meta)
        ds.save_meta()
        if ds_id not in self._order():
            self.set_order([d.id for d in self.list() if d.id != ds_id] + [ds_id])
        return ds

    def rename(self, ds_id, label):
        ds = self.get(ds_id)
        if ds:
            ds.meta["bezeichnung"] = label.strip() or default_label(ds.path, ds.crop_um)
            ds.save_meta()

    def set_compare(self, ds_id, on):
        ds = self.get(ds_id)
        if ds and ds.in_compare != bool(on):
            ds.meta["im_vergleich"] = bool(on)
            ds.save_meta()

    def remove(self, ds_id):
        shutil.rmtree(self.root / ds_id, ignore_errors=True)
        self.set_order([i for i in self._order() if i != ds_id])


# ======================================================================== Vergleich
def subgrains_on(ds):
    """Subkornerkennung des Datensatzes: True/False; Datensätze vor 1.3 hatten sie immer an"""
    k = ds.metrics
    if "subkoerner" in k:
        return bool(k["subkoerner"])
    return True if "einstellungen" in k else None


METRIC_ROWS = [  # (Beschriftung, Beschriftung engl., Schlüssel, Nachkommastellen)
    ("Schrittweite (µm)", "Step size (µm)", "schrittweite_um", 3),
    ("Subkornerkennung", "Subgrain detection", "_subkoerner", None),
    ("Kratzer automatisch erkannt", "Scratches detected automatically", "kratzer_n", 0),
    ("manuelle Kratzermarkierungen", "Manual scratch marks", "_kratzer_manuell", 0),
    ("Kratzerfläche (%)", "Scratch area (%)", "kratzer_flaeche_pct", 2),
    ("Nulllösungen Rohdaten (%)", "Zero solutions raw data (%)", "nullloesung_roh_pct", 2),
    ("Fehlindizierung Rohdaten (%)", "Misindexing raw data (%)", "fehlindizierung_roh_pct", 2),
    ("Nulllösungs-Cluster (Mg2Si) (%)", "Zero-solution clusters (Mg2Si) (%)", "cluster_pct", 2),
    ("HAGB-Anteil der Grenzlänge (%)", "HAGB fraction of boundary length (%)", "hagb_anteil_pct", 1),
    ("LAGB-Anteil der Grenzlänge (%)", "LAGB fraction of boundary length (%)", "lagb_anteil_pct", 1),
    ("Sehnenlänge HAGB RD (µm)", "Mean intercept HAGB RD (µm)", "sehne_hagb_RD_um", 2),
    ("Sehnenlänge HAGB CD (µm)", "Mean intercept HAGB CD (µm)", "sehne_hagb_CD_um", 2),
    ("Sehnenlänge Grenzen >= LAGB RD (µm)", "Mean intercept boundaries >= LAGB RD (µm)", "sehne_alle_RD_um", 2),
    ("Sehnenlänge Grenzen >= LAGB CD (µm)", "Mean intercept boundaries >= LAGB CD (µm)", "sehne_alle_CD_um", 2),
    ("ASTM-Korngrößenkennzahl G", "ASTM grain size number G", "astm_G", 1),
    ("Korn-ECD Mittelwert (µm)", "Grain ECD mean (µm)", "korn_ecd_mittel_um", 2),
    ("Korn-ECD Median (µm)", "Grain ECD median (µm)", "korn_ecd_median_um", 2),
    ("Subkorn-ECD Mittelwert (µm)", "Subgrain ECD mean (µm)", "subkorn_ecd_mittel_um", 2),
    ("Subkorn-ECD Median (µm)", "Subgrain ECD median (µm)", "subkorn_ecd_median_um", 2),
    ("KAM Mittelwert (°)", "KAM mean (°)", "kam_mittel_deg", 3),
    ("rekristallisiert (%)", "recrystallized (%)", "rx_pct", 2),
    ("nicht rekristallisiert (%)", "non-recrystallized (%)", "_nicht_rx", 2),
    ("substrukturiert (%)", "substructured (%)", "sub_pct", 2),
    ("verformt (%)", "deformed (%)", "def_pct", 2),
    ("rekristallisiert 0–1 µm von Partikeln (%)", "recrystallized 0–1 µm from particles (%)", ("rx_abstand_pct", 0), 1),
    ("rekristallisierte Körner ECD Median (µm)", "recrystallized grains ECD median (µm)", "rx_korn_ecd_median_um", 2),
    ("<110> || CD (15°) (%)", "<110> || CD (15°) (%)", "faser110_CD_pct", 1),
]


def metric(ds, key):
    k = ds.metrics
    if isinstance(key, tuple):
        v = k.get(key[0])
        return v[key[1]] if isinstance(v, list) and len(v) > key[1] else None
    if key == "_subkoerner":
        return subgrains_on(ds)
    if key == "_kratzer_manuell":
        return len(ds.strokes) if "kratzer_manuell" in k else None
    if key == "_nicht_rx":
        if k.get("nicht_rx_pct") is not None:
            return k["nicht_rx_pct"]
        parts = [k.get("sub_pct"), k.get("def_pct")]
        return sum(p for p in parts if p is not None) if any(p is not None for p in parts) else None
    return k.get(key)


def _cell(v, nd):
    if isinstance(v, bool):
        return L("an", "on") if v else L("aus", "off")
    return dez(v, nd)


def metrics_table(datasets):
    """Zeilen (Beschriftung, [Werte als Text]) für Tabelle und CSV"""
    return [(L(lab, lab_en), [_cell(metric(d, key), nd) for d in datasets]) for lab, lab_en, key, nd in METRIC_ROWS]


def metrics_csv(datasets):
    lines = [csv_join([L("Kennwert", "Quantity")] + [d.label for d in datasets]),
             csv_join([L("Datei", "File")] + [d.metrics.get("datei", "") for d in datasets])]
    lines += [csv_join([lab] + ["" if v == "–" else v for v in vals]) for lab, vals in metrics_table(datasets)]
    centers, curves = misori_curves(datasets)
    if curves:
        lines += ["", L("MISORIENTIERUNGSVERTEILUNG (Anteil der Grenzen >= LAGB, %, Klassenmitte in °)",
                        "MISORIENTATION DISTRIBUTION (fraction of boundaries >= LAGB, %, class centre in °)"),
                  csv_join([L("Klassenmitte (°)", "Class centre (°)")] + [d.label for d in datasets])]
        for i, c in enumerate(centers):
            lines.append(csv_join([dez(c, 0)] + [dez(curves[d.id][i], 2) if d.id in curves else "" for d in datasets]))
    return "\n".join(lines) + "\n"


def colors_for(datasets):
    return {d.id: PALETTE[i % len(PALETTE)] for i, d in enumerate(datasets)}


def _common(datasets, key):
    vals = {str(d.metrics.get(key)) for d in datasets}
    return datasets[0].metrics.get(key) if len(vals) == 1 else None


def _title(datasets, what):
    t = _common(datasets, "titel")
    return f"{t} – {what}" if t else what


def misori_curves(datasets, width=2.0):
    """2°-Klassen ab der LAGB-Grenze; Anteil an allen Grenzschnitten >= LAGB in %"""
    lagb = min((d.metrics.get("lagb") or 2.0) for d in datasets) if datasets else 2.0
    edges = np.arange(lagb, 62.0 + 1e-9, width)
    if edges[-1] < 62.8:
        edges = np.append(edges, 62.8)
    centers = 0.5 * (edges[:-1] + edges[1:])
    curves = {}
    for d in datasets:
        h = d.metrics.get("misori_hist")
        if not h:
            continue
        e = np.asarray(h["kanten"], float)
        n = np.asarray(h["anzahl"], float)
        mids = 0.5 * (e[:-1] + e[1:])
        tot = n[mids >= (d.metrics.get("lagb") or lagb)].sum()
        if tot <= 0:
            continue
        idx = np.digitize(mids, edges) - 1
        cnt = np.bincount(idx[(idx >= 0) & (idx < len(centers))], weights=n[(idx >= 0) & (idx < len(centers))],
                          minlength=len(centers))
        curves[d.id] = cnt / tot * 100
    return centers, curves


def tick_labels(datasets):
    """Achsbeschriftungen: gemeinsamen ersten Teil (z. B. Legierung) weglassen, Rest zeilenweise umbrechen"""
    parts = [[p.strip() for p in d.label.split(", ")] for d in datasets]
    common = parts[0][0] if parts and len(parts[0]) > 1 and all(len(p) > 1 and p[0] == parts[0][0] for p in parts) else None
    out = []
    for p in parts:
        rest = p[1:] if common else p
        lines = []
        for chunk in rest:
            for sub in chunk.split(" / "):
                lines += textwrap.wrap(sub, 13, break_on_hyphens=False) or [sub]
        out.append("\n".join(lines))
    return out, common


def fig_compare_kennwerte(datasets):
    n = len(datasets)
    col = colors_for(datasets)
    cols = [col[d.id] for d in datasets]
    w = max(18.0, min(32.0, 4 * (2.0 + 0.7 * n)))
    fig, axs = plt.subplots(1, 4, figsize=(w, 5.6), dpi=100)
    x = np.arange(n)
    labels, common = tick_labels(datasets)
    fs = 10 if n <= 4 else (9 if n <= 6 else 8)
    rot = 0 if n <= 6 else 35

    def bars(ax, vals, title, ylabel, nd=1, offset=0.0, width=0.8, **kw):
        v = np.array([_nan(a) for a in vals])
        b = ax.bar(x + offset, np.nan_to_num(v), width, color=kw.pop("color", cols), **kw)
        for xi, vi in zip(x + offset, v):
            ax.text(xi, (0 if np.isnan(vi) else vi), dez(vi, nd), ha="center", va="bottom", fontsize=9)
        ax.set_title(title, fontsize=12)
        ax.set_ylabel(ylabel)
        return b

    lagb, thc = _common(datasets, "lagb"), _common(datasets, "theta_c")
    bars(axs[0], [metric(d, "subkorn_ecd_median_um") for d in datasets],
         L("Subkorngröße", "Subgrain size") + f" ({L('Median-ECD', 'median ECD')}, ≥ {dez(lagb, 0) if lagb else 'LAGB'}°)", "µm")
    if any(subgrains_on(d) is False for d in datasets):
        axs[0].text(0.5, 0.97, L("„–“: ohne Subkornerkennung", "“–”: without subgrain detection"), transform=axs[0].transAxes,
                    ha="center", va="top", fontsize=9, color="0.4")
    bars(axs[1], [metric(d, "hagb_anteil_pct") for d in datasets], L("HAGB-Anteil der Grenzlänge", "HAGB fraction of boundary length"), "%")
    rx_title = L("Rekristallisiert", "Recrystallized") + f" (GOS ≤ {dez(thc, 0) if thc else 'θc'}°)"
    bars(axs[2], [metric(d, "rx_pct") for d in datasets], rx_title, "%", offset=-0.2, width=0.4)
    bars(axs[2], [metric(d, ("rx_abstand_pct", 0)) for d in datasets], "", "%", offset=0.2, width=0.4,
         alpha=0.45, hatch="//", edgecolor="k", linewidth=0)
    axs[2].set_title(rx_title, fontsize=12)
    axs[2].legend(handles=[Patch(fc="0.55", label=L("gesamt", "total")),
                           Patch(fc="0.85", hatch="//", ec="k", lw=0, label=L("0–1 µm von Partikeln", "0–1 µm from particles"))],
                  frameon=False, fontsize=9, loc="upper right")
    bars(axs[3], [metric(d, "faser110_CD_pct") for d in datasets], L("<110> ∥ Stauchrichtung (±15°)", "<110> ∥ compression direction (±15°)"),
         L("% (regellos ≈ 20)", "% (random ≈ 20)"))
    axs[3].axhline(20, color="k", ls=":", lw=1)
    for ax in axs:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=fs, rotation=rot, ha="right" if rot else "center")
        ax.set_xlim(-0.6, n - 0.4)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_formatter(comma)
        lo, hi = ax.get_ylim()
        ax.set_ylim(0, hi * 1.08 if hi > 0 else 1)
    fig.suptitle(_title(datasets, L("EBSD-Vergleich", "EBSD comparison") + (f" {common}" if common else "")), fontsize=14)
    fig.tight_layout()
    return fig


def fig_compare_misori(datasets):
    centers, curves = misori_curves(datasets)
    col = colors_for(datasets)
    fig, ax = plt.subplots(figsize=(10, 6.2), dpi=100)
    for d in datasets:
        if d.id in curves:
            ax.plot(centers, curves[d.id], "-o", ms=3.5, lw=1.8, color=col[d.id], label=d.label)
    hagb = _common(datasets, "hagb")
    for h in ({hagb} if hagb else {d.metrics.get("hagb") for d in datasets if d.metrics.get("hagb")}):
        ax.axvline(h, color="0.4", ls="--", lw=1.2)
    lagb = _common(datasets, "lagb")
    ax.set_xlabel(L("Misorientierung (°)", "Misorientation (°)"))
    ax.set_ylabel(L("Anteil der Grenzen", "Fraction of boundaries") + f" ≥ {dez(lagb, 0) if lagb else 'LAGB'}° (%)")
    ax.set_xlim(0, 64)
    ax.set_ylim(0, None)
    ax.grid(alpha=0.3)
    ax.xaxis.set_major_formatter(comma)
    ax.yaxis.set_major_formatter(comma)
    ax.set_title(_title(datasets, L("Misorientierungsverteilung", "Misorientation distribution"))
                 + L("\nLinienschnitt RD + CD, Längenanteil, 2°-Klassen", "\nLine intercept RD + CD, length fraction, 2° classes"),
                 fontsize=12)
    if curves:
        ax.legend(frameon=False, fontsize=10)
    else:
        ax.text(0.5, 0.5, L("Keine Korngrenzauswertung in den gewählten Datensätzen",
                            "No grain boundary analysis in the selected datasets"), transform=ax.transAxes,
                ha="center", va="center", color="0.4")
    fig.tight_layout()
    return fig


COMPARE_BUILDERS = {"V1": fig_compare_kennwerte, "V2": fig_compare_misori}


def compare_png(fid, datasets, dpi=100):
    fig = COMPARE_BUILDERS[fid](datasets)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    plt.close(fig)
    return buf.getvalue()


def compare_csv_name():
    return L("Vergleich_Kennwerte.csv", "Comparison_key_figures.csv")


def export_compare(datasets, folder, formats, dpi=300, fids=None):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for fid, _, fname in compare_figures():
        if fids and fid not in fids:
            continue
        fig = COMPARE_BUILDERS[fid](datasets)
        for fmt in formats:
            path = folder / f"{fname}.{'tif' if fmt == 'tiff' else fmt}"
            kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
            fig.savefig(path, format=fmt, dpi=dpi, **kw)
            written.append(path)
        plt.close(fig)
    path = folder / compare_csv_name()
    path.write_text(metrics_csv(datasets), encoding="utf-8-sig")
    written.append(path)
    return written


def export_stored(ds, folder):
    """Kennwerte, Kornlisten und Vorschaubilder eines gespeicherten Datensatzes kopieren (ohne Neuberechnung)"""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    written = []
    for fid, _, _ in ds.figures:
        b = ds.preview(fid)
        if b:
            p = folder / f"{fig_fname(fid)}{L('_Vorschau', '_preview')}.png"
            p.write_bytes(b)
            written.append(p)
    p = folder / kennwerte_name()
    p.write_text(ds.kennwerte(), encoding="utf-8")
    written.append(p)
    for name, text in ds.tables().items():
        p = folder / name
        p.write_text(text, encoding="utf-8-sig")
        written.append(p)
    return written
