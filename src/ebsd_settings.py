# -*- coding: utf-8 -*-
"""Einstellungen, Parameterliste und Abbildungsliste der EBSD-Auswertung (ohne schwere Importe,
damit die Oberfläche sofort startet)."""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields

from ebsd_i18n import L

APP_NAME = "EBSD-Auswertung"
APP_VERSION = "1.3"
ROLES = ("Matrix", "Partikel", "Fehlindizierung")          # interne Kennungen (so auch in gespeicherten Datensätzen)
ROLE_NAMES = {"Matrix": ("Matrix", "Matrix"), "Partikel": ("Partikel", "Particle"),
              "Fehlindizierung": ("Fehlindizierung", "Misindexing")}


def role_name(role):
    """angezeigter Name einer Phasenrolle"""
    de, en = ROLE_NAMES.get(role, (role, role))
    return L(de, en)


def app_title():
    return L("EBSD-Auswertung", "EBSD Analysis")


class Abgebrochen(Exception):
    """Auswertung wurde vom Benutzer abgebrochen"""


# ======================================================================== Einstellungen
@dataclass
class Settings:
    title: str = "AW-6082 (Legierung 750-3)"
    # Module
    scratch: bool = True          # Kratzer automatisch erkennen und entfernen
    mg2si: bool = True            # große Nulllösungs-Cluster als Mg2Si werten
    ipf: bool = True              # IPF-Karten
    kam: bool = True              # KAM-Karte
    gb: bool = True               # Korngrenz-Misorientierung (Linienschnitt)
    texture: bool = True          # inverse Polfiguren, Polfiguren
    rx: bool = True               # Rekristallisation (GOS)
    grainsize: bool = True        # Korngröße, Subkorngröße, Aspektverhältnis
    subgrains: bool = True        # Subkörner erkennen (aus, wenn die Schrittweite für Subkörner zu groß ist)
    # Parameter
    lagb: float = 2.0
    hagb: float = 15.0
    kam_max: float = 5.0
    kam_vmax: float = 2.0
    gap_bridge: int = 5
    mg2si_min_area: float = 1.4
    mg2si_max_bc: float = 100.0
    particle_min_area: float = 0.28
    pf_halfwidth: float = 5.0
    rx_theta_c: float = 2.0
    rx_min_px: int = 10
    psn_min_ecd: float = 1.0
    gs_min_px: int = 10
    scratch_min_len: float = 40.0
    scratch_rel_bc: float = 0.75
    scratch_window: float = 9.0

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        names = {f.name: f.type for f in fields(cls)}
        kw = {}
        for k, v in (d or {}).items():
            if k in names:
                kw[k] = v
        s = cls(**kw)
        for f in fields(cls):                              # Typen absichern (JSON/QSettings)
            v = getattr(s, f.name)
            default = getattr(cls(), f.name)
            try:
                if isinstance(default, bool):
                    v = v if isinstance(v, bool) else str(v).lower() in ("1", "true", "yes")
                elif isinstance(default, int):
                    v = int(v)
                elif isinstance(default, float):
                    v = float(v)
                else:
                    v = str(v)
            except (TypeError, ValueError):
                v = default
            setattr(s, f.name, v)
        return s


# (Feld, Beschriftung, Einheit, min, max, Schritt, Nachkommastellen, Hinweis, Beschriftung engl., Hinweis engl.)
PARAM_INFO = [
    ("lagb", "LAGB, untere Grenze", "°", 0.5, 10, 0.5, 1, "Kleinwinkelkorngrenzen ab diesem Winkel",
     "LAGB, lower limit", "Low-angle grain boundaries from this angle"),
    ("hagb", "HAGB, Grenze", "°", 5, 30, 1, 1, "Großwinkelkorngrenzen ab diesem Winkel; Kornrekonstruktion",
     "HAGB, threshold", "High-angle grain boundaries from this angle; grain reconstruction"),
    ("kam_max", "KAM-Schwellwert", "°", 1, 15, 0.5, 1, "Nachbarpaare darüber werden in der KAM ignoriert",
     "KAM threshold", "Neighbour pairs above this angle are ignored in the KAM"),
    ("kam_vmax", "KAM-Farbskala bis", "°", 0.5, 10, 0.25, 2, "Obergrenze der Farbskala der KAM-Karte",
     "KAM colour scale up to", "Upper limit of the colour scale of the KAM map"),
    ("gap_bridge", "Linienschnitt: max. Lücke", "Punkte", 0, 20, 1, 0, "überbrückte nicht indizierte Punkte",
     "Line intercept: max. gap", "bridged non-indexed points"),
    ("mg2si_min_area", "Nulllösungs-Cluster, Mindestfläche", "µm²", 0.1, 50, 0.1, 2, "kleinere Nulllösungen werden ergänzt",
     "Zero-solution cluster, min. area", "smaller zero solutions are filled in from the neighbours"),
    ("mg2si_max_bc", "Nulllösungs-Cluster, max. mittl. BC", "", 10, 255, 5, 0, "Cluster mit höherem Bandkontrast werden ergänzt",
     "Zero-solution cluster, max. mean BC", "clusters with higher band contrast are filled in"),
    ("particle_min_area", "Partikel, Mindestfläche", "µm²", 0.0, 10, 0.02, 2, "kleinere Partikelcluster gelten als Fehlindizierung",
     "Particles, min. area", "smaller particle clusters are treated as misindexing"),
    ("pf_halfwidth", "ODF-Halbwertsbreite", "°", 2, 20, 0.5, 1, "MTEX-Kernel (HWHM); 5° ≙ FWHM 10°",
     "ODF half-width", "MTEX kernel (HWHM); 5° ≙ FWHM 10°"),
    ("rx_theta_c", "θc (GOS-Grenze, Subkorngrenze)", "°", 0.5, 10, 0.25, 2, "rekristallisiert: GOS ≤ θc",
     "θc (GOS threshold, subgrain boundary)", "recrystallized: GOS ≤ θc"),
    ("rx_min_px", "min. Punkte je (Sub-)Korn", "Punkte", 3, 200, 1, 0, "kleinere Körner: nicht klassifizierbar",
     "min. points per (sub)grain", "smaller grains: not classifiable"),
    ("gs_min_px", "Korngröße: min. Punkte je (Sub-)Korn", "Punkte", 1, 500, 1, 0,
     "kleinere (Sub-)Körner und Randkörner gehen nicht in die Korngrößenstatistik ein",
     "Grain size: min. points per (sub)grain", "smaller (sub)grains and edge grains are excluded from the grain size statistics"),
    ("psn_min_ecd", "PSN: Partikel ab", "µm", 0.2, 10, 0.1, 1, "Mindestgröße der Partikel für den PSN-Indikator",
     "PSN: particles from", "Minimum particle size for the PSN indicator"),
    ("scratch_min_len", "Kratzer: Mindestlänge", "µm", 5, 500, 5, 0, "geschlossener Abschnitt",
     "Scratches: min. length", "continuous section"),
    ("scratch_rel_bc", "Kratzer: BC-Schwelle relativ", "", 0.3, 0.95, 0.05, 2, "BC < Wert × Umgebung = Kandidat",
     "Scratches: relative BC threshold", "BC < value × surroundings = candidate"),
    ("scratch_window", "Kratzer: Fensterlänge", "µm", 2, 50, 0.5, 1, "gleitende Dichte entlang der Linie",
     "Scratches: window length", "sliding density along the line"),
]


def param_text(entry):
    """(Beschriftung, Einheit, Hinweis) eines PARAM_INFO-Eintrags in der eingestellten Sprache"""
    key, lab, unit, lo, hi, step, nd, tip, lab_en, tip_en = entry
    unit = L(unit, "points") if unit == "Punkte" else unit
    return L(lab, lab_en), unit, L(tip, tip_en)


# (Kennung, Titel, Dateiname, Modul, Titel engl., Dateiname engl.)
FIGURES = [
    ("00", "Kratzererkennung", "00_Kratzer_Erkennung", "scratch", "Scratch detection", "00_Scratch_detection"),
    ("01", "IPF-Karte ‖ CD", "01_IPF-Karte_CD", "ipf", "IPF map ‖ CD", "01_IPF_map_CD"),
    ("02", "IPF-Karte ‖ RD", "02_IPF-Karte_RD", "ipf", "IPF map ‖ RD", "02_IPF_map_RD"),
    ("03", "IPF-Karte ‖ TD", "03_IPF-Karte_TD", "ipf", "IPF map ‖ TD", "03_IPF_map_TD"),
    ("04", "Inverse Polfiguren", "04_Inverse_Polfiguren", "texture", "Inverse pole figures", "04_Inverse_pole_figures"),
    ("05", "Polfiguren", "05_Polfiguren", "texture", "Pole figures", "05_Pole_figures"),
    ("06", "KAM-Karte", "06_KAM-Karte", "kam", "KAM map", "06_KAM_map"),
    ("07", "Korngrenz-Misorientierung", "07_KG-Misorientierung_Histogramm", "gb",
     "Grain boundary misorientation", "07_GB_misorientation_histogram"),
    ("08", "Phasenverteilung", "08_Phasenverteilung", None, "Phase distribution", "08_Phase_distribution"),
    ("09", "Rekristallisation – Karte", "09_Rekristallisation_Karte", "rx", "Recrystallization – map", "09_Recrystallization_map"),
    ("10", "Rekristallisation – Statistik", "10_Rekristallisation_Statistik", "rx",
     "Recrystallization – statistics", "10_Recrystallization_statistics"),
    ("11", "Korngröße – Karte", "11_Korngroesse_Karte", "grainsize", "Grain size – map", "11_Grain_size_map"),
    ("12", "Subkorngröße – Karte", "12_Subkorngroesse_Karte", "grainsize", "Subgrain size – map", "12_Subgrain_size_map"),
    ("13", "Kornaspektverhältnis – Karte", "13_Kornaspektverhaeltnis_Karte", "grainsize",
     "Grain aspect ratio – map", "13_Grain_aspect_ratio_map"),
    ("14", "Korngröße und Kornform – Statistik", "14_Korngroesse_Statistik", "grainsize",
     "Grain size and shape – statistics", "14_Grain_size_statistics"),
]
FORMATS = {"png": "PNG", "pdf": "PDF", "svg": "SVG", "tiff": "TIFF"}

# (Kennung, Titel, Dateiname, Titel engl., Dateiname engl.) – Vergleich mehrerer gespeicherter Datensätze
COMPARE_FIGURES = [
    ("V1", "Kennwerte-Vergleich", "V1_Vergleich_Kennwerte", "Comparison of key figures", "V1_Comparison_key_figures"),
    ("V2", "Misorientierungsverteilung", "V2_Vergleich_Misorientierung", "Misorientation distribution",
     "V2_Comparison_misorientation"),
]


def figures():
    """[(Kennung, Titel, Dateiname, Modul)] in der eingestellten Sprache"""
    return [(f[0], L(f[1], f[4]), L(f[2], f[5]), f[3]) for f in FIGURES]


def compare_figures():
    """[(Kennung, Titel, Dateiname)] in der eingestellten Sprache"""
    return [(f[0], L(f[1], f[3]), L(f[2], f[4])) for f in COMPARE_FIGURES]


def fig_title(fid, fallback=None):
    for f in FIGURES:
        if f[0] == fid:
            return L(f[1], f[4])
    return fallback or fid


def fig_fname(fid):
    for f in FIGURES:
        if f[0] == fid:
            return L(f[2], f[5])
    return fid


def kennwerte_name():
    return L("Kennwerte_und_Parameter.txt", "Results_and_parameters.txt")


def table_names():
    """Dateinamen der Kornlisten (Körner, Subkörner)"""
    return L("Kornliste_Koerner.csv", "Grain_list_grains.csv"), L("Kornliste_Subkoerner.csv", "Grain_list_subgrains.csv")
