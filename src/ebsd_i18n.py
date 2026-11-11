# -*- coding: utf-8 -*-
"""Sprachumschaltung Deutsch / Englisch der EBSD-Auswertung (ohne schwere Importe).

    from ebsd_i18n import L, dez
    L("Korngröße", "Grain size")      -> Text in der eingestellten Sprache
    dez(3.14159, 2)                   -> "3,14" (Deutsch) bzw. "3.14" (Englisch)

Die Sprache gilt für Oberfläche, Abbildungen, Kennwerte und Kornlisten. Abbildungen und Kennwerte werden in der
Sprache erzeugt, die beim Auswerten bzw. Erstellen der Abbildung eingestellt ist.
"""
from __future__ import annotations

import math

LANGUAGES = {"de": "Deutsch", "en": "English"}
_LANG = "de"


def set_language(code):
    global _LANG
    _LANG = "en" if str(code or "").lower().startswith("en") else "de"
    return _LANG


def get_language():
    return _LANG


def is_en():
    return _LANG == "en"


def L(de, en):
    """Text in der eingestellten Sprache"""
    return en if _LANG == "en" else de


def _finite(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def dez(x, nd=1):
    """Zahl mit nd Nachkommastellen, Dezimalkomma (Deutsch) bzw. Dezimalpunkt (Englisch); ungültig -> „–“"""
    v = _finite(x)
    if v is None:
        return "–"
    s = f"{v:.{nd}f}"
    return s if _LANG == "en" else s.replace(".", ",")


def num_g(x):
    """kompakte Zahl ohne überflüssige Nullen (Achsbeschriftungen)"""
    v = _finite(x)
    if v is None:
        return ""
    s = f"{v:g}"
    return s if _LANG == "en" else s.replace(".", ",")


def thou(n):
    """ganze Zahl mit Tausendertrennzeichen: 12.345 (Deutsch) bzw. 12,345 (Englisch)"""
    try:
        s = f"{int(n):,}"
    except (TypeError, ValueError):
        return "–"
    return s if _LANG == "en" else s.replace(",", ".")


def csv_sep():
    """Trennzeichen der CSV-Dateien: Semikolon (deutsches Excel) bzw. Komma (englisches Excel)"""
    return "," if _LANG == "en" else ";"


def csv_num(x, nd=3):
    v = _finite(x)
    if v is None:
        return ""
    s = f"{v:.{nd}f}"
    return s if _LANG == "en" else s.replace(".", ",")


def csv_join(fields):
    """eine CSV-Zeile; Felder mit Trennzeichen oder Anführungszeichen werden in Anführungszeichen gesetzt"""
    sep = csv_sep()
    out = []
    for f in fields:
        f = "" if f is None else str(f)
        if sep in f or '"' in f or "\n" in f:
            f = '"' + f.replace('"', '""') + '"'
        out.append(f)
    return sep.join(out)


def csv_note():
    return L("Semikolon, Dezimalkomma", "comma-separated, decimal point")
