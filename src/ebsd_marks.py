# -*- coding: utf-8 -*-
"""Manuelle Kratzermarkierungen (nur NumPy, damit die Übersicht sie ohne MTEX anzeigen kann).

Eine Markierung ist eine gerade Linie mit Breite in Kartenkoordinaten (µm, x = RD nach rechts, y = CD nach unten):
    {"modus": "add" | "remove", "x0": …, "y0": …, "x1": …, "y1": …, "breite": …}
"add"     markiert die überstrichenen Punkte als Kratzer,
"remove"  hebt dort jede Kratzermarkierung auf (automatisch erkannte und manuell gesetzte).
Die Linien werden in der gezeichneten Reihenfolge angewendet; eine spätere Linie überschreibt eine frühere.
Kratzer gesamt = (automatisch erkannt UND NICHT aufgehoben) ODER manuell markiert.
"""
from __future__ import annotations

import numpy as np

MODES = ("add", "remove")


def clean_strokes(strokes):
    """Liste aus JSON/QSettings prüfen; ungültige Einträge werden verworfen"""
    out = []
    for s in strokes or []:
        try:
            d = {"modus": str(s.get("modus", "add")), "x0": float(s["x0"]), "y0": float(s["y0"]),
                 "x1": float(s["x1"]), "y1": float(s["y1"]), "breite": float(s.get("breite", 1.0))}
        except (AttributeError, KeyError, TypeError, ValueError):
            continue
        if d["modus"] in MODES and d["breite"] > 0 and all(np.isfinite(v) for k, v in d.items() if k != "modus"):
            out.append(d)
    return out


def segment_mask(shape, x0, y0, x1, y1, half_width):
    """Punkte (Pixelmitten, Spalte x / Zeile y in Pixeln) im Abstand ≤ half_width zur Strecke.
    Liefert (Zeilen-Slice, Spalten-Slice, Maske) für den betroffenen Kartenausschnitt oder None."""
    ny, nx = shape
    hw = max(float(half_width), 0.5)
    pad = hw + 1.0
    c0 = max(0, int(np.floor(min(x0, x1) - pad)))
    c1 = min(nx, int(np.ceil(max(x0, x1) + pad)) + 1)
    r0 = max(0, int(np.floor(min(y0, y1) - pad)))
    r1 = min(ny, int(np.ceil(max(y0, y1) + pad)) + 1)
    if c0 >= c1 or r0 >= r1:
        return None
    yy, xx = np.mgrid[r0:r1, c0:c1].astype(float)
    dx, dy = x1 - x0, y1 - y0
    l2 = dx * dx + dy * dy
    t = np.clip(((xx - x0) * dx + (yy - y0) * dy) / l2, 0, 1) if l2 > 0 else 0.0
    d = np.hypot(xx - (x0 + t * dx), yy - (y0 + t * dy))
    return slice(r0, r1), slice(c0, c1), d <= hw


def stroke_px(stroke, step):
    """Linie in µm -> Pixelkoordinaten (Pixel c reicht von c·Schrittweite bis (c+1)·Schrittweite)"""
    return (stroke["x0"] / step - 0.5, stroke["y0"] / step - 0.5, stroke["x1"] / step - 0.5, stroke["y1"] / step - 0.5,
            0.5 * stroke["breite"] / step)


def masks(strokes, shape, step):
    """(add, remove): Punkte, die manuell als Kratzer markiert bzw. deren Markierung aufgehoben wurde"""
    add = np.zeros(shape, bool)
    rem = np.zeros(shape, bool)
    for s in clean_strokes(strokes):
        out = segment_mask(shape, *stroke_px(s, step))
        if out is None:
            continue
        rs, cs, m = out
        if s["modus"] == "add":
            add[rs, cs] |= m
            rem[rs, cs] &= ~m
        else:
            rem[rs, cs] |= m
            add[rs, cs] &= ~m
    return add, rem


def combine(auto, add, rem):
    """Kratzermaske aus automatischer Erkennung und manuellen Markierungen"""
    return (auto & ~rem) | add
