# -*- coding: utf-8 -*-
"""
Auswertungskern der EBSD-Auswertung (MTEX für Python) – ohne Oberfläche.

    proj = EBSDProject.load("datei.ctf")
    roles = proj.default_roles()                       # {Phasenindex: "Matrix" | "Partikel" | "Fehlindizierung"}
    res = analyze(proj, Settings(), roles, crop_um=None, progress=print_progress)
    fig = res.build("01")                              # matplotlib-Figure, z. B. IPF-Karte ‖ CD
    res.export(["01", "06"], "ordner", ["png", "pdf"], dpi=300)

Die Berechnungen entsprechen ebsd_auswertung_mtex.py; zusätzlich sind die Phasenrollen frei wählbar
(eine Matrixphase, beliebig viele Partikelphasen, Fehlindizierungen) und die Module einzeln schaltbar.
"""
from __future__ import annotations

import io
import os
import time
import warnings
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.collections import Collection
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from matplotlib.image import AxesImage
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from ebsd_settings import (APP_NAME, APP_VERSION, Abgebrochen, Settings, app_title, fig_fname, figures,
                           kennwerte_name, role_name, table_names)
from ebsd_i18n import L as tx, csv_join, csv_num, dez, num_g, thou
import ebsd_marks
from mtex import (EBSD, Miller, calcAngleDistribution, calcDensity, calcGrains, degree, fibre, ipfTSLKey,
                  orientation, plot, plotIPF, plotPF, rotation, specimenFrame, xvector, yvector, zvector)


plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "savefig.facecolor": "white",
                     "font.family": "DejaVu Sans", "mathtext.default": "regular"})
comma = FuncFormatter(lambda v, p: num_g(v))           # Dezimaltrennzeichen der eingestellten Sprache
S8 = np.ones((3, 3), bool)


def chunks(n, size=500_000):
    for i in range(0, n, size):
        yield slice(i, min(i + size, n))


GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "theta": "θ", "eta": "η", "kappa": "κ",
         "sigma": "σ", "tau": "τ", "omega": "ω", "chi": "χ", "epsilon": "ε", "mu": "μ", "pi": "π"}


def pretty(name):
    """Phasenname für Abbildungen: 'alpha_c_AlFeMnSi' -> 'α-AlFeMnSi', 'Al5FeSi' -> 'Al$_5$FeSi'"""
    parts = name.split("_")
    if parts and parts[0].lower() in GREEK:
        rest = [p for p in parts[1:] if len(p) > 1] or parts[1:]
        s = GREEK[parts[0].lower()] + "-" + "".join(rest)
    else:
        s = name.replace("_", " ")
    out, i = "", 0
    while i < len(s):                                   # Ziffern hinter Elementsymbolen tiefstellen
        if s[i].isdigit() and i > 0 and s[i - 1].isalpha():
            j = i
            while j < len(s) and s[j].isdigit():
                j += 1
            out += f"$_{{{s[i:j]}}}$"
            i = j
        else:
            out += s[i]
            i += 1
    return out


def plain(name):
    """Phasenname für Textdateien"""
    return pretty(name).replace("$_{", "").replace("}$", "")


# ======================================================================== Projekt (geladene Datei)
class EBSDProject:
    """geladene ctf-Datei mit Phasenübersicht"""

    def __init__(self, path, ebsd, step):
        self.path = Path(path)
        self.ebsd = ebsd
        self.step = step
        self.ny, self.nx = ebsd.shape
        self.names = list(ebsd.mineralList)
        self.ph = np.asarray(ebsd.phaseId).reshape(self.ny, self.nx).astype(np.int16)
        self.bc = np.asarray(ebsd.prop["bc"], float).reshape(self.ny, self.nx) if "bc" in ebsd.prop else np.full((self.ny, self.nx), 128.0)
        self.mad = np.asarray(ebsd.prop["mad"], float).reshape(self.ny, self.nx) if "mad" in ebsd.prop else np.full((self.ny, self.nx), np.nan)
        self.phases = []
        for i, nm in enumerate(self.names):
            m = self.ph == i
            sym = ""
            if i:
                try:
                    sym = str(ebsd.CSList[i].pointGroup)
                except Exception:
                    sym = ""
            self.phases.append(dict(index=i, name=nm if i else tx("Nulllösung", "Zero solution"), symmetry=sym, frac=m.mean() * 100,
                                    mad=float(np.nanmedian(self.mad[m])) if m.any() else np.nan,
                                    bc=float(np.median(self.bc[m])) if m.any() else np.nan))

    @property
    def width_um(self):
        return self.nx * self.step

    @property
    def height_um(self):
        return self.ny * self.step

    @classmethod
    def load(cls, path, log=print):
        path = Path(path)
        t0 = time.time()
        frame = specimenFrame("Stauchprobe", axesNames=("RD", "CD", "TD"), how2plot="y↓→x")
        frame.makeDefault()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ebsd = EBSD.load(str(path), frame=frame)
        if getattr(ebsd, "ndim", 2) != 2:
            raise ValueError(tx("Die Karte liegt nicht als Quadratraster vor (z. B. Hexagonalraster) – nicht unterstützt.",
                                "The map is not a square grid (e.g. hexagonal grid) – not supported."))
        hdr = ebsd.meta.get("header", {}) if hasattr(ebsd, "meta") else {}
        step = float(hdr.get("XStep", np.nan))
        if not np.isfinite(step) or step <= 0:
            step = float(np.median(np.abs(np.diff(np.asarray(ebsd.x)[0]))))
        p = cls(path, ebsd, step)
        log(tx(f"geladen: {path.name} – {p.nx} × {p.ny} Punkte, Schrittweite {dez(step, 4)} µm ({time.time() - t0:.1f} s)",
               f"loaded: {path.name} – {p.nx} × {p.ny} points, step size {dez(step, 4)} µm ({time.time() - t0:.1f} s)"))
        return p

    def default_roles(self):
        """häufigste Phase = Matrix; übrige Phasen mit MAD-Median ≤ 0,8° = Partikel, sonst Fehlindizierung"""
        idx = [p for p in self.phases if p["index"] > 0]
        if not idx:
            return {}
        matrix = max(idx, key=lambda p: p["frac"])["index"]
        roles = {}
        for p in idx:
            if p["index"] == matrix:
                roles[p["index"]] = "Matrix"
            elif np.isfinite(p["mad"]) and p["mad"] <= 0.8:
                roles[p["index"]] = "Partikel"
            else:
                roles[p["index"]] = "Fehlindizierung"
        return roles


# ======================================================================== Kratzererkennung
def detect_scratches(ph, bc, step, misidx, min_len_um, win_um, rel_bc, dtheta=0.5, max_checks=500, cancel=None):
    """Kratzer = lange, gerade, schmale Linien aus Nulllösungen, Fehlindizierungen und lokal
    niedrigem Bandkontrast (Hough-Transformation + Dichteprofil entlang der Linie)."""
    NY, NX = ph.shape
    bg = ndi.percentile_filter(np.where(ph > 0, bc, 0).astype(float), 75, size=9)
    rel = np.where(bg > 0, bc / np.maximum(bg, 1), 0.0)
    cand = (ph == 0) | misidx | (rel < rel_bc)
    blob = ndi.binary_opening(cand, structure=np.ones((5, 5), bool))
    thin = cand & ~ndi.binary_dilation(blob, structure=S8)
    w = max(11, int(round(win_um / step)) | 1)
    Lcore = max(3 * w, int(round(min(min_len_um / step, 0.6 * min(NX, NY)))))
    th = np.deg2rad(np.arange(0, 180, dtheta))
    cs, sn = np.cos(th), np.sin(th)
    diag = int(np.ceil(np.hypot(NY, NX)))
    nr = 2 * diag + 2

    def votes(ys, xs, sign=1.0, acc=None):
        acc = np.zeros((len(th), nr)) if acc is None else acc
        for i in range(len(th)):
            r = xs * cs[i] + ys * sn[i] + diag
            r0 = np.floor(r).astype(np.int64)
            f = r - r0
            acc[i] += sign * (np.bincount(r0, 1 - f, minlength=nr) + np.bincount(r0 + 1, f, minlength=nr))
        return acc

    ys, xs = np.nonzero(thin)
    ys = ys.astype(float)
    xs = xs.astype(float)
    if len(xs) == 0 or min(NX, NY) < 3 * w:
        return np.zeros((NY, NX), bool), [], dict(cand_frac=0, p0=0, thr=0, window_px=w, min_len_px=Lcore)
    acc = votes(ys, xs)
    exc = acc - ndi.median_filter(acc, size=(1, 41), mode="nearest")
    tt_all = np.arange(-diag - 2, diag + 3)

    def hits_along(n, rho, band):
        d = np.array([-n[1], n[0]])
        sel = np.abs(xs * n[0] + ys * n[1] - rho) <= band
        px = rho * n[0] + tt_all * d[0]
        py = rho * n[1] + tt_all * d[1]
        inside = (px >= 0) & (px <= NX - 1) & (py >= 0) & (py <= NY - 1)
        if inside.sum() < 3 * w:
            return None
        tmin, tmax = tt_all[inside].min(), tt_all[inside].max()
        hit = np.zeros(tmax - tmin + 1, bool)
        ti = np.round(xs[sel] * d[0] + ys[sel] * d[1]).astype(np.int64) - tmin
        hit[ti[(ti >= 0) & (ti < len(hit))]] = True
        return d, tmin, hit

    rng = np.random.default_rng(0)
    p0s = []
    for _ in range(200):
        a = rng.uniform(0, np.pi)
        n = np.array([np.cos(a), np.sin(a)])
        lo = min(0, n[0] * (NX - 1)) + min(0, n[1] * (NY - 1))
        hi = max(0, n[0] * (NX - 1)) + max(0, n[1] * (NY - 1))
        out = hits_along(n, rng.uniform(lo + 0.2 * (hi - lo), hi - 0.2 * (hi - lo)), 1.0)
        if out is not None:
            p0s.append(out[2].mean())
    p0 = float(np.median(p0s)) if p0s else 0.3
    thr = float(np.clip(p0 + 0.3, 0.6, 0.85))
    thr_sub = float(np.clip(p0 + 0.2, 0.45, thr))
    info = dict(cand_frac=thin.mean() * 100, p0=p0, thr=thr, window_px=w, min_len_px=Lcore)

    supp = np.zeros_like(exc, bool)
    mask = np.zeros((NY, NX), bool)
    lines = []
    YY, XX = np.mgrid[0:NY, 0:NX]
    for _ in range(max_checks):
        if cancel is not None:
            cancel()
        e = np.where(supp, -np.inf, exc)
        i, j = np.unravel_index(np.argmax(e), e.shape)
        if e[i, j] < 0.35 * Lcore:
            break
        supp[max(0, i - 2):i + 3, max(0, j - 4):j + 5] = True
        n = np.array([cs[i], sn[i]])
        rho = j - diag
        sel = np.abs(xs * n[0] + ys * n[1] - rho) <= 3
        if sel.sum() < 0.5 * Lcore:
            continue
        P = np.stack([xs[sel], ys[sel]], 1)
        c0 = P.mean(0)
        dvec = np.linalg.svd(P - c0, full_matrices=False)[2][0]
        n = np.array([-dvec[1], dvec[0]])
        rho = float(c0 @ n)
        out = hits_along(n, rho, 1.0)
        if out is None:
            continue
        d, tmin, hit = out
        dens = np.convolve(hit.astype(float), np.ones(w) / w, mode="same")
        lab, _ = ndi.label(dens >= thr)
        core = [(s[0].start, s[0].stop) for s in ndi.find_objects(lab)]
        if not any(b - a >= Lcore for a, b in core):
            continue
        lab, _ = ndi.label(dens >= thr_sub)
        segs = [(s[0].start + tmin - w // 4, s[0].stop + tmin + w // 4) for s in ndi.find_objects(lab)
                if s[0].stop - s[0].start >= w]
        near = np.abs(xs * n[0] + ys * n[1] - rho) <= 8
        tt = xs[near] * d[0] + ys[near] * d[1]
        ins = np.zeros(len(tt), bool)
        for a, b in segs:
            ins |= (tt >= a) & (tt < b)
        dd = np.round(xs[near] * n[0] + ys[near] * n[1] - rho).astype(int)[ins]
        seglen = sum(b - a for a, b in segs)
        prof = np.array([(dd == o).sum() / max(seglen, 1) for o in range(-8, 9)])
        bgp = np.median(prof[np.r_[0:3, 14:17]])
        on = prof >= bgp + 0.4 * (prof[6:11].max() - bgp)
        lo = hi = 0
        while lo - 1 >= -8 and on[lo - 1 + 8]:
            lo -= 1
        while hi + 1 <= 8 and on[hi + 1 + 8]:
            hi += 1
        hw_lo, hw_hi = min(abs(lo) + 1.0, 6.0), min(hi + 1.0, 6.0)
        dist = XX * n[0] + YY * n[1] - rho
        tpos = XX * d[0] + YY * d[1]
        band = (dist >= -hw_lo - 0.5) & (dist <= hw_hi + 0.5)
        m = np.zeros((NY, NX), bool)
        for a, b in segs:
            m |= band & (tpos >= a) & (tpos < b)
        mask |= m
        ends = [((rho * n[0] + a * d[0], rho * n[1] + a * d[1]), (rho * n[0] + b * d[0], rho * n[1] + b * d[1])) for a, b in segs]
        lines.append(dict(angle=float(np.degrees(np.arctan2(-d[1], d[0])) % 180), length_um=seglen * step,
                          width_um=(hw_lo + hw_hi + 1) * step, n_segs=len(segs), ends=ends))
        rem = m[ys.astype(int), xs.astype(int)]
        if rem.any():
            votes(ys[rem], xs[rem], -1.0, acc)
            ys, xs = ys[~rem], xs[~rem]
            exc = acc - ndi.median_filter(acc, size=(1, 41), mode="nearest")
    return mask, lines, info


# ======================================================================== Abbildungsbausteine
SIDE_W = 3.3
MAP_W_MAX = 9.5
MAP_H_MAX = 7.0
LEFT_MIN = 5.6
M_L, M_R, M_B = 0.3, 0.3, 0.3
GAP = 0.4
TITLE_H = 1.0
BLOCK_GAP = 0.22


def nice_length(target):
    e = 10 ** np.floor(np.log10(target))
    for m in (5, 2, 1):
        if m * e <= target:
            return m * e
    return e


class Block:
    def __init__(self, h, draw):
        self.h, self.draw = h, draw


def r_h(fig, r):
    return r[3] * fig.get_size_inches()[1]


def fit_text(fig, txt, max_w_in):
    r = fig.canvas.get_renderer()
    bb = txt.get_window_extent(renderer=r)
    w_in = bb.width / fig.dpi
    if w_in > max_w_in:
        txt.set_fontsize(txt.get_fontsize() * max_w_in / w_in * 0.98)


def text_block(lines, fs=10.5, head=None, mono=False):
    n_head = head.count("\n") + 1 if head else 0
    h = 0.06 + (len(lines) + 1.25 * n_head + (0.3 if head else 0)) * fs * 1.32 / 72

    def draw(fig, r):
        ax = fig.add_axes(r)
        ax.axis("off")
        w_in = r[2] * fig.get_size_inches()[0]
        y = 1.0
        if head:
            t = ax.text(0, y, head, fontsize=fs + 1.5, weight="bold", va="top", linespacing=1.25, transform=ax.transAxes)
            fit_text(fig, t, w_in)
            y -= (1.25 * n_head + 0.3) * fs * 1.32 / 72 / r_h(fig, r)
        t = ax.text(0, y, "\n".join(lines), fontsize=fs, va="top", linespacing=1.32, transform=ax.transAxes,
                    family="DejaVu Sans Mono" if mono else None)
        fit_text(fig, t, w_in)
    return Block(h, draw)


def legend_block(handles, fs=10.5, head=None):
    h = 0.08 + len(handles) * 0.255 + (0.3 if head else 0)

    def draw(fig, r):
        ax = fig.add_axes(r)
        ax.axis("off")
        if head:
            ax.text(0, 1, head, fontsize=fs + 1, weight="bold", va="top", transform=ax.transAxes)
        ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=fs, borderaxespad=0,
                  bbox_to_anchor=(0, 1 - (0.3 / r_h(fig, r) if head else 0)), handlelength=1.8, labelspacing=0.45)
    return Block(h, draw)


def glyph_block(h=1.45):
    def draw(fig, r):
        ax = fig.add_axes(r)
        W, H = r[2] * fig.get_size_inches()[0], r[3] * fig.get_size_inches()[1]
        ax.set_xlim(0, W)
        ax.set_ylim(0, H)
        ax.axis("off")
        ox, oy, L = 0.35, H - 0.3, min(1.0, H - 0.55)
        kw = dict(arrowprops=dict(arrowstyle="-|>", lw=1.8, color="k", mutation_scale=16))
        ax.annotate("", xy=(ox + L, oy), xytext=(ox + 0.1, oy), **kw)
        ax.annotate("", xy=(ox, oy - L), xytext=(ox, oy - 0.1), **kw)
        ax.add_patch(Circle((ox, oy), 0.09, fc="white", ec="k", lw=1.4))
        ax.plot([ox - 0.055, ox + 0.055], [oy - 0.055, oy + 0.055], "k-", lw=1.2)
        ax.plot([ox - 0.055, ox + 0.055], [oy + 0.055, oy - 0.055], "k-", lw=1.2)
        ax.text(ox + L + 0.08, oy, "RD", va="center", fontsize=11.5)
        ax.text(ox, oy - L - 0.08, tx("CD (Stauchrichtung)", "CD (compression direction)"), ha="left", va="top", fontsize=11.5)
        ax.text(ox + 0.16, oy + 0.1, "TD", ha="left", va="bottom", fontsize=11.5)
    return Block(h, draw)


def ipfkey_block(key, label, matrix_name, h=2.35):
    def draw(fig, r):
        H = r_h(fig, r)
        head_h = 0.5
        axt = fig.add_axes([r[0], r[1] + r[3] * (1 - head_h / H), r[2], r[3] * head_h / H])
        axt.axis("off")
        axt.text(0, 1, tx(f"IPF-Farbschlüssel ‖ {label}", f"IPF colour key ‖ {label}"), fontsize=12, weight="bold", va="top", transform=axt.transAxes)
        axt.text(0, 0.35, f"{matrix_name}", fontsize=10.5, va="top", transform=axt.transAxes)
        axk = fig.add_axes([r[0] + 0.05 * r[2], r[1], 0.85 * r[2], r[3] * (1 - head_h / H)])
        plot(key, parent=axk)
        axk._ebsd_map = True                                 # Farbschlüssel in Vektorexporten als Raster
    return Block(h, draw)


def cbar_block(cmap, norm, label, ticks=None, extend="max", h=0.95):
    def draw(fig, r):
        H = r_h(fig, r)
        cax = fig.add_axes([r[0] + 0.05 * r[2], r[1] + r[3] * (0.55 / H), 0.85 * r[2], r[3] * (0.22 / H)])
        cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax, orientation="horizontal",
                          extend=extend, ticks=ticks)
        cb.ax.xaxis.set_major_formatter(comma)
        cb.ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, p: ""))
        cb.ax.set_title(label, fontsize=11, loc="left")
    return Block(h, draw)


def axes_block(h, fn, pad=(0.62, 0.55, 0.12, 0.3)):
    def draw(fig, r):
        FW, FH = fig.get_size_inches()
        x, y, w_, h_ = r[0] * FW, r[1] * FH, r[2] * FW, r[3] * FH
        ax = fig.add_axes([(x + pad[0]) / FW, (y + pad[1]) / FH, (w_ - pad[0] - pad[2]) / FW, (h_ - pad[1] - pad[3]) / FH])
        fn(ax)
    return Block(h, draw)


class MapLayout:
    """Karten mit Seitenleiste in festen Zoll-Maßen (unabhängig von Kartengröße und Seitenverhältnis)"""

    def __init__(self, ebsd, step):
        self.ebsd = ebsd
        self.NY, self.NX = ebsd.shape
        self.W_UM, self.H_UM = self.NX * step, self.NY * step
        s = min(MAP_W_MAX / self.W_UM, MAP_H_MAX / self.H_UM)
        self.mw, self.mh = self.W_UM * s, self.H_UM * s

    def lw_for(self, base):
        px_pt = self.mw * 72 / self.NX
        return float(np.clip(base * px_pt, 0.25, 1.3))

    def add_scalebar(self, ax, mw_in):
        L = nice_length(0.2 * self.W_UM)
        thick = 3.5 / 72 * self.W_UM / mw_in
        sb = AnchoredSizeBar(ax.transData, L, f"{num_g(L)} µm", "lower left", pad=0.35, borderpad=0.6, sep=4, frameon=True,
                             size_vertical=thick, fontproperties=fm.FontProperties(size=11))
        sb.patch.set_alpha(0.85)
        sb.patch.set_edgecolor("none")
        ax.add_artist(sb)

    def draw_map(self, ax, rgb, boundaries, mw_in):
        plot(self.ebsd, np.asarray(rgb, float).reshape(-1, 3), parent=ax, micronbar=False, refFrame=False)
        for gb_, col, lw in boundaries:
            if gb_ is not None and len(gb_):
                plot(gb_, parent=ax, lineColor=col, lineWidth=lw, micronbar=False, refFrame=False)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_linewidth(0.8)
        ax._ebsd_map = True                                  # Kennzeichnung für Raster in Vektorexporten
        self.add_scalebar(ax, mw_in)

    def figure(self, rgb, title, subtitle, blocks, boundaries=()):
        left_w = max(self.mw, LEFT_MIN)
        side_h = sum(b.h for b in blocks) + BLOCK_GAP * (len(blocks) - 1)
        body_h = max(self.mh, side_h)
        FW = M_L + left_w + GAP + SIDE_W + M_R
        FH = TITLE_H + body_h + M_B
        fig = plt.figure(figsize=(FW, FH), dpi=100)

        def rect(x, ytop, w_, h_):
            return [x / FW, (FH - ytop - h_) / FH, w_ / FW, h_ / FH]

        ax = fig.add_axes(rect(M_L + (left_w - self.mw) / 2, TITLE_H, self.mw, self.mh))
        self.draw_map(ax, rgb, boundaries, self.mw)
        t1 = fig.text(M_L / FW, 1 - 0.36 / FH, title, fontsize=16, weight="bold", va="center")
        t2 = fig.text(M_L / FW, 1 - 0.72 / FH, subtitle, fontsize=11, va="center", color="0.25")
        for t in (t1, t2):
            fit_text(fig, t, FW - 2 * M_L)
        y = TITLE_H
        x = M_L + left_w + GAP
        for b in blocks:
            b.draw(fig, rect(x, y, SIDE_W, b.h))
            y += b.h + BLOCK_GAP
        return fig


def mrd_levels(vmax):
    top = max(np.ceil(vmax * 2) / 2, 1.5)
    step = 0.25 if top <= 4 else (0.5 if top <= 10 else 1.0)
    return np.arange(0, top + step / 2, step)


def level_cbar(fig, cax, levels, cmap_name="turbo"):
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(levels[0], levels[-1])
    cols = ListedColormap([cmap(norm(l)) for l in levels[:-1]])
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cols, norm=BoundaryNorm(levels, cols.N)), cax=cax,
                      ticks=levels[::max(1, len(levels) // 8)])
    cb.set_label("MRD")
    cb.ax.yaxis.set_major_formatter(comma)


# ======================================================================== Abbildungsliste
class Result:
    """Ergebnis einer Auswertung: Abbildungen auf Anforderung, Kennwerte, Export"""

    def __init__(self, builders, kennwerte, settings, project, crop_um, summary, tables=None, metrics=None,
                 scratch_auto=None, crop_px=(0, 0), strokes=None):
        self.builders = builders
        self.metrics = metrics or {}                     # Kennwerte für den Datensatzvergleich (JSON-fähig)
        self.tables = tables or {}                       # {Dateiname: CSV-Text}, z. B. Kornlisten
        self.kennwerte = kennwerte
        self.settings = settings
        self.project = project
        self.crop_um = crop_um
        self.summary = summary
        self.scratch_auto = scratch_auto                 # automatisch erkannte Kratzer (Ausschnitt), für die Übersicht
        self.crop_px = crop_px                           # (Zeile, Spalte) der linken oberen Ecke des Ausschnitts
        self.strokes = strokes or []                     # verwendete manuelle Kratzermarkierungen

    @property
    def available(self):
        return [(fid, title, fname) for fid, title, fname, _ in figures() if fid in self.builders]

    def build(self, fid):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return self.builders[fid]()

    def preview_png(self, fid, dpi=110):
        fig = self.build(fid)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi)
        plt.close(fig)
        return buf.getvalue()

    def default_folder(self):
        p = self.project.path
        name = p.stem.replace(" ", "_") + tx("_Auswertung", "_analysis")
        if self.crop_um is not None:
            name += tx("_Ausschnitt_", "_crop_") + "_".join(f"{v:g}" for v in self.crop_um)
        return p.with_name(name)

    def save_kennwerte(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / kennwerte_name()
        path.write_text(self.kennwerte, encoding="utf-8")
        return path

    def save_tables(self, folder):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        out = []
        for name, text in self.tables.items():
            path = folder / name
            path.write_text(text, encoding="utf-8-sig")       # mit BOM, damit Excel die Umlaute erkennt
            out.append(path)
        return out

    def export(self, fids, folder, formats, dpi=300, rasterize_vector=True, progress=None, cancel=None):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        written = []
        n = max(1, len(fids))
        for k, fid in enumerate(fids):
            if cancel is not None:
                cancel()
            name = fig_fname(fid)
            if progress:
                progress(int(100 * k / n), f"Export {name}")
            fig = self.build(fid)
            for fmt in formats:
                path = folder / f"{name}.{'tif' if fmt == 'tiff' else fmt}"
                vector = fmt in ("pdf", "svg")
                for ax in fig.axes:
                    if getattr(ax, "_ebsd_map", False):
                        for ch in ax.get_children():
                            if isinstance(ch, (AxesImage, Collection)):
                                ch.set_rasterized(bool(vector and rasterize_vector))
                kw = {}
                if fmt == "tiff":
                    kw["pil_kwargs"] = {"compression": "tiff_lzw"}
                fig.savefig(path, format=fmt, dpi=dpi, **kw)
                written.append(path)
            plt.close(fig)
        written.append(self.save_kennwerte(folder))
        written += self.save_tables(folder)
        if progress:
            progress(100, tx("Export abgeschlossen", "Export finished"))
        return written


def crop_pixels(project, crop_um):
    """Ausschnitt in µm -> (Zeile0, Zeile1, Spalte0, Spalte1) der Gesamtkarte"""
    if crop_um is None:
        return 0, project.ny, 0, project.nx
    x0, y0, x1, y1 = crop_um
    st = project.step
    c0, r0 = max(0, int(round(x0 / st))), max(0, int(round(y0 / st)))
    c1, r1 = min(project.nx, int(round(x1 / st))), min(project.ny, int(round(y1 / st)))
    return r0, r1, c0, c1


def preview_scratches(project, S: Settings, roles: dict, crop_um=None, cancel=None):
    """nur die automatische Kratzererkennung (für die Anzeige in der Übersicht); Maske der Gesamtkarte"""
    r0, r1, c0, c1 = crop_pixels(project, crop_um)
    ph = project.ph[r0:r1, c0:c1]
    bc = project.bc[r0:r1, c0:c1]
    iMIS = [i for i, r in roles.items() if r == "Fehlindizierung"]
    m, lines, _ = detect_scratches(ph, bc, project.step, np.isin(ph, iMIS), S.scratch_min_len, S.scratch_window,
                                   S.scratch_rel_bc, cancel=cancel)
    full = np.zeros((project.ny, project.nx), bool)
    full[r0:r1, c0:c1] = m
    return full, len(lines)


# ======================================================================== Auswertung
def analyze(project: EBSDProject, S: Settings, roles: dict, crop_um=None, progress=None, cancel=None, log=None,
            manual_scratch=None):
    """führt die Auswertung aus und liefert ein Result. progress(prozent, text); cancel() wirft Abgebrochen.
    manual_scratch: Liste manueller Kratzermarkierungen (siehe ebsd_marks) in Koordinaten der Gesamtkarte."""
    T0 = time.time()

    def say(msg):
        if log:
            log(f"[{time.time() - T0:6.1f}s] {msg}")

    def prog(p, msg):
        if cancel is not None:
            cancel()
        if progress:
            progress(int(p), msg)
        say(msg)

    # ------------------------------------------------------------------ Rollen
    iAL = [i for i, r in roles.items() if r == "Matrix"]
    if len(iAL) != 1:
        raise ValueError(tx("Es muss genau eine Phase als Matrix gewählt sein.", "Exactly one phase must be set as matrix."))
    iAL = iAL[0]
    iPARTS = sorted(i for i, r in roles.items() if r == "Partikel")
    iMIS = sorted(i for i, r in roles.items() if r == "Fehlindizierung")
    names = project.names
    PH_MATRIX = names[iAL]
    MATRIX_LABEL = pretty(PH_MATRIX)
    SUBG = bool(S.subgrains)

    # ------------------------------------------------------------------ Ausschnitt
    prog(1, tx("Daten vorbereiten", "Preparing data"))
    STEP = project.step
    ebsd = project.ebsd
    r0, r1, c0, c1 = crop_pixels(project, crop_um)
    if crop_um is not None:
        if c1 - c0 < 20 or r1 - r0 < 20:
            raise ValueError(tx("Der Ausschnitt ist zu klein (mindestens 20 × 20 Messpunkte).",
                                "The selected area is too small (at least 20 × 20 points)."))
        ebsd = ebsd[r0:r1, c0:c1]
    NY, NX = ebsd.shape
    NPIX = NY * NX
    W_UM, H_UM = NX * STEP, NY * STEP
    csAl = ebsd.CSList[iAL]
    ph = np.asarray(ebsd.phaseId).reshape(NY, NX).astype(np.int16)
    bc = np.asarray(ebsd.prop["bc"], float).reshape(NY, NX) if "bc" in ebsd.prop else np.full((NY, NX), 128.0)
    mad = np.asarray(ebsd.prop["mad"], float).reshape(NY, NX) if "mad" in ebsd.prop else np.full((NY, NX), np.nan)
    Q = np.ascontiguousarray(ebsd.rotations.data, float).reshape(NY, NX, 4)
    raw_frac = {i: (ph == i).mean() * 100 for i in range(len(names))}
    phase_stats = {i: (np.nanmedian(mad[ph == i]) if (ph == i).any() else np.nan,
                       np.median(bc[ph == i]) if (ph == i).any() else np.nan) for i in range(1, len(names))}
    PX_AREA = STEP ** 2
    MG2SI_MIN_PX = max(4, int(round(S.mg2si_min_area / PX_AREA)))
    ALPHA_MIN_PX = max(1, int(round(S.particle_min_area / PX_AREA)))
    PSN_MIN_PX = max(1, int(np.ceil(np.pi / 4 * S.psn_min_ecd ** 2 / PX_AREA)))
    AL, CLUSTER = 1, 5
    PCODE = {ip: 10 + k for k, ip in enumerate(iPARTS)}          # Code der Partikelphasen in L
    CLUSTER_LABEL = (tx("Mg$_2$Si (Nulllösungs-Cluster)", "Mg$_2$Si (zero-solution cluster)") if S.mg2si
                     else tx("nicht indiziert (Cluster)", "non-indexed (cluster)"))
    CLUSTER_PLAIN = (tx("Mg2Si (Nulllösungs-Cluster)", "Mg2Si (zero-solution cluster)") if S.mg2si
                     else tx("nicht indizierte Cluster", "non-indexed clusters"))

    def mis_angle(qa, qb):
        out = np.empty(len(qa))
        for s in chunks(len(qa)):
            out[s] = orientation(rotation(qa[s]), csAl).angle(orientation(rotation(qb[s]), csAl)) / degree
        return out

    # ------------------------------------------------------------------ Kratzer (automatisch + manuell)
    misidx = np.isin(ph, iMIS)
    strokes = ebsd_marks.clean_strokes(manual_scratch)
    if S.scratch:
        prog(3, tx("Kratzer erkennen", "Detecting scratches"))
        scratch_auto, scratch_lines, scratch_info = detect_scratches(ph, bc, STEP, misidx, S.scratch_min_len,
                                                                     S.scratch_window, S.scratch_rel_bc, cancel=cancel)
    else:
        scratch_auto, scratch_lines, scratch_info = np.zeros((NY, NX), bool), [], {}
    if strokes:
        m_add, m_rem = ebsd_marks.masks(strokes, (project.ny, project.nx), STEP)
        m_add, m_rem = m_add[r0:r1, c0:c1], m_rem[r0:r1, c0:c1]
    else:
        m_add = m_rem = np.zeros((NY, NX), bool)
    scratch = ebsd_marks.combine(scratch_auto, m_add, m_rem)
    scr_removed = scratch_auto & m_rem                          # automatische Markierung manuell aufgehoben
    scr_added = scratch & ~scratch_auto                         # nur manuell markiert
    if m_rem.any() and scratch_lines:                           # weitgehend aufgehobene Linien nicht mehr zählen
        kept = []
        for ln in scratch_lines:
            tot = left = 0
            for (xa, ya), (xb, yb) in ln["ends"]:
                out = ebsd_marks.segment_mask((NY, NX), xa, ya, xb, yb, 0.5 * ln["width_um"] / STEP)
                if out is not None:
                    rs, cs, mm = out
                    tot += int((mm & scratch_auto[rs, cs]).sum())
                    left += int((mm & scratch_auto[rs, cs] & scratch[rs, cs]).sum())
            if tot == 0 or left >= 0.5 * tot:
                kept.append(ln)
        n_dropped = len(scratch_lines) - len(kept)
        scratch_lines = kept
    else:
        n_dropped = 0
    n_add = sum(1 for s in strokes if s["modus"] == "add")
    n_rem = len(strokes) - n_add
    HAS_SCRATCH = bool(S.scratch or strokes)
    if HAS_SCRATCH:
        say(tx(f"Kratzer: {len(scratch_lines)} automatisch erkannt, {len(strokes)} manuelle Markierung(en), "
               f"{dez(scratch.mean() * 100, 2)} % der Fläche",
               f"Scratches: {len(scratch_lines)} detected automatically, {len(strokes)} manual mark(s), "
               f"{dez(scratch.mean() * 100, 2)} % of the area"))

    # ------------------------------------------------------------------ Phasenbereinigung
    prog(15, tx("Phasen bereinigen und ergänzen", "Cleaning and filling phases"))
    vr = (ph == iAL) & ~scratch
    zero = (ph == 0) & ~scratch
    op = ndi.binary_opening(zero, structure=S8)
    lab, n = ndi.label(op, structure=S8)
    sz = np.bincount(lab.ravel())[1:]
    mbc = ndi.mean(bc, lab, np.arange(1, n + 1)) if n else np.zeros(0)
    keep = np.zeros(n + 1, bool)
    keep[1:] = (sz >= MG2SI_MIN_PX) & (mbc < S.mg2si_max_bc)
    mg = ndi.binary_fill_holes(keep[lab])
    mg |= ndi.binary_dilation(mg, S8) & zero & (bc < S.mg2si_max_bc)
    if scratch.any():
        mg |= ndi.binary_closing(mg, structure=np.ones((7, 7), bool)) & scratch
    n_mg = ndi.label(mg, structure=S8)[1]

    part_keep = {}
    for ip in iPARTS:
        a0 = (ph == ip) & ~scratch
        la, _ = ndi.label(a0, structure=S8)
        sza = np.bincount(la.ravel())
        sza[0] = 0
        part_keep[ip] = sza[la] >= ALPHA_MIN_PX

    L = np.zeros((NY, NX), np.int16)
    L[vr] = AL
    for ip in iPARTS:
        L[part_keep[ip]] = PCODE[ip]
    L[mg] = CLUSTER
    idx = np.arange(NPIX).reshape(NY, NX)
    codes = [AL] + [PCODE[ip] for ip in iPARTS]
    K8 = np.ones((3, 3), np.uint8)
    K8[1, 1] = 0
    for kmin in range(7, 0, -1):
        if cancel is not None:
            cancel()
        while True:
            cnt = np.stack([ndi.correlate((L == c).astype(np.uint8), K8, mode="constant") for c in codes])
            tot = cnt.sum(0)
            cy, cx = np.nonzero((L == 0) & (tot >= kmin))
            if not len(cy):
                break
            L[cy, cx] = np.asarray(codes)[np.argmax(cnt[:, cy, cx], axis=0)]
    src = np.full((NY, NX), -1, np.int64)
    origs = [(AL, vr)] + [(PCODE[ip], part_keep[ip]) for ip in iPARTS]
    for code, orig in origs:
        if orig.any():
            _, (iy, ix) = ndi.distance_transform_edt(~orig, return_indices=True)
            m = L == code
            src[m] = (iy * NX + ix)[m]
    L[L == 0] = CLUSTER
    clean_frac = {"Al": (L == AL).mean() * 100, "Cluster": (L == CLUSTER).mean() * 100}
    for ip in iPARTS:
        clean_frac[ip] = (L == PCODE[ip]).mean() * 100
    PART = np.isin(L, list(PCODE.values())) if iPARTS else np.zeros((NY, NX), bool)

    Qf = Q.reshape(-1, 4)
    Qc = Qf.copy()
    fill = (src.ravel() >= 0) & (src.ravel() != idx.ravel())
    Qc[fill] = Qf[src.ravel()[fill]]
    ebsdC = ebsd.copy()
    pc = np.zeros((NY, NX), np.int8)
    pc[L == AL] = iAL
    for ip in iPARTS:
        pc[L == PCODE[ip]] = ip
    ebsdC.phaseId = pc
    ebsdC.rotations = rotation(Qc.reshape(NY, NX, 4))

    # ------------------------------------------------------------------ KAM
    kam, kv = None, np.zeros(0)
    if S.kam or S.rx:
        prog(28, tx("KAM berechnen", "Calculating KAM"))
        ebsdK = ebsd.copy()
        ebsdK.phaseId = np.where(vr, iAL, 0).astype(np.int8)
        kam = np.asarray(ebsdK.KAM(radius=1.01 * np.sqrt(2) * STEP, threshold=S.kam_max * degree), float).reshape(NY, NX) / degree
        kam[~vr] = np.nan
        kv = kam[np.isfinite(kam)]
        del ebsdK
        say(tx(f"KAM: Mittel {dez(kv.mean(), 3)}°, Median {dez(np.median(kv), 3)}°",
               f"KAM: mean {dez(kv.mean(), 3)}°, median {dez(np.median(kv), 3)}°"))

    # ------------------------------------------------------------------ Körner
    need_grains = S.ipf or S.rx or S.grainsize
    grainsL = GL = None
    nGL = 0
    if need_grains:
        prog(33, tx("Körner rekonstruieren (MTEX calcGrains)", "Reconstructing grains (MTEX calcGrains)"))

        def grains_at(theta):
            e = ebsdC.copy()
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                g = calcGrains(e, angle=theta * degree)
            G = np.asarray(e.grainId).reshape(NY, NX).astype(np.int64)
            pos = np.full(int(max(G.max(), g.id.max())) + 1, -1, np.int64)
            pos[g.id] = np.arange(len(g.id))
            return g, np.where(G > 0, pos[np.maximum(G, 0)], -1)

        grainsH, GH = grains_at(S.hagb)
        if cancel is not None:
            cancel()
        grainsS, GS = grains_at(S.rx_theta_c) if (S.rx and SUBG) else (None, None)
        if SUBG:
            if S.rx and abs(S.rx_theta_c - S.lagb) < 1e-9:
                grainsL, GL = grainsS, GS
            else:
                grainsL, GL = grains_at(S.lagb)
            nGL = len(grainsL.id)
        nGH = len(grainsH.id)
        nGS = len(grainsS.id) if grainsS is not None else 0
        if SUBG:
            say(tx(f"Körner: {nGH} (≥ {dez(S.hagb, 0)}°), Subkörner: {nGL} (≥ {dez(S.lagb, 1)}°)",
                   f"Grains: {nGH} (≥ {dez(S.hagb, 0)}°), subgrains: {nGL} (≥ {dez(S.lagb, 1)}°)")
                + (tx(f", Subkörner für GOS: {nGS} (≥ {dez(S.rx_theta_c, 1)}°)",
                      f", subgrains for GOS: {nGS} (≥ {dez(S.rx_theta_c, 1)}°)") if S.rx and grainsS is not grainsL else ""))
        else:
            say(tx(f"Körner: {nGH} (≥ {dez(S.hagb, 0)}°); Subkornerkennung ausgeschaltet",
                   f"Grains: {nGH} (≥ {dez(S.hagb, 0)}°); subgrain detection switched off"))

        def al_segments(gB):
            try:
                g = gB[PH_MATRIX, PH_MATRIX]
            except Exception:
                return None, np.zeros(0), np.zeros(0)
            if not len(g):
                return None, np.zeros(0), np.zeros(0)
            return g, g.misorientation.angle() / degree, g.segLength

        gbH, aH, _ = al_segments(grainsH.boundary)
        gbHi, aHi, _ = al_segments(grainsH.innerBoundary)
        gbHi = gbHi[aHi >= S.hagb] if gbHi is not None and (aHi >= S.hagb).any() else None
        if SUBG:
            gbL1, aL1, sL1 = al_segments(grainsL.boundary)
            gbL2, aL2, sL2 = al_segments(grainsL.innerBoundary)
            gbL1 = gbL1[(aL1 >= S.lagb) & (aL1 < S.hagb)] if gbL1 is not None and ((aL1 >= S.lagb) & (aL1 < S.hagb)).any() else None
            gbL2 = gbL2[(aL2 >= S.lagb) & (aL2 < S.hagb)] if gbL2 is not None and ((aL2 >= S.lagb) & (aL2 < S.hagb)).any() else None
            a_all, s_all = np.concatenate([aL1, aL2]), np.concatenate([sL1, sL2])
            hagb_mtex = s_all[a_all >= S.hagb].sum() / max(s_all[a_all >= S.lagb].sum(), 1e-9) * 100
        else:
            gbL1 = gbL2 = None
            hagb_mtex = np.nan

    # ------------------------------------------------------------------ Linienschnitt
    if S.gb or S.grainsize:
        prog(50, tx("Linienschnitt (Korngrenzen, Sehnenlängen)", "Line intercept (grain boundaries, intercept lengths)"))

        def intercept_angles(G):
            out = {}
            for name, V, LL, QQ in [("RD", vr, L, Q), ("CD", vr.T, L.T, Q.transpose(1, 0, 2))]:
                yy, xx = np.nonzero(V)
                sel = (yy[1:] == yy[:-1]) & ((xx[1:] - xx[:-1] - 1) <= G)
                ay, ax_, bx = yy[:-1][sel], xx[:-1][sel], xx[1:][sel]
                g = bx - ax_ - 1
                ok = np.ones(len(ay), bool)
                for k in range(1, G + 1):
                    m = g >= k
                    ok[m] &= LL[ay[m], ax_[m] + k] == AL
                ay, ax_, bx, g = ay[ok], ax_[ok], bx[ok], g[ok]
                ang = mis_angle(QQ[ay, ax_], QQ[ay, bx]) if len(ay) else np.zeros(0)
                out[name] = (ang, (g + 1).sum() * STEP)
            return out

        ic = intercept_angles(int(S.gap_bridge))
        mli = {k: ic[k][1] / max((ic[k][0] >= S.hagb).sum(), 1) for k in ("RD", "CD")}
        mli2 = {k: ic[k][1] / max((ic[k][0] >= S.lagb).sum(), 1) for k in ("RD", "CD")}
        tot_len = ic["RD"][1] + ic["CD"][1]
        n_h = sum(int((ic[k][0] >= S.hagb).sum()) for k in ("RD", "CD"))
        n_l = sum(int((ic[k][0] >= S.lagb).sum()) for k in ("RD", "CD"))
        mli["ges"] = tot_len / max(n_h, 1)
        mli2["ges"] = tot_len / max(n_l, 1)
        mli["n"], mli2["n"] = n_h, n_l
        if not SUBG:                                     # Sehnenlängen der Subkörner nur mit Subkornerkennung
            mli2 = {"RD": np.nan, "CD": np.nan, "ges": np.nan, "n": 0}

    if S.gb:
        gb_ang = np.concatenate([ic["RD"][0], ic["CD"][0]])
        gb = gb_ang[gb_ang >= S.lagb]
        f_lagb = (gb < S.hagb).mean() * 100 if gb.size else np.nan
        f_hagb = 100 - f_lagb
        cls = [((gb >= a) & (gb < b)).mean() * 100 if gb.size else np.nan for a, b in [(S.lagb, 5), (5, S.hagb), (S.hagb, 70)]]
        g0 = intercept_angles(0)
        a0 = np.concatenate([g0["RD"][0], g0["CD"][0]])
        a0 = a0[a0 >= S.lagb]
        hagb_nogap = (a0 >= S.hagb).mean() * 100 if a0.size else np.nan
        dens_mk, w_mk = calcAngleDistribution(csAl, csAl)
        say(f"LAGB {dez(f_lagb)} %, HAGB {dez(f_hagb)} %")

    # ------------------------------------------------------------------ Textur
    dirs = {"CD": yvector, "RD": xvector, "TD": zvector}
    dir_long = {"CD": tx("Stauchrichtung CD", "compression direction CD"),
                "RD": tx("Radialrichtung RD", "radial direction RD"),
                "TD": tx("Schliffnormale TD", "section normal TD")}
    fib = fibre(Miller(1, 1, 0, csAl), yvector)
    if S.texture:
        prog(58, tx("Textur (ODF, Polfiguren)", "Texture (ODF, pole figures)"))
        oriAl = ebsdC[PH_MATRIX].orientations
        odf = calcDensity(oriAl, halfwidth=S.pf_halfwidth * degree)
        ipf_max = {k: float(odf.radon(r=v).max()[0]) for k, v in dirs.items()}
        hfam = {"{100}": Miller(1, 0, 0, csAl), "{110}": Miller(1, 1, 0, csAl), "{111}": Miller(1, 1, 1, csAl)}
        pf_max = {k: float(odf.radon(h).max()[0]) for k, h in hfam.items()}
        fib_all = (fib.angle(oriAl) / degree < 15).mean() * 100

    # ------------------------------------------------------------------ Rekristallisation
    RX, SUB, DEF, NC = 1, 2, 3, 4
    if S.rx:
        prog(66, tx("Rekristallisation (GOS)", "Recrystallization (GOS)"))

        def gos_raw(Gp, grains):
            n = len(grains.id)
            m = vr & (Gp >= 0)
            g = Gp[m]
            mean = np.asarray(grains.meanRotation.data, float).reshape(-1, 4)
            ang = mis_angle(Q[m], mean[g])
            nidx = np.bincount(g, minlength=n)
            return np.bincount(g, weights=ang, minlength=n) / np.maximum(nidx, 1), nidx

        gos_g, nidx_g = gos_raw(GH, grainsH)
        m_al = (L == AL) & (GH >= 0)
        if SUBG:
            gos_s, nidx_s = gos_raw(GS, grainsS)
            m_al &= GS >= 0
        area_g = np.bincount(GH[m_al], minlength=nGH)
        A_al = max(area_g.sum(), 1)
        valid_g = nidx_g >= S.rx_min_px
        cls_g = np.full(nGH, NC, np.int8)
        cls_g[valid_g & (gos_g <= S.rx_theta_c)] = RX
        if SUBG:
            area_s = np.bincount(GS[m_al], minlength=nGS)
            sub_to_grain = np.full(nGS, -1)
            sub_to_grain[GS[m_al]] = GH[m_al]
            sub_bad = (gos_s > S.rx_theta_c) & (nidx_s >= S.rx_min_px) & (sub_to_grain >= 0)
            grain_has_bad = np.bincount(sub_to_grain[sub_bad], minlength=nGH) > 0
            cls_g[valid_g & (gos_g > S.rx_theta_c) & ~grain_has_bad] = SUB
            cls_g[valid_g & (gos_g > S.rx_theta_c) & grain_has_bad] = DEF
            c5_frac = {k: area_g[cls_g == v].sum() / A_al * 100 for k, v in [("rx", RX), ("sub", SUB), ("def", DEF), ("nc", NC)]}
            sub_cls = np.where(nidx_s < S.rx_min_px, NC, np.where(gos_s <= S.rx_theta_c, SUB, DEF)).astype(np.int8)
            CL = np.zeros((NY, NX), np.int8)
            CL[m_al] = sub_cls[GS[m_al]]
            gcls_pix = np.zeros((NY, NX), np.int8)
            gcls_pix[m_al] = cls_g[GH[m_al]]
            CL[m_al & (gcls_pix == RX)] = RX
            CL[m_al & (gcls_pix == NC) & (CL == SUB)] = NC
        else:                                            # ohne Subkörner: nur Kornebene, rekristallisiert / nicht
            cls_g[valid_g & (gos_g > S.rx_theta_c)] = DEF
            CL = np.zeros((NY, NX), np.int8)
            CL[m_al] = cls_g[GH[m_al]]
        big3 = (np.sort(area_g)[::-1][:3] / A_al * 100).tolist()
        rx_frac = {k: (CL == v).sum() / A_al * 100 for k, v in [("rx", RX), ("sub", SUB), ("def", DEF), ("nc", NC)]}
        if not SUBG:
            rx_frac["sub"] = np.nan
        nonrx_frac = (np.nansum([rx_frac["sub"], rx_frac["def"]]))
        if cancel is not None:
            cancel()

        if SUBG:
            perim = np.zeros(nGS)
            hperim = np.zeros(nGS)
            for dy, dx in ((0, 1), (1, 0)):
                A1, B1 = GS[:NY - dy, :NX - dx], GS[dy:, dx:]
                ok = (A1 >= 0) & (B1 >= 0) & (A1 != B1) & (L[:NY - dy, :NX - dx] == AL) & (L[dy:, dx:] == AL)
                ia = idx[:NY - dy, :NX - dx][ok]
                ib = idx[dy:, dx:][ok]
                hi = mis_angle(Qc[ia], Qc[ib]) >= S.hagb
                for lab_ in (A1[ok], B1[ok]):
                    perim += np.bincount(lab_, minlength=nGS)
                    hperim += np.bincount(lab_, weights=hi.astype(float), minlength=nGS)
            f_hagb_perim = hperim / np.maximum(perim, 1)
            low_sub = (gos_s <= S.rx_theta_c) & (nidx_s >= S.rx_min_px)
            fh_scan = np.arange(0.2, 1.001, 0.05)
            rx_sub_scan = np.array([area_s[low_sub & (f_hagb_perim >= f_ - 1e-9)].sum() / A_al * 100 for f_ in fh_scan])
            rx_sub = {f_: area_s[low_sub & (f_hagb_perim >= f_ - 1e-9)].sum() / A_al * 100 for f_ in (0.5, 0.75, 0.9)}
            low_sub_frac = area_s[low_sub].sum() / A_al * 100
        theta_scan = np.arange(0.5, 5.01, 0.25)
        rx_scan = np.array([area_g[valid_g & (gos_g <= t_)].sum() / A_al * 100 for t_ in theta_scan])

        yy, xx = np.nonzero(m_al)
        gg = GH[m_al]
        nn = np.maximum(area_g, 1)

        def gsum(w):
            return np.bincount(gg, weights=w, minlength=nGH)

        mx, my = gsum(xx) / nn, gsum(yy) / nn
        cxx = gsum(xx * xx.astype(float)) / nn - mx ** 2 + 1 / 12
        cyy = gsum(yy * yy.astype(float)) / nn - my ** 2 + 1 / 12
        cxy = gsum(xx * yy.astype(float)) / nn - mx * my
        tr_, det_ = cxx + cyy, cxx * cyy - cxy ** 2
        dsc = np.sqrt(np.maximum(tr_ ** 2 / 4 - det_, 0))
        aspect_g = np.sqrt((tr_ / 2 + dsc) / np.maximum(tr_ / 2 - dsc, 1e-6))
        ecd_g = 2 * np.sqrt(area_g * PX_AREA / np.pi)
        rxg = cls_g == RX
        touch = np.zeros(nGH, bool)
        for edge in (GH[0], GH[-1], GH[:, 0], GH[:, -1]):
            touch[edge[edge >= 0]] = True
        n_rx = int(rxg.sum())
        rx_ecd_mean = ecd_g[rxg].mean() if n_rx else np.nan
        rx_ecd_med = np.median(ecd_g[rxg]) if n_rx else np.nan
        rx_ar_med = np.median(aspect_g[rxg]) if n_rx else np.nan
        rx_equiax = (aspect_g[rxg] < 2).mean() * 100 if n_rx else np.nan
        kam_cls = {k: (np.nanmean(kam[(CL == v) & vr]) if ((CL == v) & vr).any() else np.nan)
                   for k, v in [("rx", RX), ("sub", SUB), ("def", DEF)]}

        def fiber_share(mask):
            m = mask & vr
            if not m.any():
                return np.nan
            return (fib.angle(orientation(rotation(Q[m]), csAl)) / degree < 15).mean() * 100

        fib_cls = {k: fiber_share(CL == v) for k, v in [("rx", RX), ("sub", SUB), ("def", DEF)]}

        part = PART | ((L == CLUSTER) if S.mg2si else False)
        plab, _ = ndi.label(part, structure=S8)
        psz = np.bincount(plab.ravel())
        psz[0] = 0
        big_part = psz[plab] >= PSN_MIN_PX
        dist_um = ndi.distance_transform_edt(~big_part) * STEP if big_part.any() else np.full((NY, NX), np.inf)
        dist_bins = [0, 1, 2, 3, 5, 10, np.inf]
        dist_lbl = ["0–1", "1–2", "2–3", "3–5", "5–10", "> 10"]
        classified = (CL > 0) & (CL != NC)
        rx_vs_dist, area_vs_dist = [], []
        for a_, b_ in zip(dist_bins[:-1], dist_bins[1:]):
            m = classified & (dist_um > a_) & (dist_um <= b_)
            rx_vs_dist.append((CL[m] == RX).mean() * 100 if m.any() else np.nan)
            area_vs_dist.append(m.mean() * 100)
        psn_names = ", ".join([pretty(names[i]) for i in iPARTS] + (["Mg$_2$Si"] if S.mg2si else [])) or "–"
        say(tx(f"rekristallisiert {dez(rx_frac['rx'], 2)} %, {n_rx} Körner",
               f"recrystallized {dez(rx_frac['rx'], 2)} %, {n_rx} grains"))

    CLS_NAME = {RX: tx("rekristallisiert", "recrystallized"), SUB: tx("substrukturiert", "substructured"),
                DEF: tx("verformt", "deformed") if SUBG else tx("nicht rekristallisiert", "non-recrystallized"),
                NC: tx("nicht klassifizierbar", "not classifiable")}

    # ------------------------------------------------------------------ Korngröße und Kornform
    if S.grainsize:
        prog(71, tx("Korngröße, Subkorngröße, Aspektverhältnis", "Grain size, subgrain size, aspect ratio")
             if SUBG else tx("Korngröße, Aspektverhältnis", "Grain size, aspect ratio"))
        m_alg = L == AL
        A_AL_UM2 = max(m_alg.sum(), 1) * PX_AREA

        def shape_stats(Gp, n, grains):
            """Fläche, ECD und Ellipse gleicher zweiter Flächenmomente je (Sub-)Korn (nur Matrixpunkte)"""
            m = m_alg & (Gp >= 0)
            yy, xx = np.nonzero(m)
            gg = Gp[m]
            xf, yf = xx.astype(float), yy.astype(float)
            area = np.bincount(gg, minlength=n).astype(float)
            nn = np.maximum(area, 1)

            def s_(w):
                return np.bincount(gg, weights=w, minlength=n)

            mx, my = s_(xf) / nn, s_(yf) / nn
            cxx = s_(xf * xf) / nn - mx ** 2 + 1 / 12
            cyy = s_(yf * yf) / nn - my ** 2 + 1 / 12
            cxy = s_(xf * yf) / nn - mx * my
            tr_, det_ = cxx + cyy, cxx * cyy - cxy ** 2
            dsc = np.sqrt(np.maximum(tr_ ** 2 / 4 - det_, 0))
            l1 = np.maximum(tr_ / 2 + dsc, 1 / 12)
            l2 = np.clip(tr_ / 2 - dsc, 1 / 12, None)
            # Hauptachsenwinkel zu RD, gegen den Uhrzeigersinn in Kartenansicht (Bildzeilen laufen nach unten)
            ang = -np.degrees(0.5 * np.arctan2(2 * cxy, cxx - cyy))
            ang = np.where(ang <= -90, ang + 180, ang)
            touch = np.zeros(n, bool)
            for edge in (Gp[0], Gp[-1], Gp[:, 0], Gp[:, -1]):
                touch[edge[edge >= 0]] = True
            mi = vr & (Gp >= 0)
            nidx = np.bincount(Gp[mi], minlength=n)
            present = area > 0
            valid = present & ~touch & (nidx >= S.gs_min_px)
            eul = np.degrees(np.column_stack(rotation(np.asarray(grains.meanRotation.data, float).reshape(-1, 4)).Euler()))
            return dict(n=n, area=area * PX_AREA, ecd=2 * np.sqrt(area * PX_AREA / np.pi), major=4 * np.sqrt(l1) * STEP,
                        minor=4 * np.sqrt(l2) * STEP, ar=np.sqrt(l1 / l2), ang=ang, cx=mx * STEP, cy=my * STEP,
                        touch=touch, nidx=nidx, present=present, valid=valid, ids=np.asarray(grains.id), euler=eul)

        SUM_KEYS = ("ecd_mean", "ecd_sd", "ecd_med", "d10", "d90", "ecd_aw", "area_mean", "area_frac", "ar_mean", "ar_med",
                    "ar_aw", "ar_lt2", "ang_mean", "ang_R", "rd_frac", "n_el")

        def size_summary(d):
            if d is None:
                out = dict(n=0, n_all=0, n_edge=0, n_small=0)
                out.update({k: np.nan for k in SUM_KEYS})
                return out
            v = d["valid"]
            e, a, ar, ang = d["ecd"][v], d["area"][v], d["ar"][v], d["ang"][v]
            n = int(v.sum())
            out = dict(n=n, n_all=int(d["present"].sum()), n_edge=int((d["present"] & d["touch"]).sum()),
                       n_small=int((d["present"] & ~d["touch"] & (d["nidx"] < S.gs_min_px)).sum()))
            out.update({k: np.nan for k in SUM_KEYS})
            if not n:
                return out
            out.update(ecd_mean=e.mean(), ecd_sd=e.std(ddof=1) if n > 1 else np.nan, ecd_med=np.median(e),
                       d10=np.percentile(e, 10), d90=np.percentile(e, 90), ecd_aw=(a * e).sum() / a.sum(),
                       area_mean=a.mean(), area_frac=a.sum() / A_AL_UM2 * 100, ar_mean=ar.mean(), ar_med=np.median(ar),
                       ar_aw=(a * ar).sum() / a.sum(), ar_lt2=(ar < 2).mean() * 100)
            el = ar >= 1.5                                     # Richtung nur für gestreckte Körner aussagekräftig
            out["n_el"] = int(el.sum())
            if el.any():
                w = a[el]
                c2, s2 = (w * np.cos(np.radians(2 * ang[el]))).sum() / w.sum(), (w * np.sin(np.radians(2 * ang[el]))).sum() / w.sum()
                out.update(ang_mean=0.5 * np.degrees(np.arctan2(s2, c2)), ang_R=np.hypot(c2, s2),
                           rd_frac=w[np.abs(ang[el]) < 45].sum() / w.sum() * 100)
            return out

        gsH = shape_stats(GH, nGH, grainsH)
        gsL = shape_stats(GL, nGL, grainsL) if SUBG else None
        sumH, sumL = size_summary(gsH), size_summary(gsL)
        astm_G = (-6.643856 * np.log10(mli["ges"] / 1000) - 3.288) if mli["n"] else np.nan
        if SUBG:
            # Zuordnung Subkorn → Korn (Subkörner liegen vollständig in einem Korn, da HAGB ⊂ Grenzen ≥ LAGB)
            sub_parent = np.full(nGL, -1, np.int64)
            mm = m_alg & (GL >= 0) & (GH >= 0)
            sub_parent[GL[mm]] = GH[mm]
            sub_per_grain = np.bincount(sub_parent[gsL["present"] & (sub_parent >= 0)], minlength=nGH)
            spg_mean = sub_per_grain[gsH["present"]].mean() if gsH["present"].any() else np.nan
            say(tx(f"Korngröße (ECD, Anzahlmittel): Körner {dez(sumH['ecd_mean'], 2)} µm (n = {sumH['n']}), "
                   f"Subkörner {dez(sumL['ecd_mean'], 2)} µm (n = {sumL['n']}); Aspektverhältnis Körner (Median) {dez(sumH['ar_med'], 2)}",
                   f"Grain size (ECD, number mean): grains {dez(sumH['ecd_mean'], 2)} µm (n = {sumH['n']}), "
                   f"subgrains {dez(sumL['ecd_mean'], 2)} µm (n = {sumL['n']}); grain aspect ratio (median) {dez(sumH['ar_med'], 2)}"))
        else:
            spg_mean = np.nan
            say(tx(f"Korngröße (ECD, Anzahlmittel): Körner {dez(sumH['ecd_mean'], 2)} µm (n = {sumH['n']}); "
                   f"Aspektverhältnis (Median) {dez(sumH['ar_med'], 2)}",
                   f"Grain size (ECD, number mean): grains {dez(sumH['ecd_mean'], 2)} µm (n = {sumH['n']}); "
                   f"aspect ratio (median) {dez(sumH['ar_med'], 2)}"))
        # Auflösung der Subkörner: Median-ECD in Schrittweiten
        sub_res = sumL["ecd_med"] / STEP if SUBG and sumL["n"] else np.nan
        sub_warn = None
        if np.isfinite(sub_res) and sub_res < 5:
            sub_warn = tx(f"Hinweis: Der Median-ECD der Subkörner entspricht nur {dez(sub_res, 1)} Schrittweiten – "
                          "die Subkörner sind mit dieser Schrittweite kaum aufgelöst. Subkornerkennung ggf. abschalten.",
                          f"Note: the median ECD of the subgrains is only {dez(sub_res, 1)} step sizes – the subgrains are "
                          "hardly resolved at this step size. Consider switching subgrain detection off.")
            say(sub_warn)

        def grain_table(d, sub):
            head = [tx("MTEX-ID", "MTEX ID"), "x / µm", "y / µm", tx("Punkte (Matrix, bereinigt)", "points (matrix, cleaned)"),
                    tx("davon indiziert", "thereof indexed"), tx("Fläche / µm²", "area / µm²"), "ECD / µm",
                    tx("Hauptachse / µm", "major axis / µm"), tx("Nebenachse / µm", "minor axis / µm"),
                    tx("Aspektverhältnis", "aspect ratio"), tx("Winkel Hauptachse zu RD / °", "angle major axis to RD / °"),
                    tx("Randkorn", "edge grain"), tx("in Statistik", "in statistics"), "phi1 / °", "Phi / °", "phi2 / °"]
            if sub:
                head.append(tx("MTEX-ID Korn", "MTEX ID grain"))
            elif S.rx:
                head += ["GOS / °", tx("Klasse", "class")]
            rows = [csv_join(head)]
            for i in np.nonzero(d["present"])[0]:
                r = [str(d["ids"][i]), csv_num(d["cx"][i], 2), csv_num(d["cy"][i], 2), str(int(round(d["area"][i] / PX_AREA))),
                     str(int(d["nidx"][i])), csv_num(d["area"][i], 3), csv_num(d["ecd"][i]), csv_num(d["major"][i]),
                     csv_num(d["minor"][i]), csv_num(d["ar"][i]), csv_num(d["ang"][i], 1), str(int(d["touch"][i])),
                     str(int(d["valid"][i])), *[csv_num(v_, 2) for v_ in d["euler"][i]]]
                if sub:
                    r.append(str(gsH["ids"][sub_parent[i]]) if sub_parent[i] >= 0 else "")
                elif S.rx:
                    r += [csv_num(gos_g[i], 3), CLS_NAME[int(cls_g[i])]]
                rows.append(csv_join(r))
            return "\n".join(rows) + "\n"

    # ======================================================================== Abbildungen
    prog(74, tx("Abbildungen vorbereiten", "Preparing figures"))
    lay = MapLayout(ebsd, STEP)
    TITLE = S.title.strip()
    T_PRE = f"{TITLE} – " if TITLE else ""
    PART_COLS = [np.array([0.08, 0.08, 0.08]), np.array([0.45, 0.22, 0.05]), np.array([0.1, 0.25, 0.55]),
                 np.array([0.4, 0.0, 0.4]), np.array([0.0, 0.35, 0.3])]
    PCOL = {ip: PART_COLS[k % len(PART_COLS)] for k, ip in enumerate(iPARTS)}
    COL_MG = np.array([0.55, 0.55, 0.55])
    COL_LAGB = (0.72, 0.72, 0.72)
    COL_CLS = {RX: (0.12, 0.35, 0.95), SUB: (1.0, 0.85, 0.1), DEF: (0.85, 0.12, 0.12), NC: (0.8, 0.8, 0.8)}
    COL_SCR_AUTO, COL_SCR_ADD, COL_SCR_REM = (0.95, 0.1, 0.1), (0.85, 0.1, 0.85), (0.25, 0.7, 1.0)
    LW_HAGB, LW_LAGB = lay.lw_for(1.0), lay.lw_for(0.6)
    sub_info = tx(f"{NX} × {NY} Punkte, Schrittweite {dez(STEP, 3)} µm, {dez(W_UM, 0)} × {dez(H_UM, 0)} µm",
                  f"{NX} × {NY} points, step size {dez(STEP, 3)} µm, {dez(W_UM, 0)} × {dez(H_UM, 0)} µm")
    if HAS_SCRATCH:
        how = (tx(f"{len(scratch_lines)} automatisch", f"{len(scratch_lines)} automatic") if S.scratch else "") \
            + (("; " if S.scratch else "") + tx("manuell bearbeitet", "edited manually") if strokes else "")
        sub_info += tx(f"  |  Kratzer entfernt ({how}, {dez(scratch.mean() * 100, 1)} % der Fläche)",
                       f"  |  scratches removed ({how}, {dez(scratch.mean() * 100, 1)} % of the area)")
    if not SUBG and need_grains:
        sub_info += tx("  |  ohne Subkornerkennung", "  |  without subgrain detection")
    b_part = [Patch(fc=PCOL[ip], label=pretty(names[ip])) for ip in iPARTS] + [Patch(fc=COL_MG, label=CLUSTER_LABEL)]
    b_gb = [Line2D([], [], color="k", lw=2, label=f"HAGB ≥ {dez(S.hagb, 0)}°")]
    if SUBG:
        b_gb.append(Line2D([], [], color=COL_LAGB, lw=2, label=f"LAGB {dez(S.lagb, 0)}–{dez(S.hagb, 0)}°"))
    bc_img = (bc / max(bc.max(), 1))[..., None] * np.ones(3) * 0.85 + 0.1
    lbl_bc = tx("Bandkontrast", "band contrast")

    def paint_particles(rgb):
        for ip in iPARTS:
            rgb[L == PCODE[ip]] = PCOL[ip]
        rgb[L == CLUSTER] = COL_MG
        return rgb

    builders = {}

    # --- 00 Kratzer
    def fig_00():
        rgb = bc_img.copy()
        rgb[scratch & ~scr_added] = COL_SCR_AUTO
        rgb[scr_added] = COL_SCR_ADD
        rgb[scr_removed] = COL_SCR_REM
        lines_txt = [f"{k:>2}  {dez(l['angle'], 0):>4}°  {dez(l['length_um'], 0):>5} µm  {dez(l['width_um'], 1):>4} µm"
                     for k, l in enumerate(scratch_lines, 1)][:20]
        if len(scratch_lines) > 20:
            lines_txt.append(tx(f"… und {len(scratch_lines) - 20} weitere", f"… and {len(scratch_lines) - 20} more"))
        if S.scratch:
            info = [tx(f"automatisch erkannt: {len(scratch_lines)} Kratzer", f"detected automatically: {len(scratch_lines)} scratch" + ("" if len(scratch_lines) == 1 else "es"))]
            if n_dropped:
                info.append(tx(f"   ({n_dropped} weitere manuell aufgehoben)", f"   ({n_dropped} more removed manually)"))
        else:
            info = [tx("automatische Erkennung ausgeschaltet", "automatic detection switched off")]
        if strokes:
            info += [tx(f"manuell: {n_add} markiert, {n_rem} aufgehoben", f"manual: {n_add} marked, {n_rem} unmarked"),
                     tx(f"   +{dez(scr_added.mean() * 100, 2)} % / −{dez(scr_removed.mean() * 100, 2)} % der Fläche",
                        f"   +{dez(scr_added.mean() * 100, 2)} % / −{dez(scr_removed.mean() * 100, 2)} % of the area")]
        info.append(tx(f"markiert gesamt: {dez(scratch.mean() * 100, 2)} % der Fläche",
                       f"marked in total: {dez(scratch.mean() * 100, 2)} % of the area"))
        if S.scratch:
            info += [tx(f"Mindestlänge {dez(S.scratch_min_len, 0)} µm, Fenster {dez(S.scratch_window, 0)} µm",
                        f"min. length {dez(S.scratch_min_len, 0)} µm, window {dez(S.scratch_window, 0)} µm"),
                     tx("Kandidaten: Nulllösung, Fehlindizierung,", "candidates: zero solution, misindexing,"),
                     tx(f"   BC < {dez(S.scratch_rel_bc * 100, 0)} % der Umgebung",
                        f"   BC < {dez(S.scratch_rel_bc * 100, 0)} % of the surroundings")]
        blocks = [text_block(info, head=tx("Kratzererkennung", "Scratch detection"))]
        if lines_txt:
            blocks.append(text_block([tx("Nr  Richt.*    Länge  Breite", "No  Dir.*    Length  Width")] + lines_txt
                                     + [tx("* Winkel zu RD, gegen den Uhrzeigersinn", "* angle to RD, counter-clockwise")],
                                     fs=9.5, mono=True))
        leg = [Patch(fc=(0.6, 0.6, 0.6), label=lbl_bc)]
        if S.scratch:
            leg.append(Patch(fc=COL_SCR_AUTO, label=tx("automatisch als Kratzer markiert", "marked as scratch automatically")))
        if scr_added.any():
            leg.append(Patch(fc=COL_SCR_ADD, label=tx("manuell als Kratzer markiert", "marked as scratch manually")))
        if scr_removed.any():
            leg.append(Patch(fc=COL_SCR_REM, label=tx("Markierung manuell aufgehoben", "marking removed manually")))
        blocks += [legend_block(leg), glyph_block()]
        return lay.figure(rgb, T_PRE + tx("Kratzererkennung", "Scratch detection"),
                          sub_info + tx("  |  markierte Punkte werden aus den Nachbarn ergänzt",
                                        "  |  marked points are filled in from the neighbours"), blocks)

    if HAS_SCRATCH:
        builders["00"] = fig_00

    # --- 01–03 IPF-Karten
    if S.ipf:
        alm = (L == AL).ravel()
        oriC = orientation(rotation(Qc[alm]), csAl)

        def make_ipf(name):
            def f():
                key = ipfTSLKey(csAl)
                key.ipfDirection = dirs[name]
                rgb = np.ones((NPIX, 3))
                rgb[alm] = np.clip(np.asarray(key.orientation2color(oriC)).reshape(-1, 3), 0, 1)
                rgb = paint_particles(rgb.reshape(NY, NX, 3))
                blocks = [ipfkey_block(key, name, f"{MATRIX_LABEL} ({project.phases[iAL]['symmetry']})"),
                          glyph_block(), legend_block(b_part + b_gb)]
                return lay.figure(rgb, T_PRE + tx(f"IPF-Karte ‖ {dir_long[name]}", f"IPF map ‖ {dir_long[name]}"), sub_info, blocks,
                                  boundaries=[(gbL1, COL_LAGB, LW_LAGB), (gbL2, COL_LAGB, LW_LAGB), (gbH, "k", LW_HAGB), (gbHi, "k", LW_HAGB)])
            return f

        builders["01"], builders["02"], builders["03"] = make_ipf("CD"), make_ipf("RD"), make_ipf("TD")

    # --- 04/05 Textur
    if S.texture:
        def fig_04():
            lev = mrd_levels(max(ipf_max.values()))
            fig = plt.figure(figsize=(14, 5.2), dpi=100)
            for k, name in enumerate(["CD", "RD", "TD"]):
                ax = fig.add_axes([0.02 + k * 0.3, 0.06, 0.27, 0.72])
                plotIPF(odf, dirs[name], parent=ax, contourf=lev, colorRange=[lev[0], lev[-1]], colormap="turbo")
                ax.set_title(f"IPF ‖ {dir_long[name]}\nmax. {dez(ipf_max[name], 2)} MRD", fontsize=12)
            level_cbar(fig, fig.add_axes([0.925, 0.15, 0.015, 0.6]), lev)
            t = fig.text(0.02, 0.95, T_PRE + tx(f"Inverse Polfiguren {MATRIX_LABEL} (MTEX, Kernel-Halbwertsbreite {dez(S.pf_halfwidth, 0)}°)",
                                                f"Inverse pole figures {MATRIX_LABEL} (MTEX, kernel half-width {dez(S.pf_halfwidth, 0)}°)"),
                         fontsize=15, weight="bold", va="center")
            fit_text(fig, t, 13.6)
            return fig

        def fig_05():
            lev = mrd_levels(max(pf_max.values()))
            fig = plt.figure(figsize=(13.5, 9.6), dpi=100)
            rows = [("z↓→x", tx("Projektion ‖ CD\n(Stauchachse im Zentrum)", "Projection ‖ CD\n(compression axis in the centre)")),
                    ("y↓→x", tx("Projektion ‖ TD\n(Schliffebene, Kartenlage)", "Projection ‖ TD\n(section plane, map orientation)"))]
            for r_, (conv, rlabel) in enumerate(rows):
                for k, (fam, h) in enumerate(hfam.items()):
                    ax = fig.add_axes([0.1 + k * 0.27, 0.5 - r_ * 0.45, 0.24, 0.36])
                    plotPF(odf, h, parent=ax, contourf=lev, colorRange=[lev[0], lev[-1]], colormap="turbo", how2plot=conv)
                    ax.set_title(f"{fam}   max. {dez(pf_max[fam], 2)} MRD", fontsize=12)
                fig.text(0.04, 0.68 - r_ * 0.45, rlabel, rotation=90, va="center", ha="center", fontsize=12, weight="bold")
            level_cbar(fig, fig.add_axes([0.93, 0.2, 0.015, 0.55]), lev)
            t = fig.text(0.02, 0.965, T_PRE + tx(f"Polfiguren {MATRIX_LABEL} (MTEX, flächentreu, Kernel-Halbwertsbreite {dez(S.pf_halfwidth, 0)}°)",
                                                 f"Pole figures {MATRIX_LABEL} (MTEX, equal-area, kernel half-width {dez(S.pf_halfwidth, 0)}°)"),
                         fontsize=15, weight="bold", va="center")
            fit_text(fig, t, 13.1)
            return fig

        builders["04"], builders["05"] = fig_04, fig_05

    # --- 06 KAM
    if S.kam:
        def fig_06():
            cmap_k = plt.get_cmap("turbo")
            rgb = cmap_k(np.clip(np.nan_to_num(kam) / S.kam_vmax, 0, 1))[..., :3]
            rgb[~np.isfinite(kam)] = 0
            rgb[PART] = 1.0
            rgb[L == CLUSTER] = 0.6
            rgb[scratch & (L == AL)] = 0.35

            def kam_hist(ax):
                hh, be = np.histogram(kv, bins=np.arange(0, S.kam_max + 0.05, 0.05))
                ax.bar(be[:-1], hh / max(hh.sum(), 1) * 100, width=0.05, align="edge",
                       color=cmap_k(np.clip(be[:-1] / S.kam_vmax, 0, 1)), ec="none")
                if kv.size:
                    ax.axvline(kv.mean(), color="k", ls="--", lw=1)
                ax.set_xlim(0, min(3, S.kam_max))
                ax.set_xlabel("KAM / °")
                ax.set_ylabel(tx("Anteil / %", "Fraction / %"))
                ax.xaxis.set_major_formatter(comma)
                ax.yaxis.set_major_formatter(comma)
                mv, md = kv.mean() if kv.size else np.nan, np.median(kv) if kv.size else np.nan
                ax.set_title(tx(f"Mittelwert {dez(mv, 2)}°, Median {dez(md, 2)}°", f"mean {dez(mv, 2)}°, median {dez(md, 2)}°"), fontsize=10)

            blocks = [text_block([tx("MTEX ebsd.KAM, 3×3-Nachbarschaft", "MTEX ebsd.KAM, 3×3 neighbourhood"),
                                  tx(f"Schwellwert {dez(S.kam_max, 0)}°, Schrittweite {dez(STEP, 3)} µm",
                                     f"threshold {dez(S.kam_max, 0)}°, step size {dez(STEP, 3)} µm"),
                                  tx(f"nur original indizierte Punkte ({plain(PH_MATRIX)})",
                                     f"originally indexed points only ({plain(PH_MATRIX)})")],
                                 head="Kernel Average Misorientation"),
                      cbar_block(cmap_k, Normalize(0, S.kam_vmax), "KAM / °"),
                      axes_block(2.1, kam_hist),
                      legend_block([Patch(fc="k", label=tx("nicht indiziert", "non-indexed")),
                                    Patch(fc="white", ec="k", label=tx("Partikelphasen", "particle phases")),
                                    Patch(fc=(0.6, 0.6, 0.6), label=CLUSTER_LABEL)]
                                   + ([Patch(fc=(0.35, 0.35, 0.35), label=tx("Kratzer (ausgeblendet)", "scratches (masked)"))]
                                      if scratch.any() else [])),
                      glyph_block(1.3)]
            return lay.figure(rgb, f"{T_PRE}Kernel Average Misorientation (KAM)", sub_info, blocks)

        builders["06"] = fig_06

    # --- 07 Misorientierung
    if S.gb:
        def fig_07():
            bins = np.arange(S.lagb, 64, 1.0)
            h, _ = np.histogram(gb, bins=bins)
            h = h / max(gb.size, 1) * 100
            cent = bins[:-1] + 0.5
            pm = np.interp(np.radians(cent), w_mk, dens_mk)
            pm = pm / pm.sum()
            hm_scaled = pm / max(pm[bins[:-1] >= S.hagb].sum(), 1e-12) * f_hagb
            fig, ax = plt.subplots(figsize=(10.5, 6), dpi=100)
            cols = np.where(bins[:-1] < S.hagb, "#4C72B0", "#C44E52")
            ax.bar(bins[:-1], h, width=1.0, align="edge", color=cols, ec="white", lw=0.3)
            ax.plot(cent, hm_scaled, "k-", lw=1.5, label=tx("Mackenzie (MTEX), normiert auf HAGB-Anteil",
                                                          "Mackenzie (MTEX), normalized to HAGB fraction"))
            ax.axvline(S.hagb, color="k", ls="--", lw=1)
            ax.set_xlabel(tx("Misorientierungswinkel / °", "Misorientation angle / °"))
            ax.set_ylabel(tx("Anteil an der Korngrenzlänge / %", "Fraction of boundary length / %"))
            ax.set_xlim(0, 65)
            ax.yaxis.set_major_formatter(comma)
            txt = (f"LAGB ({dez(S.lagb, 0)}–{dez(S.hagb, 0)}°):  {dez(f_lagb)} %\nHAGB (≥ {dez(S.hagb, 0)}°):  {dez(f_hagb)} %\n\n"
                   f"{dez(S.lagb, 0)}–5°:  {dez(cls[0])} %\n5–{dez(S.hagb, 0)}°: {dez(cls[1])} %\n≥ {dez(S.hagb, 0)}°:  {dez(cls[2])} %")
            ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha="right", va="top", fontsize=12, family="DejaVu Sans Mono",
                    bbox=dict(fc="white", ec="0.6"))
            ax.legend(handles=[Patch(fc="#4C72B0", label="LAGB"), Patch(fc="#C44E52", label="HAGB"), ax.lines[0]],
                      loc="upper center", frameon=False)
            ax.set_title(T_PRE + tx(f"Korngrenz-Misorientierungsverteilung ({plain(PH_MATRIX)}, Längenanteil)\n"
                                    f"Linienschnitt in RD und CD, Nulllösungen ≤ {int(S.gap_bridge)} Punkte überbrückt, "
                                    f"n = {thou(gb.size)} Grenzschnitte",
                                    f"Grain boundary misorientation distribution ({plain(PH_MATRIX)}, length fraction)\n"
                                    f"Line intercept in RD and CD, zero solutions ≤ {int(S.gap_bridge)} points bridged, "
                                    f"n = {thou(gb.size)} boundary intercepts"), fontsize=12)
            fig.tight_layout()
            return fig

        builders["07"] = fig_07

    # --- 08 Phasenverteilung
    RAW_PCOLS = [(0.85, 0.1, 0.1), (1.0, 0.55, 0.0), (0.75, 0.2, 0.6), (0.55, 0.35, 0.1)]
    MIS_COLS = [(0.1, 0.35, 0.95), (0.1, 0.8, 0.2), (0.6, 0.2, 0.8), (0.0, 0.6, 0.6), (0.95, 0.85, 0.1)]
    if scratch.any():
        _, (iy_, ix_) = ndi.distance_transform_edt(scratch, return_indices=True)
        bc_fill = bc[iy_, ix_]
    else:
        bc_fill = bc
    bc_img_clean = (bc_fill / max(bc.max(), 1))[..., None] * np.ones(3) * 0.85 + 0.1

    def fig_08():
        def phase_rgb(raw):
            rgb = bc_img.copy() if raw else bc_img_clean.copy()
            if raw:
                rgb[ph == 0] = 0
                for k_, ip in enumerate(iPARTS):
                    rgb[ph == ip] = RAW_PCOLS[k_ % len(RAW_PCOLS)]
                for k_, i_ in enumerate(iMIS):
                    rgb[ph == i_] = MIS_COLS[k_ % len(MIS_COLS)]
            else:
                for k_, ip in enumerate(iPARTS):
                    rgb[L == PCODE[ip]] = RAW_PCOLS[k_ % len(RAW_PCOLS)]
                rgb[L == CLUSTER] = [0.0, 0.75, 1.0]
            return rgb

        s2 = min(6.3 / W_UM, 5.2 / H_UM)
        mw2, mh2 = W_UM * s2, H_UM * s2
        pw = max(mw2, 5.0)
        FW = 0.3 + 2 * pw + 0.4 + 0.3
        leg_h, bar_h = 0.55, 2.9
        FH = 0.95 + mh2 + 0.35 + leg_h + bar_h + 0.25
        fig = plt.figure(figsize=(FW, FH), dpi=100)

        def frect(x, ytop, w_, h_):
            return [x / FW, (FH - ytop - h_) / FH, w_ / FW, h_ / FH]

        for k, (raw, ttl) in enumerate([(True, tx("Rohdaten (AZtec-Indizierung)", "raw data (AZtec indexing)")),
                                        (False, tx("bereinigt", "cleaned") + (tx(" (Kratzer aus Nachbarn ergänzt)",
                                                                                " (scratches filled from neighbours)")
                                                                             if scratch.any() else ""))]):
            ax = fig.add_axes(frect(0.3 + k * (pw + 0.4) + (pw - mw2) / 2, 1.3, mw2, mh2))
            lay.draw_map(ax, phase_rgb(raw), (), mw2)
            ax.set_title(ttl, fontsize=12)
        lbl_zero = tx("Nulllösung", "zero solution")
        handles = ([Patch(fc=(0.7,) * 3, label=f"{MATRIX_LABEL} ({lbl_bc})"), Patch(fc="k", label=lbl_zero)]
                   + [Patch(fc=RAW_PCOLS[k_ % len(RAW_PCOLS)], label=pretty(names[ip])) for k_, ip in enumerate(iPARTS)]
                   + [Patch(fc=MIS_COLS[k_ % len(MIS_COLS)], label=f"{pretty(names[i])} ({tx('Fehlindiz.', 'misindexed')})")
                      for k_, i in enumerate(iMIS)]
                   + [Patch(fc=[0, 0.75, 1], label=CLUSTER_LABEL)])
        axl = fig.add_axes(frect(0.3, 1.3 + mh2 + 0.25, FW - 0.6, leg_h))
        axl.axis("off")
        axl.legend(handles=handles, loc="center", ncol=min(len(handles), 4), frameon=False, fontsize=10.5)
        bw = min(7.0, FW * 0.55)
        axb = fig.add_axes(frect(1.0, 1.3 + mh2 + 0.35 + leg_h + 0.1, bw, bar_h - 0.6))
        labels = ([pretty(names[ip]) for ip in iPARTS] + [pretty(names[i]) for i in iMIS]
                  + [lbl_zero, tx("Kratzer", "scratches"), "Mg$_2$Si" if S.mg2si else "Cluster"])
        rawv = [raw_frac[ip] for ip in iPARTS] + [raw_frac[i] for i in iMIS] + [raw_frac[0], scratch.mean() * 100, 0]
        clv = [clean_frac[ip] for ip in iPARTS] + [0] * len(iMIS) + [0, 0, clean_frac["Cluster"]]
        xx_ = np.arange(len(labels))
        b1 = axb.bar(xx_ - 0.2, rawv, 0.4, color="0.65", label=tx("Rohdaten", "raw data"))
        b2 = axb.bar(xx_ + 0.2, clv, 0.4, color="#C44E52", label=tx("bereinigt", "cleaned"))
        top = max(rawv + clv + [0.1])
        for bars in (b1, b2):
            for b in bars:
                if b.get_height() > 0:
                    axb.text(b.get_x() + b.get_width() / 2, b.get_height() + top * 0.012, dez(b.get_height(), 2), ha="center", fontsize=9)
        axb.set_xticks(xx_, labels, fontsize=9.5, rotation=15)
        axb.set_ylabel(tx("Flächenanteil / %", "Area fraction / %"))
        axb.yaxis.set_major_formatter(comma)
        axb.legend(frameon=False)
        axb.set_ylim(0, top * 1.2)
        axt = fig.add_axes(frect(1.0 + bw + 0.5, 1.3 + mh2 + 0.35 + leg_h + 0.1, max(FW - (1.5 + bw) - 0.3, 2.0), bar_h - 0.4))
        axt.axis("off")
        mis_txt = ", ".join(pretty(names[i]) for i in iMIS) or "–"
        part_txt = "\n".join(tx(f"{pretty(names[ip])} (bereinigt): {dez(clean_frac[ip], 2)} %",
                                f"{pretty(names[ip])} (cleaned): {dez(clean_frac[ip], 2)} %") for ip in iPARTS)
        n_scr = len(scratch_lines) + n_add
        axt.text(0, 1, tx("Bereinigung\n", "Cleaning\n")
                 + (tx(f"• Kratzer ({n_scr}) → aus Nachbarn ergänzt\n", f"• scratches ({n_scr}) → filled from neighbours\n")
                    if scratch.any() else "")
                 + tx(f"• Fehlindizierungen ({mis_txt}) → nicht indiziert\n", f"• misindexing ({mis_txt}) → non-indexed\n")
                 + tx(f"• Partikelcluster < {ALPHA_MIN_PX} Punkte → nicht indiziert\n",
                      f"• particle clusters < {ALPHA_MIN_PX} points → non-indexed\n")
                 + tx(f"• Nulllösungs-Cluster: Opening 3×3,\n   ≥ {MG2SI_MIN_PX} Punkte, mittl. BC < {dez(S.mg2si_max_bc, 0)}\n",
                      f"• zero-solution clusters: opening 3×3,\n   ≥ {MG2SI_MIN_PX} points, mean BC < {dez(S.mg2si_max_bc, 0)}\n")
                 + tx("• übrige Punkte: Phase per Nachbar-Mehrheit,\n   Orientierung vom nächsten Messpunkt\n\n",
                      "• remaining points: phase by neighbour majority,\n   orientation from the nearest measured point\n\n")
                 + (part_txt + "\n" if part_txt else "")
                 + tx(f"{CLUSTER_LABEL}: {dez(clean_frac['Cluster'], 2)} %  ({n_mg} Cluster)",
                      f"{CLUSTER_LABEL}: {dez(clean_frac['Cluster'], 2)} %  ({n_mg} clusters)"),
                 va="top", fontsize=10.5, transform=axt.transAxes)
        t = fig.text(0.3 / FW, 1 - 0.4 / FH, T_PRE + tx("Phasenverteilung", "Phase distribution"), fontsize=16, weight="bold", va="center")
        fig.text(0.3 / FW, 1 - 0.75 / FH, tx(f"Bezug: gesamte Messfläche {dez(W_UM, 0)} µm × {dez(H_UM, 0)} µm",
                                             f"Reference: total measured area {dez(W_UM, 0)} µm × {dez(H_UM, 0)} µm"), fontsize=11,
                 va="center", color="0.25")
        fit_text(fig, t, FW - 0.6)
        return fig

    builders["08"] = fig_08

    # --- 09/10 Rekristallisation
    if S.rx:
        RX_KEYS = ["rx", "sub", "def"] if SUBG else ["rx", "def"]
        RX_CODE = {"rx": RX, "sub": SUB, "def": DEF}
        RX_SHORT = {"rx": tx("rekrist.", "recryst."), "sub": tx("substrukt.", "substruct."),
                    "def": tx("verformt", "deformed") if SUBG else tx("nicht rekr.", "non-recr."), "nc": tx("n. klass.", "n. class.")}

        def fig_09():
            rgb = np.ones((NY, NX, 3))
            for k_, c_ in COL_CLS.items():
                rgb[CL == k_] = c_
            rgb = paint_particles(rgb)
            keys = RX_KEYS + ["nc"]

            def rx_bars(ax):
                vals = [rx_frac[k] for k in keys]
                bars = ax.bar(range(len(keys)), vals, color=[COL_CLS[RX_CODE.get(k, NC)] for k in keys], ec="k", lw=0.5)
                for b_, v_ in zip(bars, vals):
                    ax.text(b_.get_x() + b_.get_width() / 2, v_ + max(vals) * 0.02, dez(v_, 1), ha="center", fontsize=9.5)
                ax.set_xticks(range(len(keys)), [RX_SHORT[k] for k in keys], rotation=30, fontsize=9.5)
                ax.set_ylabel(tx(f"Anteil {plain(PH_MATRIX)} / %", f"Fraction {plain(PH_MATRIX)} / %"))
                ax.set_ylim(0, max(vals) * 1.18 + 0.1)
                ax.yaxis.set_major_formatter(comma)

            if SUBG:
                expl = [tx(f"rekristallisiert: Korn (HAGB ≥ {dez(S.hagb, 0)}°)", f"recrystallized: grain (HAGB ≥ {dez(S.hagb, 0)}°)"),
                        tx("   mit GOS ≤ θc", "   with GOS ≤ θc"),
                        tx("substrukturiert / verformt:", "substructured / deformed:"),
                        tx("   Subkorn (≥ θc) mit GOS ≤ / > θc", "   subgrain (≥ θc) with GOS ≤ / > θc"),
                        tx(f"θc = {dez(S.rx_theta_c, 1)}°, min. {S.rx_min_px} indizierte Punkte",
                           f"θc = {dez(S.rx_theta_c, 1)}°, min. {S.rx_min_px} indexed points")]
                leg_cls = [Patch(fc=COL_CLS[RX], label=tx(f"rekristallisiert (GOS ≤ {dez(S.rx_theta_c, 1)}°)",
                                                          f"recrystallized (GOS ≤ {dez(S.rx_theta_c, 1)}°)")),
                           Patch(fc=COL_CLS[SUB], label=tx("substrukturiert", "substructured")),
                           Patch(fc=COL_CLS[DEF], label=tx("verformt", "deformed"))]
            else:
                expl = [tx(f"Korn (HAGB ≥ {dez(S.hagb, 0)}°) mit", f"grain (HAGB ≥ {dez(S.hagb, 0)}°) with"),
                        tx("   GOS ≤ θc: rekristallisiert", "   GOS ≤ θc: recrystallized"),
                        tx("   GOS > θc: nicht rekristallisiert", "   GOS > θc: non-recrystallized"),
                        tx("ohne Subkornerkennung: keine Trennung", "without subgrain detection: no separation"),
                        tx("   substrukturiert / verformt", "   substructured / deformed"),
                        tx(f"θc = {dez(S.rx_theta_c, 1)}°, min. {S.rx_min_px} indizierte Punkte",
                           f"θc = {dez(S.rx_theta_c, 1)}°, min. {S.rx_min_px} indexed points")]
                leg_cls = [Patch(fc=COL_CLS[RX], label=tx(f"rekristallisiert (GOS ≤ {dez(S.rx_theta_c, 1)}°)",
                                                          f"recrystallized (GOS ≤ {dez(S.rx_theta_c, 1)}°)")),
                           Patch(fc=COL_CLS[DEF], label=tx(f"nicht rekristallisiert (GOS > {dez(S.rx_theta_c, 1)}°)",
                                                           f"non-recrystallized (GOS > {dez(S.rx_theta_c, 1)}°)"))]
            leg_cls.append(Patch(fc=COL_CLS[NC], label=tx(f"nicht klassifizierbar (< {S.rx_min_px} Punkte)",
                                                          f"not classifiable (< {S.rx_min_px} points)")))
            blocks = [text_block(expl, head=tx("Klassifizierung\n(in Anlehnung an Channel 5)", "Classification\n(following Channel 5)")
                                 if SUBG else tx("Klassifizierung\n(Kornebene)", "Classification\n(grain level)"), fs=10),
                      axes_block(2.2, rx_bars, pad=(0.62, 0.75, 0.12, 0.15)),
                      legend_block(leg_cls + b_part + b_gb[:1]),
                      glyph_block(1.3)]
            return lay.figure(rgb, T_PRE + tx(f"Rekristallisierter Anteil {dez(rx_frac['rx'], 1)} % ({plain(PH_MATRIX)})",
                                              f"Recrystallized fraction {dez(rx_frac['rx'], 1)} % ({plain(PH_MATRIX)})"), sub_info, blocks,
                              boundaries=[(gbH, "k", LW_HAGB), (gbHi, "k", LW_HAGB)])

        def fig_10():
            fig, axs = plt.subplots(2, 3, figsize=(20, 10.5), dpi=100)
            lbl_rxfrac = tx("rekristallisierter Anteil / %", "recrystallized fraction / %")
            ax = axs[0, 0]
            b_gos = np.arange(0, 10.01, 0.25)
            hg, _ = np.histogram(gos_g[valid_g], bins=b_gos, weights=area_g[valid_g])
            hg = hg / A_al * 100
            ax.bar(b_gos[:-1], hg, width=0.25, align="edge", ec="white", lw=0.3,
                   color=[COL_CLS[RX] if b_ < S.rx_theta_c else "0.55" for b_ in b_gos[:-1]])
            ax.axvline(S.rx_theta_c, color="k", ls="--", lw=1)
            ax.set_xlabel("Grain Orientation Spread (GOS) / °")
            ax.set_ylabel(tx(f"Flächenanteil {plain(PH_MATRIX)} / %", f"Area fraction {plain(PH_MATRIX)} / %"))
            ax.set_title(tx("(a) GOS-Verteilung, flächengewichtet", "(a) GOS distribution, area-weighted"))
            g10 = dez(area_g[valid_g & (gos_g > 10)].sum() / A_al * 100, 1)
            ax.text(0.97, 0.95, tx(f"GOS > 10°: {g10} % der Fläche\n(nicht dargestellt)", f"GOS > 10°: {g10} % of the area\n(not shown)"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=10)
            ax = axs[0, 1]
            ax.plot(theta_scan, rx_scan, "o-", color=COL_CLS[RX], ms=4)
            ax.axvline(S.rx_theta_c, color="k", ls="--", lw=1)
            for th_ in sorted({1.0, S.rx_theta_c}):
                v_ = area_g[valid_g & (gos_g <= th_)].sum() / A_al * 100
                ax.annotate(f"{dez(v_, 1)} % (≤ {dez(th_, 1)}°)", (th_, v_), xytext=(10, -18), textcoords="offset points", fontsize=10)
            ax.set_xlabel(tx("GOS-Grenzwert θc / °", "GOS threshold θc / °"))
            ax.set_ylabel(lbl_rxfrac)
            ax.set_title(tx("(b) Sensitivität gegenüber dem Grenzwert", "(b) Sensitivity to the threshold"))
            ax.set_xlim(0, 5.2)
            ax.set_ylim(0, None)
            ax.grid(alpha=0.3)
            ax = axs[0, 2]
            ax.set_title(tx("(c) Gegenprüfung auf Subkornebene (CDRX-Übergang)", "(c) Cross-check at subgrain level (CDRX transition)"))
            if SUBG:
                ax.plot(fh_scan * 100, rx_sub_scan, "s-", color="0.25", ms=4,
                        label=tx(f"Subkörner mit GOS ≤ {dez(S.rx_theta_c, 1)}°", f"subgrains with GOS ≤ {dez(S.rx_theta_c, 1)}°"))
                ax.axhline(rx_frac["rx"], color=COL_CLS[RX], ls="--", lw=1.5, label=tx("Kornebene (Hauptergebnis)", "grain level (main result)"))
                for f_ in (0.5, 0.75, 0.9):
                    ax.annotate(f"{dez(rx_sub[f_], 1)} %", (f_ * 100, rx_sub[f_]), xytext=(6, 6), textcoords="offset points", fontsize=10)
                ax.set_xlabel(tx("min. HAGB-Anteil am Subkornumfang / %", "min. HAGB fraction of the subgrain perimeter / %"))
                ax.set_ylabel(lbl_rxfrac)
                ax.set_ylim(0, None)
                ax.grid(alpha=0.3)
                ax.legend(frameon=False, fontsize=10)
                ax.xaxis.set_major_formatter(comma)
                ax.yaxis.set_major_formatter(comma)
            else:
                ax.axis("off")
                ax.text(0.5, 0.5, tx("Subkornerkennung ausgeschaltet –\nkeine Gegenprüfung auf Subkornebene",
                                     "Subgrain detection switched off –\nno cross-check at subgrain level"),
                        transform=ax.transAxes, ha="center", va="center", fontsize=12, color="0.4")
            ax = axs[1, 0]
            if n_rx:
                b_ecd = np.arange(0, np.ceil(np.percentile(ecd_g[rxg], 99)) + 0.5, 0.5)
                ax.hist(ecd_g[rxg], bins=b_ecd, color=COL_CLS[RX], ec="white", lw=0.3)
                ax.axvline(rx_ecd_mean, color="k", ls="--", lw=1)
            ax.set_xlabel(tx("Äquivalentkreisdurchmesser / µm", "Equivalent circle diameter / µm"))
            ax.set_ylabel(tx("Anzahl Körner", "Number of grains"))
            ax.set_title(tx("(d) Korngröße der rekristallisierten Körner", "(d) Size of the recrystallized grains"))
            ax.text(0.97, 0.95, tx(f"n = {n_rx}\nMittelwert {dez(rx_ecd_mean, 2)} µm\nMedian {dez(rx_ecd_med, 2)} µm\n"
                                   f"Streckungsgrad (Median) {dez(rx_ar_med, 2)}\nAnteil mit Streckungsgrad < 2: {dez(rx_equiax, 0)} %",
                                   f"n = {n_rx}\nmean {dez(rx_ecd_mean, 2)} µm\nmedian {dez(rx_ecd_med, 2)} µm\n"
                                   f"aspect ratio (median) {dez(rx_ar_med, 2)}\nfraction with aspect ratio < 2: {dez(rx_equiax, 0)} %"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=10.5)
            for a_ in (axs[0, 0], axs[0, 1], axs[1, 0]):
                a_.xaxis.set_major_formatter(comma)
                a_.yaxis.set_major_formatter(comma)
            ax = axs[1, 1]
            rv = np.nan_to_num(rx_vs_dist)
            bars = ax.bar(range(len(dist_lbl)), rv, color=COL_CLS[RX], ec="k", lw=0.5)
            ax.axhline(rx_frac["rx"] / max(100 - rx_frac["nc"], 1e-9) * 100, color="k", ls="--", lw=1,
                       label=tx("Mittelwert Gesamtfläche", "mean of total area"))
            for b_, v_, a_ in zip(bars, rx_vs_dist, area_vs_dist):
                ax.text(b_.get_x() + b_.get_width() / 2, (0 if not np.isfinite(v_) else v_) + 0.3, f"{dez(v_, 1)} %", ha="center", fontsize=10)
                ax.text(b_.get_x() + b_.get_width() / 2, 0.3, f"{dez(a_, 0)} %\n" + tx("Fläche", "area"), ha="center", fontsize=8.5, color="white")
            ax.set_xticks(range(len(dist_lbl)), dist_lbl)
            ax.set_xlabel(tx(f"Abstand zum nächsten Partikel ≥ {dez(S.psn_min_ecd, 1)} µm ({psn_names}) / µm",
                             f"Distance to the nearest particle ≥ {dez(S.psn_min_ecd, 1)} µm ({psn_names}) / µm"))
            ax.set_ylabel(lbl_rxfrac)
            ax.set_title(tx("(e) Rekristallisation in Partikelnähe", "(e) Recrystallization near particles"))
            ax.yaxis.set_major_formatter(comma)
            ax.legend(frameon=False, loc="upper right")
            ax.set_ylim(0, max(np.nanmax(rv) if rv.size else 1, 1) * 1.25)
            ax = axs[1, 2]
            kk = RX_KEYS
            bars = ax.bar(range(len(kk)), [np.nan_to_num(fib_cls[k]) for k in kk], color=[COL_CLS[RX_CODE[k]] for k in kk], ec="k", lw=0.5)
            for b_, k in zip(bars, kk):
                ax.text(b_.get_x() + b_.get_width() / 2, np.nan_to_num(fib_cls[k]) + 1, f"{dez(fib_cls[k], 1)} %\nKAM {dez(kam_cls[k], 2)}°",
                        ha="center", fontsize=10)
            ax.axhline(20.4, color="k", ls=":", lw=1.2, label=tx("regellos (≈ 20 %)", "random (≈ 20 %)"))
            ax.set_xticks(range(len(kk)), [CLS_NAME[RX_CODE[k]] for k in kk])
            ax.set_ylabel(tx("Flächenanteil mit <110> ‖ CD (± 15°) / %", "Area fraction with <110> ‖ CD (± 15°) / %"))
            ax.set_ylim(0, 100)
            ax.set_title(tx("(f) Faseranteil <110> ‖ CD und mittlere KAM je Klasse", "(f) Fibre fraction <110> ‖ CD and mean KAM per class"))
            ax.legend(frameon=False, loc="upper left")
            ax.yaxis.set_major_formatter(comma)
            fig.suptitle(T_PRE + tx(f"Rekristallisationsanalyse (θc = {dez(S.rx_theta_c, 1)}°)",
                                    f"Recrystallization analysis (θc = {dez(S.rx_theta_c, 1)}°)")
                         + ("" if SUBG else tx(" – ohne Subkornerkennung", " – without subgrain detection")), fontsize=15, weight="bold")
            fig.tight_layout()
            return fig

        builders["09"], builders["10"] = fig_09, fig_10

    # --- 11–14 Korngröße und Kornform
    if S.grainsize:
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import NullFormatter
        COL_G, COL_S = "#1F4E9C", "#E08A1E"
        COL_EXCL = (0.86, 0.86, 0.86)
        b_excl = Patch(fc=COL_EXCL, label=tx(f"Randkorn / < {S.gs_min_px} Punkte (nicht gewertet)",
                                             f"edge grain / < {S.gs_min_px} points (excluded)"))
        lbl_ecd = tx("Äquivalentkreisdurchmesser ECD / µm", "Equivalent circle diameter ECD / µm")
        lbl_num = tx("Anzahlanteil / %", "Number fraction / %")
        lbl_area = tx("Flächenanteil / %", "Area fraction / %")

        def log_ticks(lo, hi):
            cand = [0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
            t = [c for c in cand if lo * 0.999 <= c <= hi * 1.001]
            return t if len(t) >= 2 else [lo, hi]

        def ecd_norm(d):
            e = d["ecd"][d["valid"]]
            if not e.size:
                return LogNorm(STEP, 10 * STEP)
            lo, hi = np.percentile(e, 2), np.percentile(e, 98)
            return LogNorm(max(lo, STEP), max(hi, 1.5 * max(lo, STEP)))

        def grain_rgb(Gp, d, val, cmap, norm):
            rgb = np.ones((NY, NX, 3))
            m = m_alg & (Gp >= 0)
            col = np.asarray(cmap(norm(np.clip(val, norm.vmin, norm.vmax))))[:, :3]
            col[~d["valid"]] = COL_EXCL
            rgb[m] = col[Gp[m]]
            return paint_particles(rgb)

        def log_axis(ax, lo, hi):
            ax.set_xscale("log")
            ax.set_xlim(lo, hi)
            ax.set_xticks(log_ticks(lo, hi))
            ax.xaxis.set_major_formatter(comma)
            ax.xaxis.set_minor_formatter(NullFormatter())

        def stat_lines(sm, sub):
            kind = tx("Subkörner", "subgrains") if sub else tx("Körner", "grains")
            return [tx(f"gewertet: {thou(sm['n'])} {kind}", f"evaluated: {thou(sm['n'])} {kind}"),
                    tx(f"   ({sm['n_edge']} Rand, {sm['n_small']} < {S.gs_min_px} Punkte ausgeschl.)",
                       f"   ({sm['n_edge']} edge, {sm['n_small']} < {S.gs_min_px} points excluded)"),
                    tx(f"   Flächenanteil gewertet: {dez(sm['area_frac'], 1)} %", f"   area fraction evaluated: {dez(sm['area_frac'], 1)} %"),
                    tx(f"ECD Mittelwert: {dez(sm['ecd_mean'], 2)} ± {dez(sm['ecd_sd'], 2)} µm",
                       f"ECD mean: {dez(sm['ecd_mean'], 2)} ± {dez(sm['ecd_sd'], 2)} µm"),
                    tx(f"ECD Median: {dez(sm['ecd_med'], 2)} µm", f"ECD median: {dez(sm['ecd_med'], 2)} µm"),
                    f"ECD d10 / d90: {dez(sm['d10'], 2)} / {dez(sm['d90'], 2)} µm",
                    tx(f"ECD flächengewichtet: {dez(sm['ecd_aw'], 2)} µm", f"ECD area-weighted: {dez(sm['ecd_aw'], 2)} µm")]

        def ecd_hist_block(d, cmap, norm):
            def fn(ax):
                e = d["ecd"][d["valid"]]
                if e.size:
                    bins = np.logspace(np.log10(e.min() * 0.95), np.log10(e.max() * 1.05), 26)
                    h, be = np.histogram(e, bins=bins)
                    cen = np.sqrt(be[:-1] * be[1:])
                    ax.bar(be[:-1], h / e.size * 100, width=np.diff(be), align="edge", ec="white", lw=0.3,
                           color=cmap(norm(np.clip(cen, norm.vmin, norm.vmax))))
                    ax.axvline(e.mean(), color="k", ls="--", lw=1)
                    log_axis(ax, be[0], be[-1])
                ax.set_xlabel("ECD / µm")
                ax.set_ylabel(lbl_num)
                ax.yaxis.set_major_formatter(comma)
            return axes_block(2.1, fn)

        def make_size_map(d, sm, Gp, sub):
            def f():
                cmap = plt.get_cmap("viridis")
                norm = ecd_norm(d)
                rgb = grain_rgb(Gp, d, d["ecd"], cmap, norm)
                if sub:
                    head = tx(f"Subkorngröße\n(Grenzen ≥ {dez(S.lagb, 1)}°)", f"Subgrain size\n(boundaries ≥ {dez(S.lagb, 1)}°)")
                    bnd = [(gbL1, COL_LAGB, LW_LAGB), (gbL2, COL_LAGB, LW_LAGB), (gbH, "k", LW_HAGB), (gbHi, "k", LW_HAGB)]
                    leg = [b_excl] + b_part + b_gb
                    ttl = T_PRE + tx(f"Subkorngröße (ECD, Subkörner ≥ {dez(S.lagb, 1)}°): Mittelwert {dez(sm['ecd_mean'], 2)} µm",
                                     f"Subgrain size (ECD, subgrains ≥ {dez(S.lagb, 1)}°): mean {dez(sm['ecd_mean'], 2)} µm")
                    extra = [tx(f"Subkörner je Korn (Mittel): {dez(spg_mean, 1)}", f"subgrains per grain (mean): {dez(spg_mean, 1)}"),
                             tx(f"mittl. Sehnenlänge (≥ {dez(S.lagb, 0)}°): {dez(mli2['ges'], 2)} µm",
                                f"mean linear intercept (≥ {dez(S.lagb, 0)}°): {dez(mli2['ges'], 2)} µm")]
                else:
                    head = tx(f"Korngröße\n(HAGB ≥ {dez(S.hagb, 0)}°)", f"Grain size\n(HAGB ≥ {dez(S.hagb, 0)}°)")
                    bnd = [(gbH, "k", LW_HAGB), (gbHi, "k", LW_HAGB)]
                    leg = [b_excl] + b_part + b_gb[:1]
                    ttl = T_PRE + tx(f"Korngröße (ECD, HAGB ≥ {dez(S.hagb, 0)}°): Mittelwert {dez(sm['ecd_mean'], 2)} µm",
                                     f"Grain size (ECD, HAGB ≥ {dez(S.hagb, 0)}°): mean {dez(sm['ecd_mean'], 2)} µm")
                    extra = [tx(f"mittl. Sehnenlänge (HAGB): {dez(mli['ges'], 2)} µm", f"mean linear intercept (HAGB): {dez(mli['ges'], 2)} µm"),
                             tx(f"Korngrößenkennzahl G (ASTM E112): {dez(astm_G, 1)}", f"grain size number G (ASTM E112): {dez(astm_G, 1)}")]
                blocks = [text_block(stat_lines(sm, sub) + extra, head=head, fs=10),
                          cbar_block(cmap, norm, tx("Äquivalentkreisdurchmesser ECD / µm", "Equivalent circle diameter ECD / µm"),
                                     ticks=log_ticks(norm.vmin, norm.vmax), extend="both"),
                          ecd_hist_block(d, cmap, norm), legend_block(leg), glyph_block(1.3)]
                return lay.figure(rgb, ttl, sub_info, blocks, boundaries=bnd)
            return f

        def fig_13():
            cmap = plt.get_cmap("plasma")
            norm = Normalize(1, 4)
            rgb = grain_rgb(GH, gsH, gsH["ar"], cmap, norm)

            def ar_hist(ax):
                ar = gsH["ar"][gsH["valid"]]
                bins = np.arange(1, 6.01, 0.2)
                if ar.size:
                    h, be = np.histogram(np.clip(ar, 1, 5.99), bins=bins)
                    ax.bar(be[:-1], h / ar.size * 100, width=0.2, align="edge", ec="white", lw=0.3,
                           color=cmap(norm(np.clip(be[:-1] + 0.1, 1, 4))))
                    ax.axvline(np.median(ar), color="k", ls="--", lw=1)
                ax.set_xlim(1, 6)
                ax.set_xlabel(tx("Aspektverhältnis (≥ 6 im letzten Balken)", "Aspect ratio (≥ 6 in the last bar)"))
                ax.set_ylabel(lbl_num)
                ax.xaxis.set_major_formatter(comma)
                ax.yaxis.set_major_formatter(comma)

            ratio = dez(mli['RD'] / max(mli['CD'], 1e-9), 2)
            lines = [tx(f"gewertet: {thou(sumH['n'])} Körner", f"evaluated: {thou(sumH['n'])} grains"),
                     tx(f"Mittelwert: {dez(sumH['ar_mean'], 2)}, Median: {dez(sumH['ar_med'], 2)}",
                        f"mean: {dez(sumH['ar_mean'], 2)}, median: {dez(sumH['ar_med'], 2)}"),
                     tx(f"flächengewichtet: {dez(sumH['ar_aw'], 2)}", f"area-weighted: {dez(sumH['ar_aw'], 2)}"),
                     tx(f"Anteil < 2 (gleichachsig): {dez(sumH['ar_lt2'], 1)} %", f"fraction < 2 (equiaxed): {dez(sumH['ar_lt2'], 1)} %"),
                     tx("Hauptachse (AV ≥ 1,5, flächengew.):", "major axis (AR ≥ 1.5, area-weighted):"),
                     tx(f"   {dez(sumH['ang_mean'], 0)}° zu RD, Ausrichtung R = {dez(sumH['ang_R'], 2)}",
                        f"   {dez(sumH['ang_mean'], 0)}° to RD, alignment R = {dez(sumH['ang_R'], 2)}"),
                     tx(f"   näher an RD als an CD: {dez(sumH['rd_frac'], 0)} %", f"   closer to RD than to CD: {dez(sumH['rd_frac'], 0)} %"),
                     tx(f"Sehnenlängen RD/CD (HAGB): {ratio}", f"intercept ratio RD/CD (HAGB): {ratio}")]
            blocks = [text_block(lines, head=tx("Aspektverhältnis\n(Ellipse gleicher Flächenmomente)",
                                                "Aspect ratio\n(ellipse of equal area moments)"), fs=10),
                      cbar_block(cmap, norm, tx("Aspektverhältnis Haupt-/Nebenachse", "Aspect ratio major/minor axis"),
                                 ticks=[1, 1.5, 2, 2.5, 3, 3.5, 4]),
                      axes_block(2.1, ar_hist), legend_block([b_excl] + b_part + b_gb[:1]), glyph_block(1.3)]
            return lay.figure(rgb, T_PRE + tx(f"Kornaspektverhältnis (HAGB ≥ {dez(S.hagb, 0)}°): Median {dez(sumH['ar_med'], 2)}",
                                              f"Grain aspect ratio (HAGB ≥ {dez(S.hagb, 0)}°): median {dez(sumH['ar_med'], 2)}"),
                              sub_info, blocks, boundaries=[(gbH, "k", LW_HAGB), (gbHi, "k", LW_HAGB)])

        def fig_14():
            fig, axs = plt.subplots(2, 3, figsize=(20, 10.5), dpi=100)
            sets = [(gsH, sumH, COL_G, tx(f"Körner (HAGB ≥ {dez(S.hagb, 0)}°)", f"grains (HAGB ≥ {dez(S.hagb, 0)}°)"))]
            if SUBG:
                sets.append((gsL, sumL, COL_S, tx(f"Subkörner (≥ {dez(S.lagb, 1)}°)", f"subgrains (≥ {dez(S.lagb, 1)}°)")))
            e_all = np.concatenate([d["ecd"][d["valid"]] for d, *_ in sets])
            if e_all.size:
                lo, hi = e_all.min() * 0.95, e_all.max() * 1.05
            else:
                lo, hi = STEP, 10 * STEP
            bins = np.logspace(np.log10(lo), np.log10(hi), 36)
            # (a) Anzahlverteilung
            ax = axs[0, 0]
            for d, sm, c, lab in sets:
                e = d["ecd"][d["valid"]]
                if e.size:
                    h, _ = np.histogram(e, bins=bins)
                    ax.stairs(h / e.size * 100, bins, color=c, lw=2, label=f"{lab}, n = {thou(sm['n'])}")
                    ax.axvline(sm["ecd_mean"], color=c, ls="--", lw=1.2)
            log_axis(ax, lo, hi)
            ax.set_xlabel(lbl_ecd)
            ax.set_ylabel(lbl_num)
            ax.set_title(tx("(a) Größenverteilung, anzahlgewichtet (-- Mittelwert)", "(a) Size distribution, number-weighted (-- mean)"))
            ax.legend(frameon=False, fontsize=10)
            # (b) Flächenverteilung
            ax = axs[0, 1]
            for d, sm, c, lab in sets:
                v = d["valid"]
                if v.any():
                    h, _ = np.histogram(d["ecd"][v], bins=bins, weights=d["area"][v])
                    ax.stairs(h / d["area"][v].sum() * 100, bins, color=c, lw=2, label=lab)
                    ax.axvline(sm["ecd_aw"], color=c, ls="--", lw=1.2)
            log_axis(ax, lo, hi)
            ax.set_xlabel(lbl_ecd)
            ax.set_ylabel(lbl_area)
            ax.set_title(tx("(b) Größenverteilung, flächengewichtet (-- Mittelwert)", "(b) Size distribution, area-weighted (-- mean)"))
            ax.legend(frameon=False, fontsize=10)
            # (c) Kennwerttabelle
            ax = axs[0, 2]
            ax.axis("off")
            rows = [("", tx("Körner", "Grains"), tx("Subkörner", "Subgrains")),
                    (tx("gewertet", "evaluated"), thou(sumH['n']), thou(sumL['n'])),
                    (tx("ECD Mittelwert / µm", "ECD mean / µm"), dez(sumH["ecd_mean"], 2), dez(sumL["ecd_mean"], 2)),
                    (tx("ECD Standardabw. / µm", "ECD std. dev. / µm"), dez(sumH["ecd_sd"], 2), dez(sumL["ecd_sd"], 2)),
                    (tx("ECD Median / µm", "ECD median / µm"), dez(sumH["ecd_med"], 2), dez(sumL["ecd_med"], 2)),
                    ("ECD d10 / d90 / µm", f"{dez(sumH['d10'], 2)} / {dez(sumH['d90'], 2)}", f"{dez(sumL['d10'], 2)} / {dez(sumL['d90'], 2)}"),
                    (tx("ECD flächengew. / µm", "ECD area-weighted / µm"), dez(sumH["ecd_aw"], 2), dez(sumL["ecd_aw"], 2)),
                    (tx("Sehnenlänge ges. / µm", "Mean intercept total / µm"), dez(mli["ges"], 2), dez(mli2["ges"], 2)),
                    (tx("Sehnenlänge RD / µm", "Mean intercept RD / µm"), dez(mli["RD"], 2), dez(mli2["RD"], 2)),
                    (tx("Sehnenlänge CD / µm", "Mean intercept CD / µm"), dez(mli["CD"], 2), dez(mli2["CD"], 2)),
                    (tx("Aspektverh. Median", "Aspect ratio median"), dez(sumH["ar_med"], 2), dez(sumL["ar_med"], 2)),
                    (tx("Aspektverh. Mittel", "Aspect ratio mean"), dez(sumH["ar_mean"], 2), dez(sumL["ar_mean"], 2)),
                    (tx("Anteil AV < 2 / %", "Fraction AR < 2 / %"), dez(sumH["ar_lt2"], 1), dez(sumL["ar_lt2"], 1)),
                    (tx("Kennzahl G (ASTM E112)", "Grain size number G (ASTM E112)"), dez(astm_G, 1), "–"),
                    (tx("Subkörner je Korn", "Subgrains per grain"), "–", dez(spg_mean, 1))]
            if SUBG:
                widths = [0.5, 0.25, 0.25]
            else:
                rows = [r[:2] for r in rows if r[0] != tx("Subkörner je Korn", "Subgrains per grain")]
                widths = [0.62, 0.38]
            tb = ax.table(cellText=[list(r) for r in rows[1:]], colLabels=list(rows[0]), cellLoc="right", colLoc="right",
                          colWidths=widths, bbox=[0, 0.07, 1, 0.93])
            tb.auto_set_font_size(False)
            tb.set_fontsize(11)
            for (r_, c_), cell in tb.get_celld().items():
                cell.set_edgecolor("0.8")
                if c_ == 0:
                    cell.set_text_props(ha="left")
                    cell._loc = "left"
                if r_ == 0:
                    cell.set_text_props(weight="bold")
            ax.set_title(tx(f"(c) Kennwerte (ohne Randkörner, ≥ {S.gs_min_px} Punkte)",
                            f"(c) Key figures (without edge grains, ≥ {S.gs_min_px} points)"))
            ax.text(0.0, 0.0, tx("Sehnenlängen: Linienschnitt in RD und CD; G aus der mittleren Sehnenlänge der HAGB",
                                 "Intercepts: line intercept in RD and CD; G from the mean linear intercept of the HAGB")
                    + ("" if SUBG else tx("; ohne Subkornerkennung", "; without subgrain detection")),
                    transform=ax.transAxes, fontsize=10, va="bottom", color="0.3")
            # (d) Aspektverhältnis
            ax = axs[1, 0]
            b_ar = np.arange(1, 6.01, 0.2)
            for d, sm, c, lab in sets:
                ar = d["ar"][d["valid"]]
                if ar.size:
                    h, _ = np.histogram(np.clip(ar, 1, 5.99), bins=b_ar)
                    ax.stairs(h / ar.size * 100, b_ar, color=c, lw=2, label=f"{lab}, " + tx("Median", "median") + f" {dez(sm['ar_med'], 2)}")
                    ax.axvline(sm["ar_med"], color=c, ls="--", lw=1.2)
            ax.set_xlim(1, 6)
            ax.set_xlabel(tx("Aspektverhältnis Haupt-/Nebenachse (≥ 6 im letzten Intervall)",
                             "Aspect ratio major/minor axis (≥ 6 in the last interval)"))
            ax.set_ylabel(lbl_num)
            ax.set_title(tx("(d) Aspektverhältnis (-- Median)", "(d) Aspect ratio (-- median)"))
            ax.legend(frameon=False, fontsize=10)
            # (e) Aspektverhältnis über Korngröße
            ax = axs[1, 1]
            v = gsH["valid"]
            ax.scatter(gsH["ecd"][v], gsH["ar"][v], s=7, alpha=0.4, color=COL_G, lw=0, rasterized=True)
            for d, sm, c, lab in sets:
                vv = d["valid"]
                if vv.sum() > 20:
                    e, ar = d["ecd"][vv], d["ar"][vv]
                    eb = np.logspace(np.log10(e.min()), np.log10(e.max() * 1.0001), 13)
                    k = np.digitize(e, eb) - 1
                    med = [np.median(ar[k == j]) if (k == j).sum() >= 5 else np.nan for j in range(len(eb) - 1)]
                    ax.plot(np.sqrt(eb[:-1] * eb[1:]), med, "o-", color=c, ms=4, lw=1.8,
                            label=tx(f"Median je Größenklasse, {lab}", f"median per size class, {lab}"))
            log_axis(ax, lo, hi)
            ax.set_ylim(1, 6)
            ax.set_xlabel(lbl_ecd)
            ax.set_ylabel(tx("Aspektverhältnis", "Aspect ratio"))
            ax.set_title(tx("(e) Aspektverhältnis über Korngröße (Punkte: Körner)", "(e) Aspect ratio vs. grain size (dots: grains)"))
            if ax.get_legend_handles_labels()[0]:
                ax.legend(frameon=False, fontsize=9.5, loc="upper left")
            # (f) Ausrichtung der Hauptachse
            ax = axs[1, 2]
            b_ang = np.arange(-90, 90.1, 10)
            for d, sm, c, lab in sets:
                vv = d["valid"] & (d["ar"] >= 1.5)
                if vv.any():
                    h, _ = np.histogram(d["ang"][vv], bins=b_ang, weights=d["area"][vv])
                    ax.stairs(h / d["area"][vv].sum() * 100, b_ang, color=c, lw=2,
                              label=f"{lab}: R = {dez(sm['ang_R'], 2)}")
            ax.axhline(100 / 18, color="0.5", ls=":", lw=1.2, label=tx("regellos", "random"))
            ax.set_xticks([-90, -45, 0, 45, 90], ["−90°\n(CD)", "−45°", "0°\n(RD)", "45°", "90°\n(CD)"])
            ax.set_xlim(-90, 90)
            ax.set_ylim(0, None)
            ax.set_xlabel(tx("Winkel der Hauptachse zu RD (gegen den Uhrzeigersinn)", "Angle of the major axis to RD (counter-clockwise)"))
            ax.set_ylabel(lbl_area)
            ax.set_title(tx("(f) Ausrichtung der Hauptachse (Aspektverhältnis ≥ 1,5)", "(f) Orientation of the major axis (aspect ratio ≥ 1.5)"))
            ax.legend(frameon=False, fontsize=10)
            for a_ in axs.ravel():
                if a_ is not axs[0, 2]:
                    a_.yaxis.set_major_formatter(comma)
            for a_ in (axs[1, 0],):
                a_.xaxis.set_major_formatter(comma)
            if SUBG:
                sup = tx(f"Korngröße und Kornform ({plain(PH_MATRIX)}; Körner HAGB ≥ {dez(S.hagb, 0)}°, Subkörner ≥ {dez(S.lagb, 1)}°)",
                         f"Grain size and shape ({plain(PH_MATRIX)}; grains HAGB ≥ {dez(S.hagb, 0)}°, subgrains ≥ {dez(S.lagb, 1)}°)")
            else:
                sup = tx(f"Korngröße und Kornform ({plain(PH_MATRIX)}; Körner HAGB ≥ {dez(S.hagb, 0)}°, ohne Subkornerkennung)",
                         f"Grain size and shape ({plain(PH_MATRIX)}; grains HAGB ≥ {dez(S.hagb, 0)}°, without subgrain detection)")
            fig.suptitle(T_PRE + sup, fontsize=15, weight="bold")
            fig.tight_layout()
            return fig

        builders["11"] = make_size_map(gsH, sumH, GH, False)
        if SUBG:
            builders["12"] = make_size_map(gsL, sumL, GL, True)
        builders["13"], builders["14"] = fig_13, fig_14

    # ======================================================================== Kennwerte
    import mtex as _mtex
    on_off = [(tx("Kratzererkennung", "scratch detection"), S.scratch), (tx("Mg2Si-Erkennung", "Mg2Si detection"), S.mg2si),
              (tx("IPF-Karten", "IPF maps"), S.ipf), ("KAM", S.kam), (tx("Korngrenzen", "grain boundaries"), S.gb),
              (tx("Textur", "texture"), S.texture), (tx("Rekristallisation", "recrystallization"), S.rx),
              (tx("Korngröße", "grain size"), S.grainsize), (tx("Subkornerkennung", "subgrain detection"), SUBG)]
    K = [tx(f"EBSD-Auswertung {TITLE}  ({APP_NAME} {APP_VERSION}, MTEX für Python {getattr(_mtex, '__version__', '?')})",
            f"EBSD analysis {TITLE}  ({app_title()} {APP_VERSION}, MTEX for Python {getattr(_mtex, '__version__', '?')})"),
         tx("Datei: ", "File: ") + project.path.name
         + (tx("  |  Ausschnitt ", "  |  area ") + f"{tuple(round(v, 1) for v in crop_um)} µm" if crop_um is not None else ""),
         tx(f"Karte: {NX} x {NY} Punkte, Schrittweite {dez(STEP, 4)} µm, Fläche {dez(W_UM, 1)} x {dez(H_UM, 1)} µm",
            f"Map: {NX} x {NY} points, step size {dez(STEP, 4)} µm, area {dez(W_UM, 1)} x {dez(H_UM, 1)} µm"),
         tx("Probenkoordinaten: RD = x (horizontal), CD (Stauchrichtung) = y (vertikal), TD = z (Schliffnormale); MTEX-Darstellung y↓→x",
            "Specimen axes: RD = x (horizontal), CD (compression direction) = y (vertical), TD = z (section normal); MTEX display y↓→x"),
         tx("Phasenrollen: ", "Phase roles: ") + "; ".join(f"{names[i]} = {role_name(roles[i])}" for i in sorted(roles)),
         tx("Module: ", "Modules: ") + ", ".join(n for n, on in on_off if on),
         "", tx("KRATZER", "SCRATCHES")]
    if S.scratch:
        K += [tx(f"automatisch erkannt: {len(scratch_lines)}; Mindestlänge {dez(S.scratch_min_len, 0)} µm, "
                 f"Fenster {dez(S.scratch_window, 0)} µm, Kandidaten-BC < {dez(S.scratch_rel_bc * 100, 0)} % lokal",
                 f"detected automatically: {len(scratch_lines)}; min. length {dez(S.scratch_min_len, 0)} µm, "
                 f"window {dez(S.scratch_window, 0)} µm, candidate BC < {dez(S.scratch_rel_bc * 100, 0)} % of local value")]
        K += [tx(f"Kratzer {k}: Richtung {dez(l['angle'], 1)}° zu RD, Länge {dez(l['length_um'], 1)} µm in {l['n_segs']} Abschnitt(en), "
                 f"Breite {dez(l['width_um'], 1)} µm",
                 f"Scratch {k}: direction {dez(l['angle'], 1)}° to RD, length {dez(l['length_um'], 1)} µm in {l['n_segs']} section(s), "
                 f"width {dez(l['width_um'], 1)} µm") for k, l in enumerate(scratch_lines, 1)]
        if n_dropped:
            K += [tx(f"{n_dropped} automatisch erkannte Kratzer manuell aufgehoben",
                     f"{n_dropped} automatically detected scratch(es) removed manually")]
    else:
        K += [tx("automatische Kratzererkennung ausgeschaltet", "automatic scratch detection switched off")]
    if strokes:
        K += [tx(f"manuelle Markierungen: {n_add} Linie(n) markiert, {n_rem} Linie(n) aufgehoben; "
                 f"zusätzlich markiert {dez(scr_added.mean() * 100, 2)} %, aufgehoben {dez(scr_removed.mean() * 100, 2)} % der Fläche",
                 f"manual marks: {n_add} line(s) marked, {n_rem} line(s) unmarked; "
                 f"additionally marked {dez(scr_added.mean() * 100, 2)} %, unmarked {dez(scr_removed.mean() * 100, 2)} % of the area")]
        K += ["  " + ("+" if s_["modus"] == "add" else "−")
              + f" ({dez(s_['x0'], 1)}; {dez(s_['y0'], 1)}) → ({dez(s_['x1'], 1)}; {dez(s_['y1'], 1)}) µm, "
              + tx("Breite ", "width ") + f"{dez(s_['breite'], 1)} µm" for s_ in strokes]
    if HAS_SCRATCH:
        K += [tx(f"Kratzerfläche gesamt: {dez(scratch.mean() * 100, 2)} %", f"Scratch area in total: {dez(scratch.mean() * 100, 2)} %")]
    K += ["", tx("PHASENANTEILE (Flächenanteil bezogen auf gesamte Messfläche, %)",
                 "PHASE FRACTIONS (area fraction of the total measured area, %)"),
          tx("Phase;Rolle;Rohdaten;bereinigt;MAD-Median (°);BC-Median", "Phase;role;raw data;cleaned;MAD median (°);BC median"),
          f"{PH_MATRIX};{role_name('Matrix')};{dez(raw_frac[iAL], 2)};{dez(clean_frac['Al'], 2)};{dez(phase_stats[iAL][0], 2)};{dez(phase_stats[iAL][1], 0)}"]
    K += [f"{names[ip]};{role_name('Partikel')};{dez(raw_frac[ip], 2)};{dez(clean_frac[ip], 2)};{dez(phase_stats[ip][0], 2)};"
          f"{dez(phase_stats[ip][1], 0)}" for ip in iPARTS]
    K += [f"{names[i]};{role_name('Fehlindizierung')};{dez(raw_frac[i], 2)};0;{dez(phase_stats[i][0], 2)};{dez(phase_stats[i][1], 0)}"
          for i in iMIS]
    K += [tx("Nulllösung", "Zero solution") + f";-;{dez(raw_frac[0], 2)};0;-;-",
          f"{CLUSTER_PLAIN};-;-;{dez(clean_frac['Cluster'], 2)};-;-",
          tx(f"Anzahl Cluster: {n_mg}; Mindestgröße {MG2SI_MIN_PX} Punkte; Partikel-Mindestgröße {ALPHA_MIN_PX} Punkte",
             f"Number of clusters: {n_mg}; min. size {MG2SI_MIN_PX} points; min. particle size {ALPHA_MIN_PX} points")]
    if S.gb:
        K += ["", tx(f"KORNGRENZEN ({PH_MATRIX}, Längenanteil, Linienschnitt RD+CD, Lückenüberbrückung <= {int(S.gap_bridge)} Punkte, "
                     "Winkel mit MTEX)",
                     f"GRAIN BOUNDARIES ({PH_MATRIX}, length fraction, line intercept RD+CD, gap bridging <= {int(S.gap_bridge)} points, "
                     "angles with MTEX)"),
              f"LAGB {dez(S.lagb, 0)}-{dez(S.hagb, 0)}°: {dez(f_lagb, 1)} %;HAGB >={dez(S.hagb, 0)}°: {dez(f_hagb, 1)} %",
              f"{dez(S.lagb, 0)}-5°: {dez(cls[0], 1)} %;5-{dez(S.hagb, 0)}°: {dez(cls[1], 1)} %;>={dez(S.hagb, 0)}°: {dez(cls[2], 1)} %",
              tx(f"Querprüfung HAGB-Anteil: ohne Lückenüberbrückung {dez(hagb_nogap, 1)} %",
                 f"Cross-check HAGB fraction: without gap bridging {dez(hagb_nogap, 1)} %")
              + (tx(f"; MTEX-Grenzsegmente der bereinigten Karte (Längenanteil) {dez(hagb_mtex, 1)} %",
                    f"; MTEX boundary segments of the cleaned map (length fraction) {dez(hagb_mtex, 1)} %") if need_grains and SUBG else ""),
              tx(f"Mittl. Sehnenlänge HAGB: RD {dez(mli['RD'], 2)} µm, CD {dez(mli['CD'], 2)} µm, "
                 f"Verhältnis RD/CD {dez(mli['RD'] / max(mli['CD'], 1e-9), 2)}",
                 f"Mean linear intercept HAGB: RD {dez(mli['RD'], 2)} µm, CD {dez(mli['CD'], 2)} µm, "
                 f"ratio RD/CD {dez(mli['RD'] / max(mli['CD'], 1e-9), 2)}")]
        if SUBG:
            K += [tx(f"Mittl. Sehnenlänge alle Grenzen >={dez(S.lagb, 0)}°: RD {dez(mli2['RD'], 2)} µm, CD {dez(mli2['CD'], 2)} µm",
                     f"Mean linear intercept all boundaries >={dez(S.lagb, 0)}°: RD {dez(mli2['RD'], 2)} µm, CD {dez(mli2['CD'], 2)} µm")]
    if kam is not None and kv.size:
        K += ["", tx(f"KAM (MTEX ebsd.KAM, 3x3-Nachbarschaft, Schwellwert {dez(S.kam_max, 0)}°, nur original indizierte Punkte "
                     f"{PH_MATRIX} ohne Kratzer)",
                     f"KAM (MTEX ebsd.KAM, 3x3 neighbourhood, threshold {dez(S.kam_max, 0)}°, originally indexed points of "
                     f"{PH_MATRIX} only, without scratches)"),
              tx(f"Mittelwert {dez(kv.mean(), 3)}°;Median {dez(np.median(kv), 3)}°;90%-Quantil {dez(np.percentile(kv, 90), 3)}°",
                 f"mean {dez(kv.mean(), 3)}°;median {dez(np.median(kv), 3)}°;90% quantile {dez(np.percentile(kv, 90), 3)}°")]
    if S.texture:
        K += ["", tx(f"TEXTUR (MTEX calcDensity, Kernel-Halbwertsbreite {dez(S.pf_halfwidth, 0)}°, alle Punkte {PH_MATRIX} der "
                     "bereinigten Karte; max. MRD)",
                     f"TEXTURE (MTEX calcDensity, kernel half-width {dez(S.pf_halfwidth, 0)}°, all points of {PH_MATRIX} in the "
                     "cleaned map; max. MRD)"),
              *[f"IPF {k}: {dez(ipf_max[k], 2)}" for k in ["CD", "RD", "TD"]],
              *[f"PF {k}: {dez(v, 2)}" for k, v in pf_max.items()],
              tx(f"Flächenanteil <110> || CD (15°): {dez(fib_all, 1)} % (regellos ca. 20 % bei kubischer Symmetrie)",
                 f"Area fraction <110> || CD (15°): {dez(fib_all, 1)} % (random approx. 20 % for cubic symmetry)")]
    if S.rx:
        K += ["", tx(f"REKRISTALLISATION (Flächenanteil {PH_MATRIX}; Körner HAGB >= {dez(S.hagb, 0)}°, theta_c = {dez(S.rx_theta_c, 1)}°, "
                     f"min. {S.rx_min_px} indizierte Punkte/Korn)",
                     f"RECRYSTALLIZATION (area fraction {PH_MATRIX}; grains HAGB >= {dez(S.hagb, 0)}°, theta_c = {dez(S.rx_theta_c, 1)}°, "
                     f"min. {S.rx_min_px} indexed points/grain)")]
        if SUBG:
            K += [tx(f"Körner (MTEX calcGrains): {nGH} (klassifizierbar {int(valid_g.sum())}), Subkörner (>= {dez(S.rx_theta_c, 1)}°): {nGS}",
                     f"Grains (MTEX calcGrains): {nGH} (classifiable {int(valid_g.sum())}), subgrains (>= {dez(S.rx_theta_c, 1)}°): {nGS}"),
                  tx(f"rekristallisiert (Korn-GOS <= theta_c): {dez(rx_frac['rx'], 2)} %  ({n_rx} Körner)",
                     f"recrystallized (grain GOS <= theta_c): {dez(rx_frac['rx'], 2)} %  ({n_rx} grains)"),
                  tx(f"substrukturiert (Subkorn-GOS <= theta_c): {dez(rx_frac['sub'], 2)} %;verformt (Subkorn-GOS > theta_c): "
                     f"{dez(rx_frac['def'], 2)} %;nicht klassifizierbar: {dez(rx_frac['nc'], 2)} %",
                     f"substructured (subgrain GOS <= theta_c): {dez(rx_frac['sub'], 2)} %;deformed (subgrain GOS > theta_c): "
                     f"{dez(rx_frac['def'], 2)} %;not classifiable: {dez(rx_frac['nc'], 2)} %"),
                  tx(f"zum Vergleich klassisch auf Kornebene (Channel-5-analog): substrukturiert {dez(c5_frac['sub'], 2)} %, "
                     f"verformt {dez(c5_frac['def'], 2)} % (3 größte Körner: {', '.join(dez(v_, 1) for v_ in big3)} % der Fläche)",
                     f"for comparison, classic grain-level classification (Channel 5 analogue): substructured {dez(c5_frac['sub'], 2)} %, "
                     f"deformed {dez(c5_frac['def'], 2)} % (3 largest grains: {', '.join(dez(v_, 1) for v_ in big3)} % of the area)"),
                  tx(f"Gegenprüfung Subkornebene (Subkornfläche mit GOS <= theta_c gesamt {dez(low_sub_frac, 1)} %): ",
                     f"Cross-check at subgrain level (subgrain area with GOS <= theta_c in total {dez(low_sub_frac, 1)} %): ")
                  + "; ".join(tx(f"HAGB-Umfang >= {int(f_ * 100)} %: {dez(v_, 1)} %", f"HAGB perimeter >= {int(f_ * 100)} %: {dez(v_, 1)} %")
                              for f_, v_ in rx_sub.items())]
        else:
            K += [tx(f"Körner (MTEX calcGrains): {nGH} (klassifizierbar {int(valid_g.sum())}); Subkornerkennung ausgeschaltet – "
                     "nur Kornebene, keine Trennung substrukturiert/verformt",
                     f"Grains (MTEX calcGrains): {nGH} (classifiable {int(valid_g.sum())}); subgrain detection switched off – "
                     "grain level only, no separation substructured/deformed"),
                  tx(f"rekristallisiert (Korn-GOS <= theta_c): {dez(rx_frac['rx'], 2)} %  ({n_rx} Körner)",
                     f"recrystallized (grain GOS <= theta_c): {dez(rx_frac['rx'], 2)} %  ({n_rx} grains)"),
                  tx(f"nicht rekristallisiert (Korn-GOS > theta_c): {dez(rx_frac['def'], 2)} %;nicht klassifizierbar: {dez(rx_frac['nc'], 2)} % "
                     f"(3 größte Körner: {', '.join(dez(v_, 1) for v_ in big3)} % der Fläche)",
                     f"non-recrystallized (grain GOS > theta_c): {dez(rx_frac['def'], 2)} %;not classifiable: {dez(rx_frac['nc'], 2)} % "
                     f"(3 largest grains: {', '.join(dez(v_, 1) for v_ in big3)} % of the area)")]
        K += [tx("Sensitivität rekristallisierter Anteil: ", "Sensitivity of the recrystallized fraction: ")
              + "; ".join(f"GOS <= {dez(t_, 2)}°: {dez(v_, 1)} %" for t_, v_ in zip(theta_scan, rx_scan) if round(t_, 2) in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0)),
              tx(f"rekristallisierte Körner: ECD Mittelwert {dez(rx_ecd_mean, 2)} µm, Median {dez(rx_ecd_med, 2)} µm; "
                 f"Streckungsgrad Median {dez(rx_ar_med, 2)}; Anteil Streckungsgrad < 2: {dez(rx_equiax, 0)} %",
                 f"recrystallized grains: ECD mean {dez(rx_ecd_mean, 2)} µm, median {dez(rx_ecd_med, 2)} µm; "
                 f"aspect ratio median {dez(rx_ar_med, 2)}; fraction with aspect ratio < 2: {dez(rx_equiax, 0)} %"),
              tx(f"rekristallisierte Körner am Kartenrand: {int((rxg & touch).sum())}",
                 f"recrystallized grains at the map edge: {int((rxg & touch).sum())}")]
        if SUBG:
            K += [tx(f"mittlere KAM: rekristallisiert {dez(kam_cls['rx'], 3)}°, substrukturiert {dez(kam_cls['sub'], 3)}°, "
                     f"verformt {dez(kam_cls['def'], 3)}°",
                     f"mean KAM: recrystallized {dez(kam_cls['rx'], 3)}°, substructured {dez(kam_cls['sub'], 3)}°, "
                     f"deformed {dez(kam_cls['def'], 3)}°"),
                  tx(f"Flächenanteil <110> || CD (15°): rekristallisiert {dez(fib_cls['rx'], 1)} %, substrukturiert {dez(fib_cls['sub'], 1)} %, "
                     f"verformt {dez(fib_cls['def'], 1)} %",
                     f"Area fraction <110> || CD (15°): recrystallized {dez(fib_cls['rx'], 1)} %, substructured {dez(fib_cls['sub'], 1)} %, "
                     f"deformed {dez(fib_cls['def'], 1)} %")]
        else:
            K += [tx(f"mittlere KAM: rekristallisiert {dez(kam_cls['rx'], 3)}°, nicht rekristallisiert {dez(kam_cls['def'], 3)}°",
                     f"mean KAM: recrystallized {dez(kam_cls['rx'], 3)}°, non-recrystallized {dez(kam_cls['def'], 3)}°"),
                  tx(f"Flächenanteil <110> || CD (15°): rekristallisiert {dez(fib_cls['rx'], 1)} %, nicht rekristallisiert {dez(fib_cls['def'], 1)} %",
                     f"Area fraction <110> || CD (15°): recrystallized {dez(fib_cls['rx'], 1)} %, non-recrystallized {dez(fib_cls['def'], 1)} %")]
        K += [tx("rekristallisierter Anteil vs. Abstand zu Partikeln (", "recrystallized fraction vs. distance to particles (")
              + plain(psn_names.replace("$_2$", "2")) + "): " + "; ".join(
                  tx(f"{l_} µm: {dez(v_, 1)} % (Flächenanteil {dez(a_, 1)} %)", f"{l_} µm: {dez(v_, 1)} % (area fraction {dez(a_, 1)} %)")
                  for l_, v_, a_ in zip(dist_lbl, rx_vs_dist, area_vs_dist))]
    tables = {}
    if S.grainsize:
        def gs_row(label, key, nd=2):
            return f"{label};{dez(sumH[key], nd)}" + (f";{dez(sumL[key], nd)}" if SUBG else "")

        def gs_raw(label, a, b):
            return f"{label};{a}" + (f";{b}" if SUBG else "")

        K += ["", tx(f"KORNGRÖSSE UND KORNFORM ({PH_MATRIX}, bereinigte Karte; Körner: HAGB >= {dez(S.hagb, 0)}°, ",
                     f"GRAIN SIZE AND SHAPE ({PH_MATRIX}, cleaned map; grains: HAGB >= {dez(S.hagb, 0)}°, ")
              + (tx(f"Subkörner: Grenzen >= {dez(S.lagb, 1)}°, MTEX calcGrains)", f"subgrains: boundaries >= {dez(S.lagb, 1)}°, MTEX calcGrains)")
                 if SUBG else tx("ohne Subkornerkennung, MTEX calcGrains)", "without subgrain detection, MTEX calcGrains)")),
              tx(f"Statistik ohne Randkörner und ohne (Sub-)Körner mit < {S.gs_min_px} indizierten Punkten; "
                 "ECD = Äquivalentkreisdurchmesser; Aspektverhältnis = Haupt-/Nebenachse der Ellipse gleicher Flächenmomente 2. Ordnung",
                 f"Statistics without edge grains and without (sub)grains with < {S.gs_min_px} indexed points; "
                 "ECD = equivalent circle diameter; aspect ratio = major/minor axis of the ellipse with equal second area moments"),
              tx("Kennwert;Körner", "Quantity;grains") + (tx(";Subkörner", ";subgrains") if SUBG else ""),
              gs_raw(tx("Anzahl gesamt", "number in total"), sumH['n_all'], sumL['n_all']),
              gs_raw(tx("davon Randkörner", "thereof edge grains"), sumH['n_edge'], sumL['n_edge']),
              gs_raw(tx(f"davon zu klein (< {S.gs_min_px} Punkte)", f"thereof too small (< {S.gs_min_px} points)"), sumH['n_small'], sumL['n_small']),
              gs_raw(tx("gewertet", "evaluated"), sumH['n'], sumL['n']),
              gs_row(tx("Flächenanteil gewertet (%)", "area fraction evaluated (%)"), "area_frac", 1),
              gs_row(tx("ECD Mittelwert anzahlgewichtet (µm)", "ECD mean number-weighted (µm)"), "ecd_mean"),
              gs_row(tx("ECD Standardabweichung (µm)", "ECD standard deviation (µm)"), "ecd_sd"),
              gs_row(tx("ECD Median (µm)", "ECD median (µm)"), "ecd_med"),
              gs_row("ECD d10 (µm)", "d10"),
              gs_row("ECD d90 (µm)", "d90"),
              gs_row(tx("ECD Mittelwert flächengewichtet (µm)", "ECD mean area-weighted (µm)"), "ecd_aw"),
              gs_row(tx("mittlere Kornfläche (µm²)", "mean grain area (µm²)"), "area_mean"),
              gs_raw(tx("mittl. Sehnenlänge gesamt RD+CD (µm)", "mean linear intercept total RD+CD (µm)"), dez(mli['ges'], 2), dez(mli2['ges'], 2)),
              gs_raw(tx("mittl. Sehnenlänge RD (µm)", "mean linear intercept RD (µm)"), dez(mli['RD'], 2), dez(mli2['RD'], 2)),
              gs_raw(tx("mittl. Sehnenlänge CD (µm)", "mean linear intercept CD (µm)"), dez(mli['CD'], 2), dez(mli2['CD'], 2)),
              gs_raw(tx("Sehnenlängenverhältnis RD/CD", "intercept ratio RD/CD"), dez(mli['RD'] / max(mli['CD'], 1e-9), 2),
                     dez(mli2['RD'] / max(mli2['CD'], 1e-9), 2)),
              gs_row(tx("Aspektverhältnis Mittelwert", "aspect ratio mean"), "ar_mean"),
              gs_row(tx("Aspektverhältnis Median", "aspect ratio median"), "ar_med"),
              gs_row(tx("Aspektverhältnis flächengewichtet", "aspect ratio area-weighted"), "ar_aw"),
              gs_row(tx("Anteil Aspektverhältnis < 2 (%)", "fraction aspect ratio < 2 (%)"), "ar_lt2", 1),
              gs_raw(tx("Hauptachse (AV >= 1,5) mittlere Richtung zu RD (°, flächengewichtet)",
                        "major axis (AR >= 1.5) mean direction to RD (°, area-weighted)"), dez(sumH['ang_mean'], 1), dez(sumL['ang_mean'], 1)),
              gs_row(tx("Ausrichtungsgrad R (0 = regellos, 1 = parallel)", "degree of alignment R (0 = random, 1 = parallel)"), "ang_R"),
              gs_row(tx("Flächenanteil Hauptachse näher an RD als an CD (%)", "area fraction with major axis closer to RD than to CD (%)"),
                     "rd_frac", 1),
              tx(f"Korngrößenkennzahl G nach ASTM E112 aus mittl. Sehnenlänge HAGB: {dez(astm_G, 1)}",
                 f"Grain size number G according to ASTM E112 from the mean linear intercept of the HAGB: {dez(astm_G, 1)}")]
        if SUBG:
            K += [tx(f"Subkörner je Korn (Mittel über alle Körner): {dez(spg_mean, 2)}",
                     f"Subgrains per grain (mean over all grains): {dez(spg_mean, 2)}")]
        if sub_warn:
            K += [sub_warn]
        names_t = table_names()
        K += [tx("Sehnenlängen aus dem Linienschnitt (RD- und CD-Zeilen, alle Matrix-Messpunkte, "
                 f"Lückenüberbrückung <= {int(S.gap_bridge)} Punkte); Einzelwerte je (Sub-)Korn: ",
                 "Intercept lengths from the line intercept (RD and CD lines, all matrix points, "
                 f"gap bridging <= {int(S.gap_bridge)} points); values per (sub)grain: ")
              + ", ".join(names_t if SUBG else names_t[:1])]
        tables[names_t[0]] = grain_table(gsH, False)
        if SUBG:
            tables[names_t[1]] = grain_table(gsL, True)
    K += ["", tx("PARAMETER", "PARAMETERS"), "; ".join(f"{k} = {v}" for k, v in S.to_dict().items() if k != "title")]

    summary = {"Punkte": f"{NX} × {NY}", "Kratzer": len(scratch_lines) if S.scratch else None}
    if S.rx:
        summary["rekristallisiert"] = rx_frac["rx"]
    if S.grainsize:
        summary["Korngröße"] = sumH["ecd_mean"]
        if SUBG:
            summary["Subkorngröße"] = sumL["ecd_mean"]
            if sub_warn:
                summary["Hinweis"] = sub_warn

    # ------------------------------------------------------------------ Kennwerte für den Datensatzvergleich
    def num(x):
        try:
            x = float(x)
        except (TypeError, ValueError):
            return None
        return x if np.isfinite(x) else None

    M = {"datei": project.path.name, "pfad": str(project.path), "titel": TITLE,
         "ausschnitt_um": [float(v) for v in crop_um] if crop_um is not None else None,
         "punkte": [int(NX), int(NY)], "schrittweite_um": num(STEP), "flaeche_um": [num(W_UM), num(H_UM)],
         "rollen": {names[i]: r for i, r in roles.items()}, "einstellungen": S.to_dict(),
         "lagb": num(S.lagb), "hagb": num(S.hagb), "theta_c": num(S.rx_theta_c), "subkoerner": SUBG,
         "kratzer_n": len(scratch_lines) if S.scratch else None,
         "kratzer_flaeche_pct": num(scratch.mean() * 100) if HAS_SCRATCH else None,
         "kratzer_manuell": strokes,
         "kratzer_manuell_plus_pct": num(scr_added.mean() * 100) if strokes else None,
         "kratzer_manuell_minus_pct": num(scr_removed.mean() * 100) if strokes else None,
         "nullloesung_roh_pct": num(raw_frac[0]), "fehlindizierung_roh_pct": num(sum(raw_frac[i] for i in iMIS)),
         "cluster_pct": num(clean_frac["Cluster"])}
    if S.gb or S.grainsize:
        M.update(sehne_hagb_RD_um=num(mli["RD"]), sehne_hagb_CD_um=num(mli["CD"]),
                 sehne_alle_RD_um=num(mli2["RD"]), sehne_alle_CD_um=num(mli2["CD"]))
    if S.gb:
        cnt, edges = np.histogram(gb, bins=np.arange(0.0, 64.0, 1.0))
        M.update(hagb_anteil_pct=num(f_hagb), lagb_anteil_pct=num(f_lagb), hagb_ohne_luecke_pct=num(hagb_nogap),
                 n_grenzschnitte=int(gb.size), misori_hist={"kanten": edges.tolist(), "anzahl": cnt.astype(int).tolist()})
    if kam is not None and kv.size:
        M["kam_mittel_deg"] = num(kv.mean())
    if S.texture:
        M["faser110_CD_pct"] = num(fib_all)
    if S.rx:
        M.update(rx_pct=num(rx_frac["rx"]), sub_pct=num(rx_frac["sub"]) if SUBG else None,
                 def_pct=num(rx_frac["def"]) if SUBG else None, nicht_rx_pct=num(nonrx_frac),
                 rx_abstand_klassen=list(dist_lbl), rx_abstand_pct=[num(v) for v in rx_vs_dist],
                 rx_korn_ecd_median_um=num(rx_ecd_med))
    if S.grainsize:
        M.update(subkorn_ecd_median_um=num(sumL["ecd_med"]), subkorn_ecd_mittel_um=num(sumL["ecd_mean"]),
                 korn_ecd_median_um=num(sumH["ecd_med"]), korn_ecd_mittel_um=num(sumH["ecd_mean"]),
                 astm_G=num(astm_G), subkorn_aufloesung_schritte=num(sub_res))
    say(tx(f"Auswertung fertig ({time.time() - T0:.0f} s)", f"Analysis finished ({time.time() - T0:.0f} s)"))
    return Result(builders, "\n".join(K), S, project, crop_um, summary, tables, metrics=M,
                  scratch_auto=scratch_auto, crop_px=(r0, c0), strokes=strokes)
