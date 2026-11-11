# -*- coding: utf-8 -*-
"""
EBSD-Auswertung – grafische Oberfläche (PySide6) für die MTEX-Auswertung von Oxford-ctf-Dateien.

Start:      python ebsd_gui.py                      (bzw. EBSD-Auswertung.exe)
Selbsttest: EBSD-Auswertung.exe --selftest datei.ctf ausgabeordner [--crop x0 y0 x1 y1] [--lang en]
Sprache:    Menü „Sprache / Language“ (wird nach einem Neustart übernommen) oder --lang de|en
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path


# ------------------------------------------------------------------ Umgebung (vor numba/matplotlib)
def app_data_dir():
    base = os.environ.get("LOCALAPPDATA") or os.path.join(os.path.expanduser("~"), ".local", "share")
    d = Path(base) / "EBSD-Auswertung"
    d.mkdir(parents=True, exist_ok=True)
    return d


def resource_dir():
    return Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))


# ohne Konsole (exe) gibt es kein stdout/stderr: Ausgaben in eine Protokolldatei umleiten
if sys.stdout is None or sys.stderr is None:
    _logf = open(app_data_dir() / "ausgabe.log", "w", encoding="utf-8", buffering=1)
    sys.stdout = sys.stdout or _logf
    sys.stderr = sys.stderr or _logf

os.environ.setdefault("NUMBA_CACHE_DIR", str(app_data_dir() / "numba-cache"))
os.environ.setdefault("MPLCONFIGDIR", str(app_data_dir() / "matplotlib"))
os.environ["QT_API"] = "pyside6"

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.widgets import RectangleSelector  # noqa: E402

from PySide6.QtCore import (QLibraryInfo, QLocale, QProcess, QSettings, Qt, QThread, QTimer, QTranslator,  # noqa: E402
                            QUrl, Signal)
from PySide6.QtGui import (QAction, QActionGroup, QDesktopServices, QFont, QIcon, QKeySequence,  # noqa: E402
                           QPixmap)
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QButtonGroup, QCheckBox, QComboBox,  # noqa: E402
                               QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog, QFormLayout, QFrame,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit,
                               QListWidget, QListWidgetItem, QMainWindow, QMessageBox, QPlainTextEdit, QProgressBar,
                               QPushButton, QScrollArea, QSpinBox, QSplitter, QStatusBar, QTableWidget,
                               QTableWidgetItem, QTabWidget, QTextBrowser, QToolButton, QVBoxLayout, QWidget)

from ebsd_i18n import L, LANGUAGES, dez, get_language, set_language  # noqa: E402
from ebsd_settings import (APP_NAME, APP_VERSION, PARAM_INFO, ROLES, Abgebrochen, Settings,  # noqa: E402
                           app_title, compare_figures, fig_title, kennwerte_name, param_text, role_name)
import ebsd_datasets as dsets  # noqa: E402
import ebsd_marks  # noqa: E402

core = None          # ebsd_core wird beim ersten Laden im Hintergrund importiert (MTEX, numba, scipy)
ORG = "TU Bergakademie Freiberg IMF"


def load_core():
    global core
    if core is None:
        import ebsd_core as _core
        core = _core
    return core


def fmt(x, nd=1):
    return dez(x, nd)


HELP_DE = """
<h2>{name} {ver} – Kurzanleitung</h2>
<ol>
<li><b>ctf-Datei öffnen</b> (Datei → Öffnen oder Schaltfläche). Die Phasen aus dem Dateikopf erscheinen in der
Tabelle, rechts die Übersicht (Bandkontrast).</li>
<li><b>Phasenrollen</b> festlegen: genau eine <i>Matrix</i>, beliebig viele <i>Partikel</i>-Phasen (bleiben erhalten),
<i>Fehlindizierung</i> (wird verworfen und aus den Nachbarn ergänzt). Vorschlag: häufigste Phase = Matrix,
Phasen mit MAD-Median ≤ 0,8° = Partikel, sonst Fehlindizierung. Fehlindizierungen erkennt man an hohem MAD und
niedrigem Bandkontrast.</li>
<li><b>Funktionen</b> an- oder abwählen. Nicht benötigte Auswertungen sparen Rechenzeit; die Kornrekonstruktion
läuft nur, wenn IPF-Karten, Rekristallisation oder Korngröße gewählt sind.</li>
<li>Optional einen <b>Ausschnitt</b> wählen: in der Übersicht (Werkzeug <i>Ausschnitt</i>) ein Rechteck mit der Maus
aufziehen oder die Koordinaten eintragen.</li>
<li>Optional <b>Kratzer manuell bearbeiten</b> (siehe unten).</li>
<li><b>Auswertung starten</b>. Danach erscheinen die Abbildungen unter <i>Ergebnisse</i> als Vorschau.</li>
<li>Jede Auswertung wird automatisch als <b>Datensatz</b> gespeichert (Reiter <i>Datensätze / Vergleich</i>) und
steht nach einem Neustart wieder zur Verfügung: Vorschauen, Kennwerte und Kornlisten ohne Neuberechnung. Dieselbe Datei
mit demselben Ausschnitt ersetzt den vorhandenen Datensatz.</li>
<li><b>Vergleich:</b> angehakte Datensätze erscheinen im Kennwerte-Vergleich (Subkorngröße, HAGB-Anteil,
rekristallisierter Anteil gesamt und nahe Partikeln, ⟨110⟩ ∥ CD) und in der überlagerten Misorientierungsverteilung;
Bezeichnung und Reihenfolge sind frei wählbar, Export als Abbildung und CSV.</li>
<li>Abbildungen ankreuzen, Formate (PNG, PDF, SVG, TIFF) und Auflösung wählen und <b>exportieren</b>.
Die Datei <i>Kennwerte_und_Parameter.txt</i> wird immer mitgeschrieben, mit der Korngrößenauswertung zusätzlich
<i>Kornliste_Koerner.csv</i> und <i>Kornliste_Subkoerner.csv</i> (Semikolon, Dezimalkomma, für Excel).</li>
</ol>
<h3>Kratzer manuell markieren oder Markierung aufheben (neu in 1.3)</h3>
<ul>
<li>In der <i>Übersicht</i> das Werkzeug <b>Kratzer markieren</b> wählen und mit gedrückter linker Maustaste eine Linie
entlang des Kratzers ziehen. Die <b>Breite</b> (µm) wird in der Werkzeugleiste eingestellt.</li>
<li>Mit <b>Markierung aufheben</b> eine Linie über einen fälschlich erkannten Kratzer ziehen: dort gilt nichts mehr als
Kratzer – weder automatisch erkannte noch manuell markierte Punkte.</li>
<li><b>Automatische Erkennung anzeigen</b> führt nur die Kratzererkennung aus (für die Gesamtkarte bzw. den gewählten
Ausschnitt) und blendet das Ergebnis rot ein, damit Fehlerkennungen gezielt aufgehoben werden können. Nach einer
Auswertung wird die automatische Erkennung ebenfalls eingeblendet.</li>
<li>Farben: rot = automatisch erkannt, magenta = manuell markiert, hellblau = Markierung aufgehoben.</li>
<li><b>Rückgängig</b> (Strg+Z) entfernt die letzte Linie, <b>Alle löschen</b> alle manuellen Linien dieser Datei.
Mit dem Lupen- und Verschiebe-Werkzeug lässt sich für genaues Zeichnen in die Karte zoomen.</li>
<li>Die Markierungen werden je ctf-Datei gespeichert, beim erneuten Öffnen wiederhergestellt, im Datensatz abgelegt und
bei „Neu berechnen“ verwendet. Sie gelten auch, wenn die automatische Kratzererkennung ausgeschaltet ist.</li>
<li>Abbildung 00 und die Kennwerte weisen automatische und manuelle Anteile getrennt aus.</li>
</ul>
<h3>Subkornerkennung an-/abwählen (neu in 1.3)</h3>
<ul>
<li>Funktion <b>Subkörner erkennen</b>: Bei großer Schrittweite sind Subkörner nicht aufgelöst; dann die Funktion
abwählen. Faustregel: ein Subkorn sollte mindestens etwa 5 Schrittweiten groß sein. Das Programm gibt einen Hinweis aus,
wenn der Median-ECD der Subkörner darunter liegt.</li>
<li>Ohne Subkornerkennung entfallen Subkornrekonstruktion, LAGB-Linien in den Karten, Abbildung 12, die
Subkorn-Kornliste und die Gegenprüfung auf Subkornebene. Die Rekristallisation wird dann nur auf Kornebene
klassifiziert: rekristallisiert (GOS ≤ θc) bzw. nicht rekristallisiert (GOS &gt; θc), ohne Trennung in
substrukturiert/verformt. Die Misorientierungsverteilung (Linienschnitt) bleibt unverändert.</li>
</ul>
<h3>Hinweise</h3>
<ul>
<li><b>Sprache:</b> Menü <i>Sprache / Language</i>. Oberfläche, Abbildungen, Kennwerte und Kornlisten erscheinen in der
gewählten Sprache (englisch: CSV mit Komma als Trennzeichen und Dezimalpunkt). Der Wechsel erfordert einen Neustart.</li>
<li><b>Koordinaten:</b> RD = horizontal (x), CD = Stauchrichtung vertikal (y), TD = Schliffnormale (z).
Die Euler-Winkel werden mit der MTEX-Standardkorrektur für Oxford-Daten auf das Kartensystem bezogen.</li>
<li><b>Mg₂Si-Erkennung:</b> große Nulllösungs-Cluster mit niedrigem Bandkontrast werden als Mg₂Si gewertet
(präparationsbedingt nicht indizierbar). Das ist eine Interpretation – stichprobenartig per EDX prüfen.
Ausgeschaltet bleiben die Cluster „nicht indiziert“ und zählen nicht als Partikel.</li>
<li><b>Kratzer:</b> lange gerade Linien aus Nulllösungen/niedrigem Bandkontrast; die Punkte werden nicht als
Messwerte verwendet und für Karten aus den Nachbarn ergänzt. Ergebnis in Abbildung 00 prüfen.</li>
<li><b>Korngröße:</b> Körner = von HAGB umschlossene Bereiche, Subkörner = von Grenzen ≥ LAGB umschlossene
Bereiche der bereinigten Karte. Randkörner und (Sub-)Körner mit weniger indizierten Punkten als eingestellt
(Reiter Parameter) gehen nicht in die Statistik ein, bleiben aber in den Kornlisten. Größe als
Äquivalentkreisdurchmesser (anzahl- und flächengewichtet) sowie als mittlere Sehnenlänge aus dem Linienschnitt
in RD und CD (daraus die Kennzahl G nach ASTM E112). Aspektverhältnis = Haupt-/Nebenachse der Ellipse mit
gleichen Flächenmomenten 2. Ordnung; Winkel der Hauptachse zu RD gegen den Uhrzeigersinn.</li>
<li><b>KAM</b> hängt von der Schrittweite ab – nur Karten gleicher Schrittweite vergleichen.</li>
<li><b>Rekristallisation:</b> EBSD nach dem Versuch unterscheidet nicht zwischen dynamischer und
metadynamischer/statischer Rekristallisation; der Anteil gilt nur bei direktem Abschrecken als DRX-Anteil.</li>
<li><b>Vektorexporte (PDF/SVG):</b> Karten werden als Rasterbild eingebettet (sonst sehr große Dateien),
Beschriftungen, Legenden und Diagramme bleiben Vektorgrafik.</li>
<li><b>Gespeicherte Datensätze exportieren:</b> Für Abbildungen in voller Auflösung wird der Datensatz aus der
ctf-Datei mit den gespeicherten Rollen, Parametern, Kratzermarkierungen und dem Ausschnitt neu berechnet; alternativ
werden die Vorschaubilder, Kennwerte und Kornlisten kopiert. Ablage: %LOCALAPPDATA%\\EBSD-Auswertung\\datensaetze.</li>
<li>Der erste Lauf nach dem Programmstart ist langsamer, weil Rechenroutinen kompiliert werden
(Zwischenspeicher unter %LOCALAPPDATA%\\EBSD-Auswertung).</li>
</ul>
"""

HELP_EN = """
<h2>{name} {ver} – Quick guide</h2>
<ol>
<li><b>Open a ctf file</b> (File → Open or the button). The phases from the file header appear in the table,
the overview (band contrast) on the right.</li>
<li>Set the <b>phase roles</b>: exactly one <i>Matrix</i>, any number of <i>Particle</i> phases (kept),
<i>Misindexing</i> (discarded and filled in from the neighbours). Suggestion: most frequent phase = matrix,
phases with MAD median ≤ 0.8° = particle, otherwise misindexing. Misindexing shows high MAD and low band contrast.</li>
<li>Switch <b>functions</b> on or off. Analyses that are not needed save computing time; grain reconstruction only
runs if IPF maps, recrystallization or grain size are selected.</li>
<li>Optionally select an <b>area</b>: drag a rectangle in the overview (tool <i>Area</i>) or enter the coordinates.</li>
<li>Optionally <b>edit scratches manually</b> (see below).</li>
<li><b>Start analysis</b>. The figures then appear as previews under <i>Results</i>.</li>
<li>Each analysis is stored automatically as a <b>dataset</b> (tab <i>Datasets / Comparison</i>) and is available
again after a restart: previews, key figures and grain lists without recalculation. The same file with the same area
replaces the existing dataset.</li>
<li><b>Comparison:</b> ticked datasets appear in the key-figure comparison (subgrain size, HAGB fraction,
recrystallized fraction in total and near particles, ⟨110⟩ ∥ CD) and in the overlaid misorientation distribution;
label and order can be chosen freely, export as figure and CSV.</li>
<li>Tick figures, choose formats (PNG, PDF, SVG, TIFF) and resolution and <b>export</b>.
<i>Results_and_parameters.txt</i> is always written, with the grain size analysis also <i>Grain_list_grains.csv</i>
and <i>Grain_list_subgrains.csv</i> (comma-separated, decimal point).</li>
</ol>
<h3>Mark scratches manually or remove a marking (new in 1.3)</h3>
<ul>
<li>In the <i>Overview</i>, choose the tool <b>Mark scratch</b> and draw a line along the scratch with the left mouse
button held down. The <b>width</b> (µm) is set in the toolbar.</li>
<li>With <b>Unmark</b>, draw a line over a falsely detected scratch: nothing counts as a scratch there any more –
neither automatically detected nor manually marked points.</li>
<li><b>Show automatic detection</b> runs only the scratch detection (for the whole map or the selected area) and shows
the result in red so that false detections can be removed specifically. After an analysis, the automatic detection
is shown as well.</li>
<li>Colours: red = detected automatically, magenta = marked manually, light blue = marking removed.</li>
<li><b>Undo</b> (Ctrl+Z) removes the last line, <b>Clear all</b> removes all manual lines of this file.
Use the magnifier and pan tools to zoom into the map for precise drawing.</li>
<li>The marks are saved per ctf file, restored when the file is opened again, stored in the dataset and used for
“Recalculate”. They also apply when automatic scratch detection is switched off.</li>
<li>Figure 00 and the key figures report automatic and manual contributions separately.</li>
</ul>
<h3>Switch subgrain detection on/off (new in 1.3)</h3>
<ul>
<li>Function <b>Detect subgrains</b>: at a large step size, subgrains are not resolved; then switch the function off.
Rule of thumb: a subgrain should be at least about 5 step sizes in size. The program issues a note if the median
ECD of the subgrains is below that.</li>
<li>Without subgrain detection, the subgrain reconstruction, LAGB lines in the maps, figure 12, the subgrain list and
the cross-check at subgrain level are omitted. Recrystallization is then classified at grain level only:
recrystallized (GOS ≤ θc) or non-recrystallized (GOS &gt; θc), without separating substructured/deformed.
The misorientation distribution (line intercept) is unchanged.</li>
</ul>
<h3>Notes</h3>
<ul>
<li><b>Language:</b> menu <i>Sprache / Language</i>. User interface, figures, key figures and grain lists appear in the
selected language (German: CSV with semicolon and decimal comma). Changing the language requires a restart.</li>
<li><b>Coordinates:</b> RD = horizontal (x), CD = compression direction vertical (y), TD = section normal (z).
The Euler angles are referred to the map frame with the MTEX standard correction for Oxford data.</li>
<li><b>Mg₂Si detection:</b> large zero-solution clusters with low band contrast are interpreted as Mg₂Si
(not indexable due to preparation). This is an interpretation – check by EDX on a sample basis.
When switched off, the clusters remain “non-indexed” and do not count as particles.</li>
<li><b>Scratches:</b> long straight lines of zero solutions/low band contrast; the points are not used as measured
values and are filled in from the neighbours for maps. Check the result in figure 00.</li>
<li><b>Grain size:</b> grains = areas enclosed by HAGB, subgrains = areas enclosed by boundaries ≥ LAGB in the
cleaned map. Edge grains and (sub)grains with fewer indexed points than set (tab Parameters) are excluded from the
statistics but remain in the grain lists. Size as equivalent circle diameter (number- and area-weighted) and as mean
linear intercept in RD and CD (giving the grain size number G according to ASTM E112). Aspect ratio = major/minor
axis of the ellipse with equal second area moments; angle of the major axis to RD counter-clockwise.</li>
<li><b>KAM</b> depends on the step size – compare only maps with the same step size.</li>
<li><b>Recrystallization:</b> post-test EBSD cannot distinguish between dynamic and metadynamic/static
recrystallization; the fraction only counts as DRX fraction for direct quenching.</li>
<li><b>Vector exports (PDF/SVG):</b> maps are embedded as raster images (otherwise very large files);
labels, legends and diagrams remain vector graphics.</li>
<li><b>Exporting stored datasets:</b> for full-resolution figures the dataset is recalculated from the ctf file with
the stored roles, parameters, scratch marks and area; alternatively, the previews, key figures and grain lists are
copied. Storage: %LOCALAPPDATA%\\EBSD-Auswertung\\datensaetze.</li>
<li>The first run after starting the program is slower because routines are compiled
(cache in %LOCALAPPDATA%\\EBSD-Auswertung).</li>
</ul>
"""


def help_html():
    return L(HELP_DE, HELP_EN).format(name=app_title(), ver=APP_VERSION)


# ------------------------------------------------------------------ manuelle Kratzermarkierungen je Datei
class MarkStore:
    """speichert die manuellen Kratzermarkierungen je ctf-Datei (auch ohne Auswertung)"""

    def __init__(self, path):
        self.path = Path(path)

    def _all(self):
        try:
            return json.loads(self.path.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            return {}

    def get(self, ctf):
        return ebsd_marks.clean_strokes(self._all().get(str(Path(ctf)).lower()))

    def has(self, ctf):
        return str(Path(ctf)).lower() in self._all()

    def set(self, ctf, strokes):
        d = self._all()
        key = str(Path(ctf)).lower()
        if strokes:
            d[key] = strokes
        else:
            d.pop(key, None)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(d, ensure_ascii=False), encoding="utf-8")
        tmp.replace(self.path)


# ------------------------------------------------------------------ Hintergrundaufgaben
class Task(QThread):
    progress = Signal(int, str)
    message = Signal(str)
    done = Signal(object)
    failed = Signal(str, bool)          # Text, abgebrochen?

    def __init__(self, fn, parent=None):
        super().__init__(parent)
        self.fn = fn
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def check(self):
        if self._cancel:
            raise Abgebrochen()

    def run(self):
        try:
            out = self.fn(self)
        except Abgebrochen:
            self.failed.emit(L("Vorgang abgebrochen.", "Operation cancelled."), True)
        except Exception as e:  # noqa: BLE001
            self.failed.emit(f"{type(e).__name__}: {e}\n\n{traceback.format_exc()}", False)
        else:
            self.done.emit(out)


# ------------------------------------------------------------------ Übersicht: Ausschnitt wählen, Kratzer markieren
OV_AUTO = (235, 30, 30, 200)           # automatisch erkannt
OV_ADD = (215, 25, 215, 225)           # manuell markiert
OV_REM_AUTO = (60, 175, 255, 215)      # automatische Markierung aufgehoben
OV_REM = (60, 175, 255, 70)            # Bereich einer Aufhebe-Linie ohne automatische Markierung


class NavBar(NavigationToolbar2QT):
    """matplotlib-Werkzeugleiste, beschränkt auf Ansicht, Verschieben und Zoomen"""

    def __init__(self, canvas, parent):
        tips = {"Home": L("Gesamte Karte anzeigen", "Show the whole map"),
                "Back": L("vorherige Ansicht", "previous view"), "Forward": L("nächste Ansicht", "next view"),
                "Pan": L("Verschieben (linke Maustaste) bzw. Zoomen (rechte Maustaste)", "Pan (left mouse button) or zoom (right mouse button)"),
                "Zoom": L("auf Rechteck zoomen", "zoom to rectangle")}
        items = [(t[0], tips[t[0]], t[2], t[3]) for t in NavigationToolbar2QT.toolitems if t[0] in tips]
        self.toolitems = items[:3] + [(None, None, None, None)] + items[3:]
        super().__init__(canvas, parent, coordinates=True)


class OverviewCanvas(FigureCanvasQTAgg):
    rect_selected = Signal(float, float, float, float)
    stroke_drawn = Signal(float, float, float, float)

    def __init__(self, parent=None):
        self.fig = Figure(figsize=(8, 6), dpi=100, facecolor="white")
        super().__init__(self.fig)
        self.setParent(parent)
        self.ax = self.fig.add_axes([0.08, 0.09, 0.9, 0.86])
        self.ax.set_axis_off()
        self.ax.text(0.5, 0.5, L("Noch keine Datei geladen.\n\nDatei → Öffnen (Strg+O)", "No file loaded yet.\n\nFile → Open (Ctrl+O)"),
                     ha="center", va="center", transform=self.ax.transAxes, fontsize=13, color="0.4")
        self.selector = None
        self._mute = False
        self.extent = None
        self.overlay = None
        self.mode = "crop"                 # "crop" | "add" | "remove"
        self.mark_width = 2.0              # µm
        self._drag = None
        self._line = None
        self.mpl_connect("button_press_event", self._press)
        self.mpl_connect("motion_notify_event", self._move)
        self.mpl_connect("button_release_event", self._release)

    def show_map(self, img, w_um, h_um, title):
        self.fig.clf()
        self.ax = self.fig.add_axes([0.08, 0.09, 0.9, 0.86])
        self.ax.imshow(img, cmap="gray", extent=(0, w_um, h_um, 0), interpolation="nearest")
        self.overlay = self.ax.imshow(np.zeros(img.shape + (4,), np.uint8), extent=(0, w_um, h_um, 0),
                                      interpolation="nearest", zorder=3)
        self.ax.set_xlabel(L("x = RD / µm", "x = RD / µm"))
        self.ax.set_ylabel(L("y = CD / µm", "y = CD / µm"))
        self.ax.set_title(title, fontsize=11)
        self.ax.grid(color="y", alpha=0.3, lw=0.5)
        self.extent = (w_um, h_um)
        self.selector = RectangleSelector(self.ax, self._on_select, useblit=True, button=[1], minspanx=2, minspany=2,
                                          spancoords="data", interactive=True,
                                          props=dict(facecolor="red", edgecolor="red", alpha=0.15, fill=True, linewidth=2))
        self.selector.set_visible(False)
        self.set_mode(self.mode)
        self.draw_idle()

    def set_overlay(self, rgba):
        if self.overlay is not None and rgba is not None:
            self.overlay.set_data(rgba)
            self.draw_idle()

    def set_mode(self, mode):
        self.mode = mode
        if self.selector is not None:
            self.selector.set_active(mode == "crop")
        self.setCursor(Qt.ArrowCursor if mode == "crop" else Qt.CrossCursor)

    def _on_select(self, e1, e2):
        if self._mute or self.mode != "crop":
            return
        x0, x1 = sorted((e1.xdata, e2.xdata))
        y0, y1 = sorted((e1.ydata, e2.ydata))
        w, h = self.extent
        self.rect_selected.emit(max(0.0, x0), max(0.0, y0), min(w, x1), min(h, y1))

    def set_rect(self, rect):
        if self.selector is None:
            return
        self._mute = True
        if rect is None:
            self.selector.set_visible(False)
        else:
            x0, y0, x1, y1 = rect
            self.selector.extents = (x0, x1, y0, y1)
            self.selector.set_visible(True)
        self.selector.set_active(self.mode == "crop")
        self._mute = False
        self.draw_idle()

    # --- Kratzerlinien zeichnen
    def _lw_points(self):
        x0, x1 = self.ax.get_xlim()
        return max(1.0, self.mark_width * self.ax.bbox.width / max(abs(x1 - x0), 1e-9) * 72 / self.fig.dpi)

    def _press(self, ev):
        if (self.mode == "crop" or self.extent is None or ev.button != 1 or ev.inaxes is not self.ax
                or ev.xdata is None or self.widgetlock.locked()):
            return
        self._drag = (ev.xdata, ev.ydata)
        col = np.array(OV_ADD if self.mode == "add" else OV_REM_AUTO) / 255
        self._line = Line2D([ev.xdata, ev.xdata], [ev.ydata, ev.ydata], color=col[:3], alpha=0.6, lw=self._lw_points(),
                            solid_capstyle="round", zorder=5)
        self.ax.add_line(self._line)
        self.draw_idle()

    def _move(self, ev):
        if self._drag is None or ev.inaxes is not self.ax or ev.xdata is None:
            return
        self._line.set_data([self._drag[0], ev.xdata], [self._drag[1], ev.ydata])
        self.draw_idle()

    def _release(self, ev):
        if self._drag is None:
            return
        x0, y0 = self._drag
        xs, ys = self._line.get_data()
        x1, y1 = (ev.xdata, ev.ydata) if (ev.inaxes is self.ax and ev.xdata is not None) else (xs[-1], ys[-1])
        self._line.remove()
        self._line = None
        self._drag = None
        self.draw_idle()
        self.stroke_drawn.emit(float(x0), float(y0), float(x1), float(y1))


# ------------------------------------------------------------------ Vorschau
class Preview(QScrollArea):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.label = QLabel(L("Noch keine Ergebnisse.\n\nNach der Auswertung erscheinen hier die Abbildungen.",
                              "No results yet.\n\nThe figures appear here after the analysis."))
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setStyleSheet("color: #666; background: white;")
        self.setWidget(self.label)
        self.setWidgetResizable(True)
        self.setAlignment(Qt.AlignCenter)
        self.pix = None
        self.fit = True
        self.scale = 1.0

    def set_pixmap(self, pix):
        self.pix = pix
        self.refresh()

    def set_fit(self):
        self.fit = True
        self.refresh()

    def zoom(self, factor=None, absolute=None):
        if self.pix is None:
            return
        if self.fit:
            self.scale = self._fit_scale()
        self.fit = False
        self.scale = absolute if absolute is not None else max(0.1, min(4.0, self.scale * factor))
        self.refresh()

    def _fit_scale(self):
        vw, vh = self.viewport().width() - 4, self.viewport().height() - 4
        return min(vw / self.pix.width(), vh / self.pix.height(), 1.0)

    def refresh(self):
        if self.pix is None:
            return
        s = self._fit_scale() if self.fit else self.scale
        self.setWidgetResizable(self.fit)
        p = self.pix.scaled(int(self.pix.width() * s), int(self.pix.height() * s), Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.label.setPixmap(p)
        self.label.setStyleSheet("background: white;")
        if not self.fit:
            self.label.resize(p.size())

    def resizeEvent(self, e):
        super().resizeEvent(e)
        if self.fit and self.pix is not None:
            self.refresh()


# ------------------------------------------------------------------ Hauptfenster
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"{app_title()} {APP_VERSION} – " + L("MTEX für Python", "MTEX for Python"))
        icon = resource_dir() / "icon.png"
        if icon.exists():
            self.setWindowIcon(QIcon(str(icon)))
        self.qs = QSettings(ORG, APP_NAME)
        self.project = None
        self.result = None
        self.previews = {}
        self.task = None
        self.export_folder_user = False
        self.autotest = None            # {"out": Ordner, "crop": (...)} – automatischer Durchlauf zum Testen
        self.store = dsets.DatasetStore(app_data_dir() / "datensaetze")
        self.marks = MarkStore(app_data_dir() / "kratzer_markierungen.json")
        self.strokes = []               # manuelle Kratzermarkierungen der geladenen Datei
        self.auto_scratch = None        # automatisch erkannte Kratzer (Gesamtkarte) zur Anzeige
        self.auto_scratch_src = ""
        self.live_ds = None             # Datensatz-ID des Ergebnisses im Speicher (self.result)
        self.stored_view = None         # gespeicherter Datensatz, dessen Ergebnisse gerade angezeigt werden
        self._pending = None            # Neuberechnung eines gespeicherten Datensatzes
        self._ds_mute = False
        self._restart = False
        self._build_ui()
        self._restore()
        self.refresh_datasets()

    # ---------------------------------------------------------------- Aufbau
    def _build_ui(self):
        # Menü
        m = self.menuBar().addMenu(L("&Datei", "&File"))
        a = QAction(L("ctf-Datei &öffnen …", "&Open ctf file …"), self, shortcut=QKeySequence.Open, triggered=self.open_file)
        m.addAction(a)
        m.addSeparator()
        m.addAction(QAction(L("Einstellungen &speichern …", "&Save settings …"), self, triggered=self.save_settings_file))
        m.addAction(QAction(L("Einstellungen &laden …", "&Load settings …"), self, triggered=self.load_settings_file))
        m.addAction(QAction(L("Parameter auf &Standardwerte", "Reset parameters to &defaults"), self, triggered=self.reset_params))
        m.addSeparator()
        m.addAction(QAction(L("Ordner der gespeicherten &Datensätze öffnen", "Open folder of stored &datasets"), self,
                            triggered=lambda: QDesktopServices.openUrl(QUrl.fromLocalFile(str(self.store.root)))))
        m.addSeparator()
        m.addAction(QAction(L("&Beenden", "E&xit"), self, shortcut=QKeySequence.Quit, triggered=self.close))
        e = self.menuBar().addMenu(L("&Bearbeiten", "&Edit"))
        self.act_undo = QAction(L("Letzte Kratzermarkierung &rückgängig", "&Undo last scratch mark"), self,
                                shortcut=QKeySequence.Undo, triggered=self.undo_stroke)
        e.addAction(self.act_undo)
        self.act_clear = QAction(L("Alle Kratzermarkierungen dieser Datei &löschen", "&Clear all scratch marks of this file"), self,
                                 triggered=self.clear_strokes)
        e.addAction(self.act_clear)
        lm = self.menuBar().addMenu("&Sprache / Language")
        grp = QActionGroup(self)
        for code, name in LANGUAGES.items():
            act = QAction(name, self, checkable=True)
            act.setChecked(code == get_language())
            act.triggered.connect(lambda _=False, c=code: self.change_language(c))
            grp.addAction(act)
            lm.addAction(act)
        h = self.menuBar().addMenu(L("&Hilfe", "&Help"))
        h.addAction(QAction(L("&Kurzanleitung", "&Quick guide"), self, shortcut=QKeySequence.HelpContents, triggered=self.show_help))
        h.addAction(QAction(L("Ü&ber", "&About"), self, triggered=self.show_about))

        split = QSplitter(Qt.Horizontal)
        self.setCentralWidget(split)

        # ---------- linke Seite
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(6, 6, 6, 6)
        self.ltabs = QTabWidget()
        lv.addWidget(self.ltabs, 1)

        # Tab Auswertung
        page = QWidget()
        pv = QVBoxLayout(page)

        g = QGroupBox(L("1  Datei", "1  File"))
        gl = QVBoxLayout(g)
        self.btn_open = QPushButton(L("ctf-Datei öffnen …", "Open ctf file …"))
        self.btn_open.clicked.connect(self.open_file)
        self.lbl_file = QLabel(L("keine Datei geladen", "no file loaded"))
        self.lbl_file.setWordWrap(True)
        self.lbl_file.setTextInteractionFlags(Qt.TextSelectableByMouse)
        gl.addWidget(self.btn_open)
        gl.addWidget(self.lbl_file)
        pv.addWidget(g)

        g = QGroupBox(L("2  Phasen und ihre Rolle", "2  Phases and their role"))
        gl = QVBoxLayout(g)
        self.tbl = QTableWidget(0, 5)
        self.tbl.setHorizontalHeaderLabels([L("Phase", "Phase"), L("Anteil / %", "Fraction / %"), L("MAD / °", "MAD / °"), "BC",
                                            L("Rolle", "Role")])
        self.tbl.verticalHeader().setVisible(False)
        self.tbl.setSelectionMode(QAbstractItemView.NoSelection)
        self.tbl.setEditTriggers(QAbstractItemView.NoEditTriggers)
        hh = self.tbl.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.Stretch)
        for c in (1, 2, 3, 4):
            hh.setSectionResizeMode(c, QHeaderView.ResizeToContents)
        self.tbl.setMinimumHeight(150)
        gl.addWidget(self.tbl)
        self.lbl_zero = QLabel("")
        hint = QLabel(L("Matrix: genau eine Phase · Partikel: bleiben erhalten · Fehlindizierung: wird verworfen und "
                        "ergänzt. Hinweis auf Fehlindizierung: hoher MAD, niedriger Bandkontrast (BC).",
                        "Matrix: exactly one phase · Particle: kept · Misindexing: discarded and filled in. "
                        "Indication of misindexing: high MAD, low band contrast (BC)."))
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #555; font-size: 9pt;")
        gl.addWidget(self.lbl_zero)
        gl.addWidget(hint)
        pv.addWidget(g)

        g = QGroupBox(L("3  Funktionen", "3  Functions"))
        gl = QGridLayout(g)
        self.chk = {}
        items = [("scratch", L("Kratzer automatisch erkennen und entfernen", "Detect and remove scratches automatically"),
                  L("Lange gerade Linien aus Nulllösungen/niedrigem Bandkontrast werden markiert, nicht als Messwerte "
                    "verwendet und aus den Nachbarn ergänzt. Ergebnis: Abbildung 00. Manuelle Markierungen (Reiter "
                    "„Übersicht“) gelten auch bei ausgeschalteter automatischer Erkennung.",
                    "Long straight lines of zero solutions/low band contrast are marked, not used as measured values and "
                    "filled in from the neighbours. Result: figure 00. Manual marks (tab “Overview”) also apply when "
                    "automatic detection is switched off.")),
                 ("mg2si", L("Mg₂Si-Erkennung (Nulllösungs-Cluster)", "Mg₂Si detection (zero-solution clusters)"),
                  L("Große zusammenhängende Nulllösungen mit niedrigem BC als Mg₂Si werten. Aus: Cluster bleiben "
                    "„nicht indiziert“ und zählen nicht als Partikel.",
                    "Interpret large connected zero solutions with low BC as Mg₂Si. Off: clusters remain "
                    "“non-indexed” and do not count as particles.")),
                 ("ipf", L("IPF-Karten (CD, RD, TD)", "IPF maps (CD, RD, TD)"),
                  L("Orientierungskarten mit Korngrenzen (benötigt Kornrekonstruktion).",
                    "Orientation maps with grain boundaries (requires grain reconstruction).")),
                 ("kam", L("KAM-Karte", "KAM map"), L("Kernel Average Misorientation, 3×3-Nachbarschaft.",
                                                       "Kernel average misorientation, 3×3 neighbourhood.")),
                 ("gb", L("Korngrenz-Misorientierung", "Grain boundary misorientation"),
                  L("Linienschnitt mit Lückenüberbrückung, LAGB/HAGB-Anteile, Mackenzie.",
                    "Line intercept with gap bridging, LAGB/HAGB fractions, Mackenzie.")),
                 ("texture", L("Textur (inverse Polfiguren, Polfiguren)", "Texture (inverse pole figures, pole figures)"),
                  L("ODF mit MTEX calcDensity.", "ODF with MTEX calcDensity.")),
                 ("rx", L("Rekristallisation (GOS)", "Recrystallization (GOS)"),
                  L("Klassifizierung rekristallisiert/substrukturiert/verformt, Gegenprüfung, PSN-Indikator "
                    "(benötigt Kornrekonstruktion und KAM).",
                    "Classification recrystallized/substructured/deformed, cross-check, PSN indicator "
                    "(requires grain reconstruction and KAM).")),
                 ("grainsize", L("Korngröße, Subkorngröße, Aspektverhältnis", "Grain size, subgrain size, aspect ratio"),
                  L("Äquivalentkreisdurchmesser der Körner (HAGB) und Subkörner (≥ LAGB), mittlere Sehnenlänge, "
                    "Aspektverhältnis und Ausrichtung der Hauptachse; Kornlisten als CSV (benötigt Kornrekonstruktion).",
                    "Equivalent circle diameter of grains (HAGB) and subgrains (≥ LAGB), mean linear intercept, "
                    "aspect ratio and orientation of the major axis; grain lists as CSV (requires grain reconstruction).")),
                 ("subgrains", L("Subkörner erkennen (Grenzen ≥ LAGB)", "Detect subgrains (boundaries ≥ LAGB)"),
                  L("Subkornrekonstruktion für LAGB-Linien, Subkorngröße (Abb. 12), Rekristallisation "
                    "(substrukturiert/verformt) und Gegenprüfung auf Subkornebene.\nAbwählen, wenn die Schrittweite für "
                    "Subkörner zu groß ist (Faustregel: Subkorn ≥ ca. 5 Schrittweiten). Ohne Subkörner wird die "
                    "Rekristallisation nur auf Kornebene klassifiziert (rekristallisiert / nicht rekristallisiert).",
                    "Subgrain reconstruction for LAGB lines, subgrain size (fig. 12), recrystallization "
                    "(substructured/deformed) and cross-check at subgrain level.\nSwitch off if the step size is too "
                    "large for subgrains (rule of thumb: subgrain ≥ approx. 5 step sizes). Without subgrains, "
                    "recrystallization is classified at grain level only (recrystallized / non-recrystallized)."))]
        lbl = QLabel(L("Bereinigung", "Cleaning"))
        lbl.setStyleSheet("font-weight: bold;")
        gl.addWidget(lbl, 0, 0)
        row = 1
        for key, text, tip in items:
            if key == "ipf":
                lbl = QLabel(L("Auswertungen", "Analyses"))
                lbl.setStyleSheet("font-weight: bold; margin-top: 4px;")
                gl.addWidget(lbl, row, 0)
                row += 1
            cb = QCheckBox(text)
            cb.setChecked(True)
            cb.setToolTip(tip)
            self.chk[key] = cb
            gl.addWidget(cb, row, 0)
            row += 1
            if key == "scratch":
                self.lbl_marks = QLabel("")
                self.lbl_marks.setStyleSheet("color: #555; font-size: 9pt; margin-left: 22px;")
                self.lbl_marks.setWordWrap(True)
                gl.addWidget(self.lbl_marks, row, 0)
                row += 1
        self.chk["subgrains"].setStyleSheet("margin-left: 18px;")
        pv.addWidget(g)

        g = QGroupBox(L("4  Ausschnitt", "4  Area"))
        gl = QGridLayout(g)
        self.chk_crop = QCheckBox(L("nur Ausschnitt auswerten", "analyse selected area only"))
        self.chk_crop.toggled.connect(self._crop_toggled)
        gl.addWidget(self.chk_crop, 0, 0, 1, 4)
        self.crop_spins = []
        for k, (lab, r, c) in enumerate([("x0", 1, 0), ("y0", 1, 2), ("x1", 2, 0), ("y1", 2, 2)]):
            sp = QDoubleSpinBox()
            sp.setDecimals(1)
            sp.setRange(0, 1e5)
            sp.setSuffix(" µm")
            sp.valueChanged.connect(self._crop_spin_changed)
            gl.addWidget(QLabel(lab), r, c)
            gl.addWidget(sp, r, c + 1)
            self.crop_spins.append(sp)
        self.btn_full = QPushButton(L("Gesamte Karte", "Whole map"))
        self.btn_full.clicked.connect(lambda: self.chk_crop.setChecked(False))
        gl.addWidget(self.btn_full, 3, 0, 1, 2)
        lab = QLabel(L("oder in der Übersicht ein Rechteck mit der Maus aufziehen", "or drag a rectangle in the overview"))
        lab.setStyleSheet("color: #555; font-size: 9pt;")
        lab.setWordWrap(True)
        gl.addWidget(lab, 3, 2, 1, 2)
        pv.addWidget(g)

        g = QGroupBox(L("5  Probenbezeichnung (Titel der Abbildungen)", "5  Sample name (title of the figures)"))
        gl = QVBoxLayout(g)
        self.ed_title = QLineEdit(Settings().title)
        gl.addWidget(self.ed_title)
        pv.addWidget(g)
        pv.addStretch(1)
        sa = QScrollArea()
        sa.setWidget(page)
        sa.setWidgetResizable(True)
        sa.setFrameShape(QFrame.NoFrame)
        self.ltabs.addTab(sa, L("Auswertung", "Analysis"))

        # Tab Parameter
        page = QWidget()
        form = QFormLayout(page)
        form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.param = {}
        for entry in PARAM_INFO:
            key, _, _, lo, hi, step, nd = entry[:7]
            label, unit, tip = param_text(entry)
            if nd == 0:
                sp = QSpinBox()
                sp.setRange(int(lo), int(hi))
                sp.setSingleStep(int(step))
            else:
                sp = QDoubleSpinBox()
                sp.setDecimals(nd)
                sp.setRange(lo, hi)
                sp.setSingleStep(step)
            if unit:
                sp.setSuffix(f" {unit}")
            sp.setToolTip(tip)
            lab = QLabel(label)
            lab.setToolTip(tip)
            form.addRow(lab, sp)
            self.param[key] = sp
        btn = QPushButton(L("Standardwerte", "Defaults"))
        btn.clicked.connect(self.reset_params)
        form.addRow("", btn)
        note = QLabel(L("Fahren Sie mit der Maus über einen Parameter für eine Erläuterung. "
                        "Die Parameter werden beim Beenden gespeichert.",
                        "Hover over a parameter for an explanation. The parameters are saved on exit."))
        note.setWordWrap(True)
        note.setStyleSheet("color: #555; font-size: 9pt;")
        form.addRow(note)
        sa = QScrollArea()
        sa.setWidget(page)
        sa.setWidgetResizable(True)
        sa.setFrameShape(QFrame.NoFrame)
        self.ltabs.addTab(sa, L("Parameter", "Parameters"))

        # Start / Fortschritt
        hb = QHBoxLayout()
        self.btn_run = QPushButton(L("▶  Auswertung starten", "▶  Start analysis"))
        self.btn_run.setMinimumHeight(38)
        f = self.btn_run.font()
        f.setBold(True)
        f.setPointSize(f.pointSize() + 1)
        self.btn_run.setFont(f)
        self.btn_run.clicked.connect(self.run_analysis)
        self.btn_cancel = QPushButton(L("Abbrechen", "Cancel"))
        self.btn_cancel.setEnabled(False)
        self.btn_cancel.clicked.connect(self.cancel_task)
        hb.addWidget(self.btn_run, 3)
        hb.addWidget(self.btn_cancel, 1)
        lv.addLayout(hb)
        self.bar = QProgressBar()
        self.bar.setRange(0, 100)
        self.bar.setValue(0)
        lv.addWidget(self.bar)
        split.addWidget(left)

        # ---------- rechte Seite
        self.rtabs = QTabWidget()
        split.addWidget(self.rtabs)
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setSizes([480, 1020])
        self.rtabs.addTab(self._build_overview_page(), L("Übersicht / Ausschnitt / Kratzer", "Overview / area / scratches"))

        page = QWidget()
        v = QVBoxLayout(page)
        sp = QSplitter(Qt.Horizontal)
        lw = QWidget()
        lwl = QVBoxLayout(lw)
        lwl.setContentsMargins(0, 0, 0, 0)
        self.lbl_shown = QLabel("")
        self.lbl_shown.setWordWrap(True)
        self.lbl_shown.setStyleSheet("color: #14466B; font-weight: bold;")
        lwl.addWidget(self.lbl_shown)
        lwl.addWidget(QLabel(L("Abbildungen (Haken = exportieren)", "Figures (tick = export)")))
        self.lst = QListWidget()
        self.lst.currentItemChanged.connect(self._show_preview)
        lwl.addWidget(self.lst, 1)
        hb = QHBoxLayout()
        b1, b2 = QPushButton(L("alle", "all")), QPushButton(L("keine", "none"))
        b1.clicked.connect(lambda: self._check_all(True))
        b2.clicked.connect(lambda: self._check_all(False))
        hb.addWidget(b1)
        hb.addWidget(b2)
        lwl.addLayout(hb)
        sp.addWidget(lw)
        pw = QWidget()
        pwl = QVBoxLayout(pw)
        pwl.setContentsMargins(0, 0, 0, 0)
        tb = QHBoxLayout()
        self.lbl_prev = QLabel(L("Vorschau", "Preview"))
        self.lbl_prev.setStyleSheet("font-weight: bold;")
        tb.addWidget(self.lbl_prev, 1)
        for text, fn in [(L("Einpassen", "Fit"), lambda: self.preview.set_fit()), ("100 %", lambda: self.preview.zoom(absolute=1.0)),
                         ("−", lambda: self.preview.zoom(1 / 1.25)), ("+", lambda: self.preview.zoom(1.25))]:
            b = QToolButton()
            b.setText(text)
            b.clicked.connect(fn)
            tb.addWidget(b)
        pwl.addLayout(tb)
        self.preview = Preview()
        pwl.addWidget(self.preview, 1)
        sp.addWidget(pw)
        sp.setSizes([250, 900])
        v.addWidget(sp, 1)

        g = QGroupBox("Export")
        gl = QGridLayout(g)
        gl.addWidget(QLabel(L("Zielordner:", "Target folder:")), 0, 0)
        self.ed_out = QLineEdit()
        self.ed_out.textEdited.connect(lambda _: setattr(self, "export_folder_user", True))
        gl.addWidget(self.ed_out, 0, 1, 1, 5)
        b = QPushButton("…")
        b.setMaximumWidth(36)
        b.clicked.connect(self.choose_folder)
        gl.addWidget(b, 0, 6)
        gl.addWidget(QLabel(L("Formate:", "Formats:")), 1, 0)
        self.fmt = {}
        for k, (key, name) in enumerate({"png": "PNG", "pdf": "PDF", "svg": "SVG", "tiff": "TIFF"}.items()):
            cb = QCheckBox(name)
            cb.setChecked(key == "png")
            self.fmt[key] = cb
            gl.addWidget(cb, 1, 1 + k)
        gl.addWidget(QLabel(L("Auflösung:", "Resolution:")), 2, 0)
        self.sp_dpi = QSpinBox()
        self.sp_dpi.setRange(72, 1200)
        self.sp_dpi.setSingleStep(50)
        self.sp_dpi.setValue(300)
        self.sp_dpi.setSuffix(" dpi")
        gl.addWidget(self.sp_dpi, 2, 1)
        self.chk_raster = QCheckBox(L("Karten in PDF/SVG als Rasterbild (kleinere Dateien)", "Maps in PDF/SVG as raster image (smaller files)"))
        self.chk_raster.setChecked(True)
        gl.addWidget(self.chk_raster, 2, 2, 1, 4)
        self.btn_export = QPushButton(L("Ausgewählte Abbildungen exportieren", "Export selected figures"))
        self.btn_export.setEnabled(False)
        self.btn_export.clicked.connect(self.export)
        self.btn_folder = QPushButton(L("Ordner öffnen", "Open folder"))
        self.btn_folder.clicked.connect(self.open_folder)
        gl.addWidget(self.btn_export, 3, 1, 1, 4)
        gl.addWidget(self.btn_folder, 3, 5, 1, 2)
        v.addWidget(g)
        self.rtabs.addTab(page, L("Ergebnisse", "Results"))
        self.rtabs.addTab(self._build_datasets_page(), L("Datensätze / Vergleich", "Datasets / Comparison"))

        page = QWidget()
        v = QVBoxLayout(page)
        self.txt_kenn = QPlainTextEdit()
        self.txt_kenn.setReadOnly(True)
        self.txt_kenn.setFont(QFont("Consolas", 9))
        self.txt_kenn.setLineWrapMode(QPlainTextEdit.NoWrap)
        v.addWidget(self.txt_kenn, 1)
        hb = QHBoxLayout()
        b = QPushButton(L("In Zwischenablage kopieren", "Copy to clipboard"))
        b.clicked.connect(lambda: QApplication.clipboard().setText(self.txt_kenn.toPlainText()))
        hb.addWidget(b)
        b = QPushButton(L("Speichern …", "Save …"))
        b.clicked.connect(self.save_kennwerte)
        hb.addWidget(b)
        hb.addStretch(1)
        v.addLayout(hb)
        self.rtabs.addTab(page, L("Kennwerte", "Key figures"))

        self.txt_log = QPlainTextEdit()
        self.txt_log.setReadOnly(True)
        self.txt_log.setFont(QFont("Consolas", 9))
        self.rtabs.addTab(self.txt_log, L("Protokoll", "Log"))

        self.setStatusBar(QStatusBar())
        self.statusBar().showMessage(L("Bereit – bitte eine ctf-Datei öffnen.", "Ready – please open a ctf file."))
        self._set_enabled_for_project(False)
        self.resize(1500, 950)

    # ---------------------------------------------------------------- Übersicht (Aufbau)
    def _build_overview_page(self):
        page = QWidget()
        v = QVBoxLayout(page)
        self.overview = OverviewCanvas(page)
        self.overview.rect_selected.connect(self._crop_from_mouse)
        self.overview.stroke_drawn.connect(self._stroke_drawn)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel(L("Werkzeug:", "Tool:")))
        self.mode_group = QButtonGroup(self)
        self.mode_btns = {}
        for mode, text, tip in [
                ("crop", L("▭ Ausschnitt", "▭ Area"), L("Rechteck aufziehen: Ausschnitt für die Auswertung wählen",
                                                       "Drag a rectangle: select the area to be analysed")),
                ("add", L("✎ Kratzer markieren", "✎ Mark scratch"),
                 L("Linie entlang eines Kratzers ziehen: die Punkte gelten als Kratzer und werden aus den Nachbarn ergänzt",
                   "Draw a line along a scratch: the points count as scratch and are filled in from the neighbours")),
                ("remove", L("⌫ Markierung aufheben", "⌫ Unmark"),
                 L("Linie über einen fälschlich erkannten Kratzer ziehen: dort wird jede Kratzermarkierung aufgehoben",
                   "Draw a line over a falsely detected scratch: every scratch marking is removed there"))]:
            b = QToolButton()
            b.setText(text)
            b.setToolTip(tip)
            b.setCheckable(True)
            b.setToolButtonStyle(Qt.ToolButtonTextOnly)
            self.mode_group.addButton(b)
            self.mode_btns[mode] = b
            b.clicked.connect(lambda _=False, m=mode: self.set_tool(m))
            row1.addWidget(b)
        self.mode_btns["crop"].setChecked(True)
        row1.addSpacing(10)
        row1.addWidget(QLabel(L("Breite:", "Width:")))
        self.sp_width = QDoubleSpinBox()
        self.sp_width.setDecimals(1)
        self.sp_width.setRange(0.1, 200)
        self.sp_width.setSingleStep(0.5)
        self.sp_width.setValue(2.0)
        self.sp_width.setSuffix(" µm")
        self.sp_width.setToolTip(L("Breite der gezeichneten Linie", "Width of the drawn line"))
        self.sp_width.valueChanged.connect(lambda val: setattr(self.overview, "mark_width", float(val)))
        row1.addWidget(self.sp_width)
        row1.addSpacing(10)
        self.btn_undo = QToolButton()
        self.btn_undo.setText(L("↶ Rückgängig", "↶ Undo"))
        self.btn_undo.setToolTip(L("letzte Linie entfernen (Strg+Z)", "remove the last line (Ctrl+Z)"))
        self.btn_undo.clicked.connect(self.undo_stroke)
        row1.addWidget(self.btn_undo)
        self.btn_clear = QToolButton()
        self.btn_clear.setText(L("Alle löschen", "Clear all"))
        self.btn_clear.setToolTip(L("alle manuellen Kratzermarkierungen dieser Datei entfernen",
                                    "remove all manual scratch marks of this file"))
        self.btn_clear.clicked.connect(self.clear_strokes)
        row1.addWidget(self.btn_clear)
        row1.addStretch(1)
        v.addLayout(row1)

        row2 = QHBoxLayout()
        self.navbar = NavBar(self.overview, page)
        self.navbar.setIconSize(self.navbar.iconSize() * 0.85)
        row2.addWidget(self.navbar, 1)
        self.btn_detect = QPushButton(L("Automatische Erkennung anzeigen", "Show automatic detection"))
        self.btn_detect.setToolTip(L("nur die automatische Kratzererkennung ausführen (Gesamtkarte bzw. gewählter Ausschnitt, "
                                     "Parameter aus dem Reiter „Parameter“) und rot einblenden",
                                     "run only the automatic scratch detection (whole map or selected area, parameters "
                                     "from the tab “Parameters”) and show it in red"))
        self.btn_detect.clicked.connect(self.detect_scratches_preview)
        row2.addWidget(self.btn_detect)
        self.chk_show_marks = QCheckBox(L("Markierungen einblenden", "Show marks"))
        self.chk_show_marks.setChecked(True)
        self.chk_show_marks.toggled.connect(lambda _: self.update_overlay())
        row2.addWidget(self.chk_show_marks)
        v.addLayout(row2)
        v.addWidget(self.overview, 1)
        self.lbl_ov = QLabel("")
        self.lbl_ov.setWordWrap(True)
        self.lbl_ov.setStyleSheet("color: #555;")
        v.addWidget(self.lbl_ov)
        self._update_ov_hint()
        return page

    def _update_ov_hint(self):
        mode = self.overview.mode
        txt = {"crop": L("Bandkontrast der Gesamtkarte. Mit gedrückter linker Maustaste ein Rechteck aufziehen, um einen "
                         "Ausschnitt zu wählen; das Rechteck lässt sich an den Ecken anpassen.",
                         "Band contrast of the whole map. Drag a rectangle with the left mouse button to select an area; "
                         "the rectangle can be adjusted at its corners."),
               "add": L("Kratzer markieren: mit gedrückter linker Maustaste eine Linie entlang des Kratzers ziehen. "
                        "Zum genauen Zeichnen mit der Lupe hineinzoomen (danach Lupe wieder abwählen).",
                        "Mark scratch: draw a line along the scratch with the left mouse button held down. "
                        "Zoom in with the magnifier for precise drawing (deselect the magnifier afterwards)."),
               "remove": L("Markierung aufheben: mit gedrückter linker Maustaste eine Linie über die falsche Markierung ziehen.",
                           "Unmark: draw a line over the wrong marking with the left mouse button held down.")}[mode]
        legend = L("  Farben: rot = automatisch erkannt, magenta = manuell markiert, hellblau = aufgehoben.",
                   "  Colours: red = detected automatically, magenta = marked manually, light blue = unmarked.")
        self.lbl_ov.setText(txt + ("" if mode == "crop" else legend))

    # ---------------------------------------------------------------- Kratzer manuell markieren (Logik)
    def set_tool(self, mode):
        if self.navbar.mode:                       # Lupe/Verschieben der matplotlib-Leiste beenden
            (self.navbar.zoom if str(self.navbar.mode) == "zoom rect" else self.navbar.pan)()
        self.overview.set_mode(mode)
        self.mode_btns[mode].setChecked(True)
        self._update_ov_hint()
        if mode != "crop" and self.project is not None and self.rtabs.currentIndex() != 0:
            self.rtabs.setCurrentIndex(0)

    def _stroke_drawn(self, x0, y0, x1, y1):
        if self.project is None or self.overview.mode not in ("add", "remove"):
            return
        self.strokes.append({"modus": self.overview.mode, "x0": round(x0, 3), "y0": round(y0, 3), "x1": round(x1, 3),
                             "y1": round(y1, 3), "breite": round(self.sp_width.value(), 3)})
        self._strokes_changed()

    def undo_stroke(self):
        if self.strokes:
            self.strokes.pop()
            self._strokes_changed()

    def clear_strokes(self):
        if not self.strokes:
            return
        if QMessageBox.question(self, app_title(), L(f"Alle {len(self.strokes)} manuellen Kratzermarkierungen dieser Datei löschen?",
                                                     f"Delete all {len(self.strokes)} manual scratch marks of this file?")) != QMessageBox.Yes:
            return
        self.strokes = []
        self._strokes_changed()

    def _strokes_changed(self, save=True):
        if save and self.project is not None:
            try:
                self.marks.set(self.project.path, self.strokes)
            except OSError as e:
                self.log(L(f"Kratzermarkierungen konnten nicht gespeichert werden: {e}", f"Scratch marks could not be saved: {e}"))
        self.update_overlay()
        self._update_marks_label()

    def _update_marks_label(self):
        n_add = sum(1 for s in self.strokes if s["modus"] == "add")
        n_rem = len(self.strokes) - n_add
        if self.strokes:
            t = L(f"manuell: {n_add} Linie(n) markiert, {n_rem} aufgehoben (Reiter „Übersicht“)",
                  f"manual: {n_add} line(s) marked, {n_rem} unmarked (tab “Overview”)")
        else:
            t = L("manuelle Markierungen: keine (Reiter „Übersicht“ → Kratzer markieren)",
                  "manual marks: none (tab “Overview” → Mark scratch)")
        self.lbl_marks.setText(t)
        for w in (self.btn_undo, self.btn_clear, self.act_undo, self.act_clear):
            w.setEnabled(bool(self.strokes))

    def update_overlay(self):
        p = self.project
        if p is None:
            return
        rgba = np.zeros((p.ny, p.nx, 4), np.uint8)
        if self.chk_show_marks.isChecked():
            auto = self.auto_scratch if self.auto_scratch is not None else np.zeros((p.ny, p.nx), bool)
            add, rem = ebsd_marks.masks(self.strokes, (p.ny, p.nx), p.step) if self.strokes else (None, None)
            rgba[auto] = OV_AUTO
            if rem is not None:
                rgba[rem & ~auto] = OV_REM
                rgba[rem & auto] = OV_REM_AUTO
                rgba[add] = OV_ADD
        self.overview.set_overlay(rgba)

    def detect_scratches_preview(self):
        if self.project is None:
            return
        proj, s, roles = self.project, self.current_settings(), self.roles()
        crop = self.crop_rect() if self.chk_crop.isChecked() else None
        self.log(L("Automatische Kratzererkennung (Vorschau) …", "Automatic scratch detection (preview) …"))

        def job(task):
            c = load_core()
            return c.preview_scratches(proj, s, roles, crop_um=crop, cancel=task.check)

        self.bar.setRange(0, 0)
        self._start(job, self._detect_done)

    def _detect_done(self, out):
        mask, n = out
        self.bar.setRange(0, 100)
        self.bar.setValue(100)
        self.auto_scratch = mask
        self.auto_scratch_src = L("Vorschau", "preview")
        self.update_overlay()
        msg = L(f"Automatische Erkennung: {n} Kratzer, {dez(mask.mean() * 100, 2)} % der Fläche (rot eingeblendet)",
                f"Automatic detection: {n} scratch(es), {dez(mask.mean() * 100, 2)} % of the area (shown in red)")
        self.log(msg)
        self.statusBar().showMessage(msg)
        self.rtabs.setCurrentIndex(0)

    def _restore_strokes(self, proj):
        """Markierungen der Datei: zuerst aus der Markierungsablage, sonst aus dem jüngsten gespeicherten Datensatz"""
        if self.marks.has(proj.path):
            return self.marks.get(proj.path), L("gespeicherte Markierungen", "saved marks")
        key = str(proj.path).lower()
        best = None
        for ds in self.store.list():
            if str(ds.path).lower() == key and ds.strokes:
                if best is None or ds.meta.get("geaendert", 0) > best.meta.get("geaendert", 0):
                    best = ds
        if best is not None:
            return ebsd_marks.clean_strokes(best.strokes), L(f"Datensatz „{best.label}“", f"dataset “{best.label}”")
        return [], ""

    # ---------------------------------------------------------------- Datensätze / Vergleich (Aufbau)
    def _build_datasets_page(self):
        page = QWidget()
        v = QVBoxLayout(page)
        sp = QSplitter(Qt.Horizontal)

        lw = QWidget()
        lv = QVBoxLayout(lw)
        lv.setContentsMargins(0, 0, 0, 0)
        lab = QLabel(L("Gespeicherte Datensätze (Haken = im Vergleich)", "Stored datasets (tick = in comparison)"))
        lab.setStyleSheet("font-weight: bold;")
        lv.addWidget(lab)
        self.ds_list = QListWidget()
        self.ds_list.itemChanged.connect(self._ds_item_changed)
        self.ds_list.currentItemChanged.connect(self._ds_current_changed)
        self.ds_list.itemDoubleClicked.connect(lambda it: self.show_dataset(it.data(Qt.UserRole)))
        lv.addWidget(self.ds_list, 1)
        self.lbl_ds_info = QLabel("")
        self.lbl_ds_info.setWordWrap(True)
        self.lbl_ds_info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.lbl_ds_info.setStyleSheet("color: #555; font-size: 9pt;")
        lv.addWidget(self.lbl_ds_info)
        grid = QGridLayout()
        btns = [(L("Ergebnisse anzeigen", "Show results"), self._ds_show, 0, 0,
                 L("Vorschauen und Kennwerte des Datensatzes in den Reitern „Ergebnisse“ und „Kennwerte“ anzeigen (Doppelklick)",
                   "Show previews and key figures of the dataset in the tabs “Results” and “Key figures” (double-click)")),
                (L("Neu berechnen", "Recalculate"), self._ds_recompute, 0, 1,
                 L("ctf-Datei mit den gespeicherten Rollen, Parametern, Kratzermarkierungen und dem Ausschnitt laden und erneut auswerten",
                   "Load the ctf file with the stored roles, parameters, scratch marks and area and analyse it again")),
                (L("Umbenennen …", "Rename …"), self._ds_rename, 1, 0, L("Bezeichnung in den Vergleichsabbildungen",
                                                                          "Label in the comparison figures")),
                (L("Entfernen", "Remove"), self._ds_remove, 1, 1, L("gespeicherten Datensatz löschen (die ctf-Datei bleibt erhalten)",
                                                                    "delete the stored dataset (the ctf file is kept)")),
                (L("▲ nach oben", "▲ up"), lambda: self._ds_move(-1), 2, 0, L("Reihenfolge im Vergleich", "Order in the comparison")),
                (L("▼ nach unten", "▼ down"), lambda: self._ds_move(1), 2, 1, L("Reihenfolge im Vergleich", "Order in the comparison")),
                (L("alle", "all"), lambda: self._ds_check_all(True), 3, 0, L("alle Datensätze in den Vergleich", "all datasets in the comparison")),
                (L("keine", "none"), lambda: self._ds_check_all(False), 3, 1, L("keinen Datensatz im Vergleich", "no dataset in the comparison"))]
        self.ds_buttons = []
        for text, fn, r, c, tip in btns:
            b = QPushButton(text)
            b.setToolTip(tip)
            b.clicked.connect(fn)
            grid.addWidget(b, r, c)
            self.ds_buttons.append(b)
        lv.addLayout(grid)
        sp.addWidget(lw)

        rw = QWidget()
        rv = QVBoxLayout(rw)
        rv.setContentsMargins(0, 0, 0, 0)
        tb = QHBoxLayout()
        tb.addWidget(QLabel(L("Vergleichsabbildung:", "Comparison figure:")))
        self.cmp_sel = QComboBox()
        for fid, title, _ in compare_figures():
            self.cmp_sel.addItem(f"{fid}  {title}", fid)
        self.cmp_sel.currentIndexChanged.connect(lambda _: self.refresh_compare())
        tb.addWidget(self.cmp_sel, 1)
        for text, fn in [(L("Einpassen", "Fit"), lambda: self.cmp_preview.set_fit()), ("100 %", lambda: self.cmp_preview.zoom(absolute=1.0)),
                         ("−", lambda: self.cmp_preview.zoom(1 / 1.25)), ("+", lambda: self.cmp_preview.zoom(1.25))]:
            b = QToolButton()
            b.setText(text)
            b.clicked.connect(fn)
            tb.addWidget(b)
        rv.addLayout(tb)
        vs = QSplitter(Qt.Vertical)
        self.cmp_preview = Preview()
        self.cmp_preview.label.setText(self._no_ds_text())
        vs.addWidget(self.cmp_preview)
        self.cmp_table = QTableWidget(0, 0)
        self.cmp_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.cmp_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.cmp_table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        vs.addWidget(self.cmp_table)
        vs.setSizes([560, 260])
        rv.addWidget(vs, 1)
        hb = QHBoxLayout()
        b = QPushButton(L("Vergleich exportieren …", "Export comparison …"))
        b.setToolTip(L("V1, V2 in den unter „Ergebnisse“ gewählten Formaten und der dort gewählten Auflösung, "
                       f"dazu {dsets.compare_csv_name()} (Semikolon, Dezimalkomma)",
                       "V1, V2 in the formats and resolution selected under “Results”, "
                       f"plus {dsets.compare_csv_name()} (comma-separated, decimal point)"))
        b.clicked.connect(self.export_compare)
        self.btn_cmp_export = b
        hb.addWidget(b, 2)
        b = QPushButton(L("Tabelle kopieren", "Copy table"))
        b.clicked.connect(lambda: QApplication.clipboard().setText(dsets.metrics_csv(self.compare_datasets())))
        hb.addWidget(b, 1)
        rv.addLayout(hb)
        sp.addWidget(rw)
        sp.setSizes([330, 900])
        v.addWidget(sp, 1)
        return page

    @staticmethod
    def _no_ds_text():
        return L("Noch keine gespeicherten Datensätze.\n\nNach jeder Auswertung wird der Datensatz hier abgelegt und kann mit "
                 "anderen verglichen werden.",
                 "No stored datasets yet.\n\nAfter each analysis the dataset is stored here and can be compared with others.")

    # ---------------------------------------------------------------- Datensätze / Vergleich (Logik)
    def refresh_datasets(self, select=None):
        cur = select or (self.ds_list.currentItem().data(Qt.UserRole) if self.ds_list.currentItem() else None)
        self._ds_mute = True
        self.ds_list.clear()
        sel_row = 0
        for k, ds in enumerate(self.store.list()):
            it = QListWidgetItem(ds.label)
            it.setData(Qt.UserRole, ds.id)
            it.setToolTip(ds.info())
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if ds.in_compare else Qt.Unchecked)
            self.ds_list.addItem(it)
            if ds.id == cur:
                sel_row = k
        self._ds_mute = False
        if self.ds_list.count():
            self.ds_list.setCurrentRow(sel_row)
        self._ds_current_changed(self.ds_list.currentItem())
        self.refresh_compare()

    def _ds_selected(self):
        it = self.ds_list.currentItem()
        return self.store.get(it.data(Qt.UserRole)) if it else None

    def compare_datasets(self):
        out = []
        for i in range(self.ds_list.count()):
            it = self.ds_list.item(i)
            if it.checkState() == Qt.Checked:
                ds = self.store.get(it.data(Qt.UserRole))
                if ds is not None:
                    out.append(ds)
        return out

    def _ds_item_changed(self, it):
        if self._ds_mute:
            return
        self.store.set_compare(it.data(Qt.UserRole), it.checkState() == Qt.Checked)
        self.refresh_compare()

    def _ds_current_changed(self, cur, _prev=None):
        ds = self.store.get(cur.data(Qt.UserRole)) if cur else None
        self.lbl_ds_info.setText(ds.info() if ds else "")
        for b in self.ds_buttons[:6]:
            b.setEnabled(ds is not None)

    def _ds_check_all(self, on):
        self._ds_mute = True
        for i in range(self.ds_list.count()):
            it = self.ds_list.item(i)
            it.setCheckState(Qt.Checked if on else Qt.Unchecked)
            self.store.set_compare(it.data(Qt.UserRole), on)
        self._ds_mute = False
        self.refresh_compare()

    def _ds_move(self, d):
        r = self.ds_list.currentRow()
        if r < 0 or not (0 <= r + d < self.ds_list.count()):
            return
        ids = [self.ds_list.item(i).data(Qt.UserRole) for i in range(self.ds_list.count())]
        ids[r], ids[r + d] = ids[r + d], ids[r]
        self.store.set_order(ids)
        self.refresh_datasets(select=ids[r + d])

    def _ds_rename(self):
        ds = self._ds_selected()
        if ds is None:
            return
        text, ok = QInputDialog.getText(self, app_title(), L("Bezeichnung des Datensatzes (erscheint in den Vergleichsabbildungen):",
                                                             "Label of the dataset (appears in the comparison figures):"),
                                        QLineEdit.Normal, ds.label)
        if ok:
            self.store.rename(ds.id, text)
            self.refresh_datasets(select=ds.id)
            if self.stored_view is not None and self.stored_view.id == ds.id:
                self.show_dataset(ds.id)

    def _ds_remove(self):
        ds = self._ds_selected()
        if ds is None:
            return
        if QMessageBox.question(self, app_title(), L(f"Gespeicherten Datensatz „{ds.label}“ entfernen?\n\n"
                                                     "Die ctf-Datei und bereits exportierte Dateien bleiben erhalten.",
                                                     f"Remove the stored dataset “{ds.label}”?\n\n"
                                                     "The ctf file and files already exported are kept.")) != QMessageBox.Yes:
            return
        self.store.remove(ds.id)
        if self.live_ds == ds.id:
            self.live_ds = None
        if self.stored_view is not None and self.stored_view.id == ds.id:
            self.stored_view = None
        self.refresh_datasets()

    def _ds_show(self):
        ds = self._ds_selected()
        if ds is not None:
            self.show_dataset(ds.id)

    def _ds_recompute(self):
        ds = self._ds_selected()
        if ds is not None:
            self.recompute(ds)

    def show_dataset(self, ds_id):
        """Ergebnisse eines gespeicherten Datensatzes anzeigen (ohne Neuberechnung)"""
        ds = self.store.get(ds_id)
        if ds is None:
            return
        if ds_id == self.live_ds and self.result is not None:
            self.stored_view = None
            self._fill_results(self.result.available, self.result.kennwerte)
        else:
            self.stored_view = ds
            self.previews = {}
            self._fill_results(ds.figures, ds.kennwerte())
            if not self.export_folder_user:
                p = ds.path
                name = p.stem.replace(" ", "_") + L("_Auswertung", "_analysis")
                if ds.crop_um:
                    name += L("_Ausschnitt_", "_crop_") + "_".join(f"{v:g}" for v in ds.crop_um)
                self.ed_out.setText(str(p.with_name(name)))
        self._update_shown_label()
        self.btn_export.setEnabled(True)
        self.rtabs.setCurrentIndex(1)

    def _update_shown_label(self):
        if self.stored_view is not None:
            self.lbl_shown.setText(L(f"Datensatz: {self.stored_view.label}  (gespeichert)", f"Dataset: {self.stored_view.label}  (stored)"))
        elif self.live_ds is not None and self.result is not None:
            ds = self.store.get(self.live_ds)
            self.lbl_shown.setText(L("Datensatz: ", "Dataset: ") + (ds.label if ds else self.result.project.path.name))
        else:
            self.lbl_shown.setText("")

    def _fill_results(self, figures, kennwerte):
        self.lst.blockSignals(True)
        self.lst.clear()
        for fid, title, _ in figures:
            it = QListWidgetItem(f"{fid}  {fig_title(fid, title)}")
            it.setData(Qt.UserRole, fid)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked)
            self.lst.addItem(it)
        self.lst.blockSignals(False)
        self.txt_kenn.setPlainText(kennwerte)
        if self.lst.count():
            self.lst.setCurrentRow(0)
            self._show_preview(self.lst.currentItem())

    def recompute(self, ds, export=None):
        path = ds.path
        if not path.exists():
            QMessageBox.warning(self, app_title(), L(f"Die ctf-Datei wurde nicht gefunden:\n{path}\n\n"
                                                     "Gespeicherte Vorschauen und Kennwerte bleiben verfügbar.",
                                                     f"The ctf file was not found:\n{path}\n\n"
                                                     "Stored previews and key figures remain available."))
            return
        k = ds.metrics
        self._pending = {"settings": Settings.from_dict(k.get("einstellungen", {})), "roles": k.get("rollen", {}),
                         "crop": k.get("ausschnitt_um"), "export": export, "label": ds.label,
                         "strokes": ebsd_marks.clean_strokes(k.get("kratzer_manuell"))}
        self.log(L(f"Neuberechnung des gespeicherten Datensatzes „{ds.label}“", f"Recalculating the stored dataset “{ds.label}”"))
        self.open_file(str(path))

    def _apply_pending(self):
        p = self._pending
        self.apply_settings(p["settings"])
        names = self.project.names
        for i, cb in self.role_boxes.items():
            if names[i] in p["roles"]:
                self._set_role(cb, p["roles"][names[i]])
        if p["crop"]:
            self._crop_from_mouse(*p["crop"])
            self.overview.set_rect(self.crop_rect())
        self.strokes = list(p["strokes"])
        self._strokes_changed()
        QTimer.singleShot(200, self.run_analysis)

    def _store_result(self, res, previews):
        try:
            label = self._pending.get("label") if self._pending else None
            ds = self.store.save_result(res, previews, label=label)
            self.live_ds = ds.id
            self.stored_view = None
            self.refresh_datasets(select=ds.id)
            self.log(L(f"Datensatz gespeichert: „{ds.label}“ ({ds.folder})", f"Dataset stored: “{ds.label}” ({ds.folder})"))
        except Exception as e:  # noqa: BLE001
            self.log(L("Datensatz konnte nicht gespeichert werden: ", "Dataset could not be stored: ")
                     + f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def refresh_compare(self):
        if not hasattr(self, "cmp_preview"):
            return
        dss = self.compare_datasets()
        rows = dsets.metrics_table(dss) if dss else []
        self.cmp_table.clear()
        self.cmp_table.setRowCount(len(rows))
        self.cmp_table.setColumnCount(len(dss))
        self.cmp_table.setHorizontalHeaderLabels([d.label for d in dss])
        self.cmp_table.setVerticalHeaderLabels([r[0] for r in rows])
        for r, (_, vals) in enumerate(rows):
            for c, t in enumerate(vals):
                it = QTableWidgetItem(t)
                it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.cmp_table.setItem(r, c, it)
        self.btn_cmp_export.setEnabled(bool(dss))
        if not dss:
            self.cmp_preview.pix = None
            self.cmp_preview.label.setPixmap(QPixmap())
            self.cmp_preview.label.setText(L("Keine Datensätze im Vergleich.\n\nLinks Datensätze anhaken.",
                                             "No datasets in the comparison.\n\nTick datasets on the left.")
                                           if self.ds_list.count() else self._no_ds_text())
            return
        try:
            b = dsets.compare_png(self.cmp_sel.currentData(), dss, dpi=100)
        except Exception as e:  # noqa: BLE001
            self.log(L("Vergleichsabbildung fehlgeschlagen: ", "Comparison figure failed: ") + f"{e}\n{traceback.format_exc()}")
            return
        pix = QPixmap()
        pix.loadFromData(b, "PNG")
        self.cmp_preview.fit = True
        self.cmp_preview.set_pixmap(pix)

    def export_compare(self):
        dss = self.compare_datasets()
        if not dss:
            return
        formats = [k for k, cb in self.fmt.items() if cb.isChecked()] or ["png"]
        start = self.qs.value("compare_dir", self.qs.value("last_dir", str(Path.home())))
        folder = QFileDialog.getExistingDirectory(self, L("Zielordner für den Vergleich", "Target folder for the comparison"), start)
        if not folder:
            return
        self.qs.setValue("compare_dir", folder)
        try:
            written = dsets.export_compare(dss, folder, formats, dpi=self.sp_dpi.value())
        except Exception as e:  # noqa: BLE001
            QMessageBox.critical(self, app_title(), L("Export fehlgeschlagen:\n", "Export failed:\n") + str(e))
            return
        self.log(L(f"Vergleich exportiert ({len(dss)} Datensätze): {len(written)} Dateien → {folder}",
                   f"Comparison exported ({len(dss)} datasets): {len(written)} files → {folder}"))
        self.statusBar().showMessage(L(f"Vergleich: {len(written)} Dateien geschrieben nach {folder}",
                                       f"Comparison: {len(written)} files written to {folder}"))
        self._written_box(len(written), folder)

    def _written_box(self, n, folder):
        box = QMessageBox(self)
        box.setWindowTitle(app_title())
        box.setText(L(f"{n} Dateien exportiert nach\n{folder}", f"{n} files exported to\n{folder}"))
        b_open = box.addButton(L("Ordner öffnen", "Open folder"), QMessageBox.ActionRole)
        box.addButton(QMessageBox.Ok)
        box.exec()
        if box.clickedButton() == b_open:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    # ---------------------------------------------------------------- Zustand
    def _set_enabled_for_project(self, on):
        self.btn_run.setEnabled(on)
        self.chk_crop.setEnabled(on)
        for b in self.mode_btns.values():
            b.setEnabled(on)
        self.btn_detect.setEnabled(on)
        self.sp_width.setEnabled(on)
        self._update_marks_label()
        if not on:
            self.btn_undo.setEnabled(False)
            self.btn_clear.setEnabled(False)

    def _busy(self, busy):
        for w in [self.btn_open, self.btn_run, self.btn_export, self.ltabs, self.btn_detect]:
            w.setEnabled(not busy)
        self.btn_cancel.setEnabled(busy)
        for b in getattr(self, "ds_buttons", [])[:4]:
            b.setEnabled(not busy)
        if not busy:
            self.btn_run.setEnabled(self.project is not None)
            self.btn_detect.setEnabled(self.project is not None)
            self.btn_export.setEnabled(self.result is not None or self.stored_view is not None)
            self._ds_current_changed(self.ds_list.currentItem())

    def log(self, msg):
        self.txt_log.appendPlainText(msg)

    def current_settings(self):
        s = Settings()
        s.title = self.ed_title.text()
        for k, cb in self.chk.items():
            setattr(s, k, cb.isChecked())
        for k, sp in self.param.items():
            setattr(s, k, sp.value())
        return s

    def apply_settings(self, s: Settings):
        self.ed_title.setText(s.title)
        for k, cb in self.chk.items():
            cb.setChecked(bool(getattr(s, k)))
        for k, sp in self.param.items():
            sp.setValue(getattr(s, k))

    def reset_params(self):
        d = Settings()
        for k, sp in self.param.items():
            sp.setValue(getattr(d, k))

    def _restore(self):
        try:
            s = Settings.from_dict(json.loads(self.qs.value("settings", "{}")))
            self.apply_settings(s)
        except Exception:  # noqa: BLE001
            self.apply_settings(Settings())
        fm = self.qs.value("formats", "")
        if fm:
            for k, cb in self.fmt.items():
                cb.setChecked(k in str(fm).split(","))
        try:
            self.sp_dpi.setValue(int(self.qs.value("dpi", 300)))
        except (TypeError, ValueError):
            pass
        try:
            self.sp_width.setValue(float(self.qs.value("mark_width", 2.0)))
        except (TypeError, ValueError):
            pass
        geo = self.qs.value("geometry")
        if geo is not None:
            self.restoreGeometry(geo)

    def closeEvent(self, e):
        if self.task is not None and self.task.isRunning():
            if QMessageBox.question(self, app_title(), L("Eine Berechnung läuft noch. Trotzdem beenden?",
                                                         "A calculation is still running. Exit anyway?")) != QMessageBox.Yes:
                self._restart = False
                e.ignore()
                return
            self.task.cancel()
            self.task.wait(3000)
        self.qs.setValue("settings", json.dumps(self.current_settings().to_dict()))
        self.qs.setValue("formats", ",".join(k for k, cb in self.fmt.items() if cb.isChecked()))
        self.qs.setValue("dpi", self.sp_dpi.value())
        self.qs.setValue("mark_width", self.sp_width.value())
        self.qs.setValue("geometry", self.saveGeometry())
        self.qs.sync()
        super().closeEvent(e)
        if self._restart:
            restart_program()

    # ---------------------------------------------------------------- Sprache
    def change_language(self, code):
        if code == get_language():
            return
        self.qs.setValue("language", code)
        self.qs.sync()
        box = QMessageBox(self)
        box.setWindowTitle(app_title())
        box.setText("Die Sprache wird nach einem Neustart des Programms übernommen.\n"
                    "The language will be applied after restarting the program.")
        b_now = box.addButton("Jetzt neu starten / Restart now", QMessageBox.AcceptRole)
        box.addButton("Später / Later", QMessageBox.RejectRole)
        box.exec()
        if box.clickedButton() == b_now:
            self._restart = True
            self.close()

    # ---------------------------------------------------------------- Datei
    def open_file(self, path=None):
        if not path:
            start = self.qs.value("last_dir", str(Path.home()))
            path, _ = QFileDialog.getOpenFileName(self, L("ctf-Datei öffnen", "Open ctf file"), start,
                                                  L("Oxford Channel Text (*.ctf);;Alle Dateien (*)", "Oxford Channel Text (*.ctf);;All files (*)"))
            if not path:
                return
        self.qs.setValue("last_dir", str(Path(path).parent))
        self.log(L(f"Lade {path} …", f"Loading {path} …"))
        self.statusBar().showMessage(L("Lade Datei (beim ersten Mal werden zusätzlich die Bibliotheken geladen) …",
                                       "Loading file (the first time, the libraries are loaded as well) …"))
        self.bar.setRange(0, 0)

        def job(task):
            c = load_core()
            return c.EBSDProject.load(path, log=task.message.emit)

        self._start(job, self._file_loaded)

    def _file_loaded(self, proj):
        self.bar.setRange(0, 100)
        self.bar.setValue(0)
        self.project = proj
        self.result = None
        self.live_ds = None
        self.stored_view = None
        self.lbl_shown.setText("")
        self.previews = {}
        self.lst.clear()
        self.preview.pix = None
        self.preview.label.setPixmap(QPixmap())
        self.preview.label.setText(L("Noch keine Ergebnisse.\n\nAuswertung starten.", "No results yet.\n\nStart the analysis."))
        self.txt_kenn.clear()
        self.export_folder_user = False
        self.ed_out.setText("")
        self.lbl_file.setText(f"<b>{proj.path.name}</b><br>"
                              + L(f"{proj.nx} × {proj.ny} Punkte, Schrittweite {fmt(proj.step, 4)} µm",
                                  f"{proj.nx} × {proj.ny} points, step size {fmt(proj.step, 4)} µm")
                              + f"<br>{fmt(proj.width_um, 1)} × {fmt(proj.height_um, 1)} µm")
        roles = proj.default_roles()
        rows = [p for p in proj.phases if p["index"] > 0]
        self.tbl.setRowCount(len(rows))
        self.role_boxes = {}
        for r, p in enumerate(rows):
            items = [p["name"], fmt(p["frac"], 2), fmt(p["mad"], 2), fmt(p["bc"], 0)]
            for c, t in enumerate(items):
                it = QTableWidgetItem(t)
                it.setToolTip(L(f"{p['name']} – Punktgruppe {p['symmetry']}, Anteil {fmt(p['frac'], 2)} %, "
                                f"MAD-Median {fmt(p['mad'], 2)}°, BC-Median {fmt(p['bc'], 0)}",
                                f"{p['name']} – point group {p['symmetry']}, fraction {fmt(p['frac'], 2)} %, "
                                f"MAD median {fmt(p['mad'], 2)}°, BC median {fmt(p['bc'], 0)}"))
                if c:
                    it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.tbl.setItem(r, c, it)
            cb = QComboBox()
            for role in ROLES:
                cb.addItem(role_name(role), role)
            self._set_role(cb, roles.get(p["index"], "Fehlindizierung"))
            self.tbl.setCellWidget(r, 4, cb)
            self.role_boxes[p["index"]] = cb
        hgt = self.tbl.horizontalHeader().height() + sum(self.tbl.rowHeight(r) for r in range(len(rows))) + 6
        self.tbl.setFixedHeight(max(90, min(hgt, 260)))
        z = proj.phases[0]
        self.lbl_zero.setText(L(f"Nulllösungen: {fmt(z['frac'], 2)} % der Fläche", f"Zero solutions: {fmt(z['frac'], 2)} % of the area"))
        for sp in self.crop_spins:
            sp.blockSignals(True)
        self.crop_spins[0].setRange(0, proj.width_um)
        self.crop_spins[2].setRange(0, proj.width_um)
        self.crop_spins[1].setRange(0, proj.height_um)
        self.crop_spins[3].setRange(0, proj.height_um)
        self.crop_spins[0].setValue(0)
        self.crop_spins[1].setValue(0)
        self.crop_spins[2].setValue(proj.width_um)
        self.crop_spins[3].setValue(proj.height_um)
        for sp in self.crop_spins:
            sp.blockSignals(False)
        self.chk_crop.setChecked(False)
        self.overview.show_map(proj.bc, proj.width_um, proj.height_um,
                               f"{proj.path.name} – " + L("Bandkontrast", "band contrast"))
        self.navbar.update()                          # Ansichtsverlauf der neuen Karte
        self.auto_scratch = None
        self.auto_scratch_src = ""
        self.strokes, src = ([], "") if self._pending else self._restore_strokes(proj)
        if self.strokes:
            self.log(L(f"{len(self.strokes)} manuelle Kratzermarkierung(en) übernommen ({src})",
                       f"{len(self.strokes)} manual scratch mark(s) restored ({src})"))
        self._strokes_changed(save=False)
        self.rtabs.setCurrentIndex(0)
        self._set_enabled_for_project(True)
        self._crop_toggled(False)
        self.statusBar().showMessage(L("Datei geladen – Phasenrollen prüfen und Auswertung starten.",
                                       "File loaded – check the phase roles and start the analysis."))
        if self._pending:
            self._apply_pending()
            return
        if self.autotest:
            if self.autotest.get("crop"):
                self._crop_from_mouse(*self.autotest["crop"])
            if self.autotest.get("strokes"):
                self.strokes = list(self.autotest["strokes"])
                self._strokes_changed(save=False)
            QTimer.singleShot(200, self.run_analysis)

    @staticmethod
    def _set_role(cb, role):
        i = cb.findData(role)
        if i >= 0:
            cb.setCurrentIndex(i)

    # ---------------------------------------------------------------- Ausschnitt
    def _crop_toggled(self, on):
        for sp in self.crop_spins:
            sp.setEnabled(on and self.project is not None)
        self.overview.set_rect(self.crop_rect() if on else None)

    def crop_rect(self):
        x0, y0, x1, y1 = (sp.value() for sp in self.crop_spins)
        return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))

    def _crop_spin_changed(self, _):
        if self.chk_crop.isChecked():
            self.overview.set_rect(self.crop_rect())

    def _crop_from_mouse(self, x0, y0, x1, y1):
        if self.project is None:
            return
        for sp, v in zip(self.crop_spins, (x0, y0, x1, y1)):
            sp.blockSignals(True)
            sp.setValue(v)
            sp.blockSignals(False)
        if not self.chk_crop.isChecked():
            self.chk_crop.blockSignals(True)
            self.chk_crop.setChecked(True)
            self.chk_crop.blockSignals(False)
            for sp in self.crop_spins:
                sp.setEnabled(True)
        self.statusBar().showMessage(L(f"Ausschnitt: x {fmt(x0)}–{fmt(x1)} µm, y {fmt(y0)}–{fmt(y1)} µm",
                                       f"Area: x {fmt(x0)}–{fmt(x1)} µm, y {fmt(y0)}–{fmt(y1)} µm"))

    # ---------------------------------------------------------------- Auswertung
    def roles(self):
        return {i: cb.currentData() for i, cb in self.role_boxes.items()}

    def run_analysis(self):
        if self.project is None:
            return
        roles = self.roles()
        n_matrix = sum(1 for r in roles.values() if r == "Matrix")
        if n_matrix != 1:
            QMessageBox.warning(self, app_title(), L("Bitte genau eine Phase als „Matrix“ festlegen.",
                                                     "Please set exactly one phase as “Matrix”."))
            return
        s = self.current_settings()
        if s.lagb >= s.hagb:
            QMessageBox.warning(self, app_title(), L("Die LAGB-Grenze muss kleiner als die HAGB-Grenze sein.",
                                                     "The LAGB limit must be smaller than the HAGB threshold."))
            return
        crop = self.crop_rect() if self.chk_crop.isChecked() else None
        if crop is not None:
            st = self.project.step
            if (crop[2] - crop[0]) / st < 20 or (crop[3] - crop[1]) / st < 20:
                QMessageBox.warning(self, app_title(), L("Der Ausschnitt ist zu klein (mindestens 20 × 20 Messpunkte).",
                                                         "The selected area is too small (at least 20 × 20 points)."))
                return
        proj = self.project
        strokes = [dict(x) for x in self.strokes]
        self.txt_log.appendPlainText("")
        self.log(L("Auswertung gestartet: ", "Analysis started: ") + proj.path.name
                 + (L(", Ausschnitt ", ", area ") + f"{tuple(round(v, 1) for v in crop)} µm" if crop else "")
                 + (L(f", {len(strokes)} manuelle Kratzermarkierung(en)", f", {len(strokes)} manual scratch mark(s)") if strokes else "")
                 + ("" if s.subgrains else L(", ohne Subkornerkennung", ", without subgrain detection")))

        def job(task):
            c = load_core()
            res = c.analyze(proj, s, roles, crop_um=crop, progress=lambda p, m: task.progress.emit(p, m),
                            cancel=task.check, log=task.message.emit, manual_scratch=strokes)
            previews = {}
            av = res.available
            for k, (fid, title, _) in enumerate(av):
                task.check()
                task.progress.emit(75 + int(25 * k / max(1, len(av))), L("Vorschau: ", "Preview: ") + title)
                previews[fid] = res.preview_png(fid)
            return res, previews

        self._start(job, self._analysis_done)

    def _analysis_done(self, out):
        res, previews = out
        self.btn_export.setEnabled(True)
        self.result = res
        self.previews = previews
        self.stored_view = None
        self._store_result(res, previews)
        self._fill_results(res.available, res.kennwerte)
        self._update_shown_label()
        if res.scratch_auto is not None and res.settings.scratch:  # automatische Erkennung in der Übersicht einblenden
            p = res.project
            full = np.zeros((p.ny, p.nx), bool)
            r0, c0 = res.crop_px
            h, w = res.scratch_auto.shape
            full[r0:r0 + h, c0:c0 + w] = res.scratch_auto
            self.auto_scratch = full
            self.auto_scratch_src = L("Auswertung", "analysis")
            self.update_overlay()
        if not self.export_folder_user or not self.ed_out.text():
            self.ed_out.setText(str(res.default_folder()))
        self.rtabs.setCurrentIndex(1)
        self.bar.setValue(100)
        pend, self._pending = self._pending, None
        if pend and pend.get("export"):
            ex = pend["export"]
            self.ed_out.setText(ex["folder"])
            for i in range(self.lst.count()):
                it = self.lst.item(i)
                it.setCheckState(Qt.Checked if it.data(Qt.UserRole) in ex["fids"] else Qt.Unchecked)
            QTimer.singleShot(300, self.export)
        if self.autotest:
            self.ed_out.setText(self.autotest["out"])
            for cb in self.fmt.values():
                cb.setChecked(True)
            self.sp_dpi.setValue(150)
            QTimer.singleShot(500, self.export)
        msg = L("Auswertung abgeschlossen", "Analysis finished")
        sm = res.summary
        if "rekristallisiert" in sm:
            msg += L(f" – rekristallisierter Anteil {fmt(sm['rekristallisiert'], 2)} %",
                     f" – recrystallized fraction {fmt(sm['rekristallisiert'], 2)} %")
        if "Korngröße" in sm:
            msg += L(f" – Korngröße {fmt(sm['Korngröße'], 2)} µm", f" – grain size {fmt(sm['Korngröße'], 2)} µm")
            if "Subkorngröße" in sm:
                msg += L(f", Subkorngröße {fmt(sm['Subkorngröße'], 2)} µm", f", subgrain size {fmt(sm['Subkorngröße'], 2)} µm")
            msg += " (ECD)"
        if "Hinweis" in sm:
            msg += "  ⚠ " + L("Subkörner kaum aufgelöst – siehe Kennwerte", "subgrains hardly resolved – see key figures")
        self.statusBar().showMessage(msg)

    def _show_preview(self, cur, _prev=None):
        if cur is None:
            return
        fid = cur.data(Qt.UserRole)
        b = self.previews.get(fid)
        if b is None and self.stored_view is not None:
            b = self.stored_view.preview(fid)
        if b is None:
            return
        pix = QPixmap()
        pix.loadFromData(b, "PNG")
        self.lbl_prev.setText(cur.text())
        self.preview.fit = True
        self.preview.set_pixmap(pix)

    def _check_all(self, on):
        for i in range(self.lst.count()):
            self.lst.item(i).setCheckState(Qt.Checked if on else Qt.Unchecked)

    # ---------------------------------------------------------------- Export
    def choose_folder(self):
        start = self.ed_out.text() or self.qs.value("last_dir", str(Path.home()))
        d = QFileDialog.getExistingDirectory(self, L("Zielordner wählen", "Choose target folder"), start)
        if d:
            self.ed_out.setText(d)
            self.export_folder_user = True

    def export(self):
        if self.stored_view is not None:
            self._export_stored()
            return
        if self.result is None:
            return
        fids = [self.lst.item(i).data(Qt.UserRole) for i in range(self.lst.count()) if self.lst.item(i).checkState() == Qt.Checked]
        formats = [k for k, cb in self.fmt.items() if cb.isChecked()]
        folder = self.ed_out.text().strip()
        if not fids:
            QMessageBox.information(self, app_title(), L("Bitte mindestens eine Abbildung ankreuzen.", "Please tick at least one figure."))
            return
        if not formats:
            QMessageBox.information(self, app_title(), L("Bitte mindestens ein Format wählen.", "Please choose at least one format."))
            return
        if not folder:
            QMessageBox.information(self, app_title(), L("Bitte einen Zielordner wählen.", "Please choose a target folder."))
            return
        res, dpi, raster = self.result, self.sp_dpi.value(), self.chk_raster.isChecked()
        self.log(L(f"Export: {len(fids)} Abbildung(en) × {', '.join(f.upper() for f in formats)}, {dpi} dpi → {folder}",
                   f"Export: {len(fids)} figure(s) × {', '.join(f.upper() for f in formats)}, {dpi} dpi → {folder}"))

        def job(task):
            return res.export(fids, folder, formats, dpi=dpi, rasterize_vector=raster,
                              progress=lambda p, m: task.progress.emit(p, m), cancel=task.check)

        self._start(job, self._export_done)

    def _export_stored(self):
        ds = self.stored_view
        fids = [self.lst.item(i).data(Qt.UserRole) for i in range(self.lst.count()) if self.lst.item(i).checkState() == Qt.Checked]
        folder = self.ed_out.text().strip()
        if not folder:
            QMessageBox.information(self, app_title(), L("Bitte einen Zielordner wählen.", "Please choose a target folder."))
            return
        box = QMessageBox(self)
        box.setWindowTitle(app_title())
        box.setText(L(f"„{ds.label}“ ist ein gespeicherter Datensatz.\n\nFür Abbildungen in voller Auflösung wird er aus der "
                      "ctf-Datei mit den gespeicherten Rollen, Parametern, Kratzermarkierungen und dem Ausschnitt neu berechnet "
                      "(Dauer wie eine Auswertung). Alternativ werden die Vorschaubilder, Kennwerte und Kornlisten kopiert.",
                      f"“{ds.label}” is a stored dataset.\n\nFor full-resolution figures it is recalculated from the ctf file "
                      "with the stored roles, parameters, scratch marks and area (takes as long as an analysis). "
                      "Alternatively, the previews, key figures and grain lists are copied."))
        b_calc = box.addButton(L("Neu berechnen und exportieren", "Recalculate and export"), QMessageBox.AcceptRole)
        b_copy = box.addButton(L("Vorschauen und Kennwerte kopieren", "Copy previews and key figures"), QMessageBox.ActionRole)
        box.addButton(QMessageBox.Cancel)
        box.exec()
        if box.clickedButton() == b_calc:
            self.recompute(ds, export={"folder": folder, "fids": fids})
        elif box.clickedButton() == b_copy:
            try:
                written = dsets.export_stored(ds, folder)
            except Exception as e:  # noqa: BLE001
                QMessageBox.critical(self, app_title(), L("Export fehlgeschlagen:\n", "Export failed:\n") + str(e))
                return
            self.log(L(f"Gespeicherter Datensatz „{ds.label}“: {len(written)} Dateien kopiert → {folder}",
                       f"Stored dataset “{ds.label}”: {len(written)} files copied → {folder}"))
            self._export_done(written)

    def _export_done(self, written):
        self.bar.setValue(100)
        folder = Path(written[0]).parent if written else Path(self.ed_out.text())
        self.statusBar().showMessage(L(f"{len(written)} Dateien geschrieben nach {folder}", f"{len(written)} files written to {folder}"))
        if self.autotest:
            self.grab().save(str(Path(self.autotest["out"]) / "autotest_fenster.png"))
            Path(self.autotest["out"], "autotest.log").write_text(
                self.txt_log.toPlainText() + f"\nAUTOTEST OK: {len(written)} Dateien\n", encoding="utf-8")
            QApplication.instance().quit()
            return
        self._written_box(len(written), folder)

    def open_folder(self):
        folder = self.ed_out.text().strip()
        if folder and Path(folder).exists():
            QDesktopServices.openUrl(QUrl.fromLocalFile(folder))
        else:
            QMessageBox.information(self, app_title(), L("Der Ordner existiert noch nicht – zuerst exportieren.",
                                                         "The folder does not exist yet – export first."))

    def save_kennwerte(self):
        if not self.txt_kenn.toPlainText():
            return
        start = str(Path(self.ed_out.text() or str(Path.home())) / kennwerte_name())
        path, _ = QFileDialog.getSaveFileName(self, L("Kennwerte speichern", "Save key figures"), start, L("Text (*.txt)", "Text (*.txt)"))
        if path:
            Path(path).write_text(self.txt_kenn.toPlainText(), encoding="utf-8")
            self.statusBar().showMessage(L("Kennwerte gespeichert: ", "Key figures saved: ") + path)

    # ---------------------------------------------------------------- Einstellungsdateien
    def save_settings_file(self):
        path, _ = QFileDialog.getSaveFileName(self, L("Einstellungen speichern", "Save settings"),
                                              str(Path(self.qs.value("last_dir", str(Path.home()))) / L("EBSD-Einstellungen.json", "EBSD-settings.json")),
                                              "JSON (*.json)")
        if path:
            d = {"settings": self.current_settings().to_dict(),
                 "export": {"formats": [k for k, cb in self.fmt.items() if cb.isChecked()], "dpi": self.sp_dpi.value()}}
            if self.project is not None:
                d["roles_by_name"] = {self.project.names[i]: r for i, r in self.roles().items()}
            Path(path).write_text(json.dumps(d, indent=2, ensure_ascii=False), encoding="utf-8")
            self.statusBar().showMessage(L("Einstellungen gespeichert: ", "Settings saved: ") + path)

    def load_settings_file(self):
        path, _ = QFileDialog.getOpenFileName(self, L("Einstellungen laden", "Load settings"), self.qs.value("last_dir", str(Path.home())),
                                              "JSON (*.json)")
        if not path:
            return
        try:
            d = json.loads(Path(path).read_text(encoding="utf-8"))
            self.apply_settings(Settings.from_dict(d.get("settings", {})))
            ex = d.get("export", {})
            if ex:
                for k, cb in self.fmt.items():
                    cb.setChecked(k in ex.get("formats", []))
                self.sp_dpi.setValue(int(ex.get("dpi", self.sp_dpi.value())))
            rn = d.get("roles_by_name", {})
            if rn and self.project is not None:
                for i, cb in self.role_boxes.items():
                    if self.project.names[i] in rn:
                        self._set_role(cb, rn[self.project.names[i]])
            self.statusBar().showMessage(L("Einstellungen geladen: ", "Settings loaded: ") + path)
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, app_title(), L("Die Datei konnte nicht gelesen werden:\n", "The file could not be read:\n") + str(e))

    # ---------------------------------------------------------------- Aufgaben
    def _start(self, fn, on_done):
        self.task = Task(fn, self)
        self.task.progress.connect(self._on_progress)
        self.task.message.connect(self.log)
        self.task.done.connect(lambda out: (self._busy(False), on_done(out)))
        self.task.failed.connect(self._on_failed)
        self._busy(True)
        self.task.start()

    def _on_progress(self, p, msg):
        if self.bar.maximum() == 0:
            self.bar.setRange(0, 100)
        self.bar.setValue(p)
        self.statusBar().showMessage(msg)

    def _on_failed(self, text, cancelled):
        self._pending = None
        self._busy(False)
        self.bar.setRange(0, 100)
        self.bar.setValue(0)
        self.log(text)
        if self.autotest:
            Path(self.autotest["out"]).mkdir(parents=True, exist_ok=True)
            Path(self.autotest["out"], "autotest.log").write_text(self.txt_log.toPlainText() + "\nAUTOTEST FEHLER\n", encoding="utf-8")
            QApplication.instance().quit()
            return
        if cancelled:
            self.statusBar().showMessage(L("Abgebrochen.", "Cancelled."))
        else:
            self.statusBar().showMessage(L("Fehler – Details im Protokoll.", "Error – details in the log."))
            QMessageBox.critical(self, app_title(), text.split("\n\n")[0] + L("\n\nDetails stehen im Reiter „Protokoll“.",
                                                                              "\n\nDetails are in the tab “Log”."))

    def cancel_task(self):
        if self.task is not None and self.task.isRunning():
            self.task.cancel()
            self.statusBar().showMessage(L("Abbruch angefordert – wird nach dem aktuellen Schritt beendet …",
                                           "Cancellation requested – stops after the current step …"))

    # ---------------------------------------------------------------- Hilfe
    def show_help(self):
        dlg = QDialog(self)
        dlg.setWindowTitle(f"{app_title()} – " + L("Kurzanleitung", "Quick guide"))
        v = QVBoxLayout(dlg)
        tb = QTextBrowser()
        tb.setHtml(help_html())
        v.addWidget(tb)
        bb = QDialogButtonBox(QDialogButtonBox.Close)
        bb.rejected.connect(dlg.reject)
        v.addWidget(bb)
        dlg.resize(780, 680)
        dlg.exec()

    def show_about(self):
        try:
            import importlib.metadata as md
            mv = md.version("mtex")
        except Exception:  # noqa: BLE001
            mv = "?"
        QMessageBox.about(self, L("Über ", "About ") + app_title(),
                          f"<b>{app_title()} {APP_VERSION}</b><br>"
                          + L("EBSD-Auswertung von Oxford-ctf-Dateien", "EBSD analysis of Oxford ctf files") + "<br><br>"
                          + L(f"Berechnungen mit MTEX für Python {mv}", f"Calculations with MTEX for Python {mv}")
                          + " (GPL-2.0, mtex-toolbox.github.io)<br>"
                          + L("Oberfläche", "User interface") + ": Qt / PySide6 (LGPL)<br><br>"
                          + L("Institut für Metallformung, TU Bergakademie Freiberg", "Institute of Metal Forming, TU Bergakademie Freiberg"))


# ------------------------------------------------------------------ Neustart (Sprachwechsel)
def restart_program():
    if getattr(sys, "frozen", False):
        prog, args = sys.executable, sys.argv[1:]
    else:
        prog, args = sys.executable, [str(Path(__file__).resolve())] + sys.argv[1:]
    args = [a for a in args if not a.startswith("--lang")]
    QProcess.startDetached(prog, args)


# ------------------------------------------------------------------ Kommandozeilen-Modi
def selftest(argv):
    """--selftest DATEI AUSGABE [--crop x0 y0 x1 y1] [--nosub]: Auswertung ohne Oberfläche, Export aller Formate"""
    i = argv.index("--selftest")
    ctf, out = argv[i + 1], Path(argv[i + 2])
    crop = None
    if "--crop" in argv:
        j = argv.index("--crop")
        crop = tuple(float(v) for v in argv[j + 1:j + 5])
    out.mkdir(parents=True, exist_ok=True)
    logf = open(out / "selbsttest.log", "w", encoding="utf-8")

    def log(m):
        logf.write(m + "\n")
        logf.flush()

    t0 = time.time()
    try:
        c = load_core()
        proj = c.EBSDProject.load(ctf, log=log)
        roles = proj.default_roles()
        log("Rollen: " + ", ".join(f"{proj.names[k]}={v}" for k, v in roles.items()))
        s = Settings()
        s.subgrains = "--nosub" not in argv
        res = c.analyze(proj, s, roles, crop_um=crop, log=log)
        fids = [f for f, _, _ in res.available]
        res.export(fids, out, ["png"], dpi=150)
        res.export(fids[:2], out, ["pdf", "svg", "tiff"], dpi=150)
        log(f"SELBSTTEST OK ({time.time() - t0:.0f} s)")
        return 0
    except Exception:  # noqa: BLE001
        log(traceback.format_exc())
        log("SELBSTTEST FEHLER")
        return 1
    finally:
        logf.close()


def screenshot(argv, app):
    """--screenshot BILD.png [--ctf DATEI] [--run] [--crop x0 y0 x1 y1] [--tab N] [--ltab N] [--marks]: Fenster als Bild"""
    out = argv[argv.index("--screenshot") + 1]
    w = MainWindow()
    w.resize(1500, 950)
    if "--ctf" in argv:
        c = load_core()
        proj = c.EBSDProject.load(argv[argv.index("--ctf") + 1], log=w.log)
        w._file_loaded(proj)
        if "--crop" in argv:
            j = argv.index("--crop")
            vals = [float(v) for v in argv[j + 1:j + 5]]
            w._crop_from_mouse(*vals)
            w.overview.set_rect(w.crop_rect())
        if "--marks" in argv:
            w.strokes = [{"modus": "remove", "x0": 50, "y0": 32.5, "x1": 80, "y1": 46, "breite": 6},
                         {"modus": "add", "x0": 100, "y0": 10, "x1": 100, "y1": 90, "breite": 1.5}]
            w._strokes_changed(save=False)
            w.set_tool("add")
        if "--run" in argv:
            s = w.current_settings()
            crop = w.crop_rect() if w.chk_crop.isChecked() else None
            res = c.analyze(proj, s, w.roles(), crop_um=crop, log=w.log, manual_scratch=w.strokes)
            w._analysis_done((res, {f: res.preview_png(f) for f, _, _ in res.available}))
            if "--fig" in argv:
                fids = [f for f, _, _ in res.available]
                fid = argv[argv.index("--fig") + 1]
                if fid in fids:
                    w.lst.setCurrentRow(fids.index(fid))
    if "--tab" in argv:
        w.rtabs.setCurrentIndex(int(argv[argv.index("--tab") + 1]))
    if "--ltab" in argv:
        w.ltabs.setCurrentIndex(int(argv[argv.index("--ltab") + 1]))
    w.show()

    def grab():
        w.grab().save(out)
        app.quit()

    QTimer.singleShot(1500, grab)
    return app.exec()


def hard_exit(code):
    """Prozess sicher beenden: Hintergrund-Threads von numba/OpenMP können das normale Beenden blockieren"""
    try:
        sys.stdout.flush()
        sys.stderr.flush()
    except Exception:  # noqa: BLE001
        pass
    os._exit(code)


def init_language(argv):
    """Sprache: --lang de|en, sonst gespeicherte Einstellung (Menü „Sprache / Language“), sonst Deutsch"""
    code = None
    for k, a in enumerate(argv):
        if a.startswith("--lang="):
            code = a.split("=", 1)[1]
        elif a == "--lang" and k + 1 < len(argv):
            code = argv[k + 1]
    if code is None:
        code = QSettings(ORG, APP_NAME).value("language", "de")
    return set_language(code)


def main():
    argv = sys.argv
    lang = init_language(argv)
    if "--selftest" in argv:
        hard_exit(selftest(argv))
    QApplication.setHighDpiScaleFactorRoundingPolicy(Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    if lang == "de":
        QLocale.setDefault(QLocale(QLocale.German, QLocale.Germany))      # Dezimalkomma in allen Eingabefeldern
    else:
        QLocale.setDefault(QLocale(QLocale.English, QLocale.UnitedKingdom))
    app = QApplication(argv)
    app.setApplicationName(APP_NAME)
    if lang == "de":
        tr = QTranslator(app)                                             # deutsche Qt-Standardtexte (Ja/Nein, Abbrechen …)
        import PySide6
        for d in (QLibraryInfo.path(QLibraryInfo.TranslationsPath), str(Path(PySide6.__file__).parent / "translations")):
            if tr.load(QLocale(QLocale.German, QLocale.Germany), "qtbase", "_", d):
                app.installTranslator(tr)
                break
    icon = resource_dir() / "icon.png"
    if icon.exists():
        app.setWindowIcon(QIcon(str(icon)))
    app.setStyle("Fusion")
    if "--screenshot" in argv:
        hard_exit(screenshot(argv, app))
    w = MainWindow()
    w.show()
    if "--autotest" in argv:                     # --autotest DATEI AUSGABE [--crop x0 y0 x1 y1]
        i = argv.index("--autotest")
        crop = None
        if "--crop" in argv:
            j = argv.index("--crop")
            crop = tuple(float(v) for v in argv[j + 1:j + 5])
        w.autotest = {"out": argv[i + 2], "crop": crop}
        QTimer.singleShot(300, lambda: w.open_file(argv[i + 1]))
        hard_exit(app.exec())
    files = [a for a in argv[1:] if a.lower().endswith(".ctf") and Path(a).exists()]
    if files:
        QTimer.singleShot(200, lambda: w.open_file(files[0]))
    hard_exit(app.exec())


if __name__ == "__main__":
    main()
