EBSD Analysis 1.3 – Windows program
====================================
(Deutsch: siehe LIESMICH.txt)

EBSD analysis of Oxford/AZtec ctf files with MTEX for Python:
IPF maps, pole figures and inverse pole figures, KAM, grain boundary misorientation, phase distribution,
recrystallization analysis (GOS), grain size, subgrain size and grain aspect ratio, scratch removal
(automatic and manual); stored datasets and comparison of several datasets.


LANGUAGE
--------
Menu "Sprache / Language" -> English. The language is applied after restarting the program (the program offers
to restart immediately). User interface, help, figures, key figures, grain lists and the comparison then appear
in English. English outputs use English file names (e.g. Results_and_parameters.txt, Grain_list_grains.csv)
and CSV files with comma separator and decimal point. Command line: --lang en.


START
-----
1. Unpack the ZIP file completely (e.g. to C:\Programs\EBSD-Auswertung or your Documents folder).
   Do NOT start the program directly from the ZIP file.
2. Double-click EBSD-Auswertung.exe. The folder "_internal" must stay next to the .exe.
Windows SmartScreen may warn on the first start because the program is not digitally signed:
"More info" -> "Run anyway". No installation or administrator rights are required.
Requirements: Windows 10 or 11, 64 bit; about 4 GB free memory for a typical map (approx. 2 million points).


OPERATION (short version – details in the program under Help -> Quick guide)
---------
1. File -> Open ctf file.
2. Check the phase roles: exactly one matrix; particle phases are kept; misindexing is discarded and filled in
   from the neighbours.
3. Switch functions on/off (scratch removal, Mg2Si detection, individual analyses, subgrain detection).
4. Optionally select an area: drag a rectangle in the overview.
   Optionally mark scratches manually or remove markings (overview, see below).
5. "Start analysis". Previews then appear under "Results".
6. Tick figures, choose formats (PNG, PDF, SVG, TIFF) and resolution, export.
   Results_and_parameters.txt is always written, with the grain size analysis also Grain_list_grains.csv and
   Grain_list_subgrains.csv.


NEW IN 1.3
----------
- German / English user interface and outputs (see above).
- Mark scratches manually / remove markings (tab "Overview / area / scratches"):
    Tool "Mark scratch": draw a line along the scratch (width in µm adjustable).
    Tool "Unmark": draw a line over a wrong marking – nothing counts as a scratch there any more.
    "Show automatic detection" shows the automatic detection in red (also after every analysis).
    Colours: red = automatic, magenta = marked manually, light blue = unmarked.
    Undo (Ctrl+Z), Clear all; magnifier/pan to zoom into the map.
  The marks are saved per ctf file (%LOCALAPPDATA%\EBSD-Auswertung\kratzer_markierungen.json), restored when the
  file is opened again, stored in the dataset and used for "Recalculate". They also apply when automatic detection
  is switched off. Figure 00 and the key figures show automatic and manual contributions separately.
- Subgrain detection can be switched off (function "Detect subgrains"), e.g. if the step size is too large for
  subgrains (rule of thumb: subgrain >= approx. 5 step sizes; the program issues a note if the median ECD of the
  subgrains is below that). Without subgrains, the LAGB lines, figure 12, the subgrain list and the cross-check at
  subgrain level are omitted; recrystallization is classified at grain level only (recrystallized /
  non-recrystallized). The misorientation distribution is unchanged. The dataset comparison shows whether a
  dataset was analysed with or without subgrains.
- With subgrain detection on and without manual marks, the results are identical to version 1.2.
  Datasets stored with 1.2 remain usable.


DATASETS AND COMPARISON
-----------------------
Every analysis is stored automatically as a dataset (tab "Datasets / Comparison") and is available again after a
restart. Double-click or "Show results" shows the stored results without recalculation; "Recalculate" loads the
ctf file and calculates again with the stored settings, roles, area and scratch marks.
Ticked datasets appear in the comparison: V1 key figures (subgrain size, HAGB fraction, recrystallized fraction in
total and 0–1 µm from particles, <110> || compression direction) and V2 misorientation distribution, plus a table.
Storage: %LOCALAPPDATA%\EBSD-Auswertung\datensaetze.


NOTES
-----
- The first run after starting the program is slower because routines are compiled.
- Program output/errors: %LOCALAPPDATA%\EBSD-Auswertung\ausgabe.log and tab "Log".
- Self-test without user interface (e.g. for IT):
    EBSD-Auswertung.exe --selftest "C:\path\file.ctf" "C:\path\output folder" [--lang en] [--nosub]


LICENCES
--------
Calculations: MTEX for Python (GPL-2.0, https://github.com/mtex-toolbox/pymtex).
User interface: Qt / PySide6 (LGPL-3.0). Further libraries: NumPy, SciPy, matplotlib, numba, llvmlite, finufft,
Pillow (BSD/MIT-type licences). As MTEX is licensed under the GPL, the GPL-2.0 also applies to this program;
the source code is in the package "EBSD-Auswertung_Quellcode".

Institute of Metal Forming, TU Bergakademie Freiberg
