# Changelog

## 1.4 – November 2026

- **Test data per dataset** (*Test data …*): alloy, Fe/Si content, position, test temperature and strain rate,
  **local temperature, local strain φ and local strain rate φ̇** (e.g. from the FE simulation at the measured
  position) and flow stress k_f. Decimal comma or point; *Suggest from label* reads alloy, temperature, strain rate
  and position from labels such as "750-3, 350 °C / 10 s⁻¹, centre".
- **Input field Q** (kJ/mol, default 142) for the Zener–Hollomon parameter Z = φ̇·exp(Q/RT); all datasets are
  evaluated with the same Q. Option *ln Z from*: local values (else test condition) or always the test condition.
  ln Z appears in the comparison table and CSV.
- **New comparison figures V3–V13:** key figures vs ln Z (overview or one selectable key figure, exponential fit per
  alloy), recrystallized fraction vs particle distance, decomposition of the difference in recrystallized fraction
  (particle proximity vs tendency to recrystallize) with what-if bars, maps side by side (optionally same µm scale),
  GOS classes, subgrain size distribution, subgrain size vs flow stress (Holt relation d·k_f/(G·b)), geometric DRX,
  KAM distribution and KAM per class, comparability check with warnings, sensitivity to the GOS threshold.
  V2 now shows the Mackenzie distribution and the cumulative curve.
- **Resampling to a common step size:** parameter *Resample step size to* (0 = off, nearest neighbour) and button
  *Match step size …*, which re-analyses all ticked datasets with a finer step on the coarsest grid; originals are kept.
- Additional key figures: KAM median, 90 % quantile and distribution, recrystallized fraction per particle distance
  with area fractions, GOS sensitivity, share of the three largest grains, grain counts, particle fraction, d10/d90.
- Datasets stored with 1.2/1.3 remain usable; additional key figures are read from the stored key figures text where
  possible. KAM distribution and GOS sensitivity need *Recalculate*.
- New module `ebsd_compare.py` (comparison figures V3–V13, test data, Z, Mackenzie distribution).
- The PDF documentation still describes version 1.3; the Markdown user guide covers 1.4.

## 1.3 – November 2026

- **English language:** user interface, help, figures, key figures, grain lists and the dataset comparison are
  available in English and German (menu *Sprache / Language*, or `--lang en|de`). English outputs use English file
  names and CSV files with comma separator and decimal point.
- **Manual scratch editing:** new tools *Mark scratch* and *Unmark* in the overview, adjustable line width, undo,
  clear all, zoom/pan toolbar, preview of the automatic detection. Marks are saved per ctf file, stored with each
  dataset and reused when recalculating; they also apply with automatic detection switched off. Figure 00 and the
  key figures report automatic and manual contributions separately.
- **Subgrain detection switch:** new function *Detect subgrains*. Without it, no subgrain reconstruction, no LAGB
  lines, no figure 12 and no subgrain list; recrystallization is classified at grain level only
  (recrystallized / non-recrystallized). A note is issued if the median subgrain ECD is below 5 step sizes.
- Comparison table: new rows *Subgrain detection*, *Manual scratch marks*, *non-recrystallized (%)*.
- With subgrain detection on and no manual marks, results are identical to version 1.2. Datasets stored with 1.2
  remain usable.
- New user and training documentation (English, Markdown and PDF).

## 1.2

- Stored datasets: every analysis is saved with previews, key figures and grain lists and can be shown, exported or
  recalculated later.
- Comparison of several datasets: figures V1 (key figures) and V2 (misorientation distribution) and a comparison table.

## 1.1

- Grain size (HAGB grains), subgrain size (boundaries ≥ LAGB) and aspect ratio; figures 11–14; grain lists as CSV.
  Edge grains and (sub)grains below the minimum number of points stay in the grain lists but are excluded from the
  statistics.

## 1.0

- First version: cleaning, scratch detection, IPF maps, texture, KAM, grain boundary misorientation, phase
  distribution, recrystallization (GOS).
