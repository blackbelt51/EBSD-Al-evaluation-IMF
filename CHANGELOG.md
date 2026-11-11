# Changelog

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
