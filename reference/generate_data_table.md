# Generate a multi-page landscape LaTeX data summary table

Produces a landscape longtable listing all collected incubation period
observations across pathogens, with one bold section header row per
pathogen, a khaki-coloured column header, and a DOI column.

## Usage

``` r
generate_data_table(
  datasets,
  caption = "Collected incubation period datasets by pathogen.",
  label = "tab:data_summary",
  digits = 1L
)
```

## Arguments

- datasets:

  A named list of pathogen datasets. Each element must be a named list
  of entry lists (e.g.
  `list(Nipah = datasets_Nipah, MVD = datasets_MVD)`). Pathogens appear
  in the order of `names(datasets)`.

- caption:

  LaTeX `\caption{}` string.

- label:

  LaTeX `\label{}` string.

- digits:

  Integer. Decimal places for numeric summaries (default 1).

## Value

A single character string containing the complete landscape `longtable`
environment, ready to paste into a manuscript.

## Details

Each dataset entry produces one table row. The Summary column reports:

- **Type A** (median + range): `"Median: X (Range: lo--hi)"`

- **Type B** (median + IQR): `"Median: X (IQR: Q1--Q3)"`

- **Type C** (mean + SD): `"Mean: X (SD: s)"`

- **Type D** (frequency table): prefixed *Freq. table.*, then weighted
  median and observed range

- **Type E** (interval-censored frequency table): prefixed *Freq. table
  (interval-censored).*, then overall range

For Type D entries where `n` is absent, \\n\\ is computed as
`sum(freq_count)`.

Required LaTeX packages (add to preamble):


      \usepackage{longtable}
      \usepackage{array}
      \usepackage[table]{xcolor}
      \usepackage{pdflscape}
      \newcolumntype{C}[1]{>{\centering\arraybackslash}p{#1}}
      \newcolumntype{R}[1]{>{\raggedleft\arraybackslash}p{#1}}

## See also

[`generate_results_table()`](https://cm401.github.io/ddsynth/reference/generate_results_table.md),
[`generate_results_table_split()`](https://cm401.github.io/ddsynth/reference/generate_results_table_split.md)
