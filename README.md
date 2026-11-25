# GARF: Group-Aware Random Forest matching for panel counterfactuals


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24303486.svg)](https://doi.org/10.5281/zenodo.24303486)

GARF builds the comparison group for a policy evaluation with panel data. It
encodes each unit's pre-treatment history as levels and trajectory descriptors,
trains a random forest on treatment, and matches every treated unit to untreated
donors on the forest's leaf proximity. It returns a matched panel and a set of
pre-treatment diagnostics.vThe effect is estimated in a separate step: the
companion command `garfte` (new in 1.2) gives a matching difference-in-differences
or the treated-minus-own-donors gap, by event time and overall, with an
event-study graph and the effect of each treated unit; or take the matched panel
to your own estimator.

This repository contains the Python engines (matching and estimation), the
Stata commands `garf` and `garfte`, the user guide and the do-files for the
design-stage diagnostics.

**Status:** version 1.2 (adds the estimation stage `garfte`; version 1.1 was the first public release). Binary, absorbing treatment with a single onset date.

## Requirements

- Python 3.10 or newer (the setup script installs pandas, numpy, scipy,
  scikit-learn, PyYAML and tqdm into a local `.venv`; nothing system-wide)
- For the Stata command: Stata 14 or newer. Stata's own Python integration is
  not needed.

## Install

Download the [latest release](https://github.com/GITHUB_USER/garf/releases/latest)
(or clone the repository), then run the setup once from the package folder:

```
./setup.sh        # macOS / Linux
setup.bat         # Windows (double-click)
```

In Stata, point to the package and add the command to the ado-path:

```stata
global GARF_HOME "/path/to/GARF-binary"
adopath + "$GARF_HOME/integrations/stata"
which garf
which garfte
help garf
help garfte
```

Put the two lines in your `profile.do` to make them permanent. See
`INSTALL.txt` if you need to install Python first.

## First run

```stata
import delimited "$GARF_HOME/demo/input/panel.csv", clear varnames(1) case(preserve)
save raw.dta, replace

garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) ///
      tvarying(x1 x2) tinvariant(region) store
return list

garfte, method(gap)                   // treated minus its own donors, in levels
```

The matching difference-in-differences needs a match without replacement:

```stata
use raw.dta, clear
garf, treatment(W) outcome(Y) id(id) time(year) t0(2000) ///
      tvarying(x1 x2) tinvariant(region) noreplacement k(3)
garfte, method(did)
```

`garf` replaces the data in memory with the matched panel, so save or
`preserve` your data first. From a terminal, without Stata:

```
cd demo && ./run_demo.sh
```

## What is in the repository

| Folder | Content |
|---|---|
| `engine/` | `garf_binary.py`, the matching engine, and `garf_estimate.py`, the estimation stage, with a config template and run scripts |
| `integrations/stata/` | `garf.ado`, `garfte.ado` and their help files, an example do-file, the guide source |
| `DIAGNOSTICS/` | one do-file per design-stage diagnostic; `00_master.do` runs them all |
| `demo/` | a small synthetic panel and its generator, to check the installation |
| `GARF_Stata_Guide_v1.2.pdf` | the user guide |
| `GARF_Reading_Guide_v1.2.pdf` | reading the output: every number, table and figure of `garf`, the diagnostics and `garfte`, and how to interpret it |

`README.txt` documents the engine and its configuration file in full;
`CHANGELOG.md` lists the changes by version.

## Citing

If you use GARF, please cite the paper and the software:

> Fantechi, F. (2026). *Matching Histories, not Just Averages: GARF, Group-Aware Random Forest for Counterfactuals Identification in Panel Data*. Working paper, University of Palermo.
> https://dx.doi.org/10.2139/ssrn.7719438.

> Fantechi, F. (2026). *GARF: Group-Aware Random Forest matching for panel
> counterfactuals* (version 1.2). Zenodo. https://doi.org/10.5281/zenodo.24303486

GitHub's "Cite this repository" button gives the same in BibTeX.

## Licence

MIT; see `LICENSE`.

## Contact

Bug reports and suggestions: open an issue, or write to
Federico Fantechi, University of Palermo, federico.fantechi@unipa.it.
