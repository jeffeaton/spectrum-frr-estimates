# spectrum-frr-estimates
Code for reproducing analysis of fertility rate ratios by age, CD4, and ART status for Spectrum 2019 described in Eaton _et al._ "Relative fertility of HIV positive women in the ART era: updated estimates from national household survey data" (in preparation).


## Data

The directory `/data` includes code for preparing the dataset used for analysis. The following scripts process survey datasets to create estimates of ASFR, TFR, current pregnancy and proportion sexaully active in the past twelve months by age group and HIV status.
* `data/prep-dhs.R`
* `data/prep-hsrc.R`
* `data/prep-kais.R`
* `data/prep-mics.R`
* `data/prep-phia.R`
* `data/prep-ug04ais.R`

These rely on the R package [_demogsurv_](http://github.com/mrc-ide/demogsurv/) for analysis of fertility rates from survey data.

The script `data/spec-pop1.R` extracts population stratification from Spectrum 2018 'pop1.xlsx' output files and calculates the distrribution of HIV positive women by age, CD4 stage, and ART status for each country.

The script `data/data.R` compiles survey estimates and Spectrum outputs into a dataset for analysis.

## Model

The file `model/model.stan` is the [Stan](https://mc-stan.org) model code.

The file `model/functions.R` contains a function for preparing model frame inputs to the model.

The file `model/fit.R` is the script used for fitting the model and saving outputs.

## Analysis

The script `analysis/analysis.R` contains code for reproducing all analysis, tables, and figures in the manuscript.
