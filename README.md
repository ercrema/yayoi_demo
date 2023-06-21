# Data and R scripts for the paper 'Regional demographic responses to the arrival of rice farming in prehistoric Japan'

This repository contains data and scripts used in the following paper:

Crema, E. R., Carrignon, S., Shoda, S., & Stevens, C. J., & Shoda, S. (202X). Regional demographic responses to the arrival of rice farming in prehistoric Japan. 

The repository is organised into four main directories: _data_, _analyses_, _results_, _figures_, and _src_.
The _data_ folder contains all relevant radiocarbon and settlement data, _analyses_ contains core R scripts for executing all analyses and estimates, _results_ contains R image files of all outputs, _figures_ contains all figures for the manuscript and the supplementary materials, and _src_ contains additional custom utility R functions. 

## Analyses Summary

### Settlement Data Analyses

Settlement data were obtained from the [「縄文・弥生集落データベース」("Jomon-Yayoi settlement database") of the National Museum of Japanese History](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/jomo/db_param), and stored as a CSV file (`data/site_raw.csv`). The file `data/data_prep_sites.R` contains data wrangling R scripts that converts relevant data into an R data.frame stored in the R image file `data/sitedata.RData`. Scripts for the composite Kernel Density estimates (see manuscript and ESM for details) are stored in the R file `figures/figure2.R`. 

### Bayesian Analyses of Radiocarbon Dates

14C data were obtained from [radiocarbon database of the National Museum of Japanese History](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). We primarily used the redacted version of the database described in [Kudo et al 2023](https://openarchaeologydata.metajnl.com/) but extracted specific fields from the original version to identify anthropogenic contexts. The file `data/data_prep_c14.R` contains pre-processing R scripts for generating the core dataset used in this paper, stored in the R image file `data/c14data.RData`. Bayesian analyses were carried out using the [nimbleCarbon](https://CRAN.R-project.org/package=nimbleCarbon) R package, which contains custom probability distributions and utility functions for the [NIMBLE](https://cran.r-project.org/web/packages/nimble/index.html) probabilistic programming language. The files `analyses/icar500.R` and `analyses/icar750.R` contains the R script for fitting the Bayesian models for the 500 and 750 yrs time-intervals. Posterior samples are stored in the R image files `results/icar_c14doubleRes500.RData` and `results/icar_c14doubleRes750.RData`. 


### Absolute Population Estimates

Absolute population estimates were calculated using a modified version of equation introduced by Koyama (Koyama, S.1978. Jomon Subsistence and Population. Senri Ethnological Studies, 2, 1–65.) with an updated dataset. The core calculations are included as an R script in the file `analyses/pop_dens_est.R`, with the outputs stored in the CSV files `results/pop_estimate_compare.csv` (comparison of the different estimates discussed in the supplementary materials) and `results/pop_estimate_region.csv` (estimates used in table 1 of the manuscript). Raw input data required for the calculations are Koyama's original data (`data/koyama_popestimate_1984.csv`), the number of archaeological sites [published by the Japanese Agency of Cultural Affairs](https://www.bunka.go.jp/seisaku/bunkazai/shokai/pdf/h29_03_maizotokei.pdf) (`data/maizobunkazai_2017.csv`), and a lookup table for matching administrative units (prefectures) to the regions used in this paper (`data/prefecture_data.csv`).

## File Structure

### data
 * `c14data.RData`
 * `c14db_1.0.0.Rds`
 * `c14raw_1.0.0.Rds`
 * `data_prep_c14.R`
 * `data_prep_sites.R`
 * `data_summary.R`
 * `koyama_popestimate_1984.csv`
 * `maizobunkazai_2017.csv`
 * `prefecture_data.csv`
 * `sitedata.RData`
 * `site_raw.csv`
   
### analyses
 * `icar500.R`
 * `icar750.R`
 * `pop_dens_est.R`

### results
 * `icar_c14doubleRes500.RData`
 * `icar_c14doubleRes750.RData`
 * `pop_estimate_compare.csv`
 * `pop_estimate_region.csv`
   
### figures
 * `figures_main.R`
 * `figures_esm.R`
 * `figure1.pdf` ~ `figure6.pdf`
 * `figureS1.pdf` ~ `figureS3.pdf`
 
### src
* `dbscanID.R`

## R Session Info

## Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema) and by a Philip Leverhulme Prize (PLP-2019-304) in archaeology awarded to Enrico Crema.

## Licence
CC-BY 3.0

