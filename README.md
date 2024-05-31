[![DOI](https://zenodo.org/badge/370637019.svg)](https://zenodo.org/doi/10.5281/zenodo.11398533)

# Data and R scripts for the paper 'Regional variations in the demographic response to the arrival of rice farming in prehistoric Japan'

This repository contains data and scripts used in the following paper:

Crema, E. R., Carrignon, S., Shoda, S., & Stevens, C. J., & Shoda, S. (In press). Regional variations in the demographic response to the arrival of rice farming in prehistoric Japan. _Antiquity_.

The repository is organised into five main directories: _data_, _analyses_, _results_, _figures_, _tables_, and _src_.
The _data_ folder contains all relevant radiocarbon and settlement data, _analyses_ contains core R scripts for executing all analyses and estimates, _results_ contains R image files of all outputs, _figures_ and _tables_ contain all figures and tables for the manuscript and the supplementary materials as well as R scripts required to generate them, and _src_ contains additional custom utility R functions. 

## Analyses Summary

### Settlement Data Analyses

Settlement data were obtained from the [「縄文・弥生集落データベース」("Jomon-Yayoi settlement database") of the National Museum of Japanese History](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/jomo/db_param), and stored as a CSV file (`data/site_raw.csv`). The file `data/data_prep_sites.R` contains data wrangling R scripts that convert relevant data into an R data.frame stored in the R image file `data/sitedata.RData`. Scripts for the composite Kernel Density estimates (see manuscript and ESM for details) are stored in the R file `figures/figure2.R`. Pipeline: `data/site_raw.csv` &rarr; `data/sitedata.RData`.

### Bayesian Analyses of Radiocarbon Dates

14C data were obtained from [radiocarbon database of the National Museum of Japanese History](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param). We joined the cleaned and translated version of the database (`c14db_1.1.0.csv`, obtained [here](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd_en/db_param)) to the full Japanese version of the database (`.xslx` files in `/data/rekihaku_downloads`, obtained [here](https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param)) to extract relevant field for  identifying anthropogenic contexts. The file `data/data_prep_c14.R` contains pre-processing R scripts for generating the core dataset used in this paper, stored in the R image file `data/c14data.RData`. Bayesian analyses were carried out using the [nimbleCarbon](https://CRAN.R-project.org/package=nimbleCarbon) R package, which contains custom probability distributions and utility functions for the [NIMBLE](https://cran.r-project.org/web/packages/nimble/index.html) probabilistic programming language. The files `analyses/icar500.R` and `analyses/icar750.R` contain R scripts for fitting the Bayesian models for the 500 and 750 yrs time intervals. Posterior samples are stored in the R image files `results/icar_c14doubleRes500.RData` and `results/icar_c14doubleRes750.RData`. 
Pipelines:
 1. `data/rekihaku_downloads/1200_100_T.xlsx` & `data/rekihaku_downloads/2000_1201_T.xlsx` & `data/rekihaku_downloads/3000_2001_T.xlsx` &  `data/rekihaku_downloads/5000_3001_T.xlsx` &rarr;  `data/rekihaku_downloads/bindCSV.R` &rarr; `data/rekihaku_downloads/binded.csv`.
 2. `data/rekihaku_downloads/binded.csv` & `data/rekihaku_downloads/c14db_1.1.0.csv` &rarr; `data/data_prep_c14.R` &rarr; `data/c14data.RData`.
 3.  `data/c14data.RData` &rarr; `analyses/icar500.R` &rarr; `results/icar_c14doubleRes500.RData`
 4.  `data/c14data.RData` &rarr; `analyses/icar750.R` &rarr; `results/icar_c14doubleRes750.RData`


### Absolute Population Estimates

Absolute population estimates were calculated using a modified version of the equation introduced by Koyama (Koyama, S.1978. Jomon Subsistence and Population. Senri Ethnological Studies, 2, 1–65.) with an updated dataset. The core calculations are included as an R script in the file `analyses/pop_dens_est.R`, with the outputs stored in the CSV files `results/pop_estimate_compare.csv` (comparison of the different estimates discussed in the supplementary materials) and `results/pop_estimate_region.csv` (estimates used in table 1 of the manuscript). Raw input data required for the calculations are Koyama's original data (`data/koyama_popestimate_1984.csv`), the number of archaeological sites [published by the Japanese Agency of Cultural Affairs](https://www.bunka.go.jp/seisaku/bunkazai/shokai/pdf/h29_03_maizotokei.pdf) (`data/maizobunkazai_2017.csv`), and a lookup table for matching administrative units (prefectures) to the regions used in this paper (`data/prefecture_data.csv`). Pipeline: `data/maizobunkazai_2017.csv` & `data/prefecture_data.csv` & `data/koyama_popestimate_1984.csv` &rarr; 
`analyses/pop_dens_est.R` &rarr; `results/pop_estimate_compare.csv` & `results/pop_estimate_region.csv`.

### Figures and Tables
Main (`figures/figure1.pdf` ~ `figures/figure6.pdf`) and supplementary (`figures/figureS1.pdf` ~ `figures/figureS3.pdf`) are generated using the Rscript in `figures/figures_main.R` and `figures/figures_esm.R`.
Pipelines:
 - Main Figures: `data/c14data.RData` & `data/sitedata.RData` & `results/icar_c14doubleRes500.RData` & `results/icar_c14doubleRes750.RData` &rarr;  `figures/figures_main.R` &rarr; `figures/figure1.pdf` ~ `figures/figure6.pdf`
 - Supplementary Figures: `results/pop_estimate_compare.csv` &rarr; `figures/figures_esm.R` &rarr; `figures/figureS1.pdf` ~ `figures/figureS3.pdf`


## File Structure

### data
 * `c14data.RData`
 * `data_prep_c14.R`
 * `data_prep_sites.R`
 * `data_summary.R`
 * `koyama_popestimate_1984.csv`
 * `maizobunkazai_2017.csv`
 * `prefecture_data.csv`
 * `sitedata.RData`
 * `site_raw.csv`

#### /data/rekihaku_downloads
 * `1200_100_T.xlsx`
 * `2000_1201_T.xlsx`
 * `3000_2001_T.xlsx`
 * `5000_3001_T.xlsx`
 * `c14db_1.1.0.csv`
 * `binded.csv`
 * `bindCSV.R`
   
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

### tables
 * `tables_main.R`
 * `tables_esm.R`
 * `table1.csv`
 * `tableS1.csv`
 * `tableS2.csv`

### src
* `dbscanID.R`

## R Session Info

```
attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dbscan_1.1-11       spdep_1.2-8         spData_2.2.2        coda_0.19-4        
 [5] sf_1.0-13           rnaturalearth_0.3.2 latex2exp_0.9.6     RColorBrewer_1.1-3 
 [9] here_1.0.1          rcarbon_1.5.1       nimbleCarbon_0.2.4  nimble_1.0.1       
[13] gridExtra_2.3       dplyr_1.1.2         ggplot2_3.4.2      

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.0       fastmap_1.1.1          spatstat.geom_3.2-2   
 [4] pracma_2.4.2           spatstat.explore_3.2-1 digest_0.6.31         
 [7] rpart_4.1.19           lifecycle_1.0.3        spatstat.data_3.0-1   
[10] magrittr_2.0.3         compiler_4.3.0         rlang_1.1.1           
[13] doSNOW_1.0.20          tools_4.3.0            igraph_1.5.0          
[16] utf8_1.2.3             yaml_2.3.7             knitr_1.43            
[19] sp_1.6-1               classInt_0.4-9         abind_1.4-5           
[22] KernSmooth_2.23-20     withr_2.5.0            purrr_1.0.1           
[25] numDeriv_2016.8-1.1    grid_4.3.0             polyclip_1.10-4       
[28] fansi_1.0.4            e1071_1.7-13           colorspace_2.1-0      
[31] progressr_0.13.0       scales_1.2.1           iterators_1.0.14      
[34] spatstat.utils_3.0-3   spatstat_3.0-6         cli_3.6.1             
[37] rmarkdown_2.21         generics_0.1.3         rstudioapi_0.14       
[40] httr_1.4.6             DBI_1.1.3              proxy_0.4-27          
[43] stringr_1.5.0          splines_4.3.0          spatstat.model_3.2-4  
[46] s2_1.1.4               vctrs_0.6.3            boot_1.3-28.1         
[49] Matrix_1.5-4           jsonlite_1.8.4         tensor_1.5            
[52] elevatr_0.4.5          foreach_1.5.2          units_0.8-2           
[55] snow_0.4-4             goftest_1.2-3          glue_1.6.2            
[58] spatstat.random_3.1-5  codetools_0.2-19       stringi_1.7.12        
[61] gtable_0.3.3           deldir_1.0-9           munsell_0.5.0         
[64] tibble_3.2.1           pillar_1.9.0           htmltools_0.5.5       
[67] R6_2.5.1               wk_0.7.3               rprojroot_2.0.3       
[70] evaluate_0.21          lattice_0.21-8         class_7.3-21          
[73] Rcpp_1.0.11            spatstat.linnet_3.1-1  nlme_3.1-162          
[76] spatstat.sparse_3.0-2  mgcv_1.8-42            xfun_0.39             
[79] pkgconfig_2.0.3 
```


## Funding
This research was funded by the ERC grant _Demography, Cultural Change, and the Diffusion of Rice and Millets during the Jomon-Yayoi transition in prehistoric Japan (ENCOUNTER)_ (Project N. 801953, PI: Enrico Crema) and by a Philip Leverhulme Prize (PLP-2019-304) in archaeology awarded to Enrico Crema.

## Licence
CC-BY 3.0

