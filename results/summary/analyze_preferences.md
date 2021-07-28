Analyze mutant effects and their epistatic turnover
================
Tyler Starr
10/27/2020

-   [Setup](#setup)
-   [Compute correlations in mutant binding constants between background
    and across
    ACE2s](#compute-correlations-in-mutant-binding-constants-between-background-and-across-ace2s)

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","egg","ggseqlogo","bio3d","viridis","ggrepel")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$preferences_dir)){
  dir.create(file.path(config$preferences_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggrepel_0.8.1     viridis_0.5.1     viridisLite_0.3.0 bio3d_2.4-0      
    ##  [5] ggseqlogo_0.1     egg_0.4.5         gridExtra_2.3     forcats_0.4.0    
    ##  [9] stringr_1.4.0     dplyr_0.8.3       purrr_0.3.3       readr_1.3.1      
    ## [13] tidyr_1.0.0       tibble_3.0.2      ggplot2_3.3.0     tidyverse_1.3.0  
    ## [17] data.table_1.12.8 yaml_2.2.0        knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        haven_2.2.0      colorspace_1.4-1
    ##  [5] vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0  rlang_0.4.7     
    ##  [9] pillar_1.4.5     glue_1.3.1       withr_2.1.2      DBI_1.1.0       
    ## [13] dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1     lifecycle_0.2.0 
    ## [17] munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5     
    ## [21] evaluate_0.14    parallel_3.6.2   fansi_0.4.0      broom_0.7.0     
    ## [25] Rcpp_1.0.3       scales_1.1.0     backports_1.1.5  jsonlite_1.6    
    ## [29] fs_1.3.1         hms_0.5.2        digest_0.6.23    stringi_1.4.3   
    ## [33] grid_3.6.2       cli_2.0.0        tools_3.6.2      magrittr_1.5    
    ## [37] crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2      
    ## [41] reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1 rmarkdown_2.0   
    ## [45] httr_1.4.1       rstudioapi_0.10  R6_2.4.1         compiler_3.6.2

Setup
-----

Read in tables of per-mutant and per-homolog phenotypes. Set some
factors for ordering backgrounds and sites.

``` r
dt_mutant <- data.table(read.csv(config$final_variant_scores_mut_file),stringsAsFactors=F)
#order target by order given in config
dt_mutant$target <- factor(dt_mutant$target,levels=config$mutated_targets_ordered)
#order mutant as a factor for grouping by rough biochemical grouping
dt_mutant$mutant <- factor(dt_mutant$mutant, levels=c("C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
#order sites as a factor variable
dt_mutant$position <- factor(dt_mutant$position,levels=c(455,486,493,494,498,501))

dt_wt <- data.table(read.csv(config$final_variant_scores_wt_file),stringsAsFactors=F)
#assign target as a factor in my desired overall plotting order
dt_wt[,target := factor(dt_wt$target,levels=config$targets_ordered)]
```

Compute correlations in mutant binding constants between background and across ACE2s
------------------------------------------------------------------------------------

First, just illustrate correlations with SARS2 for mutations in each
other RBD, for binding to huACE2. Want to report the RMSE, not just the
r/r-squared – because we are interested in the absolute magnitude of
deviation between backgrounds, not just the relative variance that’s
normalized in an r-squared value!

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,huACE2]
    y <- dt_mutant[target==RBD2 & position==site,huACE2]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_huACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_huACE2.pdf",sep="")))
```

RBDs versus SARS1, huACE2

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-1_Urbani_HP03L"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-1_Urbani_HP03L" & position==site,huACE2]
    y <- dt_mutant[target==RBD2 & position==site,huACE2]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-1_Urbani_HP03L" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS1_huACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS1_huACE2.pdf",sep="")))
```

And, vs. SARS2 for the rest of the ACE2s

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,cvACE2]
    y <- dt_mutant[target==RBD2 & position==site,cvACE2]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_cvACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_cvACE2.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,pgACE2]
    y <- dt_mutant[target==RBD2 & position==site,pgACE2]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_pgACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_pgACE2.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,mACE2]
    y <- dt_mutant[target==RBD2 & position==site,mACE2]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_mACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_mACE2.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,RaACE2.9479]
    y <- dt_mutant[target==RBD2 & position==site,RaACE2.9479]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_RaACE2.9479-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_RaACE2.9479.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,RaACE2.787]
    y <- dt_mutant[target==RBD2 & position==site,RaACE2.787]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_RaACE2.787-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_RaACE2.787.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,RsACE2.1434]
    y <- dt_mutant[target==RBD2 & position==site,RsACE2.1434]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_RsACE2.1434-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_RsACE2.1434.pdf",sep="")))
```

``` r
par(mfrow=c(13, 6))
for(RBD2 in levels(dt_mutant$target)[levels(dt_mutant$target)!="SARS-CoV-2"]){
  for(site in c(455, 486, 493, 494, 498, 501)){
    x <- dt_mutant[target=="SARS-CoV-2" & position==site,RsACE2.3364]
    y <- dt_mutant[target==RBD2 & position==site,RsACE2.3364]
    plot(x,y, xlim=c(5,12),ylim=c(5,12), pch=as.character(dt_mutant[target=="SARS-CoV-2" & position==site, mutant]), xlab="", ylab=RBD2, main="")
    legend("topleft",bty="n",cex=1,legend=format(mean(abs(lm(y~x)$residuals),na.rm=T),digits=3))
  }
}
```

<img src="analyze_preferences_files/figure-gfm/correlation_plots_v_SARS2_RsACE2.3364-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/correlations-by-site_SARS2_RsACE2.3364.pdf",sep="")))
```

To summarize these plots, I want to make ‘reaction coordinate’ plots
that visualize R-squared in mut binding affinites as a function of
pairwise sequence divergence. Want to only use RBD-ACE2s if the wildtype
binds wiht &gt;7 binding constant, so that there aren’t non-correlated
samples simply because nothing is binding (no variation to correlate)

Set up tables for storing correlations and identities for all ACE2- and
RBD-pairs

``` r
#data table for storing, for each RBD, its correlation in binding for pairs of ACE2s
#generate table with all combinations of ACE2_1 and ACE2_2 for each RBD+site
diffs_ACE2 <- data.table(expand.grid(site=c(455,486,493,494,498,501),ACE2_2=c("huACE2","cvACE2","pgACE2","mACE2","RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434"),ACE2_1=c("huACE2","cvACE2","pgACE2","mACE2","RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434"),RBD=as.character(unique(dt_mutant$target))))

#remove duplicates -- either ACE2_1 and _2 the same, or combinations where the _2 _1 combo is already present in the _1 _2 orientation
diffs_ACE2 <- diffs_ACE2[ACE2_1!=ACE2_2,]
diffs_ACE2 <- diffs_ACE2[ACE2_2!="huACE2",]
diffs_ACE2 <- diffs_ACE2[ACE2_1=="huACE2" | (ACE2_1=="cvACE2" & ACE2_2 %in% c("pgACE2","mACE2","RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434")) | (ACE2_1=="pgACE2" & ACE2_2 %in% c("mACE2","RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434")) | (ACE2_1=="mACE2" & ACE2_2 %in% c("RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434")) | (ACE2_1=="RaACE2.787" & ACE2_2 %in% c("RaACE2.9479","RsACE2.3364","RsACE2.1434")) | (ACE2_1=="RaACE2.9479" & ACE2_2 %in% c("RsACE2.3364","RsACE2.1434")) | (ACE2_1=="RsACE2.3364" & ACE2_2 %in% c("RsACE2.1434")),]

#remove RBD-sites where one of the ACE2s isn't bound with minimum Kd for the wildtype RBD
min_Kd <- 7
diffs_ACE2[,keep:=as.logical(NA)]
for(i in 1:nrow(diffs_ACE2)){
  ACE2_1 <- as.character(diffs_ACE2[i,ACE2_1])
  ACE2_2 <- as.character(diffs_ACE2[i,ACE2_2])
  bind_1 <- dt_wt[target==as.character(diffs_ACE2[i,RBD]),get(ACE2_1)]
  bind_2 <- dt_wt[target==as.character(diffs_ACE2[i,RBD]),get(ACE2_2)]
  if(bind_1 < min_Kd | bind_2 < min_Kd){
    diffs_ACE2[i,keep := F]
  }else if(bind_1 >= min_Kd & bind_2 >= min_Kd){
    diffs_ACE2[i,keep := T]
  }
}
diffs_ACE2 <- diffs_ACE2[keep==T,.(RBD, ACE2_1, ACE2_2, site)]

#data table for storing, for each ACE2, its difference in correlation in profiles between RBD pairs at each site. 
#generate table with all combinations of RBD_1 and RBD_2 for each ACE2+site
diffs_RBD <- data.table(expand.grid(site=c(455,486,493,494,498,501),RBD_2=as.character(unique(dt_mutant$target)),RBD_1=as.character(unique(dt_mutant$target)),ACE2=c("huACE2","cvACE2","pgACE2","mACE2","RaACE2.787","RaACE2.9479","RsACE2.3364","RsACE2.1434")))

#remove duplicates -- either RBD_1 and _2 the same, or combinations where the _2 _1 combo is already present in the _1 _2 orientation
diffs_RBD <- diffs_RBD[RBD_1!=RBD_2,]
diffs_RBD <- diffs_RBD[RBD_2!="AncSarbecovirus_MAP",]
diffs_RBD <- diffs_RBD[RBD_1=="AncSarbecovirus_MAP" | (RBD_1=="BM48-31" & RBD_2 %in% c("BtKY72","AncAsia_MAP","AncSARS2a_MAP","AncSARS2c_MAP","SARS-CoV-2","RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="BtKY72" & RBD_2 %in% c("AncAsia_MAP","AncSARS2a_MAP","AncSARS2c_MAP","SARS-CoV-2","RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="AncAsia_MAP" & RBD_2 %in% c("AncSARS2a_MAP","AncSARS2c_MAP","SARS-CoV-2","RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="AncSARS2a_MAP" & RBD_2 %in% c("AncSARS2c_MAP","SARS-CoV-2","RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="AncSARS2c_MAP" & RBD_2 %in% c("SARS-CoV-2","RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="SARS-CoV-2" & RBD_2 %in% c("RaTG13","GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="RaTG13" & RBD_2 %in% c("GD-Pangolin","AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="GD-Pangolin" & RBD_2 %in% c("AncSARS1a_MAP","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="AncSARS1a_MAP" & RBD_2 %in% c("SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="SARS-CoV-1_Urbani_HP03L" & RBD_2 %in% c("SARS-CoV-1_PC4-137_PC04","Rs7327","AncClade2_MAP")) | (RBD_1=="SARS-CoV-1_PC4-137_PC04" & RBD_2 %in% c("Rs7327","AncClade2_MAP")) | (RBD_1=="Rs7327" & RBD_2 %in% c("AncClade2_MAP")),]


#remove RBD-sites where one of the ACE2s isn't bound with minimum Kd for the wildtype RBD
min_Kd <- 7
diffs_RBD[,keep:=as.logical(NA)]
for(i in 1:nrow(diffs_RBD)){
  RBD_1 <- as.character(diffs_RBD[i,RBD_1])
  RBD_2 <- as.character(diffs_RBD[i,RBD_2])
  bind_1 <- dt_wt[target==RBD_1,get(as.character(diffs_RBD[i,ACE2]))]
  bind_2 <- dt_wt[target==RBD_2,get(as.character(diffs_RBD[i,ACE2]))]
  if(bind_1 < min_Kd | bind_2 < min_Kd){
    diffs_RBD[i,keep := F]
  }else if(bind_1 >= min_Kd & bind_2 >= min_Kd){
    diffs_RBD[i,keep := T]
  }
}
diffs_RBD <- diffs_RBD[keep==T,.(ACE2, RBD_1, RBD_2, site)]
```

Loop through to compute residual mean square error, mean absolute error,
and r in correlations of affinities for each state at each site in each
RBD-pair.

``` r
diffs_ACE2$cor <- as.numeric(NA)
for(i in 1:nrow(diffs_ACE2)){
  x <- dt_mutant[target==diffs_ACE2[i,RBD] & position==diffs_ACE2[i,site],get(as.character(diffs_ACE2[i,ACE2_1]))]
  y <- dt_mutant[target==diffs_ACE2[i,RBD] & position==diffs_ACE2[i,site],get(as.character(diffs_ACE2[i,ACE2_2]))]
  diffs_ACE2[i,cor:=cor(x,y,use="complete.obs", method="pearson")]
}

diffs_RBD$cor <- as.numeric(NA)
diffs_RBD$RMSE <- as.numeric(NA)
diffs_ACE2$MAE <- as.numeric(NA)
for(i in 1:nrow(diffs_RBD)){
  x <- dt_mutant[target==diffs_RBD[i,RBD_1] & position==diffs_RBD[i,site],get(as.character(diffs_RBD[i,ACE2]))]
  y <- dt_mutant[target==diffs_RBD[i,RBD_2] & position==diffs_RBD[i,site],get(as.character(diffs_RBD[i,ACE2]))]
  diffs_RBD[i,cor:=cor(x,y,use="complete.obs", method="pearson")]
  diffs_RBD[i,RMSE:=sqrt(mean(lm(y~x)$residuals^2,na.rm=T))]
  diffs_RBD[i,MAE:=mean(abs(lm(y~x)$residuals),na.rm=T)]
}
```

For RBD pairs, compare MAE versus pairwise RBD sequence identity across
comparisons.

First, read in alignment and add new column to table indicating %
identity

``` r
alignment <- bio3d::read.fasta(file="data/RBD_align_SSM-backgrounds.fasta")
ids <- seqidentity(alignment)

diffs_RBD[,percent_ID:=ids[as.character(RBD_1),as.character(RBD_2)],by=c("RBD_1","RBD_2")]
```

And do plots – for just human ACE2 binding, and then including all ACE2s

``` r
p_hu<- ggplot(data=diffs_RBD[ACE2=="huACE2",],aes(x=percent_ID, y=MAE))+
  geom_point(pch=16,alpha=0.25)+
  geom_smooth(method="loess",span=1)+
  facet_wrap(~site,nrow=2)+
  theme_classic()+
  xlim(1,0.7)+
  xlab("RBD pairwise sequence identity")+
  ylab("Mean absolute error")
  
p_hu
```

    ## `geom_smooth()` using formula 'y ~ x'

<img src="analyze_preferences_files/figure-gfm/RBD-MAE-versus-percent-ID_huACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/mae-v-percent-id_by-site_all-huACE2.pdf",sep=""),useDingbats=F))
```

``` r
p_all <- ggplot(data=diffs_RBD,aes(x=percent_ID, y=MAE))+
  geom_point(pch=16,alpha=0.25)+
  geom_smooth(method="loess",span=1)+
  facet_wrap(~site,nrow=2)+
  theme_classic()+xlim(1,0.7)+
  xlab("RBD pairwise sequence identity")+
  ylab("Mean absolute error")
  
p_all
```

    ## `geom_smooth()` using formula 'y ~ x'

<img src="analyze_preferences_files/figure-gfm/RBD-MAE-versus-percent-ID_all-ACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/mae-v-percent-id_by-site_all-ACE2.pdf",sep=""),useDingbats=F))
```

And, change in actual correlation constants:

``` r
p_hu<- ggplot(data=diffs_RBD[ACE2=="huACE2",],aes(x=percent_ID, y=cor))+
  geom_point(pch=16,alpha=0.25)+
  geom_smooth(method="loess",span=1.3)+
  facet_wrap(~site,nrow=2)+
  theme_classic()+
  xlim(1,0.7)+
  xlab("RBD pairwise sequence identity")+
  ylab("Pearson's r")
  
p_hu
```

    ## `geom_smooth()` using formula 'y ~ x'

<img src="analyze_preferences_files/figure-gfm/RBD-cor-versus-percent-ID_huACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/r-v-percent-id_by-site_all-huACE2.pdf",sep=""),useDingbats=F))
```

``` r
p_all <- ggplot(data=diffs_RBD,aes(x=percent_ID, y=cor))+
  geom_point(pch=16,alpha=0.25)+
  geom_smooth(method="loess",span=1)+
  facet_wrap(~site,nrow=2)+
  theme_classic()+xlim(1,0.7)+
  xlab("RBD pairwise sequence identity")+
  ylab("Pearson's r")
  
p_all
```

    ## `geom_smooth()` using formula 'y ~ x'

<img src="analyze_preferences_files/figure-gfm/RBD-cor-versus-percent-ID_all-ACE2-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/r-v-percent-id_by-site_all-ACE2.pdf",sep=""),useDingbats=F))
```

For ACE2 pairs, let’s take the median scorrelation coefficient of each
RBDs ACE2-pair correlation in mut effects per site

``` r
diffs_ACE2[,median_cor:=median(cor),by=c("ACE2_1","ACE2_2","site")]
diffs_ACE2_collapse <- unique(diffs_ACE2[,.(ACE2_1,ACE2_2,site,median_cor)])

pACE2 <- ggplot(diffs_ACE2,aes(ACE2_2,ACE2_1))+geom_tile(aes(fill=median_cor),color="black",lwd=0.1)+
  #scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3/14,5/14,7/14,9/14,11/14,14/14),na.value="yellow")+
  scale_fill_viridis(na.value="white")+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~site,nrow=2)
  #guides(y.sec=guide_axis_label_trans())+

pACE2
```

<img src="analyze_preferences_files/figure-gfm/ACE2_correlations-1.png" style="display: block; margin: auto;" />

Can we identify mouse-adaptive mutations that are specific to mouse but
don’t enhance human ACE2 binding?

``` r
ggplot(data=dt_mutant[target %in% c("BM48-31","BtKY72","GD-Pangolin","SARS-CoV-2","RaTG13","Rs7327","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04"),],aes(x=huACE2,y=mACE2))+
#ggplot(data=dt_mutant,aes(x=huACE2,y=mACE2))+
  geom_point(pch=16,alpha=0.5)+
  facet_wrap(~target,nrow=2)+
  geom_text_repel(aes(label=ifelse(((mACE2 > huACE2+1.5 & huACE2 < 7) | (mACE2 > 6 & huACE2 < 6)),as.character(paste(wildtype,position,mutant,sep="")),'')),size=3,color="gray40")+
  #geom_text_repel(aes(label=ifelse(((mACE2 > huACE2-2 & mACE2>6)),as.character(paste(wildtype,position,mutant,sep="")),'')),size=3,color="gray40")+
  geom_point(data=dt_mutant[target %in% c("BM48-31","BtKY72","GD-Pangolin","SARS-CoV-2","RaTG13","Rs7327","SARS-CoV-1_Urbani_HP03L","SARS-CoV-1_PC4-137_PC04") & as.character(wildtype)==as.character(mutant) & position==493,],pch=18,col="red",size=4)+
  theme_classic()
```

    ## Warning: Removed 7 rows containing missing values (geom_point).

    ## Warning: Removed 7 rows containing missing values (geom_text_repel).

<img src="analyze_preferences_files/figure-gfm/correlation_mACE2_SARS2-muts-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$preferences_dir,"/mACE2_v_huACE2_by_bg.pdf",sep=""),useDingbats=F))
```