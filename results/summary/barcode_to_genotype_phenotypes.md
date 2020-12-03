Compute per-barcode RBD expression score
================
Tyler Starr
10/11/2020

-   [Setup](#setup)
-   [Calculate per-variant mean
    scores](#calculate-per-variant-mean-scores)
-   [Heatmaps! Output heatmaps illustrating all wildtype variants with
    separate columns for each
    ACE2.](#heatmaps-output-heatmaps-illustrating-all-wildtype-variants-with-separate-columns-for-each-ace2.)

This notebook reads in the per-barcode titration Kds and expression
measurements from the `compute_binding_Kd` and
`compute_expression_meanF` scripts. It synthesizes these two sets of
results and calculates the final ‘mean’ phenotypes for each variant, and
generates some coverage and QC analyses.

    require("knitr")
    knitr::opts_chunk$set(echo = T)
    knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

    #list of packages to install/load
    packages = c("yaml","data.table","tidyverse","gridExtra","egg")
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
    if(!file.exists(config$final_variant_scores_dir)){
      dir.create(file.path(config$final_variant_scores_dir))
    }

Session info for reproducing environment:

    sessionInfo()

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
    ##  [1] egg_0.4.5         gridExtra_2.3     forcats_0.4.0     stringr_1.4.0    
    ##  [5] dplyr_0.8.3       purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
    ##  [9] tibble_3.0.2      ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8
    ## [13] yaml_2.2.0        knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        haven_2.2.0      colorspace_1.4-1
    ##  [5] vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0  rlang_0.4.7     
    ##  [9] pillar_1.4.5     glue_1.3.1       withr_2.1.2      DBI_1.1.0       
    ## [13] dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1     lifecycle_0.2.0 
    ## [17] munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5     
    ## [21] evaluate_0.14    fansi_0.4.0      broom_0.7.0      Rcpp_1.0.3      
    ## [25] scales_1.1.0     backports_1.1.5  jsonlite_1.6     fs_1.3.1        
    ## [29] hms_0.5.2        digest_0.6.23    stringi_1.4.3    grid_3.6.2      
    ## [33] cli_2.0.0        tools_3.6.2      magrittr_1.5     crayon_1.3.4    
    ## [37] pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2       reprex_0.3.0    
    ## [41] lubridate_1.7.4  assertthat_0.2.1 rmarkdown_2.0    httr_1.4.1      
    ## [45] rstudioapi_0.10  R6_2.4.1         compiler_3.6.2

Setup
-----

Read in tables of per-barcode expression and binding Kd measurements and
combine.

    dt_bind <- data.table(read.csv(config$Titeseq_Kds_file),stringsAsFactors=F)
    dt_expr <- data.table(read.csv(config$expression_sortseq_file),stringsAsFactors=F)

    dt <- merge(dt_bind,dt_expr)

    #assign target as a factor in my desired overall plotting order
    dt[,target := factor(dt$target,levels=config$targets_ordered)]

Calculate per-variant mean scores
---------------------------------

Unfiltered, look at distribution of binding and expression

    p0 <- ggplot(dt[!is.na(ML_meanF) & variant_class=="wildtype",],aes(x=target,y=ML_meanF))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("expression")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

    p1 <- ggplot(dt[!is.na(log10Ka_huACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_huACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("huACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p1.6 <- ggplot(dt[!is.na(log10Ka_huACE2.pool6) & variant_class=="wildtype",],aes(x=target,y=log10Ka_huACE2.pool6))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("huACE2.pool6")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p2 <- ggplot(dt[!is.na(log10Ka_cvACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_cvACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("cvACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p3 <- ggplot(dt[!is.na(log10Ka_RaACE2.787) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RaACE2.787))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RaACE2.787")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p3.6 <- ggplot(dt[!is.na(log10Ka_RaACE2.787.pool6) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RaACE2.787.pool6))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RaACE2.787.pool6")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p4 <- ggplot(dt[!is.na(log10Ka_RsACE2.3364) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RsACE2.3364))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RsACE2.3364")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p4.0 <- ggplot(dt[!is.na(log10Ka_RsACE2old) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RsACE2old))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RsACE2old")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.99,12))

    p5 <- ggplot(dt[!is.na(log10Ka_RpACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RpACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RpACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(5.99,12))

    grid.arrange(p0,p1,p1.6,p2,p3,p3.6,p4,p4.0,p5,ncol=1)

    ## Warning: Removed 3 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 3 rows containing non-finite values (stat_summary).

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/unfiltered_expression_binding-1.png" style="display: block; margin: auto;" />

Let’s add a variable that flags the top and bottom 5% of expression
scores for each variant/mutant, and see how violin plots look when
censoring these top and bottom 5% of expressed barcodes

    dt[,expression_censor_lower:=quantile(ML_meanF,0.05,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,expression_censor_upper:=quantile(ML_meanF,0.95,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    p1 <- ggplot(dt[!is.na(ML_meanF) & variant_class=="wildtype" & ML_meanF>=expression_censor_lower & ML_meanF<=expression_censor_upper,],aes(x=target,y=ML_meanF))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("expression")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

    grid.arrange(p1,ncol=1)

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/censor_5_expression-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/vioplots_wt_expression.pdf",sep="")))

Let’s do the same for binding, with 2.5% cutoffs.

    dt[,huACE2_censor_lower:=quantile(log10Ka_huACE2,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,huACE2_censor_upper:=quantile(log10Ka_huACE2,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,huACE2.pool6_censor_lower:=quantile(log10Ka_huACE2.pool6,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,huACE2.pool6_censor_upper:=quantile(log10Ka_huACE2.pool6,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,cvACE2_censor_lower:=quantile(log10Ka_cvACE2,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,cvACE2_censor_upper:=quantile(log10Ka_cvACE2,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,RaACE2.787_censor_lower:=quantile(log10Ka_RaACE2.787,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,RaACE2.787_censor_upper:=quantile(log10Ka_RaACE2.787,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,RaACE2.787.pool6_censor_lower:=quantile(log10Ka_RaACE2.787.pool6,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,RaACE2.787.pool6_censor_upper:=quantile(log10Ka_RaACE2.787.pool6,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,RsACE2.3364_censor_lower:=quantile(log10Ka_RsACE2.3364,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,RsACE2.3364_censor_upper:=quantile(log10Ka_RsACE2.3364,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,RsACE2old_censor_lower:=quantile(log10Ka_RsACE2old,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,RsACE2old_censor_upper:=quantile(log10Ka_RsACE2old,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]

    dt[,RpACE2_censor_lower:=quantile(log10Ka_RpACE2,0.025,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]
    dt[,RpACE2_censor_upper:=quantile(log10Ka_RpACE2,0.975,na.rm=T,type=7),by=c("library","target","variant_class","position","mutant")]


    p1 <- ggplot(dt[!is.na(log10Ka_huACE2) & variant_class=="wildtype" & log10Ka_huACE2>=huACE2_censor_lower & log10Ka_huACE2<=huACE2_censor_upper,],aes(x=target,y=log10Ka_huACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("huACE2 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p1.6 <- ggplot(dt[!is.na(log10Ka_huACE2.pool6) & variant_class=="wildtype" & log10Ka_huACE2.pool6>=huACE2.pool6_censor_lower & log10Ka_huACE2.pool6<=huACE2.pool6_censor_upper,],aes(x=target,y=log10Ka_huACE2.pool6))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("huACE2.pool6 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p2 <- ggplot(dt[!is.na(log10Ka_cvACE2) & variant_class=="wildtype" & log10Ka_cvACE2>=cvACE2_censor_lower & log10Ka_cvACE2<=cvACE2_censor_upper,],aes(x=target,y=log10Ka_cvACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("cvACE2 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p3 <- ggplot(dt[!is.na(log10Ka_RaACE2.787) & variant_class=="wildtype" & log10Ka_RaACE2.787>=RaACE2.787_censor_lower & log10Ka_RaACE2.787<=RaACE2.787_censor_upper,],aes(x=target,y=log10Ka_RaACE2.787))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RaACE2.787 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p3.6 <- ggplot(dt[!is.na(log10Ka_RaACE2.787.pool6) & variant_class=="wildtype" & log10Ka_RaACE2.787.pool6>=RaACE2.787.pool6_censor_lower & log10Ka_RaACE2.787.pool6<=RaACE2.787.pool6_censor_upper,],aes(x=target,y=log10Ka_RaACE2.787.pool6))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RaACE2.787.pool6 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p4 <- ggplot(dt[!is.na(log10Ka_RsACE2.3364) & variant_class=="wildtype" & log10Ka_RsACE2.3364>=RsACE2.3364_censor_lower & log10Ka_RsACE2.3364<=RsACE2.3364_censor_upper,],aes(x=target,y=log10Ka_RsACE2.3364))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RsACE2.3364 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p4.0 <- ggplot(dt[!is.na(log10Ka_RsACE2old) & variant_class=="wildtype" & log10Ka_RsACE2old>=RsACE2old_censor_lower & log10Ka_RsACE2old<=RsACE2old_censor_upper,],aes(x=target,y=log10Ka_RsACE2old))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RsACE2old affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(4.9,12))

    p5 <- ggplot(dt[!is.na(log10Ka_RpACE2) & variant_class=="wildtype" & log10Ka_RpACE2>=RpACE2_censor_lower & log10Ka_RpACE2<=RpACE2_censor_upper,],aes(x=target,y=log10Ka_RpACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RpACE2 affinity")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(5.9,12))

    grid.arrange(p1,p1.6,p2,p3,p3.6,p4,p4.0,p5,ncol=1)

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/censor_2.5_binding-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/vioplots_wt_log10Ka.pdf",sep="")))

Calculate the mean per variant, the standard deviation, and the number
of (post-filter) barcodes on which a variant score was determined

    #apply the censors to NA out the phenotypes outside the range
    dt[log10Ka_huACE2 < huACE2_censor_lower | log10Ka_huACE2 > huACE2_censor_upper,log10Ka_huACE2:=NA]
    dt[log10Ka_huACE2.pool6 < huACE2.pool6_censor_lower | log10Ka_huACE2.pool6 > huACE2.pool6_censor_upper,log10Ka_huACE2.pool6:=NA]
    dt[log10Ka_cvACE2 < cvACE2_censor_lower | log10Ka_cvACE2 > cvACE2_censor_upper,log10Ka_cvACE2:=NA]
    dt[log10Ka_RaACE2.787 < RaACE2.787_censor_lower | log10Ka_RaACE2.787 > RaACE2.787_censor_upper,log10Ka_RaACE2.787:=NA]
    dt[log10Ka_RaACE2.787.pool6 < RaACE2.787.pool6_censor_lower | log10Ka_RaACE2.787.pool6 > RaACE2.787.pool6_censor_upper,log10Ka_RaACE2.787.pool6:=NA]
    dt[log10Ka_RsACE2.3364 < RsACE2.3364_censor_lower | log10Ka_RsACE2.3364 > RsACE2.3364_censor_upper,log10Ka_RsACE2.3364:=NA]
    dt[log10Ka_RsACE2old < RsACE2old_censor_lower | log10Ka_RsACE2old > RsACE2old_censor_upper,log10Ka_RsACE2old:=NA]
    dt[log10Ka_RpACE2 < RpACE2_censor_lower | log10Ka_RpACE2 > RpACE2_censor_upper,log10Ka_RpACE2:=NA]

    dt[,mean_huACE2_log10Ka:=mean(log10Ka_huACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_huACE2_log10Ka:=sd(log10Ka_huACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_huACE2_log10Ka:=sum(!is.na(log10Ka_huACE2)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_huACE2.pool6_log10Ka:=mean(log10Ka_huACE2.pool6,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_huACE2.pool6_log10Ka:=sd(log10Ka_huACE2.pool6,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_huACE2.pool6_log10Ka:=sum(!is.na(log10Ka_huACE2.pool6)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_cvACE2_log10Ka:=mean(log10Ka_cvACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_cvACE2_log10Ka:=sd(log10Ka_cvACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_cvACE2_log10Ka:=sum(!is.na(log10Ka_cvACE2)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_RaACE2.787_log10Ka:=mean(log10Ka_RaACE2.787,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_RaACE2.787_log10Ka:=sd(log10Ka_RaACE2.787,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_RaACE2.787_log10Ka:=sum(!is.na(log10Ka_RaACE2.787)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_RaACE2.787.pool6_log10Ka:=mean(log10Ka_RaACE2.787.pool6,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_RaACE2.787.pool6_log10Ka:=sd(log10Ka_RaACE2.787.pool6,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_RaACE2.787.pool6_log10Ka:=sum(!is.na(log10Ka_RaACE2.787.pool6)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_RsACE2.3364_log10Ka:=mean(log10Ka_RsACE2.3364,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_RsACE2.3364_log10Ka:=sd(log10Ka_RsACE2.3364,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_RsACE2.3364_log10Ka:=sum(!is.na(log10Ka_RsACE2.3364)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_RsACE2old_log10Ka:=mean(log10Ka_RsACE2old,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_RsACE2old_log10Ka:=sd(log10Ka_RsACE2old,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_RsACE2old_log10Ka:=sum(!is.na(log10Ka_RsACE2old)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_RpACE2_log10Ka:=mean(log10Ka_RpACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_RpACE2_log10Ka:=sd(log10Ka_RpACE2,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_RpACE2_log10Ka:=sum(!is.na(log10Ka_RpACE2)),by=c("library","target","variant_class","position","mutant")]

    dt[,mean_expression:=mean(ML_meanF,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,sd_expression:=sd(ML_meanF,na.rm=T),by=c("library","target","variant_class","position","mutant")]
    dt[,n_expression:=sum(!is.na(ML_meanF)),by=c("library","target","variant_class","position","mutant")]

Collapse down to tables reporting just the summary statistics for each
genotype.

    dt_wildtype <- dt[variant_class=="wildtype",.(library,target,
                                                  mean_expression,sd_expression,n_expression,
                                                  mean_huACE2_log10Ka,sd_huACE2_log10Ka,n_huACE2_log10Ka,
                                                  mean_huACE2.pool6_log10Ka,sd_huACE2.pool6_log10Ka,n_huACE2.pool6_log10Ka,
                                                  mean_cvACE2_log10Ka,sd_cvACE2_log10Ka,n_cvACE2_log10Ka,
                                                  mean_RaACE2.787_log10Ka,sd_RaACE2.787_log10Ka,n_RaACE2.787_log10Ka,
                                                  mean_RaACE2.787.pool6_log10Ka,sd_RaACE2.787.pool6_log10Ka,n_RaACE2.787.pool6_log10Ka,
                                                  mean_RsACE2.3364_log10Ka,sd_RsACE2.3364_log10Ka,n_RsACE2.3364_log10Ka,
                                                  mean_RsACE2old_log10Ka,sd_RsACE2old_log10Ka,n_RsACE2old_log10Ka,
                                                  mean_RpACE2_log10Ka,sd_RpACE2_log10Ka,n_RpACE2_log10Ka)]
    dt_mutant <- dt[variant_class=="mutant",.(library,target,wildtype,position,mutant,
                                              mean_expression,sd_expression,n_expression,
                                              mean_huACE2_log10Ka,sd_huACE2_log10Ka,n_huACE2_log10Ka,
                                              mean_huACE2.pool6_log10Ka,sd_huACE2.pool6_log10Ka,n_huACE2.pool6_log10Ka,
                                              mean_cvACE2_log10Ka,sd_cvACE2_log10Ka,n_cvACE2_log10Ka,
                                              mean_RaACE2.787_log10Ka,sd_RaACE2.787_log10Ka,n_RaACE2.787_log10Ka,
                                              mean_RaACE2.787.pool6_log10Ka,sd_RaACE2.787.pool6_log10Ka,n_RaACE2.787.pool6_log10Ka,
                                              mean_RsACE2.3364_log10Ka,sd_RsACE2.3364_log10Ka,n_RsACE2.3364_log10Ka,
                                              mean_RsACE2old_log10Ka,sd_RsACE2old_log10Ka,n_RsACE2old_log10Ka,
                                              mean_RpACE2_log10Ka,sd_RpACE2_log10Ka,n_RpACE2_log10Ka)]


    dt_wildtype <- unique(dt_wildtype); setkey(dt_wildtype,library, target)
    dt_mutant <- unique(dt_mutant)

    aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
    #fill out missing values in table with a hideous loop. If you are somebody who is reading this code, I apologize.
    for(lib in c("lib1","lib2")){
      for(bg in unique(dt_mutant$target)){
        for(pos in unique(dt_mutant$position)){
          for(aa in aas){
            if(!(aa %in% as.character(dt_mutant[library==lib & target==bg & position==pos,mutant]))){
              dt_mutant <- rbind(dt_mutant,list(lib, bg, dt_mutant[library==lib & target==bg & position==pos,wildtype][1],pos,aa),fill=T)
            }
          }
        }
      }
    }
    setkey(dt_mutant,library,target,position,mutant)

We have duplicates for expression, huACE2, cvACE2, RaACE2.787,
RsACE2.3364. Let’s look at correlations!

    par(mfrow=c(5,2))
    x <- dt_wildtype[library=="lib1",mean_expression]; y <- dt_wildtype[library=="lib2",mean_expression]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="expression, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    x <- dt_mutant[library=="lib1",mean_expression]; y <- dt_mutant[library=="lib2",mean_expression]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="expression, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    #huACE2
    x <- dt_wildtype[library=="lib1" & mean_expression > 7, mean_huACE2_log10Ka]; y <- dt_wildtype[library=="lib2"& mean_expression > 7,mean_huACE2_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="huACE2 log10Ka, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n"); points(dt_wildtype[library=="lib1" & mean_expression < 7, mean_huACE2_log10Ka],dt_wildtype[library=="lib2" & mean_expression < 7, mean_huACE2_log10Ka],pch=19,col="green")

    x <- dt_mutant[library=="lib1",mean_huACE2_log10Ka]; y <- dt_mutant[library=="lib2",mean_huACE2_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="huACE2 log10Ka, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    #cvACE2
    x <- dt_wildtype[library=="lib1" & mean_expression > 7, mean_cvACE2_log10Ka]; y <- dt_wildtype[library=="lib2"& mean_expression > 7,mean_cvACE2_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="cvACE2 log10Ka, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n"); points(dt_wildtype[library=="lib1" & mean_expression < 7, mean_cvACE2_log10Ka],dt_wildtype[library=="lib2" & mean_expression < 7, mean_cvACE2_log10Ka],pch=19,col="green")

    x <- dt_mutant[library=="lib1",mean_cvACE2_log10Ka]; y <- dt_mutant[library=="lib2",mean_cvACE2_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="cvACE2 log10Ka, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    #RaACE2.787
    x <- dt_wildtype[library=="lib1" & mean_expression > 7, mean_RaACE2.787_log10Ka]; y <- dt_wildtype[library=="lib2"& mean_expression > 7,mean_RaACE2.787_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="RaACE2.787 log10Ka, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n"); points(dt_wildtype[library=="lib1" & mean_expression < 7, mean_RaACE2.787_log10Ka],dt_wildtype[library=="lib2" & mean_expression < 7, mean_RaACE2.787_log10Ka],pch=19,col="green")

    x <- dt_mutant[library=="lib1",mean_RaACE2.787_log10Ka]; y <- dt_mutant[library=="lib2",mean_RaACE2.787_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="RaACE2.787 log10Ka, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    #RsACE2.3364
    x <- dt_wildtype[library=="lib1" & mean_expression > 7, mean_RsACE2.3364_log10Ka]; y <- dt_wildtype[library=="lib2"& mean_expression > 7,mean_RsACE2.3364_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="RsACE2.3364 log10Ka, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n"); points(dt_wildtype[library=="lib1" & mean_expression < 7, mean_RsACE2.3364_log10Ka],dt_wildtype[library=="lib2" & mean_expression < 7, mean_RsACE2.3364_log10Ka],pch=19,col="green")

    x <- dt_mutant[library=="lib1",mean_RsACE2.3364_log10Ka]; y <- dt_mutant[library=="lib2",mean_RsACE2.3364_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="replicate 1",ylab="replicate 2",main="RsACE2.3364 log10Ka, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/plot_correlations-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/replicate_correlations.pdf",sep="")))

See how the “old” RsACE2 compares to new, can drop if seems well
correlated

    par(mfrow=c(1,2))
    x <- dt_wildtype[mean_expression > 7, mean_RsACE2.3364_log10Ka]; y <- dt_wildtype[mean_expression > 7,mean_RsACE2old_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="new RsACE2",ylab="old ACE2",main="RsACE2 log10Ka old v new, wildtype backgrounds");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

    x <- dt_mutant[,mean_RsACE2.3364_log10Ka]; y <- dt_mutant[,mean_RsACE2old_log10Ka]; plot(x,y,pch=19,col="#00000090",xlab="new RsACE2",ylab="old RsACE2",main="RsACE2 new v old log10Ka, miniMS muts");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/compare_RsACE2s-1.png" style="display: block; margin: auto;" />

Collapse down to mean from both replicates, and total n barcodes between
the two replicates.

    dt_wildtype[,huACE2:=mean(mean_huACE2_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,huACE2_n:=sum(n_huACE2_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,huACE2.pool6:=mean(mean_huACE2.pool6_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,huACE2.pool6_n:=sum(n_huACE2.pool6_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,cvACE2:=mean(mean_cvACE2_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,cvACE2_n:=sum(n_cvACE2_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,RaACE2.787:=mean(mean_RaACE2.787_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,RaACE2.787_n:=sum(n_RaACE2.787_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,RaACE2.787.pool6:=mean(mean_RaACE2.787.pool6_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,RaACE2.787.pool6_n:=sum(n_RaACE2.787.pool6_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,RsACE2.3364:=mean(mean_RsACE2.3364_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,RsACE2.3364_n:=sum(n_RsACE2.3364_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,RpACE2:=mean(mean_RpACE2_log10Ka,na.rm=T),by=c("target")]
    dt_wildtype[,RpACE2_n:=sum(n_RpACE2_log10Ka,na.rm=T),by=c("target")]

    dt_wildtype[,expression:=mean(mean_expression,na.rm=T),by=c("target")]
    dt_wildtype[,expression_n:=sum(n_expression,na.rm=T),by=c("target")]

    dt_wildtype <- unique(dt_wildtype[,.(target,huACE2,huACE2.pool6,cvACE2,RaACE2.787,RaACE2.787.pool6,RsACE2.3364,RpACE2,expression,huACE2_n,huACE2.pool6_n,cvACE2_n,RaACE2.787_n,RaACE2.787.pool6_n,RsACE2.3364_n,RpACE2_n,expression_n)])

    dt_mutant[,huACE2:=mean(mean_huACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,huACE2_n:=sum(n_huACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,huACE2.pool6:=mean(mean_huACE2.pool6_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,huACE2.pool6_n:=sum(n_huACE2.pool6_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,cvACE2:=mean(mean_cvACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,cvACE2_n:=sum(n_cvACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,RaACE2.787:=mean(mean_RaACE2.787_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,RaACE2.787_n:=sum(n_RaACE2.787_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,RaACE2.787.pool6:=mean(mean_RaACE2.787.pool6_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,RaACE2.787.pool6_n:=sum(n_RaACE2.787.pool6_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,RsACE2.3364:=mean(mean_RsACE2.3364_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,RsACE2.3364_n:=sum(n_RsACE2.3364_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,RpACE2:=mean(mean_RpACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,RpACE2_n:=sum(n_RpACE2_log10Ka,na.rm=T),by=c("target","position","mutant")]

    dt_mutant[,expression:=mean(mean_expression,na.rm=T),by=c("target","position","mutant")]
    dt_mutant[,expression_n:=sum(n_expression,na.rm=T),by=c("target","position","mutant")]

    dt_mutant <- unique(dt_mutant[,.(target,wildtype,position,mutant,huACE2,huACE2.pool6,cvACE2,RaACE2.787,RaACE2.787.pool6,RsACE2.3364,RpACE2,expression,huACE2_n,huACE2.pool6_n,cvACE2_n,RaACE2.787_n,RaACE2.787.pool6_n,RsACE2.3364_n,RpACE2_n,expression_n)])
    dt_mutant[as.character(wildtype)==as.character(mutant),c("huACE2_n","huACE2.pool6_n","cvACE2_n","RaACE2.787_n","RaACE2.787.pool6_n","RsACE2.3364_n","RpACE2_n","expression_n"):=NA]

Filter out the two backgrounds that were completely non-expressing. Most
barcodes were purged before the affinity measurements for these
backgrounds, so the affinities are determined from few barcodes and are
just generally unreliable because these are poorly folded/expressing
variants. (E.g. could see very high standard deviations)

    dt_wildtype[expression<7,c("huACE2","huACE2.pool6","cvACE2","RaACE2.787","RaACE2.787.pool6","RsACE2.3364","RpACE2","expression","huACE2_n","huACE2.pool6_n","cvACE2_n","RaACE2.787_n","RaACE2.787.pool6_n","RsACE2.3364_n","RpACE2_n","expression_n"):=NA]

Coverage stats on n\_barcodes for different measurements in the final
pooled measurements.

    par(mfrow=c(2,8))
    hist(dt_wildtype$huACE2_n,col="gray50",main=paste("WTs huACE2 affinity,\nmedian ",median(dt_wildtype$huACE2_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype[target %in% c("GD-Pangolin","RaTG13","BtKY72"),huACE2.pool6_n],col="gray50",main=paste("WTs huACE2.pool6 affinity,\nmedian ",median(dt_wildtype[target %in% c("GD-Pangolin","RaTG13","BtKY72"),huACE2.pool6_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype$cvACE2_n,col="gray50",main=paste("WTs cvACE2 affinity,\nmedian ",median(dt_wildtype$cvACE2_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype$RaACE2.787_n,col="gray50",main=paste("WTs RaACE2.787 affinity,\nmedian ",median(dt_wildtype$RaACE2.787_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype[target %in% c("GD-Pangolin","RaTG13","BtKY72"),RaACE2.787.pool6_n],col="gray50",main=paste("WTs RaACE2.787.pool6 affinity,\nmedian ",median(dt_wildtype[target %in% c("GD-Pangolin","RaTG13","BtKY72"),RaACE2.787.pool6_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype$RsACE2.3364_n,col="gray50",main=paste("WTs RsACE2.3364 affinity,\nmedian ",median(dt_wildtype$RsACE2.3364_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype$RpACE2_n,col="gray50",main=paste("WTs RpACE2 affinity,\nmedian ",median(dt_wildtype$RpACE2_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_wildtype$expression_n,col="gray50",main=paste("WTs expression,\nmedian ",median(dt_wildtype$expression_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

    hist(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),huACE2_n],col="gray50",main=paste("muts huACE2 affinity,\nmedian ",median(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),huACE2_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant[(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) | (target=="BtKY72"),huACE2.pool6_n],col="gray50",main=paste("muts huACE2.pool6 affinity,\nmedian ",median(dt_mutant[(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) | (target=="BtKY72"),huACE2.pool6_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant$cvACE2_n,col="gray50",main=paste("muts cvACE2 affinity,\nmedian ",median(dt_mutant$cvACE2_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),RaACE2.787_n],col="gray50",main=paste("muts RaACE2.787 affinity,\nmedian ",median(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),RaACE2.787_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant[(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) | (target=="BtKY72"),RaACE2.787.pool6_n],col="gray50",main=paste("muts RaACE2.787.pool6 affinity,\nmedian ",median(dt_mutant[(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) | (target=="BtKY72"),RaACE2.787.pool6_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant$RsACE2.3364_n,col="gray50",main=paste("muts RsACE2.3364 affinity,\nmedian ",median(dt_mutant$RsACE2.3364_n,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),RpACE2_n],col="gray50",main=paste("muts RpACE2 affinity,\nmedian ",median(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),RpACE2_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
    hist(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),expression_n],col="gray50",main=paste("muts expression,\nmedian ",median(dt_mutant[!(target %in% c("GD-Pangolin","RaTG13","BtKY72") & position==455) & !(target=="BtKY72"),expression_n],na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/n_barcode_plots-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/hist_n_barcodes.pdf",sep="")))

    #remove n info from table
    dt_wildtype <- dt_wildtype[,.(target,huACE2,huACE2.pool6,cvACE2,RaACE2.787,RaACE2.787.pool6,RsACE2.3364,RpACE2,expression)]
    dt_mutant <- dt_mutant[,.(target,wildtype,position,mutant,huACE2,huACE2.pool6,cvACE2,RaACE2.787,RaACE2.787.pool6,RsACE2.3364,RpACE2,expression)]

Reconcile pool6 and “main” pool assays. There is not perfect concordance
in Kas for the three wildtypes that are determined in both the “main”
huACE2 and RaACE2.787 assays and in the “pool6 only” assays. I show this
in the table below:

    kable(dt_wildtype[target%in% c("BtKY72","RaTG13","GD-Pangolin"),.(target,huACE2,huACE2.pool6,RaACE2.787,RaACE2.787.pool6)])

| target      |    huACE2 | huACE2.pool6 | RaACE2.787 | RaACE2.787.pool6 |
|:------------|----------:|-------------:|-----------:|-----------------:|
| BtKY72      |  5.162047 |     5.000000 |   7.036625 |         8.081915 |
| RaTG13      |  8.441850 |     8.502368 |   5.052615 |         5.182202 |
| GD-Pangolin | 10.982683 |    10.648554 |   5.993478 |         6.677448 |

I presume this is becasue I’m not perfectly dealing with ligand
depletion, and these assays would be conducted in the context of
differing competition and perhaps ligand depeltion effects. The rank
order of course matches between the two, but some of the Kds differ
between assays (moreso Ra than hu, though there is still subtle
difference in hu – Ra had the most prevalent and tightst binding of
ligands, so actually that would make sense if ligand depletion were most
exacerbated for this ligand). For the wildtype assays, I don’t need the
WT values. However, for the mutants, I should be comparing each mutant
only to the wildtype from its own internal library experiment. So,
below, I reconcile the pool6 and main assay values, and reference the
approporiate ‘wt’ value for the pool6 SSM positions.

    #first, fill in the "wildtype" state in the mutants table
    for(i in 1:nrow(dt_mutant)){
      if(as.character(dt_mutant[i,wildtype])==as.character(dt_mutant[i,mutant])){
        dt_mutant[i,c("huACE2","huACE2.pool6","cvACE2","RaACE2.787","RaACE2.787.pool6","RsACE2.3364","RpACE2","expression"):=dt_wildtype[target==dt_mutant[i,target],.(huACE2,huACE2.pool6,cvACE2,RaACE2.787,RaACE2.787.pool6,RsACE2.3364,RpACE2,expression)]]
      }
    }
    #we can remove pool6 columns now from the wildtype dt
    dt_wildtype[,c("huACE2.pool6","RaACE2.787.pool6"):=NULL]

    #reconcile pool6 with main mut measurements: for the measurements done in pool6, replace the "main" column, and remove the pool6 columns
    dt_mutant[(target=="BtKY72") | (target %in% c("RaTG13","GD-Pangolin") & position==455),huACE2:=huACE2.pool6]
    dt_mutant[(target=="BtKY72") | (target %in% c("RaTG13","GD-Pangolin") & position==455),RaACE2.787:=RaACE2.787.pool6]
    dt_mutant[,c("huACE2.pool6","RaACE2.787.pool6"):=NULL]

    #add delta_log10Kas
    for(i in 1:nrow(dt_mutant)){
      dt_mutant[i,huACE2_delta := huACE2 - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),huACE2][1]]
      dt_mutant[i,cvACE2_delta := cvACE2 - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),cvACE2][1]]
      dt_mutant[i,RaACE2.787_delta := RaACE2.787 - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),RaACE2.787][1]]
      dt_mutant[i,RsACE2.3364_delta := RsACE2.3364 - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),RsACE2.3364][1]]
      dt_mutant[i,RpACE2_delta := RpACE2 - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),RpACE2][1]]
      dt_mutant[i,expression_delta := expression - dt_mutant[target==dt_mutant[i,target] & as.character(wildtype)==as.character(mutant),expression][1]]
    }

Order factor variables for plotting

    #order target by order given in config
    dt_mutant$target <- factor(dt_mutant$target,levels=config$mutated_targets_ordered)
    #order mutant as a factor for grouping by rough biochemical grouping
    dt_mutant$mutant <- factor(dt_mutant$mutant, levels=c("C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
    #order sites as a factor variable
    dt_mutant$position <- factor(dt_mutant$position,levels=c(455,486,493,494,498,501))
    #add character vector indicating wildtype to use as plotting symbols for wt
    dt_mutant[,wildtype_indicator := ""]
    dt_mutant[as.character(mutant)==as.character(wildtype),wildtype_indicator := "x"]

Heatmaps! Output heatmaps illustrating all wildtype variants with separate columns for each ACE2.
-------------------------------------------------------------------------------------------------

    #make temp long-form data frame
    temp1 <- data.table::melt(dt_wildtype[,.(target,huACE2,cvACE2,RaACE2.787,RsACE2.3364,RpACE2)],id.vars=c("target"),measure.vars=c("huACE2","cvACE2","RaACE2.787","RsACE2.3364","RpACE2"),variable.name="ACE2",value.name="log10Ka")

    temp2 <- dt_wildtype[,.(target,expression)]

    p1 <- ggplot(temp1,aes(target,ACE2))+geom_tile(aes(fill=log10Ka),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+
      #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

    p2 <- ggplot(temp2,aes(target,y=1))+geom_tile(aes(fill=expression),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#06C528"),limits=c(5,11),values=c(0,1),na.value="yellow")+
      labs(x="",y="expression")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y = element_blank())

    ggarrange(p2,p1,nrow=2)

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_wildtypes_all-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_all_wildtypes.pdf",sep="")))

    #make temp long-form data frame
    extant <- c(config$EurAf_extant,config$SARS2_extant,config$SARS1_extant,config$Clade2_extant)

    temp1 <- temp1[target %in% extant,];temp1$target <- factor(temp1$target,levels=extant)

    p1 <- ggplot(temp1,aes(target,ACE2))+geom_tile(aes(fill=log10Ka),color="black",lwd=0.2)+
      scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+
      #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
      labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_wildtypes_phylogeny-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_extant_wildtypes.pdf",sep="")))

    #make temp long-form data frame
    ancestors <- c(config$ancestors_MAP)

    temp1 <- data.table::melt(dt_wildtype[target %in% ancestors, .(target,huACE2,cvACE2,RaACE2.787,RsACE2.3364,RpACE2)],id.vars=c("target"),measure.vars=c("huACE2","cvACE2","RaACE2.787","RsACE2.3364","RpACE2"),variable.name="ACE2",value.name="log10Ka")
    temp1$target <- factor(temp1$target,levels=ancestors)

    p1 <- ggplot(temp1,aes(target,ACE2))+geom_tile(aes(fill=log10Ka),color="black",lwd=0.2)+
      scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+
      #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
      labs(x="RBD homolog",y="measurement")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_wildtypes_MAP-ancestors-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_MAP-ancestors_wildtypes.pdf",sep="")))

    #make temp long-form data frame
    ancestors <- c(config$ancestors_MAP_v_alt)

    temp1 <- data.table::melt(dt_wildtype[target %in% ancestors, .(target,huACE2,cvACE2,RaACE2.787,RsACE2.3364,RpACE2)],id.vars=c("target"),measure.vars=c("huACE2","cvACE2","RaACE2.787","RsACE2.3364","RpACE2"),variable.name="ACE2",value.name="log10Ka")
    temp1$target <- factor(temp1$target,levels=ancestors)

    p1 <- ggplot(temp1,aes(target,ACE2))+geom_tile(aes(fill=log10Ka),color="black",lwd=0.2)+
      scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+
      #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
      labs(x="RBD homolog",y="measurement")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_wildtypes_MAP-and-alt-ancestors-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_MAP-and-alt-ancestors_wildtypes.pdf",sep="")))

Next, for the SSM libraries. Two ways of looking at these:

First, separately for each of the 14 backgrounds: 6 side-by-side plots
(for each position), each being 4\*20 (20 amino acids, 4 ACE2 targets
\[can probably skip RpACE2?\]). Metric is raw log10Ka.

Second, for each of the 6 positions: 14 side-by-side plots (for each
target), each being 4\*20 (20 amino acids, 4 ACE2 targets \[can probalby
skip RpACE2?\]). Metric is delta-log10Ka.

First, illustrating log10Ka grouped by RBD background.

    #make temp long-form data frame
    temp <- data.table::melt(dt_mutant[, .(target,position,mutant,huACE2,cvACE2,RaACE2.787,RsACE2.3364,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("huACE2","cvACE2","RaACE2.787","RsACE2.3364"),variable.name="ACE2",value.name="log10Ka")

    #for method to duplicate aa labels on right side of plot https://github.com/tidyverse/ggplot2/issues/3171
    guide_axis_label_trans <- function(label_trans = identity, ...) {
      axis_guide <- guide_axis(...)
      axis_guide$label_trans <- rlang::as_function(label_trans)
      class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
      axis_guide
    }

    guide_train.guide_axis_trans <- function(x, ...) {
      trained <- NextMethod()
      trained$key$.label <- x$label_trans(trained$key$.label)
      trained
    }

    p1 <- ggplot(temp,aes(ACE2,mutant))+geom_tile(aes(fill=log10Ka),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
      #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(target~position,nrow=7)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_log10Ka-by-target-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_log10Ka-by-target.pdf",sep="")))

Second, illustrating delta\_log10Ka grouped by SSM position.

    #make temp long-form data frame
    temp <- data.table::melt(dt_mutant[, .(target,position,mutant,huACE2_delta,cvACE2_delta,RaACE2.787_delta,RsACE2.3364_delta,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("huACE2_delta","cvACE2_delta","RaACE2.787_delta","RsACE2.3364_delta"),variable.name="ACE2",value.name="delta_log10Ka")

    p1 <- ggplot(temp,aes(ACE2,mutant))+geom_tile(aes(fill=delta_log10Ka),color="black",lwd=0.1)+
      #scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3/14,5/14,7/14,9/14,11/14,14/14),na.value="yellow")+
      scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3.5/14,7/14,10.5/14,14/14),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(position~target,nrow=6)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_delta-log10Ka-by-position-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_deltalog10Ka-by-position.pdf",sep="")))

Last data-dump heatmap to show – normalize each heatmap so its
illustrating delta-Ka relative to the *same* reference state in each
sub-panel. Easiest is probably either the ancestral state in
AncSarbecovirus, or the SARS-CoV-2 state… First, refernece state
AncSarbecovirus state

    #calculate delta_log10Ka_ancSarb, normalizes delta to be against AncSarbecovirus amino acid reference state
    for(i in 1:nrow(dt_mutant)){
      dt_mutant[i,huACE2_delta_ancSarb := huACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="AncSarbecovirus_MAP" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),huACE2]]
      dt_mutant[i,cvACE2_delta_ancSarb := cvACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="AncSarbecovirus_MAP" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),cvACE2]]
      dt_mutant[i,RaACE2.787_delta_ancSarb := RaACE2.787 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="AncSarbecovirus_MAP" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RaACE2.787]]
      dt_mutant[i,RsACE2.3364_delta_ancSarb := RsACE2.3364 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="AncSarbecovirus_MAP" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RsACE2.3364]]
      dt_mutant[i,RpACE2_delta_ancSarb := RpACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="AncSarbecovirus_MAP" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RpACE2]]
    }

    #make temp long-form data frame
    temp <- data.table::melt(dt_mutant[, .(target,position,mutant,huACE2_delta_ancSarb,cvACE2_delta_ancSarb,RaACE2.787_delta_ancSarb,RsACE2.3364_delta_ancSarb,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("huACE2_delta_ancSarb","cvACE2_delta_ancSarb","RaACE2.787_delta_ancSarb","RsACE2.3364_delta_ancSarb"),variable.name="ACE2",value.name="delta_log10Ka_ancSarb")

    p1 <- ggplot(temp,aes(ACE2,mutant))+geom_tile(aes(fill=delta_log10Ka_ancSarb),color="black",lwd=0.1)+
      #scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3/14,5/14,7/14,9/14,11/14,14/14),na.value="yellow")+
      scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3.5/14,7/14,10.5/14,14/14),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(position~target,nrow=6)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_delta-log10Ka-by-position_norm_ancestral-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_deltalog10Ka-by-position_norm_ancestral.pdf",sep="")))

Second, SARS2 state

    #calculate delta_log10Ka_SARS2, normalizes delta to be against SARS-CoV-2 amino acid reference state
    for(i in 1:nrow(dt_mutant)){
      dt_mutant[i,huACE2_delta_SARS2 := huACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="SARS-CoV-2" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),huACE2]]
      dt_mutant[i,cvACE2_delta_SARS2 := cvACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="SARS-CoV-2" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),cvACE2]]
      dt_mutant[i,RaACE2.787_delta_SARS2 := RaACE2.787 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="SARS-CoV-2" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RaACE2.787]]
      dt_mutant[i,RsACE2.3364_delta_SARS2 := RsACE2.3364 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="SARS-CoV-2" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RsACE2.3364]]
      dt_mutant[i,RpACE2_delta_SARS2 := RpACE2 - dt_mutant[target==dt_mutant[i,target] & position==dt_mutant[i,position] & as.character(mutant)==as.character(dt_mutant[target=="SARS-CoV-2" & position==dt_mutant[i,position] & as.character(wildtype)==as.character(mutant),wildtype]),RpACE2]]
    }

    #make temp long-form data frame
    temp <- data.table::melt(dt_mutant[, .(target,position,mutant,huACE2_delta_SARS2,cvACE2_delta_SARS2,RaACE2.787_delta_SARS2,RsACE2.3364_delta_SARS2,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("huACE2_delta_SARS2","cvACE2_delta_SARS2","RaACE2.787_delta_SARS2","RsACE2.3364_delta_SARS2"),variable.name="ACE2",value.name="delta_log10Ka_SARS2")

    p1 <- ggplot(temp,aes(ACE2,mutant))+geom_tile(aes(fill=delta_log10Ka_SARS2),color="black",lwd=0.1)+
      #scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3/14,5/14,7/14,9/14,11/14,14/14),na.value="yellow")+
      scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-7,7),values=c(0,3.5/14,7/14,10.5/14,14/14),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(position~target,nrow=6)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_delta-log10Ka-by-position_norm_SARS2-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_deltalog10Ka-by-position_norm_SARS2.pdf",sep="")))

Also plot expression effects of mutations in heatmaps

    temp <- dt_mutant[,.(target,position,mutant,expression,wildtype_indicator)]

    p1 <- ggplot(temp,aes(position,mutant))+geom_tile(aes(fill=expression),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#06C528"),limits=c(5,11),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_expression-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_expression-by-target.pdf",sep="")))

    temp <- dt_mutant[,.(target,position,mutant,expression_delta,wildtype_indicator)]

    p1 <- ggplot(temp,aes(position,mutant))+geom_tile(aes(fill=expression_delta),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-2,1),values=c(0,1/3,2/3,2.5/3,3/3),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_delta_expression-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_delta-expression-by-target.pdf",sep="")))

This is actually perhaps a useful zoomed in summary of the overwhelming
heatmaps above. Let’s do ligand-specific heatmaps in this arrangement.
Which metric to show? I think this is about raw affinities, so I’ll do
the white to purple heatmap showing affinities.

    p1 <- ggplot(dt_mutant,aes(position,mutant))+geom_tile(aes(fill=huACE2),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_huACE2-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_huACE2-by-target.pdf",sep="")))

    p1 <- ggplot(dt_mutant,aes(position,mutant))+geom_tile(aes(fill=cvACE2),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_cvACE2-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_cvACE2-by-target.pdf",sep="")))

    p1 <- ggplot(dt_mutant,aes(position,mutant))+geom_tile(aes(fill=RaACE2.787),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_RaACE2.787-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_RaACE2.787-by-target.pdf",sep="")))

    p1 <- ggplot(dt_mutant,aes(position,mutant))+geom_tile(aes(fill=RsACE2.3364),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_RsACE2.3364-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_RsACE2.3364-by-target.pdf",sep="")))

    p1 <- ggplot(dt_mutant,aes(position,mutant))+geom_tile(aes(fill=RpACE2),color="black",lwd=0.1)+
      scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(6,12),values=c(0,1),na.value="yellow")+
      labs(x="",y="")+theme_classic(base_size=9)+
      coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
      facet_wrap(~target,nrow=1)+
      guides(y.sec=guide_axis_label_trans())+
      geom_text(aes(label=wildtype_indicator),size=2,color="gray10")

    p1

<img src="barcode_to_genotype_phenotypes_files/figure-gfm/heatmap_SSM_RpACE2-1.png" style="display: block; margin: auto;" />

    invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_RpACE2-by-target.pdf",sep="")))

That’s my data! Now, to think…

Save output files.

    dt_wildtype %>%
      mutate_if(is.numeric, round, digits=4) %>%
      write.csv(file=config$final_variant_scores_wt_file, row.names=F,quote=F)

    dt_mutant %>%
      mutate_if(is.numeric, round, digits=4) %>%
      write.csv(file=config$final_variant_scores_mut_file, row.names=F,quote=F)
