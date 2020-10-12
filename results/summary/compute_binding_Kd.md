Compute per-barcode binding constants
================
Tyler Starr
10/8/2020

-   [Setup](#setup)
-   [Calculating mean bin for each barcode at each sample
    concentration](#calculating-mean-bin-for-each-barcode-at-each-sample-concentration)
-   [Fit titration curves](#fit-titration-curves)
-   [QC and sanity checks](#qc-and-sanity-checks)
-   [Data filtering by fit quality](#data-filtering-by-fit-quality)
-   [Data Output](#data-output)

This notebook computes and summarizes per-variant binding constants for
ACE2 variants.

    require("knitr")
    knitr::opts_chunk$set(echo = T)
    knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

    #list of packages to install/load
    packages = c("yaml","data.table","tidyverse","gridExtra")
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
    if(!file.exists(config$Titeseq_Kds_dir)){
      dir.create(file.path(config$Titeseq_Kds_dir))
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
    ##  [1] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [5] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0       tibble_3.0.2     
    ##  [9] ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0       
    ## [13] knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.5     compiler_3.6.2  
    ##  [5] dbplyr_1.4.2     tools_3.6.2      digest_0.6.23    lubridate_1.7.4 
    ##  [9] jsonlite_1.6     evaluate_0.14    lifecycle_0.2.0  gtable_0.3.0    
    ## [13] pkgconfig_2.0.3  rlang_0.4.7      reprex_0.3.0     cli_2.0.0       
    ## [17] rstudioapi_0.10  DBI_1.1.0        haven_2.2.0      xfun_0.11       
    ## [21] withr_2.1.2      xml2_1.2.2       httr_1.4.1       fs_1.3.1        
    ## [25] hms_0.5.2        generics_0.0.2   vctrs_0.3.1      grid_3.6.2      
    ## [29] tidyselect_1.1.0 glue_1.3.1       R6_2.4.1         fansi_0.4.0     
    ## [33] readxl_1.3.1     rmarkdown_2.0    modelr_0.1.5     magrittr_1.5    
    ## [37] backports_1.1.5  scales_1.1.0     ellipsis_0.3.0   htmltools_0.4.0 
    ## [41] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3   
    ## [45] munsell_0.5.0    broom_0.7.0      crayon_1.3.4

Setup
-----

Read in table of variant genotypes and barcode counts. Remove samples
corresponding to expression Sort-seq experiments, analyzed in the
accompanying notebook.

    dt <- data.table(read.csv(file=config$merged_sequencing_file,stringsAsFactors = F))

    #eliminate columns of Sortseq counts
    dt[,c("SortSeq_bin1","SortSeq_bin2","SortSeq_bin3","SortSeq_bin4"):=NULL]

    #read dataframe with list of barcode runs
    barcode_runs <- read.csv(file=config$barcode_runs,stringsAsFactors=F); barcode_runs <- subset(barcode_runs, select=-c(R1))

    #eliminate rows from barcode_runs that are not from a binding Tite-seq experiment
    barcode_runs <- barcode_runs[barcode_runs$sample_type != "SortSeq",]

    #make tables giving names of Titeseq samples and the corresponding ACE2 incubation concentrations
    samples_huACE2 <- data.frame(sample=unique(paste(barcode_runs[barcode_runs$sample_type=="huACE2","sample_type"],formatC(barcode_runs[barcode_runs$sample_type=="huACE2","concentration"], width=2,flag="0"),sep="_")),conc=c(10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12, 10^-13,0))

    samples_RaACE2 <- data.frame(sample=unique(paste(barcode_runs[barcode_runs$sample_type=="RaACE2","sample_type"],formatC(barcode_runs[barcode_runs$sample_type=="RaACE2","concentration"], width=2,flag="0"),sep="_")),conc=c(10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12, 10^-13,0))

    samples_RpACE2 <- data.frame(sample=unique(paste(barcode_runs[barcode_runs$sample_type=="RpACE2","sample_type"],formatC(barcode_runs[barcode_runs$sample_type=="RpACE2","concentration"], width=2,flag="0"),sep="_")),conc=c(10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12,0))

    samples_RsACE2 <- data.frame(sample=unique(paste(barcode_runs[barcode_runs$sample_type=="RsACE2","sample_type"],formatC(barcode_runs[barcode_runs$sample_type=="RsACE2","concentration"], width=2,flag="0"),sep="_")),conc=c(10^-6, 10^-7, 10^-8, 10^-9, 10^-10, 10^-11, 10^-12, 10^-13,0))

Calculating mean bin for each barcode at each sample concentration
------------------------------------------------------------------

Next, for each barcode at each of the ACE2 concentrations, calculate the
“mean bin” response variable. This is calculated as a simple mean, where
the value of each bin is the integer value of the bin (bin1=unbound,
bin4=highly bound) – because of how bins are defined, the mean
fluorescence of cells in each bin are equally spaced on a log-normal
scale, so mean bin correlates with simple mean fluorescence. A filtering
step is included that discards mean bin estimates for
sample/concentrations where splits of cell counts across non-consecutive
bins (1/3, 2/4, 1/4) exceed certain filter cutoffs to eliminate
egregious bimodalilty. A more robust maximum likelihood approach could
be implemented to calculate these mean bins, which requires knowledge of
the fluorescence boundaries of the sort bins. We therefore provide them
for posterity’s sake below. For the huACE2 titration sorts, the
fluorescence boundaries for bins 1-4 are as follows:

    (-288, 184), (185, 2251), (2252, 27722), (27723, 262143)

For the library 1 RaACE2, RpACE2, and RsACE2 titration sorts, the
fluorescence boundaries for bins 1-4 are as follows:

    (-288, 65), (66, 677), (678, 7054), (7055, 262143)

    #function that returns mean bin and sum of counts for four bins cell counts. Includes cutoffs for bimodal sample splits to filter out
    calc.meanbin <- function(vec, split13filter=0.4, split24filter=0.4, split14filter=0.2){
      total <- sum(vec)
      if(is.na(total) | (vec[1] > split13filter*total & vec[3] > split13filter*total) | (vec[2] > split24filter*total & vec[4] > split24filter*total) | (vec[1] > split14filter*total & vec[4] > split14filter*total)){
        return(list(NA,NA))
      }else{
        return( list((vec[1]*1+vec[2]*2+vec[3]*3+vec[4]*4)/(vec[1]+vec[2]+vec[3]+vec[4]), total) )
      }
    }
      

    #iterate through Titeseq samples, compute mean_bin and total_count for each barcode variant
    for(i in 1:nrow(samples_huACE2)){ #iterate through titeseq sample (concentration)
      meanbin_out <- paste(samples_huACE2[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
      totalcount_out <- paste(samples_huACE2[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
      bin1_in <- paste(samples_huACE2[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
      bin2_in <- paste(samples_huACE2[i,"sample"],"_bin2",sep="")
      bin3_in <- paste(samples_huACE2[i,"sample"],"_bin3",sep="")
      bin4_in <- paste(samples_huACE2[i,"sample"],"_bin4",sep="")
      dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
    }

    for(i in 1:nrow(samples_RaACE2)){ #iterate through titeseq sample (concentration)
      meanbin_out <- paste(samples_RaACE2[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
      totalcount_out <- paste(samples_RaACE2[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
      bin1_in <- paste(samples_RaACE2[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
      bin2_in <- paste(samples_RaACE2[i,"sample"],"_bin2",sep="")
      bin3_in <- paste(samples_RaACE2[i,"sample"],"_bin3",sep="")
      bin4_in <- paste(samples_RaACE2[i,"sample"],"_bin4",sep="")
      dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
    }

    for(i in 1:nrow(samples_RpACE2)){ #iterate through titeseq sample (concentration)
      meanbin_out <- paste(samples_RpACE2[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
      totalcount_out <- paste(samples_RpACE2[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
      bin1_in <- paste(samples_RpACE2[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
      bin2_in <- paste(samples_RpACE2[i,"sample"],"_bin2",sep="")
      bin3_in <- paste(samples_RpACE2[i,"sample"],"_bin3",sep="")
      bin4_in <- paste(samples_RpACE2[i,"sample"],"_bin4",sep="")
      dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
    }

    for(i in 1:nrow(samples_RsACE2)){ #iterate through titeseq sample (concentration)
      meanbin_out <- paste(samples_RsACE2[i,"sample"],"_meanbin",sep="") #define the header name for the meanbin output for the given concentration sample
      totalcount_out <- paste(samples_RsACE2[i,"sample"],"_totalcount",sep="") #define the header name for the total cell count output for the given concentration sample
      bin1_in <- paste(samples_RsACE2[i,"sample"],"_bin1",sep="") #define the header names for the input cell counts for bins1-4 of the given concnetration sample
      bin2_in <- paste(samples_RsACE2[i,"sample"],"_bin2",sep="")
      bin3_in <- paste(samples_RsACE2[i,"sample"],"_bin3",sep="")
      bin4_in <- paste(samples_RsACE2[i,"sample"],"_bin4",sep="")
      dt[,c(meanbin_out,totalcount_out) := calc.meanbin(c(get(bin1_in),get(bin2_in),get(bin3_in),get(bin4_in))),by=c("barcode","library")]
    }

Fit titration curves
--------------------

We will use nonlinear least squares regression to fit curves to each
barcode’s titration series. We will also include a minimum cell count
that is required for a meanbin estimate to be used in the titration fit,
and a minimum number of concentrations with determined meanbin that is
required for a titration to be reported.

    #For QC and filtering, output columns giving the average number of cells that were sampled for a barcode across the 9 sample concentrations, and a value for the number of meanbin estimates that were removed for being below the # of cells cutoff
    cutoff <- 10
    dt[,huACE2_avgcount := mean(c(huACE2_01_totalcount,huACE2_02_totalcount,huACE2_03_totalcount,huACE2_04_totalcount,
                                    huACE2_05_totalcount,huACE2_06_totalcount,huACE2_07_totalcount,huACE2_08_totalcount,
                                    huACE2_09_totalcount),na.rm=T),by=c("library","barcode")]
    dt[,RaACE2_avgcount := mean(c(RaACE2_01_totalcount,RaACE2_02_totalcount,RaACE2_03_totalcount,RaACE2_04_totalcount,
                                    RaACE2_05_totalcount,RaACE2_06_totalcount,RaACE2_07_totalcount,RaACE2_08_totalcount,
                                    RaACE2_09_totalcount),na.rm=T),by=c("library","barcode")]
    dt[,RpACE2_avgcount := mean(c(RpACE2_02_totalcount,RpACE2_03_totalcount,RpACE2_04_totalcount,
                                    RpACE2_05_totalcount,RpACE2_06_totalcount,RpACE2_07_totalcount,
                                    RpACE2_09_totalcount),na.rm=T),by=c("library","barcode")]
    dt[,RsACE2_avgcount := mean(c(RsACE2_01_totalcount,RsACE2_02_totalcount,RsACE2_03_totalcount,RsACE2_04_totalcount,
                                    RsACE2_05_totalcount,RsACE2_06_totalcount,RsACE2_07_totalcount,RsACE2_08_totalcount,
                                    RsACE2_09_totalcount),na.rm=T),by=c("library","barcode")]

    #number of concentrations at which meanbin is calculated from < cutoff cells or is missing b/c filtered for bimodality
    dt[,huACE2_min_cell_filtered := sum(c(c(huACE2_01_totalcount,huACE2_02_totalcount,huACE2_03_totalcount,huACE2_04_totalcount,
                                            huACE2_05_totalcount,huACE2_06_totalcount,huACE2_07_totalcount,huACE2_08_totalcount,
                                            huACE2_09_totalcount)<cutoff,is.na(c(huACE2_01_totalcount,huACE2_02_totalcount,huACE2_03_totalcount,huACE2_04_totalcount,
                                                                                 huACE2_05_totalcount,huACE2_06_totalcount,huACE2_07_totalcount,huACE2_08_totalcount,
                                                                                 huACE2_09_totalcount))),na.rm=T),by=c("library","barcode")]
    dt[,RaACE2_min_cell_filtered := sum(c(c(RaACE2_01_totalcount,RaACE2_02_totalcount,RaACE2_03_totalcount,RaACE2_04_totalcount,
                                            RaACE2_05_totalcount,RaACE2_06_totalcount,RaACE2_07_totalcount,RaACE2_08_totalcount,
                                            RaACE2_09_totalcount)<cutoff,is.na(c(RaACE2_01_totalcount,RaACE2_02_totalcount,RaACE2_03_totalcount,RaACE2_04_totalcount,
                                                                                 RaACE2_05_totalcount,RaACE2_06_totalcount,RaACE2_07_totalcount,RaACE2_08_totalcount,
                                                                                 RaACE2_09_totalcount))),na.rm=T),by=c("library","barcode")]
    dt[,RpACE2_min_cell_filtered := sum(c(c(RpACE2_02_totalcount,RpACE2_03_totalcount,RpACE2_04_totalcount,
                                            RpACE2_05_totalcount,RpACE2_06_totalcount,RpACE2_07_totalcount,
                                            RpACE2_09_totalcount)<cutoff,is.na(c(RpACE2_02_totalcount,RpACE2_03_totalcount,RpACE2_04_totalcount,
                                                                                 RpACE2_05_totalcount,RpACE2_06_totalcount,RpACE2_07_totalcount,
                                                                                 RpACE2_09_totalcount))),na.rm=T),by=c("library","barcode")]
    dt[,RsACE2_min_cell_filtered := sum(c(c(RsACE2_01_totalcount,RsACE2_02_totalcount,RsACE2_03_totalcount,RsACE2_04_totalcount,
                                            RsACE2_05_totalcount,RsACE2_06_totalcount,RsACE2_07_totalcount,RsACE2_08_totalcount,
                                            RsACE2_09_totalcount)<cutoff,is.na(c(RsACE2_01_totalcount,RsACE2_02_totalcount,RsACE2_03_totalcount,RsACE2_04_totalcount,
                                                                                 RsACE2_05_totalcount,RsACE2_06_totalcount,RsACE2_07_totalcount,RsACE2_08_totalcount,
                                                                                 RsACE2_09_totalcount))),na.rm=T),by=c("library","barcode")]


    #function that fits a nls regression to the titration series, including an option to filter below certain thresholds for average cells across all samples, and number of samples below a cutoff of cells
    fit.titration <- function(y.vals,x.vals,count.vals,min.cfu=cutoff,
                              min.means=0.8,min.average=10,Kd.start=2e-11,
                              a.start=3,a.lower=2,a.upper=3,
                              b.start=1,b.lower=1,b.upper=1.5){
      indices <- count.vals>min.cfu & !is.na(y.vals)
      y <- y.vals[indices]
      x <- x.vals[indices]
      if((length(y) < min.means*length(y.vals)) | (mean(count.vals,na.rm=T) < min.average)){ #return NAs if < min.means fraction of concentrations have above min.cfu counts or if the average count across all concentrations is below min.average
        return(list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA)))
      }else{
        fit <- nls(y ~ a*(x/(x+Kd))+b,
                   start=list(a=a.start,b=b.start,Kd=Kd.start),
                   lower=list(a=a.lower,b=b.lower,Kd=min(x.vals[x.vals>0])/100), #constrain Kd to be no lower than 1/100x the lowest concentration value
                   upper=list(a=a.upper,b=b.upper,Kd=max(x.vals[x.vals>0])*10), #constrain Kd to be no higher than the 10x highest concentration value
                   algorithm="port")  
        return(list(as.numeric(summary(fit)$coefficients["Kd","Estimate"]),
                    as.numeric(summary(fit)$coefficients["Kd","Std. Error"]),
                    as.numeric(summary(fit)$coefficients["a","Estimate"]),
                    as.numeric(summary(fit)$coefficients["b","Estimate"]),
                    as.numeric(summary(fit)$sigma),
                    list(fit)))
      }
    }

    #fit titration to huACE2 Titeseq data for each barcode
    dt[,c("Kd_huACE2","Kd_SE_huACE2","response_huACE2","baseline_huACE2","RSE_huACE2","fit_huACE2") :=
         tryCatch(fit.titration(y.vals=c(huACE2_01_meanbin,huACE2_02_meanbin,huACE2_03_meanbin,huACE2_04_meanbin,
                                         huACE2_05_meanbin,huACE2_06_meanbin,huACE2_07_meanbin,huACE2_08_meanbin,
                                         huACE2_09_meanbin),
                                x.vals=samples_huACE2$conc,
                                count.vals=c(huACE2_01_totalcount,huACE2_02_totalcount,huACE2_03_totalcount,huACE2_04_totalcount,
                                             huACE2_05_totalcount,huACE2_06_totalcount,huACE2_07_totalcount,huACE2_08_totalcount,huACE2_09_totalcount)),
                  error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))}),by=c("library","barcode")]

    #fit titration to RaACE2 Titeseq data for each barcode
    dt[,c("Kd_RaACE2","Kd_SE_RaACE2","response_RaACE2","baseline_RaACE2","RSE_RaACE2","fit_RaACE2") :=
         tryCatch(fit.titration(y.vals=c(RaACE2_01_meanbin,RaACE2_02_meanbin,RaACE2_03_meanbin,RaACE2_04_meanbin,
                                         RaACE2_05_meanbin,RaACE2_06_meanbin,RaACE2_07_meanbin,RaACE2_08_meanbin,
                                         RaACE2_09_meanbin),
                                x.vals=samples_RaACE2$conc,
                                count.vals=c(RaACE2_01_totalcount,RaACE2_02_totalcount,RaACE2_03_totalcount,RaACE2_04_totalcount,
                                             RaACE2_05_totalcount,RaACE2_06_totalcount,RaACE2_07_totalcount,RaACE2_08_totalcount,RaACE2_09_totalcount)),
                  error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))}),by=c("library","barcode")]

    #fit titration to RpACE2 Titeseq data for each barcode
    dt[,c("Kd_RpACE2","Kd_SE_RpACE2","response_RpACE2","baseline_RpACE2","RSE_RpACE2","fit_RpACE2") :=
         tryCatch(fit.titration(y.vals=c(RpACE2_02_meanbin,RpACE2_03_meanbin,RpACE2_04_meanbin,
                                         RpACE2_05_meanbin,RpACE2_06_meanbin,RpACE2_07_meanbin,
                                         RpACE2_09_meanbin),
                                x.vals=samples_RpACE2$conc,
                                count.vals=c(RpACE2_02_totalcount,RpACE2_03_totalcount,RpACE2_04_totalcount,
                                             RpACE2_05_totalcount,RpACE2_06_totalcount,RpACE2_07_totalcount,RpACE2_09_totalcount),
                                Kd.start=1e-6),
                  error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))}),by=c("library","barcode")]

    #fit titration to RsACE2 Titeseq data for each barcode
    dt[,c("Kd_RsACE2","Kd_SE_RsACE2","response_RsACE2","baseline_RsACE2","RSE_RsACE2","fit_RsACE2") :=
         tryCatch(fit.titration(y.vals=c(RsACE2_01_meanbin,RsACE2_02_meanbin,RsACE2_03_meanbin,RsACE2_04_meanbin,
                                         RsACE2_05_meanbin,RsACE2_06_meanbin,RsACE2_07_meanbin,RsACE2_08_meanbin,
                                         RsACE2_09_meanbin),
                                x.vals=samples_RsACE2$conc,
                                count.vals=c(RsACE2_01_totalcount,RsACE2_02_totalcount,RsACE2_03_totalcount,RsACE2_04_totalcount,
                                             RsACE2_05_totalcount,RsACE2_06_totalcount,RsACE2_07_totalcount,RsACE2_08_totalcount,RsACE2_09_totalcount)),
                  error=function(e){list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))}),by=c("library","barcode")]

    #the dt object is very large and interferes with session reloading if working in interactive mode
    #comment out for the pipeline run of this, but for interactive work, run this split and separate, and save of the different objects containing the fits
    # dt_models_hu <- dt[,.(library,barcode,target,variant_class,wildtype,position,mutant,fit_huACE2)]
    # dt_models_Ra <- dt[,.(library,barcode,target,variant_class,wildtype,position,mutant,fit_RaACE2)]
    # dt_models_Rp <- dt[,.(library,barcode,target,variant_class,wildtype,position,mutant,fit_RpACE2)]
    # dt_models_Rs <- dt[,.(library,barcode,target,variant_class,wildtype,position,mutant,fit_RsACE2)]
    # dt[,c("fit_huACE2","fit_RaACE2","fit_RpACE2","fit_RsACE2"):=NULL]

    # #save temp dt objects for reloading (will filter some stuff later on and might want to revert)
    # save(dt,file=paste(config$Titeseq_Kds_dir,"/dt.temp.Rda",sep=""))
    # save(dt_models_hu,file=paste(config$Titeseq_Kds_dir,"/dt_models_hu.temp.Rda",sep=""))
    # save(dt_models_Ra,file=paste(config$Titeseq_Kds_dir,"/dt_models_Ra.temp.Rda",sep=""))
    # save(dt_models_Rp,file=paste(config$Titeseq_Kds_dir,"/dt_models_Rp.temp.Rda",sep=""))
    # save(dt_models_Rs,file=paste(config$Titeseq_Kds_dir,"/dt_models_Rs.temp.Rda",sep=""))
    # #rm(dt_models_hu) #remove to aid in restart of session in interactive work
    # rm(dt_models_Ra) #remove to aid in restart of session
    # rm(dt_models_Rp) #remove to aid in restart of session
    # rm(dt_models_Rs) #remove to aid in restart of session

QC and sanity checks
--------------------

We will do some QC to make sure we got good titration curves for most of
our library barcodes. We will also spot check titration curves from
across our measurement range, and spot check curves whose fit parameters
hit the different boundary conditions of the fit variables.

We successfully generated *K*<sub>D,app</sub> estimates for 113310 of
our lib1 barcodes (78.68%) of our lib1+lib2 huACE2 titrations, 82.44%)
of the RaACE2 titrations, 85.21%) of the RpACE2 titrations, and 83.4%)
of the RsACE2 titrations.

To allow manual checks of what the data looks like for different curve
fits, I define functions that take a row from the dt table and the
corresponding table of fits, and plots the meanbin estimates and the fit
titration curve (if converged). This allows for quick and easy
troubleshooting and spot-checking of curves.

    #make functions that allow me to plot a titration for any given row from the counts data frames, for spot checking curves
    plot.titration.huACE2 <- function(row,output.text=F){
      y.vals <- c();for(sample in samples_huACE2$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
      x.vals <- samples_huACE2$conc
      count.vals <- c();for(sample in samples_huACE2$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
      if(dt[row,variant_class] =="mutant"){
        title <- paste(dt[row,target],paste(dt[row,wildtype],dt[row,position],dt[row,mutant]," huACE2",sep=""))
      }else{
        title <- paste(dt[row,target],dt[row,variant_class],"huACE2")
      }
      plot(x.vals[count.vals>cutoff],y.vals[count.vals>cutoff],xlab="[huACE2] (M)",
           ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=title)
      fit <- dt[row,fit_huACE2[[1]]]
      if(!is.na(fit)[1]){
        lines(10^c(seq(-13,-6,0.25)),predict(fit,newdata=list(x=10^c(seq(-13,-6,0.25)))))
        legend("topleft",bty="n",cex=1,legend=paste("Kd",format(dt[row,Kd_huACE2],digits=3),"M"))
      }
      if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
        dt[row,.(barcode,variant_class,wildtype,position,mutant,huACE2_avgcount,huACE2_min_cell_filtered,Kd_huACE2,Kd_SE_huACE2,baseline_huACE2,response_huACE2,RSE_huACE2)]
      }
    }

    plot.titration.RaACE2 <- function(row,output.text=F){
      y.vals <- c();for(sample in samples_RaACE2$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
      x.vals <- samples_RaACE2$conc
      count.vals <- c();for(sample in samples_RaACE2$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
      if(dt[row,variant_class] =="mutant"){
        title <- paste(dt[row,target],paste(dt[row,wildtype],dt[row,position],dt[row,mutant]," RaACE2",sep=""))
      }else{
        title <- paste(dt[row,target],dt[row,variant_class],"RaACE2")
      }
      plot(x.vals[count.vals>cutoff],y.vals[count.vals>cutoff],xlab="[RaACE2] (M)",
           ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=title)
      fit <- dt[row,fit_RaACE2[[1]]]
      if(!is.na(fit)[1]){
        lines(10^c(seq(-13,-6,0.25)),predict(fit,newdata=list(x=10^c(seq(-13,-6,0.25)))))
        legend("topleft",bty="n",cex=1,legend=paste("Kd",format(dt[row,Kd_RaACE2],digits=3),"M"))
      }
      if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
        dt[row,.(barcode,variant_class,wildtype,position,mutant,RaACE2_avgcount,RaACE2_min_cell_filtered,Kd_RaACE2,Kd_SE_RaACE2,baseline_RaACE2,response_RaACE2,RSE_RaACE2)]
      }
    }

    plot.titration.RpACE2 <- function(row,output.text=F){
      y.vals <- c();for(sample in samples_RpACE2$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
      x.vals <- samples_RpACE2$conc
      count.vals <- c();for(sample in samples_RpACE2$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
      if(dt[row,variant_class] =="mutant"){
        title <- paste(dt[row,target],paste(dt[row,wildtype],dt[row,position],dt[row,mutant]," RpACE2",sep=""))
      }else{
        title <- paste(dt[row,target],dt[row,variant_class], "RpACE2")
      }
      plot(x.vals[count.vals>cutoff],y.vals[count.vals>cutoff],xlab="[RpACE2] (M)",
           ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=title)
      fit <- dt[row,fit_RpACE2[[1]]]
      if(!is.na(fit)[1]){
        lines(10^c(seq(-13,-6,0.25)),predict(fit,newdata=list(x=10^c(seq(-13,-6,0.25)))))
        legend("topleft",bty="n",cex=1,legend=paste("Kd",format(dt[row,Kd_RpACE2],digits=3),"M"))
      }
      if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
        dt[row,.(barcode,variant_class,wildtype,position,mutant,RpACE2_avgcount,RpACE2_min_cell_filtered,Kd_RpACE2,Kd_SE_RpACE2,baseline_RpACE2,response_RpACE2,RSE_RpACE2)]
      }
    }

    plot.titration.RsACE2 <- function(row,output.text=F){
      y.vals <- c();for(sample in samples_RsACE2$sample){y.vals <- c(y.vals,paste(sample,"_meanbin",sep=""))};y.vals <- unlist(dt[row,y.vals,with=F])
      x.vals <- samples_RsACE2$conc
      count.vals <- c();for(sample in samples_RsACE2$sample){count.vals <- c(count.vals,paste(sample,"_totalcount",sep=""))};count.vals <- unlist(dt[row,count.vals,with=F])
      if(dt[row,variant_class] =="mutant"){
        title <- paste(dt[row,target],paste(dt[row,wildtype],dt[row,position],dt[row,mutant]," RsACE2",sep=""))
      }else{
        title <- paste(dt[row,target],dt[row,variant_class], "RsACE2")
      }
      plot(x.vals[count.vals>cutoff],y.vals[count.vals>cutoff],xlab="[RsACE2] (M)",
           ylab="mean bin",log="x",ylim=c(1,4),xlim=c(1e-13,1e-6),pch=19,main=title)
      fit <- dt[row,fit_RsACE2[[1]]]
      if(!is.na(fit)[1]){
        lines(10^c(seq(-13,-6,0.25)),predict(fit,newdata=list(x=10^c(seq(-13,-6,0.25)))))
        legend("topleft",bty="n",cex=1,legend=paste("Kd",format(dt[row,Kd_RsACE2],digits=3),"M"))
      }
      if(output.text==T){ #for troubleshooting and interactive work, output some info from the counts table for the given row
        dt[row,.(barcode,variant_class,wildtype,position,mutant,RsACE2_avgcount,RsACE2_min_cell_filtered,Kd_RsACE2,Kd_SE_RsACE2,baseline_RsACE2,response_RsACE2,RSE_RsACE2)]
      }
    }

Let’s look at our distribution of *K*<sub>D,app</sub> estimates.

    par(mfrow=c(4,1))
    hist(log10(dt$Kd_huACE2),col="gray40",breaks=60,xlab="log10(K_D,app), huACE2 (M)",main="",xlim=log10(range(c(dt$Kd_huACE2,dt$Kd_RaACE2,dt$Kd_RpACE2,dt$Kd_RsACE2),na.rm=T)))
    hist(log10(dt$Kd_RaACE2),col="gray40",breaks=60,xlab="log10(K_D,app), RaACE2 (M)",main="",xlim=log10(range(c(dt$Kd_huACE2,dt$Kd_RaACE2,dt$Kd_RpACE2,dt$Kd_RsACE2),na.rm=T)))
    hist(log10(dt$Kd_RpACE2),col="gray40",breaks=60,xlab="log10(K_D,app), RpACE2 (M)",main="",xlim=log10(range(c(dt$Kd_huACE2,dt$Kd_RaACE2,dt$Kd_RpACE2,dt$Kd_RsACE2),na.rm=T)))
    hist(log10(dt$Kd_RsACE2),col="gray40",breaks=60,xlab="log10(K_D,app), RsACE2 (M)",main="",xlim=log10(range(c(dt$Kd_huACE2,dt$Kd_RaACE2,dt$Kd_RpACE2,dt$Kd_RsACE2),na.rm=T)))

<img src="compute_binding_Kd_files/figure-gfm/Kd_distribution-1.png" style="display: block; margin: auto;" />

    #save pdf
    invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/hist_Kd-per-barcode.pdf",sep="")))

Let’s take a look at some of the curves with *K*<sub>D,app</sub> values
across this distribution to get a broad sense of how things look.

First, curves with *K*<sub>D,app</sub> fixed at the 10<sup>-5</sup>
maximum (RpACE2 at 10<sup>-6</sup>). We can see these are all flat-lined
curves with no response.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2==max(dt$Kd_huACE2,na.rm=T))[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2==max(dt$Kd_RaACE2,na.rm=T))[1])
    plot.titration.RpACE2(which(dt$Kd_RpACE2==max(dt$Kd_RpACE2,na.rm=T))[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2==max(dt$Kd_RsACE2,na.rm=T))[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-5_Kd-1.png" style="display: block; margin: auto;" />

Next, with *K*<sub>D,app</sub> around 10<sup>-6</sup>. These
10<sup>-6</sup> curves are similar to the 10<sup>-5</sup>, though
perhaps a bit more belieivable that the curve is bending at the highest
concentration. Overall, I think that differences in *K*<sub>D,app</sub>
within this large mode in the distribution of all variant
*K*<sub>D,app</sub>s probably isn’t reporting on much meaningful
variation.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-6 & dt$Kd_huACE2 < 1.2e-6)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-6 & dt$Kd_RaACE2 < 1.2e-6)[1])
    plot.titration.RpACE2(which(dt$Kd_RpACE2 > 9e-7 & dt$Kd_RpACE2 < 1.2e-6)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-6 & dt$Kd_RsACE2 < 1.2e-6)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-6_Kd-1.png" style="display: block; margin: auto;" />

With *K*<sub>D,app</sub> around 10<sup>-7</sup>, we seem to be picking
up more consistent signals. Many of these curves still are showing
noise, but not as consistently as the prior curves.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-7 & dt$Kd_huACE2 < 1.2e-7)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-7 & dt$Kd_RaACE2 < 1.2e-7)[1])
    plot.titration.RpACE2(which(dt$Kd_RpACE2 > 9e-8 & dt$Kd_RpACE2 < 2e-7)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-7 & dt$Kd_RsACE2 < 1.2e-7)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-7_Kd-1.png" style="display: block; margin: auto;" />

At *K*<sub>D,app</sub> of 10<sup>-8</sup>, we are likewise picking up
some signal but also some curves with lower plateaus, etc.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-8 & dt$Kd_huACE2 < 1.2e-8)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-8 & dt$Kd_RaACE2 < 1.2e-8)[1])
    plot.titration.RpACE2(which(dt$Kd_RpACE2 > 5e-9 & dt$Kd_RpACE2 < 4e-8)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-8 & dt$Kd_RsACE2 < 1.2e-8)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-8_Kd-1.png" style="display: block; margin: auto;" />

At *K*<sub>D,app</sub> of 10<sup>-9</sup>, curves look great. (No
RpACE2s titrate at this low of affinity)

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-9 & dt$Kd_huACE2 < 1.2e-9)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-9 & dt$Kd_RaACE2 < 1.2e-9)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-9 & dt$Kd_RsACE2 < 1.2e-9)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-9_Kd-1.png" style="display: block; margin: auto;" />

At *K*<sub>D,app</sub> of 10<sup>-10</sup>, curves are looking
beautiful.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-10 & dt$Kd_huACE2 < 1.2e-10)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-10 & dt$Kd_RaACE2 < 1.2e-10)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-10 & dt$Kd_RsACE2 < 1.2e-10)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-10_Kd-1.png" style="display: block; margin: auto;" />

As do curves with *K*<sub>D,app</sub> \~ 10<sup>-11</sup>.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-11 & dt$Kd_huACE2 < 1.2e-11)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-11 & dt$Kd_RaACE2 < 1.2e-11)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-11 & dt$Kd_RsACE2 < 1.2e-11)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-11_Kd-1.png" style="display: block; margin: auto;" />

And *K*<sub>D,app</sub> \~ 10<sup>-12</sup>. In contrast to the
SARS-CoV-2 DMS where this was enriched for noise once we got to this
tail, in this case, it does seem we may have many genuine mutants at
this high affinity range for some of these ligands.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$Kd_huACE2 > 1e-12 & dt$Kd_huACE2 < 1.2e-12)[1])
    plot.titration.RaACE2(which(dt$Kd_RaACE2 > 1e-12 & dt$Kd_RaACE2 < 1.2e-12)[1])
    plot.titration.RsACE2(which(dt$Kd_RsACE2 > 1e-12 & dt$Kd_RsACE2 < 2.2e-12)[1])

<img src="compute_binding_Kd_files/figure-gfm/1e-12_Kd-1.png" style="display: block; margin: auto;" />

Data filtering by fit quality
-----------------------------

Next, let’s compute a quality metric for each curve fit, and filter out
poor fits from our final dataset for downstream analyses. For each curve
fit, we will compute a *normalized* mean square residual (nMSR). This
metric computes the residual between the observed response variable and
that predicted from the titration fit, normalizes this residual by the
response range of the titration fit (which is allowed to vary between
1.5 and 3 per the titration fits above), and computes the mean-square of
these normalized residuals.

    #function to calculate mean squared residual normalized to response range
    calc.nMSR <- function(y.obs,x.vals,count.vals,response,fit,cfu.cutoff=cutoff){
      indices <- count.vals>cfu.cutoff
      y.obs <- y.obs[indices]
      x.vals <- x.vals[indices]
      y.pred <- predict(fit,newdata=list(x=x.vals))
      resid <- y.obs - y.pred
      resid.norm <- resid/response
      return(mean((resid.norm)^2,na.rm=T))
    }

    #calculate normalized MSR for each fit
    #huACE2
    dt[!is.na(Kd_huACE2), nMSR_huACE2 := calc.nMSR(y.obs=c(huACE2_01_meanbin,huACE2_02_meanbin,huACE2_03_meanbin,huACE2_04_meanbin,
                                                           huACE2_05_meanbin,huACE2_06_meanbin,huACE2_07_meanbin,huACE2_08_meanbin,huACE2_09_meanbin),
                                                   x.vals=samples_huACE2$conc,
                                                   count.vals=c(huACE2_01_totalcount,huACE2_02_totalcount,huACE2_03_totalcount,huACE2_04_totalcount,
                                                               huACE2_05_totalcount,huACE2_06_totalcount,huACE2_07_totalcount,huACE2_08_totalcount,huACE2_09_totalcount),
                                                   response=response_huACE2,
                                                   fit=fit_huACE2[[1]]),by=c("barcode","library")]
    #RaACE2
    dt[!is.na(Kd_RaACE2), nMSR_RaACE2 := calc.nMSR(y.obs=c(RaACE2_01_meanbin,RaACE2_02_meanbin,RaACE2_03_meanbin,RaACE2_04_meanbin,
                                                           RaACE2_05_meanbin,RaACE2_06_meanbin,RaACE2_07_meanbin,RaACE2_08_meanbin,RaACE2_09_meanbin),
                                                   x.vals=samples_RaACE2$conc,
                                                   count.vals=c(RaACE2_01_totalcount,RaACE2_02_totalcount,RaACE2_03_totalcount,RaACE2_04_totalcount,
                                                               RaACE2_05_totalcount,RaACE2_06_totalcount,RaACE2_07_totalcount,RaACE2_08_totalcount,RaACE2_09_totalcount),
                                                   response=response_RaACE2,
                                                   fit=fit_RaACE2[[1]]),by=c("barcode","library")]
    #RpACE2
    dt[!is.na(Kd_RpACE2), nMSR_RpACE2 := calc.nMSR(y.obs=c(RpACE2_02_meanbin,RpACE2_03_meanbin,RpACE2_04_meanbin,
                                                           RpACE2_05_meanbin,RpACE2_06_meanbin,RpACE2_07_meanbin,RpACE2_09_meanbin),
                                                   x.vals=samples_RpACE2$conc,
                                                   count.vals=c(RpACE2_02_totalcount,RpACE2_03_totalcount,RpACE2_04_totalcount,
                                                               RpACE2_05_totalcount,RpACE2_06_totalcount,RpACE2_07_totalcount,RpACE2_09_totalcount),
                                                   response=response_RpACE2,
                                                   fit=fit_RpACE2[[1]]),by=c("barcode","library")]
    #RsACE2
    dt[!is.na(Kd_RsACE2), nMSR_RsACE2 := calc.nMSR(y.obs=c(RsACE2_01_meanbin,RsACE2_02_meanbin,RsACE2_03_meanbin,RsACE2_04_meanbin,
                                                           RsACE2_05_meanbin,RsACE2_06_meanbin,RsACE2_07_meanbin,RsACE2_08_meanbin,RsACE2_09_meanbin),
                                                   x.vals=samples_RsACE2$conc,
                                                   count.vals=c(RsACE2_01_totalcount,RsACE2_02_totalcount,RsACE2_03_totalcount,RsACE2_04_totalcount,
                                                               RsACE2_05_totalcount,RsACE2_06_totalcount,RsACE2_07_totalcount,RsACE2_08_totalcount,RsACE2_09_totalcount),
                                                   response=response_RsACE2,
                                                   fit=fit_RsACE2[[1]]),by=c("barcode","library")]

Distribution of the nMSR metric in each set of fits

    par(mfrow=c(2,2))
    hist(dt$nMSR_huACE2,main="huACE2",xlab="Response-normalized mean squared residual",col="gray50",breaks=40)
    hist(dt$nMSR_RaACE2,main="RaACE2",xlab="Response-normalized mean squared residual",col="gray50",breaks=40)
    hist(dt$nMSR_RpACE2,main="RpACE2",xlab="Response-normalized mean squared residual",col="gray50",breaks=40)
    hist(dt$nMSR_RsACE2,main="RsACE2",xlab="Response-normalized mean squared residual",col="gray50",breaks=40)

<img src="compute_binding_Kd_files/figure-gfm/nMSR_distribution-1.png" style="display: block; margin: auto;" />

As we would expect, the MSR stat decreases with cell count, indicating
that higher cell counts leads to better curve fits

    plot(log10(dt$huACE2_avgcount),dt$nMSR_huACE2,pch=19,col="#00000010",xlab="average cell count (log10)",ylab="nMSR",xlim=c(1,3))

<img src="compute_binding_Kd_files/figure-gfm/nMSR_v_cell_count-1.png" style="display: block; margin: auto;" />
Let’s see what titration curves look like in different nMSR regimes.
First, let’s see titrations from each experiment at the median nMSR.
(Iterate the sample to illustrate an actual responsive curve, since
nonresponseive curves with little noise are hard to judge) Suggests that
the ‘typical’ curve looks quite nice.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$nMSR_huACE2 > quantile(dt$nMSR_huACE2,0.49,na.rm=T) & dt$nMSR_huACE2 < quantile(dt$nMSR_huACE2,0.51,na.rm=T))[1])
    plot.titration.RaACE2(which(dt$nMSR_RaACE2 > quantile(dt$nMSR_RaACE2,0.49,na.rm=T) & dt$nMSR_RaACE2 < quantile(dt$nMSR_RaACE2,0.51,na.rm=T))[4])
    plot.titration.RpACE2(which(dt$nMSR_RpACE2 > quantile(dt$nMSR_RpACE2,0.49,na.rm=T) & dt$nMSR_RpACE2 < quantile(dt$nMSR_RpACE2,0.51,na.rm=T))[4])
    plot.titration.RsACE2(which(dt$nMSR_RsACE2 > quantile(dt$nMSR_RsACE2,0.49,na.rm=T) & dt$nMSR_RsACE2 < quantile(dt$nMSR_RsACE2,0.51,na.rm=T))[4])

<img src="compute_binding_Kd_files/figure-gfm/example_curves_median_nMSR-1.png" style="display: block; margin: auto;" />
If we set a cutoff of filtering out the worst 2.5% of curves based on
nMSR, what would the borderline cases look like that we would be
retaining? They aren’t even that bad, and these would be the *worst*
curves we are keeping. So, let’s use this 2.5%ile cutoff.

    par(mfrow=c(2,2))
    plot.titration.huACE2(which(dt$nMSR_huACE2 > quantile(dt$nMSR_huACE2,0.9725,na.rm=T) & dt$nMSR_huACE2 < quantile(dt$nMSR_huACE2,0.9775,na.rm=T))[4])
    plot.titration.RaACE2(which(dt$nMSR_RaACE2 > quantile(dt$nMSR_RaACE2,0.9725,na.rm=T) & dt$nMSR_RaACE2 < quantile(dt$nMSR_RaACE2,0.9775,na.rm=T))[2])
    plot.titration.RpACE2(which(dt$nMSR_RpACE2 > quantile(dt$nMSR_RpACE2,0.9725,na.rm=T) & dt$nMSR_RpACE2 < quantile(dt$nMSR_RpACE2,0.9775,na.rm=T))[1])
    plot.titration.RsACE2(which(dt$nMSR_RsACE2 > quantile(dt$nMSR_RsACE2,0.9725,na.rm=T) & dt$nMSR_RsACE2 < quantile(dt$nMSR_RsACE2,0.9775,na.rm=T))[1])

<img src="compute_binding_Kd_files/figure-gfm/example_curves_cutoff_nMSR-1.png" style="display: block; margin: auto;" />

Next, we will apply this filtering step on normalized MSR, removing the
worst 2.5% of curves on this metric.

    dt[nMSR_huACE2 > quantile(dt$nMSR_huACE2,0.95,na.rm=T),c("Kd_huACE2","Kd_SE_huACE2","response_huACE2","baseline_huACE2","RSE_huACE2","fit_huACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))]
    dt[nMSR_RaACE2 > quantile(dt$nMSR_RaACE2,0.95,na.rm=T),c("Kd_RaACE2","Kd_SE_RaACE2","response_RaACE2","baseline_RaACE2","RSE_RaACE2","fit_RaACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))]
    dt[nMSR_RpACE2 > quantile(dt$nMSR_RpACE2,0.95,na.rm=T),c("Kd_RpACE2","Kd_SE_RpACE2","response_RpACE2","baseline_RpACE2","RSE_RpACE2","fit_RpACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))]
    dt[nMSR_RsACE2 > quantile(dt$nMSR_RsACE2,0.95,na.rm=T),c("Kd_RsACE2","Kd_SE_RsACE2","response_RsACE2","baseline_RsACE2","RSE_RsACE2","fit_RsACE2") := list(as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.list(NA))]

This leaves us with filtered *K*<sub>D,app</sub> estimates for 74.74% of
our lib1+lib2 huACE2 titrations, 78.31%) of the RaACE2 titrations,
80.95%) of the RpACE2 titrations, and 79.23%) of the RsACE2 titrations.

Last, let’s convert our *K*<sub>D,app</sub> to 1) a
log<sub>10</sub>-scale, and 2) *K*<sub>A,app</sub>, the inverse of
*K*<sub>D,app</sub>, such that higher values are associated with tighter
binding, as is more intuitive. (If we want to continue to discuss in
terms of *K*<sub>D,app</sub>, since people are often more familiar with
*K*<sub>D</sub>, we can refer to the
log<sub>10</sub>(*K*<sub>A,app</sub>) as
-log<sub>10</sub>(*K*<sub>D,app</sub>), which are identical.

    dt[,log10Kd_huACE2 := log10(Kd_huACE2),by=c("barcode","library")]
    dt[,log10Kd_RaACE2 := log10(Kd_RaACE2),by=c("barcode","library")]
    dt[,log10Kd_RpACE2 := log10(Kd_RpACE2),by=c("barcode","library")]
    dt[,log10Kd_RsACE2 := log10(Kd_RsACE2),by=c("barcode","library")]


    dt[,log10Ka_huACE2 := -log10Kd_huACE2,by=c("barcode","library")]
    dt[,log10Ka_RaACE2 := -log10Kd_RaACE2,by=c("barcode","library")]
    dt[,log10Ka_RpACE2 := -log10Kd_RpACE2,by=c("barcode","library")]
    dt[,log10Ka_RsACE2 := -log10Kd_RsACE2,by=c("barcode","library")]

Let’s visualize the final binding measurements as violin plots for the
different wildtype targets

    p1 <- ggplot(dt[!is.na(log10Ka_huACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_huACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("huACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(6,12))

    p2 <- ggplot(dt[!is.na(log10Ka_RaACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RaACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RaACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(6,12))

    p3 <- ggplot(dt[!is.na(log10Ka_RpACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RpACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RpACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(6,12))

    p4 <- ggplot(dt[!is.na(log10Ka_RsACE2) & variant_class=="wildtype",],aes(x=target,y=log10Ka_RsACE2))+
      geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
      ggtitle("RsACE2")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))+
      scale_y_continuous(limits=c(6,12))

    grid.arrange(p1,p2,p3,p4,ncol=1)

    ## Warning: Removed 18242 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 18242 rows containing non-finite values (stat_summary).

    ## Warning: Removed 12812 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 12812 rows containing non-finite values (stat_summary).

    ## Warning: Removed 14805 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 14805 rows containing non-finite values (stat_summary).

<img src="compute_binding_Kd_files/figure-gfm/final_pheno_DFE-1.png" style="display: block; margin: auto;" />

    #save pdf
    invisible(dev.print(pdf, paste(config$Titeseq_Kds_dir,"/violin-plot_log10Ka-by-target.pdf",sep="")))

Data Output
-----------

Finally, let’s output our measurements for downstream analyses.

    dt[,.(library,barcode,target,variant_class,wildtype,position,mutant,
         huACE2_avgcount,log10Kd_huACE2,log10Ka_huACE2,nMSR_huACE2,
         RaACE2_avgcount,log10Kd_RaACE2,log10Ka_RaACE2,nMSR_RaACE2,
         RpACE2_avgcount,log10Kd_RpACE2,log10Ka_RpACE2,nMSR_RpACE2,
         RsACE2_avgcount,log10Kd_RsACE2,log10Ka_RsACE2,nMSR_RsACE2)] %>%
      mutate_if(is.numeric, round, digits=5) %>%
      write.csv(file=config$Titeseq_Kds_file, row.names=F)