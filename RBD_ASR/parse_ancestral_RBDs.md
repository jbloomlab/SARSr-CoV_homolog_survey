Parse ancestral RBDs
================
Tyler Starr
7/9/2020

  - [ASR on the ML phylogeny](#asr-on-the-ml-phylogeny)

After the series of bioinformatics steps outlined in
`./RBD_ASR/README.md`, I have FastML outputs for ancestrally
reconstructed RBD sequences. These data consist of a table giving the
posterior probability of each amino acid at each position at each
internal node, followed by a second table that gives the likelihood that
each internal node has each position deleted or not given a
likelihood-based indel reconstruction. Here, we combine these two pieces
of information to generate our MAP ancestors. We also consider some
alternative ancestral sequences, including incorporating ambiguously
reconstructed positions in the context of the ML phylogeny, as well as
reconstructions under additional reasonable phylogenetic hypotheses.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in file giving concordance between RBD numbering and SARS-CoV-2 Spike numbering, along with other possibly interesting annotations
RBD_sites <- read.csv(file="../data/RBD_sites.csv",stringsAsFactors=F)

#make output directory
if(!file.exists("./parsed_sequences")){
 dir.create(file.path("./parsed_sequences"))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_haswellp-r0.3.1.so
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
    ##  [1] forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3      
    ##  [4] purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
    ##  [7] tibble_3.0.2      ggplot2_3.2.1     tidyverse_1.2.1  
    ## [10] data.table_1.12.6 yaml_2.2.0        knitr_1.25       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.2       cellranger_1.1.0 pillar_1.4.5     compiler_3.6.1  
    ##  [5] tools_3.6.1      digest_0.6.22    lubridate_1.7.4  jsonlite_1.6    
    ##  [9] evaluate_0.14    lifecycle_0.2.0  gtable_0.3.0     pkgconfig_2.0.3 
    ## [13] rlang_0.4.7      cli_1.1.0        rstudioapi_0.10  haven_2.1.1     
    ## [17] xfun_0.10        withr_2.1.2      xml2_1.2.2       httr_1.4.1      
    ## [21] generics_0.0.2   vctrs_0.3.1      hms_0.5.2        grid_3.6.1      
    ## [25] tidyselect_1.1.0 glue_1.3.1       R6_2.4.0         readxl_1.3.1    
    ## [29] rmarkdown_1.16   modelr_0.1.5     magrittr_1.5     backports_1.1.5 
    ## [33] scales_1.0.0     ellipsis_0.3.0   htmltools_0.4.0  rvest_0.3.4     
    ## [37] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    lazyeval_0.2.2  
    ## [41] munsell_0.5.0    broom_0.7.0      crayon_1.3.4

## ASR on the ML phylogeny

We’ll start by generating the best ancestral reconstructions at nodes of
interest on the ML phylogeny. We will also generate alternative
reconstructions that incorporate ambiguously reconstructed amino acid
states, so we can ultimately characterize the robustness of our
functional measurments to uncertainty in the reconstructions.

To begin, we read in the tables giving the posterior probabilities of
each amino acid at each internal node on the phylogeny, and a table that
gives the probability that a residue is a deletion at each node on the
phylogeny given the indel modeling performed by FastML.

``` r
PP_ML <- data.table(read.csv(file="RBD_ML_reconstruction/prob.marginal.csv"))
indels_ML <- read.table(file="RBD_ML_reconstruction/IndelsMarginalProb.txt",sep="\t",header=T)

#let's add the indels category as a column in PP_ML
for(i in 1:nrow(PP_ML)){
  if(PP_ML$Pos[i] %in% indels_ML$Pos){
    PP_ML$gap[i] <- indels_ML[indels_ML$Pos==PP_ML$Pos[i] & indels_ML$Node==PP_ML$Ancestral.Node[i],"Prob_Of_Indel"]
  }else{
    PP_ML$gap[i] <- 0
  }
}

#recalibrate PP in amino acid states incorporating the probability of a gap character -- that is, such that rowSum is always 1
for(i in 1:nrow(PP_ML)){
  PP_ML[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- PP_ML[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")]*(1-PP_ML[i,gap])
}
```

We want ancestral sequences for the ancestral sarbecovirus, the
ancestral non-Asian RBD, the ancestral SE-Asian RBD, the ancestor of the
three SE Asian clades, and then the lineage of ancestors leading to
SARS-CoV-1 and SARS-CoV-2. A table giving these node labels, shorthand
names, and some “Notes” information is given in the file
`data/ancestors_input.csv`. We load that table in, and append the MAP
sequence, along with some statistics about average statistical support
for ancestral states, and a list of ambiguously reconstructed states (PP
\< 0.8), and the plausible alternative states (PP \> 0.2) at that node.

``` r
ancestors <- read.csv(file="../data/ancestors.csv")

get.MAP <- function(node, table,altall=T){
  dt <- table[Ancestral.Node==node,]
  seq <- vector(mode="character",length=nrow(dt))
  PP.MAP <- vector(mode="numeric",length=nrow(dt))
  ambig.states <- vector(mode="character")
  alt.states <- vector(mode="character")
  if(altall==T){
    altall.seq <- vector(mode="character",length=nrow(dt))
  }
  for(i in dt$Pos){
    MAP.state <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap")[which(dt[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap")]==max(dt[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap")]))]
    if(MAP.state=="gap"){seq[i] <- "-"}else{seq[i] <- MAP.state}
    if(altall==T){
      altall.seq[i] <- seq[i]
    }
    PP.MAP[i] <- dt[i,get(MAP.state)]
    if(PP.MAP[i] < 0.8){
      ambig.states <- c(ambig.states,paste(MAP.state,i,sep=""))
      for(aa in c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap")){
        if(aa != MAP.state){
          if(dt[i,get(aa)]>0.2){
            alt.states <- c(alt.states,paste(aa,i,sep=""))
          }
        }
      }
      if(altall==T){
        aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap")[c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","gap") != MAP.state]
        if(max(dt[i,..aas]) >= 0.2){
          altall.state <- aas[which(dt[i,..aas]==max(dt[i,..aas]))]
          if(altall.state=="gap"){altall.seq[i] <- "-"}else{altall.seq[i] <- altall.state}
        }
      }
    }
  }
  if(altall==T){
    if(sum(altall.seq != seq)>0){
      return(list(list(seq),list(PP.MAP),list(ambig.states),list(alt.states),list(altall.seq)))
    }else{
      return(list(list(seq),list(PP.MAP),list(ambig.states),list(alt.states),list(NA)))
    }
  }else{
    return(list(list(seq),list(PP.MAP),list(ambig.states),list(alt.states)))
  }
}

for(i in 1:nrow(ancestors)){
  ancestors[i,c("MAP_seq","MAP_PP","ambig_states","alt_states","altALL_seq")] <- get.MAP(node=ancestors[i,"node_ML"],table=PP_ML)
}
```

We wanted to address how uncertainty in the tree impacts our ancestral
reconstructions. Given that the separations between the clade of
SARS-CoV-1, SARS-CoV-2, and clade 2 are essentially a polytomy, my
prediction is that the ancestral reconstructions will not differ (at
least, substantially). We aren’t even bothering to reconstruct the
ancestor of the SARS-CoV-1 clade and clade 2, given how uncertain we are
in this sister relationship on the ML tree – this is the only node that
“doesn’t exist” in the same way on these alternatively constrained
phylogenies.

To look at this, let’s load in these reconstructions. In “alt1”, we
constrained SARS-CoV-2 clade and clade2 to be sister; in “alt2”, we
constrained the clades of SARS-CoV-1 and SARS-CoV-2 to be sister. Let’s
see if the reconstructions on any of these trees change the MAP
ancestors.

``` r
PP_tree1 <- data.table(read.csv(file="RBD_codon_tree_constrained/alt1/prob.marginal.csv"))
indels_tree1 <- read.table(file="RBD_codon_tree_constrained/alt1/IndelsMarginalProb.txt",sep="\t",header=T)

#let's add the indels category as a column in PP_tree1
for(i in 1:nrow(PP_tree1)){
  if(PP_tree1$Pos[i] %in% indels_tree1$Pos){
    PP_tree1$gap[i] <- indels_tree1[indels_tree1$Pos==PP_tree1$Pos[i] & indels_tree1$Node==PP_tree1$Ancestral.Node[i],"Prob_Of_Indel"]
  }else{
    PP_tree1$gap[i] <- 0
  }
}

#recalibrate PP in amino acid states incorporating the probability of a gap character -- that is, such that rowSum is always 1
for(i in 1:nrow(PP_tree1)){
  PP_tree1[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- PP_tree1[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")]*(1-PP_tree1[i,gap])
}

#repeat for tree2
PP_tree2 <- data.table(read.csv(file="RBD_codon_tree_constrained/alt2/prob.marginal.csv"))
indels_tree2 <- read.table(file="RBD_codon_tree_constrained/alt2/IndelsMarginalProb.txt",sep="\t",header=T)

#let's add the indels category as a column in PP_tree2
for(i in 1:nrow(PP_tree2)){
  if(PP_tree2$Pos[i] %in% indels_tree2$Pos){
    PP_tree2$gap[i] <- indels_tree2[indels_tree2$Pos==PP_tree2$Pos[i] & indels_tree2$Node==PP_tree2$Ancestral.Node[i],"Prob_Of_Indel"]
  }else{
    PP_tree2$gap[i] <- 0
  }
}

#recalibrate PP in amino acid states incorporating the probability of a gap character -- that is, such that rowSum is always 1
for(i in 1:nrow(PP_tree2)){
  PP_tree2[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")] <- PP_tree2[i,c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")]*(1-PP_tree2[i,gap])
}

#for each node, add ancestors to the table only if it differs from the MAP ancestor
for(i in 1:nrow(ancestors)){
  tree1 <- get.MAP(node=ancestors[i,"node_tree1"],table=PP_tree1,altall=F)
  if(sum(tree1[[1]][[1]] != ancestors[i,"MAP_seq"][[1]])>0){
    ancestors[i,c("tree1_seq","tree1_PP","tree1_ambig_states","tree1_alt_states")] <- tree1
  }else{
    ancestors[i,c("tree1_seq","tree1_PP","tree1_ambig_states","tree1_alt_states")] <- list(list(NA),list(NA),list(NA),list(NA))
  }
  tree2 <- get.MAP(node=ancestors[i,"node_tree2"],table=PP_tree2,altall=F)
  if(sum(tree2[[1]][[1]] != ancestors[i,"MAP_seq"][[1]])>0){
    ancestors[i,c("tree2_seq","tree2_PP","tree2_ambig_states","tree2_alt_states")] <- tree2
  }else{
    ancestors[i,c("tree2_seq","tree2_PP","tree2_ambig_states","tree2_alt_states")] <- list(list(NA),list(NA),list(NA),list(NA))
  }
}
```

Let’s enumerate the substitutions that link each (MAP) ancestral
sequence to its prior ancestor, as that can be helpful to have listed
out.

``` r
#function that compares two sequences and outputs differences
enumerate_subs <- function(seq1, seq2, SARS2_numbering=F){
  diffs <- vector(mode="character")
  for(i in 1:length(seq1)){
    if(seq1[i] != seq2[i]){
      if(SARS2_numbering==T){
        diffs <- c(diffs,paste(seq1[i],RBD_sites[RBD_sites$site_alignment==i,"site_SARS2"],seq2[i],sep=""))
      }else{
        diffs <- c(diffs,paste(seq1[i],i,seq2[i],sep=""))
      }
      
    }
  }
  if(length(diffs)>0){return(paste(diffs,collapse=";"))}else{return("")}
}

#need to manually do for each node since ancestors are not always directly above one another in the table
ancestors[ancestors$node_ML=="N2","subs_alignpos"] <- ""
ancestors[ancestors$node_ML=="N2","subs_SARS2pos"] <- ""

ancestors[ancestors$node_ML=="N3","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N2","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N3","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N3","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N2","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N3","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N4","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N2","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N4","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N2","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N6","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N6","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N6","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N6","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N49","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N49","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N49","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N49","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N50","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N49","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N50","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N50","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N49","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N50","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N51","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N50","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N51","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N51","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N50","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N51","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N28","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N28","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N28","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N4","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N28","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N29","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N28","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N29","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N29","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N28","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N29","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N30","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N29","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N30","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N30","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N29","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N30","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N31","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N30","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N31","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N31","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N30","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N31","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N32","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N31","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N32","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N31","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N33","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N33","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N33","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N33","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N34","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N33","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N34","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N34","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N33","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N34","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N35","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N34","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N35","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N35","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N34","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N35","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N39","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N39","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N39","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N32","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N39","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N42","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N39","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N42","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N39","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N43","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N43","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N43","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N43","MAP_seq"][[1]],SARS2_numbering=T)

ancestors[ancestors$node_ML=="N45","subs_alignpos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N45","MAP_seq"][[1]])
ancestors[ancestors$node_ML=="N45","subs_SARS2pos"] <- enumerate_subs(seq1=ancestors[ancestors$node_ML=="N42","MAP_seq"][[1]],seq2=ancestors[ancestors$node_ML=="N45","MAP_seq"][[1]],SARS2_numbering=T)
```

The ancestors table is useful if doing anything comparative in `R`
because it keeps sequence strings split. However, for compiling
sequences, manual browsing of csv, etc., easiest to collapse some
information into a better output.

``` r
ancestors_out <- subset(ancestors,select = c(node_ML,node_tree1,node_tree2,name,notes,subs_alignpos,subs_SARS2pos))
for(i in 1:nrow(ancestors_out)){
  ancestors_out$MAP_seq[i] <- paste(ancestors$MAP_seq[i][[1]],collapse="")
  ancestors_out$MAP_PP[i] <- mean(ancestors$MAP_PP[i][[1]])
  ancestors_out$ambig_states[i] <- paste(ancestors$ambig_states[i][[1]],collapse=";")
  ancestors_out$alt_states[i] <- paste(ancestors$alt_states[i][[1]],collapse=";")
  ancestors_out$altALL_seq[i] <- paste(ancestors$altALL_seq[i][[1]],collapse="")
  ancestors_out$tree1_seq[i] <- paste(ancestors$tree1_seq[i][[1]],collapse="")
  ancestors_out$tree2_seq[i] <- paste(ancestors$tree2_seq[i][[1]],collapse="")
}
write.csv(ancestors_out,file="./parsed_sequences/ancestors_table_collapsed.csv", quote=F, row.names=F)
```
