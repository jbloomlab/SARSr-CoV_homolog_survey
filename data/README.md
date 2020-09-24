# Input data
This directory contains input data for the analysis.

## Basic information about sequences and alignments

These files are used for the basic processing of the deep sequencing data to call variants by barcode and count barcodes:

   - [PacBio_amplicons.gb](PacBio_amplicons.gb): the amplicons being sequenced by PacBio.
     Note that there are a variety of possible amplicons based on the different unmutated RBD sequences.

   - [feature_parse_specs.yaml](feature_parse_specs.yaml): how to parse the amplicon when handling the PacBio data.

   - [PacBio_runs.csv](PacBio_runs.csv): list of the PacBio runs used to call the variants.

   - [barcode_runs.csv](barcode_runs.csv): list of the Illumina runs used to count the barcodes for different samples. This file is built from [barcode_runs_orig-names.csv](barcode_runs_orig-names.csv) by the Jupyter notebook [build_barcode_runs.ipynb](build_barcode_runs.ipynb).

   - [RBD_sites.csv](RBD_sites.csv): gives site and residue information for SARS-CoV-2, including alignment of the RBD integer numbering with the Spike numbering for SARS-CoV-2 RBD, alignment to SARS-CoV, and structural annotations.
   
   - [./plasmid_maps](plasmid_maps): gives our base SARS-CoV-2 yeast display vector sequence, including the modifications made for a barcode landing pad per our library generation scheme. 2649 illustrates what the plasmid looks like after the insertion of a mutagenized amplicon with an appended N16 barcode. Note, the insert sequences vary with RBD sequence.

   - [./primers/](primers/) contains important primer sequences used in the study, including those used in amplicon barcoding, and Illumina sequencing preparation.
   
