# Summary

Analysis run by [Snakefile](../../Snakefile)
using [this config file](../../config.yaml).
See the [README in the top directory](../../README.md)
for details.

Here is the DAG of the computational workflow:
![dag.svg](dag.svg)

Here is the Markdown output of each Jupyter notebook in the
workflow:

1. [Process PacBio CCSs](process_ccs.md). Creates a [table](../variants/nucleotide_variant_table.csv)
   linking barcodes to the mutations in the variants.

2. [Count variants by barcode](count_variants.md).
   Creates a [variant counts file](../counts/variant_counts.csv)
   giving counts of each barcoded variant in each condition.

3. [Parse amino acid mutants and merge PacBio and Illumina sequencing data](merge_sequencing.md).

4. [Fit titration curves](compute_binding_Kd.md) to calculate per-barcode K<sub>D</sub>, recorded in [this file](../binding_Kds/binding_Kds.csv).

5. [Analyze Sort-seq](compute_expression_meanF.md) to calculate per-barcode RBD expression, recorded in [this file](../expression_meanFs/expression_meanFs.csv).

6. [Derive final genotype-level phenotypes from replicate barcoded sequences](barcode_to_genotype_phenotypes.md).
   Generates final phenotypes, recorded in [this file for wildtype backgrounds](../final_variant_scores/wt_variant_scores.csv) and [this file for mutants](../final_variant_scores/mut_variant_scores.csv).