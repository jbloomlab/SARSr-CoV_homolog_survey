# Phylogenetic analysis of SARS-related CoV RBD sequences

This README details the series of bioinformatic analyses that were conducted to infer ancestral sequences on the RBD phylogeny. Many steps are only pseudo-automated, so this analysis is not amenable to streamlining in a pipeline. (E.g. converting alignment files, outputting rooted phylogenies, downloads of some niche bioinformatics software, etc.)

## Sequence set

I used the base set of CoV RBD sequences used in the phylogeny presented in our DMS manuscript (sequences [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/blob/master/data/alignments/unaligned-sequences/RBD_nt.fasta)), which consiste of the curated set of SARS-related CoV RBD sequences from [Letko et al. 2020](https://www.nature.com/articles/s41564-020-0688-y), supplemented with newly described SARS-CoV-2-clade sequences from bat and pangolin.

Beyond that, I want to use this analysis to better define the relationship among SARS-CoV-1 epidemic isolates. To do this, I downloaded all sequences given in Supplementary Table 1 from [Song et al. 2005](https://www.pnas.org/content/102/7/2430). I parsed these sequences to their RBD nt sequences, and used the `cd-hit` program to eliminate any identical RBD nt sequences from within this SARS-CoV-1 set. (I have this program downloaded on my local laptop, so this was all done there). This yielded 15 unique SARS-CoV-1 sequecnes between the civet and human isolates.

I then updated the header names to conform with the nomenclature used for the SARS-related isolates (`name_accession`), adding to the name the helpful nomenclature used by Song et al. to distinguish civet (`PC`) and human (`HP`) sequences from the main 2002-2003 (`03`) or sporadic 2003-2004 (`04`) outbreaks, and indicating "Early" (`E`), "Middle" (`M`), or "Late" (`L`) epidemic phase for sequences within the primary 2002-2003 epidemic.

I concatenated this expanded SARS-CoV-1 set with the SARS-related CoV sequences from before. This yields the final set of 51 RBD nt sequences included in the `./unaligned-sequences/RBD_nt.fasta` file. I also translated each nt sequence to its amino acid sequence, which is saved in the file `.unaligned-sequences/RBD_aa.fasta`.

Finally, we may want to try using non-sarbecovirus outgroup sequences. I initially tried including adding the non-sarbecovirus sequence from the Hibecovirus strain Bat Hp-BCoV Zhejiang 2013 (GenBank KF636752) to use as an outgroup for the sarbecovirus clade. This is the closest outgroup sequence I've seen in larger phylogenies of betaCoVs.

We could also try more distant betaCoV outgroups: Rousettus CoV HKU9 (Genbank HKU9); Bat CoV GCCDC1 (MT350598); OC43 (KX344031); HKU1 (KF686346); MERS-CoV (NC_019843). A recent [preprint](https://www.biorxiv.org/content/10.1101/2020.07.07.190546v1) illustrates that the positioning of the root may be very important (and I think they may actually get it wrong for the RBD tree? It differs from what e.g. Boni et al. report, and their RBD tree differs in other respects from what I've seen in other papers.) So, instead of doing just at the level of RBD, we may want to assess robustness of rooting by comparing both spike and RBD. These do not align well when using just RBD -- we may want to extend our consideration of rooting and alignment to be the whole spike before parsing back down to RBD.

## Rooting

We need to have confidence in our rooting of our RBD tree. Inclusion of known outgroups in an RBD only alignment looked terrible -- but perhaps if we begin the alignment with whole spike protein, the alignments will look better? We can then either infer the whole spike tree (for more phylogenetic signal), or parse this down to RBD only, to try to glean insight into the proper rooting of the sarbecovirus ingroup sequences.

Started a folder `./outgroup-rooting`, which has its own `unaligned-sequences` subdirectory containing full spike sequences for the outgroups discussed above, along with our more restricted sarbecovirus set (only one SARS-CoV-1 sequence, Urbani).

First, align amino acid sequences:

```
mafft --reorder --op 4.5 ./outgroup-rooting/unaligned-sequences/outgroups_spike_aa.fasta > ./outgroup-rooting/outgroups_spike_aa_aligned.fasta
```

Also make the nt alignment from amino acid sequences using PAL2NAL:

```
pal2nal.pl ./outgroup-rooting/outgroups_spike_aa_aligned.fasta ./outgroup-rooting/unaligned-sequences/outgroups_spike_nt.fasta > ./outgroup-rooting/outgroups_spike_nt_aligned.clustal
```

Manually cleaned up alignments, deleted some super gappy regions. Saved as `./outgroup-rooting/outgroups_spike_aa_aligned_cleaned.fasta` and `./outgroup-rooting/outgroups_spike_nt_aligned_cleaned.fasta`

Try out a tree inferences using RAxML:

```
mkdir ./outgroup-rooting/spike_aa_tree
cd ./outgroup-rooting/spike_aa_tree
nohup raxmlHPC-PTHREADS -s ../outgroups_spike_aa_aligned_cleaned.fasta -n spike_aa_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10 -T 8 &
```

```
mkdir ./outgroup-rooting/spike_nt_tree
cd ./outgroup-rooting/spike_nt_tree
nohup raxmlHPC-PTHREADS -s ../outgroups_spike_nt_aligned_cleaned.fasta -n spike_nt_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q codon_partitions.txt &
```

In either of these full spike phylogenies, rooting on these outgroups puts the Europe/Africa sequences as first to diverge. Let's see if we can get these incorporated now into the RBD only trees? Both to further support this rooting, but also to enable reconstruction of the ancestral sarbecovirus. 

First, truncated the uncleaned/full alignment to just the RBD sequences. in the alignments `./outgroup-rooting/outgroups_RBD_aa_aligned.fasta` and `_nt` equivalent, and then removed outgroup-specific inserted sequences. I then removed any gaps (that is, unaligned the sequences), added back in the SARS-CoV-1 sequences that I curated above, and saved the files in `./outgroup-rooting/unaligned-sequences/`

Do the same drill -- align by amino acid sequence, then align DNA sequence using PAL2NAL

```
mafft --reorder --op 4.5 ./outgroup-rooting/unaligned-sequences/outgroups_RBD_aa.fasta > ./outgroup-rooting/outgroups_RBD_aa_aligned.fasta
```

```
pal2nal.pl ./outgroup-rooting/outgroups_RBD_aa_aligned.fasta ./outgroup-rooting/unaligned-sequences/outgroups_RBD_nt.fasta > ./outgroup-rooting/outgroups_RBD_nt_aligned.clustal
```


And then do phylogenetic inference:
```
mkdir ./outgroup-rooting/RBD_aa_tree
cd ./outgroup-rooting/RBD_aa_tree
nohup raxmlHPC-PTHREADS -s ../outgroups_RBD_aa_aligned.fasta -n RBDo_aa_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10 -T 8 &
```

```
mkdir ./outgroup-rooting/RBD_nt_tree
cd ./outgroup-rooting/RBD_nt_tree
nohup raxmlHPC-PTHREADS -s ../outgroups_RBD_nt_aligned.fasta -n RBDo_nt_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q codon_partitions.txt &
```

These trees also put the first split between the non-Asian and SE-Asian sarbecoviruses, as I had been operating under! So, that's good to see (and differs from the recent paper on bioRxiv that we had been discussing). The outgroup branches though are *incredibly* long, and it would probably look better if we pared back to just the Hibecovirus outgroup? Therefore, I add only the Hibecovirus sequence into our "main" unaligned sequences, to be used below.


## Alignment with `mafft` and `PAL2NAL`

I used `mafft` to align RBD aa sequences:

```
mafft --reorder --op 4.5 ./unaligned-sequences/RBD_aa.fasta > ./RBD_aa_aligned.fasta
```

Followed by `PAL2NAL` to align nt sequences given codon structure and the given amino acid alignment.

```
pal2nal.pl RBD_aa_aligned.fasta ./unaligned-sequences/RBD_nt.fasta > RBD_nt-codon_aligned.clustal
```

I converted our `.clustal` format nt alignment to `.fasta` using [this tool](https://www.hiv.lanl.gov/content/sequence/FORMAT_CONVERSION/form.html).

## Phylogenetic inference with `RAxML`

I next used `RAxML` to infer the ML tree. We want to use a codon-aware nt substitution model -- there is not enough amino acid diversity for some of the fine-grained resolution (i.e. the SARS-CoV-1 sequences and their nearest bat relatives), so we want to use a nt-aware model which gives us the extra information available in synonymous codon positions to gain additional resolution into relationships. So, we fit GTRGAMMA nt models with spearate data partitions (i.e. separate model parameters) fit to the first, second, and third codon positions -- this is how the RAxML manual suggests inference on codon-aligned protein-coding sequences.


```
mkdir RBD_codon_tree
cd ./RBD_codon_tree
nohup raxmlHPC-PTHREADS -s ../RBD_nt-codon_aligned.fasta -n RBD_codon_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &
```

Open up the tree `./RBD_codon_tree/RAxML_bestTree.RBD_codon_tree.txt` in FigTree, root on the Hp-BCoV\_Zhejiang outgrpus sequence, output tree as `RAxML_besttree_RBD_codon_rooted.txt` in Newick format. This will polarize the direction of nodes for ASR (not necessary, but some of the annotations that come out from ASR software are easier to interpret if polarized in the proper direction).

## Recombination analysis with `3SEQ`

We may be concerned about doing ASR on this tree if there is clear evidence for recombination, especially if it is between the four major clades of lineages. To test for recombination, I downloaded the `3seq` package from [here](https://mol.ax/software/3seq/) and installed per the instructions. I generated a p-value table within the ``../programs/3seq` directory by executing:

```
./3seq -g my3seqTable700 700
```

Next, to test for recombination within our nt/codon alignment, within the `RBD_nt-codon_3seq` subdirectory I executed:

```
../../programs/3seq/3seq -f ../RBD_nt-codon_aligned.fasta -L50
```

The main thing that crops up looking through the results of this analysis, is that there is putative recombination involving Rs4231, LYRa11, and the 02-03 SARS-CoV-1 isolates (though not the 03-04 isolates, interesting?). In looking at the stats from the run, there are ~30 unique shared states between SARS-CoV-1 sequences and Rs4231, and 83 informative characters linking LYRa11 and SARS-CoV-1 -- so it's lopsided, and then these are organized such that the maximum run suggestive of recombinant origins (the 'k' statistic) is 66. This leads to predicted breakpoints 208-519. (It seems the LYRa11 similarity is projected as being more similar to SARS-CoV-1 within this 210-519 stretch) -- though even outside of this window, it doesn't seem glaringly obvious to me that these SARS-CoV-1 isolates are closer to Rs4231 than LYRa11? And the places where the SARS-CoV-1 does share the Rs4231 state seems to be places where LYRa11 has lineage-specific mutations. It also seems to me that if there is true recombination, wouldn't we expect the hits to be significant with multiple genotypes within the originating clade? LYRa11 is a rather unique sequence, but Rs4231 has several closely related sequences that could further hammer down recombinant origins, if real. So, I think this might not be genuine?

To see if there is anything particular going on here, I took the RBD\_nt-codon alignment, and split it into two partitions: `RBD_nt-codon_aligned_rec1.fasta` contains nucleotides 1-207 and 520-606, while `RBD_nt-codon_aligned_rec2.fasta` contains nucleotides 208-519.

I then inferred phylogenies separate for these two putative recombinant regions:

```
mkdir RBD_nt-codon_recomb-test
cd RBD_nt-codon_recomb-test
nohup raxmlHPC-PTHREADS -s ./RBD_nt-codon_aligned_rec1.fasta -n RBD_codon_tree_rec1.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &

nohup raxmlHPC-PTHREADS -s ./RBD_nt-codon_aligned_rec2.fasta -n RBD_codon_tree_rec2.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &> nohup2.out &

```

Looking at the trees built from these two different sections -- there are differences from the full-length ML tree, but any differences have incredibly low bootstrap support. I do *not* take this as indication of strong support for a specific recombination history here. I will keep these analyses here in case we want to revisit -- might still be worth reconstructing e.g. the SARS-CoV-1 or SARS-CoV-1-clade ancestral sequence under the full length versus chunked ML trees, to ensure there is no major difference in sequence or function. But I am not compelled to use this chunked sequence at face value as I don't think it points to any clear signal of recombination.

## Ancestral sequence reconstruction with `FastML`

Next, we reconstruct ancestral sequences on the ML phylogeny using [`FastML`](http://fastml.tau.ac.il), which I downloaded and installed in the `../programs` directory per the instructions [here](http://fastml.tau.ac.il/source.php#run).

I want to calculate posterior probabilities of ancestral amino acids, and likelihood-based placement of indels, using the phylogeny inferred from nucleotide-level sequences (`./RBD_codon_tree`), the amino acid alignment (`./RBD_aa_aligned.fasta`), the LG model of amino acid substitutions, and optimized branch lengths. To do this, within the directory `./RBD_ML_reconstruction`, I executed the following command:

```
nohup perl ../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned.fasta --seqType AA --Tree ../RBD_codon_tree/RAxML_besttree_RBD_codon_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_ML_reconstruction &
```

The data table giving the probabilities of ancestral states is `./prob.marginal.csv`, which we can supplement with the file `./IndelsMarginalProb.txt` to indicate which of these ancestral states are censored by the reconstruction of deleted sequence at that position. The _maximum a posteriori_ (MAP) sequences are given in `./seq.marginal_IndelAndChars.txt`, but we want to supplement this with the full probability table to indicate where the MAP state has considerable uncertainty in its presence/absence or amino acid identity that we want to further characterize. The link between ancestor names and their location on the phylogeny is labeled in the file `./tree.newick.txt`. We will parse these files in the `.Rmd` analysis at `../parse_ancestors.Rmd`.

## NT level ASR

For the finer grained reconstructions within the SARS-CoV-1 and SARS-CoV-2 trees, we might do better reconstructions using the nt gene sequences. nt signals are likely saturated in the deeper nodes (and even the AA reconstructions have noticeable ambiguities in these deeper nodes), but nt might give better signal in the fine-grained clade reconstructions. To do this, within the directory `./RBD_ML_codon_reconstruction`, I ***wanted** to execute the following command:

```
nohup perl ../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_nt-codon_aligned.fasta --seqType codon --Tree ../RBD_codon_tree/RAxML_besttree_RBD_codon_rooted.txt --intersectTreeAndSeq --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_ML_codon_reconstruction &
```

But I was hitting a bug I couldn't for the life of me track, in which it was saying the sequence names in the alignment did not match the phylogeny -- which as far as I can tell is *not* a me problem. These files worked just fine when uploaded to the [FastML server](https://fastml.tau.ac.il), so I just ran this analysis there and downloaded the results to the corresponding folder. I just took the MAP ancestors from this analysis, so no need to parse.

Took the file `outFile_seq_marginal.fasta`, tranlsated to amino acid sequences, and took a look. The one sequence I was most concerned about on the AA tree (AncSARS2b, being the same sequecne as GD-Pangolin) doesn't change as dramatically as I thought it might on the codon tree -- the only change is that K417R becomes a change on the branch to GD-Pangolin. Overall, not too concerned about this amino acid versus codon discrepency, so I am not going to include these in current iteration.

## Robustness of ASR to phylogenetic uncertainty

There are some poorly supported nodes within this phylogeny. Many (e.g. within the Clade 2 sequences) we are not particularly concerned by. However, we might want to perform ancestral reconstructions on trees that differ in the global relationships between the SARS-CoV-1, SARS-CoV-2, and "clade 2" sequences, particularly as it pertains to the question of ACE2 binding in the ancestor of the SE-Asian sarbecoviruses.

The ML tree has the SE-Asia ingroup sequences with SARS-CoV-1 clade sister to "clade2". We inferred phylogenies under constraints in which we instead enforce sister relationships between the SARS-CoV-2 clade and clade2 (alt1), or between the SARS-CoV-1 and SARS-CoV-2 clades (alt2).

These constraints are set in the files `./RBD_codon_tree_constrained/alt1/constraint_alt1.txt` and `./RBD_codon_tree_constrained/alt2/constraint_alt2.txt` files.

Then, within each `altx` subdirectory, ran RAxML:

```
nohup raxmlHPC-PTHREADS -s ../../RBD_nt-codon_aligned.fasta -n RBD_codon_tree_alt1.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../../codon_partitions.txt -g constraint_alt1.txt &

nohup raxmlHPC-PTHREADS -s ../../RBD_nt-codon_aligned.fasta -n RBD_codon_tree_alt2.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../../codon_partitions.txt -g constraint_alt2.txt &
```

Opened each tree in Figtree, rooted on the outgroup, and saved the rooted tree back to this directory.

And then finally, within each, we ran FastML ancestral reconstructions:

```
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../../RBD_aa_aligned.fasta --seqType AA --Tree ./RAxML_bestTree_RBD_codon_rooted_alt1.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_codon_tree_constrained/alt1/ &
```

and 

```
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../../RBD_aa_aligned.fasta --seqType AA --Tree ./RAxML_bestTree_RBD_codon_rooted_alt2.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_codon_tree_constrained/alt2/ &
```

I then ran the script `./parse_ancestral_RBDs.Rmd` to generate a csv with the relevant ancestors, in `./parsed_sequences/ancestors_table_collapsed.csv`. Formatted these into an alignment.

For gene synthesis and downstream -- I also made `_unique` versions of the aligned ancesotrs and the original `.csv`, in which I removed redundant sequences from the alignment, and in the table, keyed which other sequence is the same as another in the relevant table entry. Finally, I used the program `codon-harmony` (currently in the BloomLab conda environemnt, not in this local repo) to design codon optimized sequences, avoiding relevant restriction enzymes and avoiding homopolymer runs of >4 nt. Within the folder `./parsed_sequences/codon-optimization/`, executed the following:

```
codon_harmony --input  ./RBD-set_aa-seq.fasta --output ./RBD-set_codon-opt.fasta --host 4932 --verbose 2 --local-homopolymer-threshold 4 --restriction-enzymes NotI NdeI XhoI SacI EcoRI --max-relax 0.2
```

## What about RaTG13 and possible recombination?

There is a suggestion that the RBM of RaTG13 was acquired via recombination with some unsampled, more distantly related strain. Without identification of this originating strain, it is hard to define the breakpoints (or know whether it's truly recombination, vs for example branch-specific positive selection focused on this interaction interface). Therefore, in lieu of trying to disentagle recombinant features, I simply want to remove RaTG13 from the alignment, and see how that impacts e.g. the GD- and GX-pangolin and deeper ancestors. If substantial, I can simply include these in the RBD panel for additional robustness considerations.

Made a directory `RBD_rm_RaTG13`, and copied the alignments renamed `RBD_nt-codon_alignment_rmRaTG13.fasta` and `RBD_aa_alignment_rmRaTG13.fasta`, and removed the `RaTG13` entry.

```
nohup raxmlHPC-PTHREADS -s ./RBD_nt-codon_aligned_rmRaTG13.fasta -n RBD_codon_tree_rmRaTG13.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &
```

Root tree and save file as `RAxML_besttree_RBD_codon_rmRaTG13_rooted.txt`.

Do ancestral reconstruction using FastML

```
mkdir ASR
cd ASR
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_rmRaTG13.fasta --seqType AA --Tree ../RAxML_besttree_RBD_codon_rmRaTG13_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_rm_RaTG13/ASR &
```

Added new sequences from AncSarbecovirus down to lineage of SARS-CoV-2 in the table in the `../parsed_sequences/ancesotrs_Table_collapsed_unique.csv` table. Nodes from AncSarbecovirus down to SARS-CoV-2:
 * AncSarbecovirus has some differences from MAP, but these are the ambiguities present in e.g. altALL and tree1 alternatives -- that is, this reflects general ambiguity in this sequence (probably in part due to the different paraphyletic relationship between Kenya and Bulgaria sequences due to the root placement, which I don't really believe anyway). I do not attribute this to the RaTG13 difference.
 * AncAsia same MAP and rmRaTG13
 * AncSARS2a same MAP and rmRaTG13
 * AncSARS2b: only difference relative to MAP is N519H. This was previously inferred to be a substitution between AncSARS2c and SARS-CoV-2. This is not expected to be a consequential mutation, and beyond that, this is not the region where we expect RaTG13 may have a recombinant portion (that is, it is outside the RBM). Therefore, I don't think it demands synthesis of another gene to test for robustness, as it is pretty easy to dismiss.
