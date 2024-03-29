# Phylogenetic analysis of SARS-related CoV RBD sequences

This README details the series of bioinformatic analyses that were conducted to infer ancestral sequences on the RBD phylogeny. Many steps are only pseudo-automated, so this analysis is not amenable to streamlining in a pipeline. (E.g. converting alignment files, outputting rooted phylogenies, downloads of some niche bioinformatics software, etc.)

## Sequence set

I used the base set of CoV RBD sequences used in the phylogeny presented in our DMS manuscript (sequences [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS/blob/master/data/alignments/unaligned-sequences/RBD_nt.fasta)), which consiste of the curated set of SARS-related CoV RBD sequences from [Letko et al. 2020](https://www.nature.com/articles/s41564-020-0688-y), supplemented with newly described SARS-CoV-2-clade sequences from bat and pangolin.

Beyond that, I want to use this analysis to better define the relationship among SARS-CoV-1 epidemic isolates. To do this, I downloaded all sequences given in Supplementary Table 1 from [Song et al. 2005](https://www.pnas.org/content/102/7/2430). I parsed these sequences to their RBD nt sequences, and used the `cd-hit` program to eliminate any identical RBD nt sequences from within this SARS-CoV-1 set. (I have this program downloaded on my local laptop, so this was all done there). This yielded 15 unique SARS-CoV-1 sequecnes between the civet and human isolates.

I then updated the header names to conform with the nomenclature used for the SARS-related isolates (`name_accession`), adding to the name the helpful nomenclature used by Song et al. to distinguish civet (`PC`) and human (`HP`) sequences from the main 2002-2003 (`03`) or sporadic 2003-2004 (`04`) outbreaks, and indicating "Early" (`E`), "Middle" (`M`), or "Late" (`L`) epidemic phase for sequences within the primary 2002-2003 epidemic.

I concatenated this expanded SARS-CoV-1 set with the SARS-related CoV sequences from before. This yields the final set of 51 RBD nt sequences included in the `./unaligned-sequences/RBD_nt.fasta` file. I also translated each nt sequence to its amino acid sequence, which is saved in the file `.unaligned-sequences/RBD_aa.fasta`.

Finally, we may want to try using non-sarbecovirus outgroup sequences. I initially tried including adding the non-sarbecovirus sequence from the Hibecovirus strain Bat Hp-BCoV Zhejiang 2013 (GenBank KF636752) to use as an outgroup for the sarbecovirus clade. This is the closest outgroup sequence I've seen in larger phylogenies of betaCoVs.

We could also try more distant betaCoV outgroups: Rousettus CoV HKU9 (Genbank HKU9); Bat CoV GCCDC1 (MT350598); OC43 (KX344031); HKU1 (KF686346); MERS-CoV (NC_019843). A recent [preprint](https://www.biorxiv.org/content/10.1101/2020.07.07.190546v1) illustrates that the positioning of the root may be very important (and their RBD tree differs in other respects from what I've seen in other papers.) So, instead of doing just at the level of RBD, we may want to assess robustness of rooting by comparing both spike and RBD. These do not align well when using just RBD -- we may want to extend our consideration of rooting and alignment to be the whole spike before parsing back down to RBD to better identify boundaries, particularly in divergent outgroups.

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

In either of these full spike phylogenies, rooting on these outgroups puts the Europe/Africa sequences as first to diverge. Let's see if this holds with the RBD only trees?

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

These trees also put the first split between the non-Asian and SE-Asian sarbecoviruses, as I had been operating under! So, that's good to see (and differs from the recent paper on bioRxiv linked above). The outgroup branches though are *incredibly* long, and it would probably look better if we pared back to just the Hibecovirus outgroup? Therefore, I add only the Hibecovirus sequence into our "main" unaligned sequences, to be used below.


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

I next used `RAxML` to infer the ML tree. We want to use a codon-aware nt substitution model -- there is not enough amino acid diversity for some of the fine-grained resolution (i.e. the SARS-CoV-1 sequences and their nearest bat relatives), so we want to use a nt-aware model which gives us the extra information available in synonymous codon positions to gain additional resolution into relationships. So, we fit GTRGAMMA nt models with separate data partitions (i.e. separate model parameters) fit to the first, second, and third codon positions -- this is how the RAxML manual suggests inference on codon-aligned protein-coding sequences.


```
mkdir RBD_codon_tree
cd ./RBD_codon_tree
nohup raxmlHPC-PTHREADS -s ../RBD_nt-codon_aligned.fasta -n RBD_codon_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &
```

Open up the tree `./RBD_codon_tree/RAxML_bestTree.RBD_codon_tree.txt` in FigTree, root on the Hp-BCoV\_Zhejiang outgrpus sequence, output tree as `RAxML_besttree_RBD_codon_rooted.txt` in Newick format. This will polarize the direction of nodes for ASR (not necessary, but some of the annotations that come out from ASR software are easier to interpret if polarized in the proper direction).

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
nohup raxmlHPC-PTHREADS -s ./RBD_nt-codon_aligned_rec1.fasta -n RBD_codon_tree_rec1.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 4 -q ../codon_partitions.txt -g constraint.txt &

nohup raxmlHPC-PTHREADS -s ./RBD_nt-codon_aligned_rec2.fasta -n RBD_codon_tree_rec2.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 4 -q ../codon_partitions.txt -g constraint.txt &> nohup2.out &

```

Looking at the trees built from these two different sections -- there are differences from the full-length ML tree, but any differences have incredibly low bootstrap support. I do *not* take this as indication of strong support for a specific recombination history here, but rather that we just nuked our phylogenetic signal by taking a small alignment and cutting it half. However, for completeness sake, we might want to do ASR on a recombination-aware tree. However, I want to try the GARD method as well for breakpoint detection. Instead of 3SEQ, which I believe requires something resembling the parental sequences to be present, the GARD approach based on determining likelihood support for alternate topologies, seems more straighforward. While I'm here I'll run the FastML, but I will likely roll with the GARD version of this below.

Output "bestTree" as rooted by outgroup hibecovirus sequence as `RAxML_besttree_RBD_codon_rec1_rooted.txt` in Newick format. This will polarize the direction of nodes for ASR.

Infer ancestors:

```
mkdir ./ASR_rec1
cd ./ASR_rec1
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_rec1.fasta --seqType AA --Tree ../RAxML_besttree_RBD_codon_rec1_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_nt-codon_recomb-test/ASR_rec1 &

mkdir ../ASR_rec2
cd ../ASR_rec2
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_rec2.fasta --seqType AA --Tree ../RAxML_besttree_RBD_codon_rec2_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_nt-codon_recomb-test/ASR_rec2 &
```


## Recombination analysis with `GARD`
The recombination analysis of 3SEQ gave results that were hard to interpret in terms of exact breakpoints and their reliability -- my by-eye checking of the breakpoints was less than convinced, though of course that's the point of an algorithm. Let's try an alternative algorithm in GARD using it's "Datamonkey" server. Ran the nt alignment through this algorithm, and saved files in `./RBD_nt-codon_recomb-test_GARD`. This analysis supports with small AIC 85.16 a breakpoint between positions 264 and 265. So, let's do the same as above, but with this single breakpoint, breaking the alignment into sub-segments: `RBD_nt-codon_aligned_GARD1.fasta` contains nucleotides 1-264, and `RBD_nt-codon_aligned_GARD2.fasta` contains nucleotides 265-606. Also translated the corresponding aa alignments for ASR.

I then inferred phylogenies separate for these two putative recombinant regions:

```
mkdir RBD_codon_GARD1_tree
cd RBD_codon_GARD1_tree
nohup raxmlHPC-PTHREADS -s ../RBD_nt-codon_aligned_GARD1.fasta -n RBD_codon_tree_GARD1.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 6 -q ../../codon_partitions.txt -g ../constraint.txt &

mkdir ../RBD_codon_GARD2_tree
cd ../RBD_codon_GARD2_tree
nohup raxmlHPC-PTHREADS -s ../RBD_nt-codon_aligned_GARD2.fasta -n RBD_codon_tree_GARD2.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 6 -q ../../codon_partitions.txt -g ../constraint.txt &
```

Same as with 3seq, the changes that are made don't really make sense -- the RsSHC014-sub-clade of SARS1-like RBDs gets moved outside as an entirely new ingroup Asia sarbecovirus clade, and GX-Pangolin becomes sister to clade to deletions/clade 2. Neither has that strong of bootstrap support. The RsSHC014 movement in particular, I feel it is likely hard to disentagle recombination from the fact that the "GARD2" segment corresponds to teh RBM, which is under elevated evolution in these Rs bats for both receptor binding evolution, and likely also antigenic evolution. Therefore, it is possible elevated evolutionary rate could end up looking similar to tree disconcordance between these two regions? The only other defining feature is that the SARS1 clade has a clade-specific 1-aa deletion relative to SARS2 clade in the "receptor binding ridge" -- and the RsSHC014 sequecnes share the exact same deletion, and other similar features within this ridge (e.g. prolines). I do therefore believe these are still "SARS1-clade" sequences, but we will see how AncSarbecovirus and AncAsia reconstructions chagne, anyway, because this will of course be a major concern of not only mine but also reviewers.

Output "bestTree" as rooted by outgroup hibecovirus sequence as `RAxML_besttree_RBD_codon_GARD1_rooted.txt` in Newick format. This will polarize the direction of nodes for ASR.

Infer ancestors:

```
mkdir ../ASR_GARD1
cd ../ASR_GARD1
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_GARD1.fasta --seqType AA --Tree ../RBD_codon_GARD1_tree/RAxML_besttree_RBD_codon_GARD1_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_nt-codon_recomb-test_GARD/ASR_GARD1 &

mkdir ../ASR_GARD2
cd ../ASR_GARD2
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_GARD2.fasta --seqType AA --Tree ../RBD_codon_GARD2_tree/RAxML_besttree_RBD_codon_GARD2_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/RBD_nt-codon_recomb-test_GARD/ASR_GARD2 &
```

Parsed sequences with the `.Rmd` script in this directory, and the `ancestors_input.csv` table which lists the matched nodes for concatenating ancestors.

How do these differ from the v1 non-recombining phylogeny?

AncSarbecovirus: differs at four positions in the GARD-segmented reconstruction: S346T, delS372a, I434L, E484S

AncAsia differs at more, consistent with the RsSHC014 clade breaking out and making things more "SARS1-like" in the receptor-binding ridge, even though i think that my by-eye parsimony (e.g. the 482 deletion) says RsSHC014 shoudl still be in the SARS1 clade. Differences include: I434L, N440S, S445T, S459G, F464Y, E471V, P479S, delG482, V483A, E484V, F486L, Y490N, Q493K, H498Y, P499T, T501S

AncSARS2a: S443A, V483Q, E484V, F486L, H498Y


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
 
## Addition of new sarbecovirus RBD sequences
As of 7 March, 2021, there have been new sarbecovirus RBD sequences reported, originally isolated from Japan, Cambodia, Thailand, Rwanda, adn Uganda. I want to incorporate these sequences into my alignment and tree, and also see if their inclusion alters the sequecnes of some key ancestors.
Update 23 April 2021, interesting new sequence just added to GISAID from Zhou et al. preprint, RsYN04, that might be a unique phylogenetic position according to Stephen Goldstein's virological post, also want to add this and perhaps functionally characterize.

Started new subdirectory, copy in working RBD alignments
```
mkdir add_new_RBDs
cd add_new_RBDs
cp ../RBD_aa_aligned.fasta ./
cp ../unaligned-sequences/RBDs_nt.fasta ./RBDs_nt_unaligned_v2.fasta
```

Downloaded new sequences, in file `RBDs_new_aa.fasta`. Added nt sequences to `RBDs_nt_unaligned_v2.fasta`: 
  - Rc-o319 (Japan, R. cornutus): Genbank LC556375
  - RacCS203 (Thailand, R. acuminatus): Genbank MW251308
  - RshSTT182 (Cambodia, R. shameli): GISAID EPI\_ISL\_852604
  - PDF-2370 (Uganda, R. spp [close to ferrumequinum], same RBD seq as PDF-2386): Genbank accession not yet available, got from supplement from Wells et al. preprint
  - PRD-0038 (Rwanda, R. clivosus): Genbank not yet available, got from Wells et al. supplement
  - RsYN04 (Yunnan, R. stheno): GISAID EPI\_ISL\_1699444

Use mafft to align new sequences in with the prior set:

```
mafft --add ./RBDs_new_aa.fasta --reorder ./RBD_aa_aligned.fasta > ./RBD_aa_aligned_v2.fasta
```

Nonparismonious double gap in the Japan sequence relative to the single-aa deletion in the SARS1 clade, I manually updated this to make it a single deletion instead of duplicate. 

Make nt alignment:

```
pal2nal.pl ./RBD_aa_aligned_v2.fasta RBDs_nt_unaligned_v2.fasta > ./RBD_nt_aligned_v2.clustal
```
Convert clustal to fasta

Infer phylogeny in RAxML:

```
mkdir RBD_codon_tree_v2
cd ./RBD_codon_tree_v2
nohup raxmlHPC-PTHREADS -s ../RBD_nt_aligned_v2.fasta -n RBD_codon_tree_v2.txt -m GTRGAMMA -f a -p 20 -N autoMRE -x 20 -T 4 -q ../codon_partitions.txt &
```

Open up the tree `./RBD_codon_tree_v2/RAxML_bestTree.RBD_codon_tree_v2.txt` in FigTree, root on the Hp-BCoV\_Zhejiang outgroup sequence, output tree as `RAxML_besttree_RBD_codon_v2_rooted.txt` in Newick format. This will polarize the direction of nodes for ASR.

Infer ancestors:

```
mkdir ../ASR_v2
cd ../ASR_v2
nohup perl ../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../RBD_aa_aligned_v2.fasta --seqType AA --Tree ../RBD_codon_tree_v2/RAxML_besttree_RBD_codon_v2_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/add_new_RBDs/ASR_v2 &
```
Parsed sequences with the `.Rmd` script in this directory. Compared the "v1" sequences to these:
  - AncSarbecovirus v2 reconstruciton has mutations I434L, F452Y, K490E, S501A compared to the v1. S501A probably slightly enhances the Ra.9479 affinity, others probably don't matter
  - AncAsia v2 has mutaitons N394S, I434L, L441Q, L452Y, G482S, Q493K, H498Y, T501A rleative to v1. From SSM data in v1 AncAsia, the three contact sites still compatible with binding (at least individually)
  - AncSARS2a v2 has mutations L441Q, S443A, K444S, G446S, L452Y, T470N, T478Q, N481S, G482S, Q493K, H498Y, T501A. Again, probably the indiviudal mutations compatible with binding, but there are many, so ti might be worth testing. (Might show binding for e.g. cvACE2 binding whihc the v1 didn't)
  - Note that the v2 AncSARS2b is closer to the v1 AncSARS1a (only differs K403R, L441Q, S443A, K444S), consistent with those having the same ingroup sequecnes (plus cambodia for v2)


## Addition of _even more_ new sarbecovirus RBD sequences (it never ends, of course)
Two new sarb sequecnes from Black Sea/Caucasus region in Russia recently reported in [this preprint](https://www.biorxiv.org/content/10.1101/2021.05.17.444362v2), dubbed Khosta-1 and Khosta-2. They also pointed out to me one more clade 3 sequence related to BM48-31 that is on [Genbank](https://www.ncbi.nlm.nih.gov/nuccore/939726458) though not published. And a [recent paper](https://www.nature.com/articles/s41598-021-94011-z) published genome of RhGB01 from Great Britain, in an R. hipposideros bat

Started new subdirectory, copy in working RBD alignments
```
mkdir add_new_RBDs/add_more_new_RBDs
cd add_new_RBDs/add_more_new_RBDs
cp ../RBDs_nt_unaligned_v2.fasta ./RBDs_nt_unaligned_v3.fasta
```

Downloaded new sequences, merged with v2 alignment (removing gaps to unalign everything) in file `RBDs_aa_unaligned_v3.fasta`. Added nt sequences to `RBDs_nt_unaligned_v3.fasta`: 
  - BB9904 (Bulgaria, R. euryale): Genbank KR559017
  - Khosta-1 (Russia, R. ferrumequinum): Genbank MZ190137
  - Khosta-2 (Russia, R. hipposideros): Genbank MZ190138
  - RhGB01 (Great Britain, R. hipposideros): Genbank MW719567
  
Use mafft to align updated amino acid sequence set 'fresh' (having undone prior alignment)

```
mafft --reorder --op 4.5 ./RBDs_aa_unaligned_v3.fasta > ./RBDs_aa_aligned_v3.fasta
```

Make nt alignment:

```
pal2nal.pl RBDs_aa_aligned_v3.fasta ./RBDs_nt_unaligned_v3.fasta > RBDs_nt_aligned_v3.clustal
```
Convert clustal to fasta, save as `RBDs_nt_aligned_v3.fasta`


Infer ML tree in RAxML based on nt and aa seqs

```
mkdir RAxML_nt_tree_v3
cd ./RAxML_nt_tree_v3
nohup raxmlHPC-PTHREADS -s ../RBDs_nt_aligned_v3.fasta -n RAxML_nt_tree_v3.txt -m GTRGAMMA -f a -p 20 -N autoMRE -x 20 -T 8 -q ../codon_partitions.txt &
```

RAxML, aa seqs:
```
mkdir ../RAxML_aa_tree_v3
cd ../RAxML_aa_tree_v3
nohup raxmlHPC-PTHREADS -s ../RBDs_aa_aligned_v3.fasta -n RAxML_aa_tree_v3.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10 -T 8 &
```

And also try in IQtree with a true codon substitution model, let IQtree find best fit model

```
mkdir ../IQtree_codon_tree_v3
cd ../IQtree_codon_tree_v3
nohup iqtree -s ../RBDs_nt_aligned_v3.fasta -st CODON -T AUTO -ntmax 8 -B 1000 -alrt 1000 --prefix IQtree_codon_tree_v3 &
```

Some variability between the three trees. All three have the root breaking up the Eur/Af sequences as paraphyletic, but first to diverge separate from Asian sarbecoviruses. The RAxML nt tree bascially pulls Khosta-2 out as a zero-branch-length extension from teh sarb root (aka, it is ancSarb), followed by RhGB01, then the remaining Af/Eur forming a clade. RsYn04 than branches from teh remaining three Asia clades as I saw previously. The RAxML aa tree pulls RhGB01 out as the zero-branch length "root" sequence, while Khosta2 then falls in with the clade of remaining Eur/Af. RsYN04 is seen to branch internal to the other three Asia clades. This is similar to what is found in the codon tree, with the only exception that RhGB01 gets some branch length away from 'the root' in its breakup of the Eur/Af monophyly.

This does still seem to support more or less our original rooting of Eur/Af / Asia divergence, except we're getting more Eur/Af paraphyly. RsYN04 position appears to differ based on nt versus codon/aa trees, though its placement is not determined with super strong node support in either scenario, though in eithe rcase it is 'it's own' group independent of other clades.

In thec ase of the RAxML aa and nt trees, ancSarb will just be inferred to be the same as either Khosta2 or RhGB01. Would be worth doing ASR in IQtree though on the codon tree, sicne there is 'branch length' to RhGB01 from the root meaning there will be sequence differences. I believe I can do this ASR *in* IQtree itself using this same codon sub model. We will try taht below.

Other than that, I wonder whether it's worth looking at constraining Eur/Af monophyly as an alternative ancestor. (Not saying I believe strongly enough we need to impose this monophyly, but rather that the ancSarb=khosta2 or RhGB01 means it's just those phenotypes. We could additionally test this alternative constraint tree ancSarb for more understanding of the robustness of its RaACE2 binding observation.) I could infer full spike trees and see if that gives any more signal for monophyletic Eur/Af sarbs?


Want to do whole spike tree with more robust outgroup inclusion to understand rooting in context of these newly added sequences.

```
mkdir ../spike
cd ../spike
```
Copied in unaligned spike sequences from `RBD_ASR/outgroup-rooting/unaligned-sequences`
Added all of the new spikes listed above, saved as `spikes_aa_unaligned_v3.fasta` and `spikes_nt_unaligned_v3.fasta`

Align spike aa seqs
```
mafft --reorder --op 4.5 ./spikes_aa_unaligned_v3.fasta > ./spikes_aa_aligned_v3.fasta
```

Also make the nt alignment from amino acid sequences using PAL2NAL, and convert to fasta with online tool linked above:

```
pal2nal.pl ./spikes_aa_aligned_v3.fasta ./spikes_nt_unaligned_v3.fasta > ./spikes_nt_aligned_v3.clustal
```

Manually cleaned up alignments, deleted regions that are gaps in sarbecovirues but insertions in outgroups. 

Saved as: `./spikes_aa_aligned_v3_cleaned.fasta` and `./spikes_nt_aligned_v3_cleaned.fasta`

Infer RAxML tree using partitioned nt model

```
mkdir ./RAxML_nt_tree
cd ./RAxML_nt_tree
nohup raxmlHPC-PTHREADS -s ../spikes_nt_aligned_v3_cleaned.fasta -n RAxML_spike_nt_tree.txt -m GTRGAMMA -f a -p 10 -N autoMRE -x 10 -T 8 -q ../codon_partitions.txt &
```

Infer RAxML tree using aa model

```
mkdir ../RAxML_aa_tree
cd ./RAxML_aa_tree
nohup raxmlHPC-PTHREADS -s ../spikes_aa_aligned_v3_cleaned.fasta -n RAxML_spike_aa_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10 -T 8 &
```

Infer IQtree tree using best-fit codon substitution model
```
mkdir ../iqtree_codon_tree
cd ../iqtree_codon_tree
nohup iqtree -s ../spikes_nt_aligned_v3_cleaned.fasta -st CODON -T AUTO -ntmax 8 -B 1000 -alrt 1000 --prefix IQtree_spike_codon_tree &

```

Interpretation:
 - Spike nt tree gives same pattern as the RBD tree: rooting induces paraphyly of Eur/Af sequences, with Khosta-2 divergence followed by RhGB01 and then monophyly of remaining Eur/Af seqs. RsYN04 placement is where I originally saw it, prior to dviergence of three other asia clades. The node support for RsYN04 placement is decent (node support 84 for monophyly of other three clades of Asia sarb), node support for paraphyly of Eur/Af modest at best (would break node supports 52 and 58 to make Eur/Af monophyletic)
 - Spike aa tree gives Eur/Af monophyly (node support decent with 82) as well as RsYN04 in the prior position I've seen (would only have to break node support 40 for clade 1a/1b/2 monophyly)
 - IQtree codon tree gives 100 support Eur/Af monophyly. It does have RsYN04 coming out as the first branch from the root, htough, which is the first time I've seen this. That has ok but not perfect node suppport.

So, my conclusions are:
 - all but IQtree, spike tree continue to support initial splits being between Eur/Af and Asia
 - 4 of the six spike+RBD trees induce Af/Eur paraphyly, though ordinarily with weak support. Two of the trees recover monophylyetic Eur/Af
 - The paraphyletic spike trees (excpet IQtree codon, RBD) have zero branch from root to a sequence I'm already ordering, implying the ancestor would just *be* that sequence. I could therefore address those putative ancestral phenotypes based on these seuqecnes that I'm already ordering, even if I don't truly believe ancestor had this same sequence and this is likely phylogentic artefact
 - It could be worth trying an IQtree ASR approach on the codon RBD tree it recovers with RhGB01 branching (wiht branch length) from ancSarb
 - It could be worth inferring the RBD tree constraining monophyletic Eur/Af as I still think that's a reasonable 'true' topology that is sometimes recovered despite the long outgroup branch

RAxML RBD, partitioned nt model with constraint for Eur/Af monophyly
```
mkdir ../RAxML_nt_tree_v3_constraint
cd ../RAxML_nt_tree_v3_constraint
nohup raxmlHPC-PTHREADS -s ../RBDs_nt_aligned_v3.fasta -n RAxML_nt_tree_v3_constraint.txt -m GTRGAMMA -f a -p 20 -N autoMRE -g ../constraint.txt -x 20 -T 8 -q ../codon_partitions.txt &
```

RAxML, aa seqs with constraint for Eur/Af monophyly:
```
mkdir ../RAxML_aa_tree_v3_constraint
cd ../RAxML_aa_tree_v3_constraint
nohup raxmlHPC-PTHREADS -s ../RBDs_aa_aligned_v3.fasta -n RAxML_aa_tree_v3_constraint.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -g ../constraint.txt -x 10 -T 8 &
```

Last, do ancestral sequence reconstruction. Will do it (1) on the constrained RBD nt tree topology, (2) using IQtree (no indel reconstruction...) with its ASR function based on its codon model and topology, and (3) on the unconstrained RBD nt tree topology (to see what results are)


Open up the following trees in FigTree, rooted, adn outputted tree with `_rooted` nomenclature as Newick tree for input to ASR: `./RAxML_nt_tree_v3/RAxML_bestTree.RAxML_nt_tree_v3.txt`, `./RAxML_nt_tree_v3_constraint/RAxML_bestTree.RAxML_nt_tree_v3_constraint.txt` and `./IQtree_codon_tree_v3/IQtree_codon_tree_v3.treefile`

Infer ancestors -- do in an `ASR` subdirectory within each of the respective phylogeny's subdirectory.

RAxML nt tree, unconstrained
```
nohup perl ../../../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../../RBDs_aa_aligned_v3.fasta --seqType AA --Tree ../RAxML_bestTree.RAxML_nt_tree_v3_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/add_new_RBDs/add_more_new_RBDs/RAxML_nt_tree_v3/ASR &
```
RAxML, constrained nt tree
```
nohup perl ../../../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../../RBDs_aa_aligned_v3.fasta --seqType AA --Tree ../RAxML_bestTree.RAxML_nt_tree_v3_constraint_rooted.txt --SubMatrix LG --OptimizeBL yes --jointReconstruction yes --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/add_new_RBDs/add_more_new_RBDs/RAxML_nt_tree_v3_constraint/ASR &
```
IQtree, codon tree (using best fit model used in the ML tree inference)
```
nohup iqtree -s ../../RBDs_nt_aligned_v3.fasta -te ../IQtree_codon_tree_v3_rooted.txt -m KOSI07+F3X4+R5 -st CODON -asr -T AUTO -ntmax 8 --prefix IQtree_codon_tree_v3_ASR &
```
Also run FastML here to infer indels? Maybe without optimizing branch length since the indel model is independent of amino acid/codon substitution modeling?

```
mkdir fastml
cd fastml
nohup perl ../../../../../../programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_file ../../../RBDs_aa_aligned_v3.fasta --seqType AA --Tree ../../IQtree_codon_tree_v3_rooted.txt --SubMatrix LG --OptimizeBL no --jointReconstruction no --indelReconstruction BOTH --outDir /fh/fast/bloom_j/computational_notebooks/tstarr/2020/SARSr-CoV_homolog_survey/RBD_ASR/add_new_RBDs/add_more_new_RBDs/IQtree_codon_tree_v3/ASR/fastml &
```

Parsed sequences with the `.Rmd` script in each ASR subdirectory.
  - AncSarbecovirus1 in the v3, nt tree, no constraint , as I expected is inferred to be identical to Khosta-2. This is because the root is put with 0 branch length to Khosta2 in this tree, as I noted above. Therefore, nothing to order here.
  - AncSarbecovirus in the v3 nt tree with constrained monophyly: how does it differ from original AncSarbecovirus?: "Q340E;T346S;D360N;A372aS;G413A;L434I;R439N;I441L;K444S;Q445S;G446-;-448G;-449N;-450N;Y452F;A475P;-482S;P483I;S484E;E490K;V501S" -- some of this is difference in placement of gaps with new sequences, together with new signal. This ancestor seems worth ordering and testing, called AncSarbecovirus\_v3
  - AncAsia2 in the v3 nt tree with constrained monophyly: how does it differ from original AncAsia? "Q340E;R403K;E406D;L434I;R439N;S440N;Q441L;A443S;Q445S;Y452L;K458H;V471E;A475P;-482G;P483V;S484E;L486F;K493Q;T498H;V501T" -- I do expect the T498 to be less of a strong (human and perhaps other) ACE2 binding state than H, and same perhaps with V501. Should probably order this along with AncAsia1 (ancestor of AncAsia1+RsYN04 -- only major difference is it's inferred to have the "region 1" deletion), as well as "AncSARS12" (ancestor of SARS1 and SARS2 clade, in case human ACE2 origin is pushed down a node in this AncAsia2).
  - IQtree reconstruction -- doesn't allow for reconstruction of the AncSarbecovirus sequecne, not sure why. This will take longer term effort to figure out how to do in codon sub model context, given lack of indel reconstruction and this issue with AncSarb. I am going to order the new extant sarbs which is the most important bit anyway, along with the "v3" ancestor above.
  


