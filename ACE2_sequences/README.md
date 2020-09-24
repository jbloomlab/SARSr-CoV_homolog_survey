# Alignment and selection of ACE2 molecules for functional study

This README details the series of bioinformatic analyses that were conducted to identify ACE2 proteins for functional study. Many steps are only pseudo-automated, so this analysis is not amenable to streamlining in a pipeline.

## Sequence set

We want to align a starting set of ACE2 sequecnes of likely interest. These include any Rhinolophus bat ACE2 sequences that we can find, in addition to human ACE2, and important species including civet, pangolin, or mouse.

We first recieved sequences of _R. sinicus_ and _R. affinis_ ACE2s described in this [preprint](https://www.biorxiv.org/content/10.1101/2020.05.13.093658v1), generously shared by Zheng-Li Shi from Wuhan Institute of Virology. This contains newly described Ra and Rs sequences from that study, in addition to 4 Rs sequences that were previously available from Genbank. We translated these ACE2 sequences to their amino acid sequences. I also added more annotations to each fasta header sequence, including provice of origin (Figure 1) if known, and Genbank number of existing.

We also downloaded all Rhinolophus ACE2 sequences available on Genbank, which added 6 additional species. Most of these are from SE Asia (Hou et al. and other previous papers from Zheng-Li Shi), with two species whose range is AFrica (R. landeri and R. alcyone), along with R. ferrumequinum, which spans from Europe to central and SE Asia (has had SARSr-CoV sampled in some of the big sampling papers).

I also added ACE2 sequences from human, mouse, civet, racoon dog, and pangolin. These are all concatenated into the file `./unaligned-sequences/ACE2_aa.fasta`.

## Alignment

I used `mafft` to align ACE2 aa sequences:

```
mafft --reorder ./unaligned-sequences/ACE2_aa.fasta > ./ACE2_aa_aligned.fasta
```

## Phylogenetic inference

Want to infer the phylogenetic relationship among these ACE2 sequences, which may be helpful in visualizing diversity for selecting sequences for protein purification.

```
mkdir ACE2_tree
cd ./ACE2_tree
nohup raxmlHPC-PTHREADS -s ../ACE2_aa_aligned.fasta -n ACE2_tree.txt -m PROTGAMMALG -f a -p 10 -N autoMRE -x 10 -T 8 &
```

