# Studying Non-Canonical G-Quadruplex repertoire in Tenacibaculum: A Case Study

PhD St. Ricardo Gomez-Reyes

### Introduction

Nucleic acid sequences rich in guanine are capable of forming four-stranded structures called G-quadruplexes, stabilized by Hoogsteen hydrogen bonding between a tetrad of guanine bases (also called Potential G-quadruplex forming sequences, PQSs). The formation of these telomeric quadruplexes has been shown to decrease the activity of the enzyme telomerase ( [8](javascript:;) ), which is responsible for elongating telomeres. Since elevated telomerase activity has been implicated in ∼85% of cancers ( [9](javascript:;) ), this has become a significant strategy for drug development ( [10](javascript:;) ) and molecules that bind to and stabilize G-quadruplexes have been identified. Recently, there has been growing interest in quadruplex-forming sequences elsewhere in the genome (Julian L. Huppert et al, 2005). 

![Figure 1](G4PromFinder_outputs/G-tetrad.jpeg)



In the genome, there exist many instances in which the PQSs possess more than four G tracks that can result in these sequences adopting dynamic structures equilibrating between multiple folds. Te ability of PQSs to fold intracellularly was first observed by immunostaining of ciliate telomeres and this study was followed by immunofuorescence of human cells to fnd folded G4s18. Te folding of PQSs to G4s has been implicated in causing strand breaks during replication in the absence of faithful helicases to resolve these roadblocks to polymerase bypass; when located in gene promoters, G4s can regulate transcription of the gene; G4s may be important in telomere biology; and they may function at origins of replication in humans. Experiments have found folded G4s are not limited to the genome as they can fold in the transcriptome ().



 Over the last few decades, computational genomics has tremendously contributed to decipher biology from genome sequences and related data. Considerable effort has been devoted to the prediction of transcription promoter and terminator sites that represent the essential "punctuation marks" for DNA transcription 



### Started

Inputs

```bash
srohit@kneipe.lavis.unam.mx/home/srohit/unamworkshop2021/raw_data
```

To do:

1. Randomize Genome Simulations for 26 genomes (2-kmer)

2. Test **G4PromFinder** for both, s

Cite: https://github.com/MarcoDiSalvo90/G4PromFinder

G4PromFinder is an algorithm for the promoter prediction in bacterial genomes. It is recommended for GC rich genomes. G4PromFinder predicts putative promoters based on AT-rich elements and G-quadruplex DNA motifs.

Input: Genome sequences file (Input file must be in fasta format)

Output: a text file containing promoter coordinates and a file containing informations about them.

Let's start

```bash
ssh srohit@kneipe.lavis.unam.mx
# zPjWZBN6aE
```

Modify the G4PromFinder.py in order to make a loop genome-to-genome (script [here](https://github.com/RJEGR/2nd-Workshop-in-Advanced-Bioinformatics/blob/main/G4PromFinder_outputs/HMG4PromFinder.py));

Then:

```bash
for i in $(ls *fa); do  python HMG4PromFinder.py $i; done
#nohup for i in $(ls *fa); do  python HMG4PromFinder.py $i; done &
```

Copy results

```bash
scp -r srohit@kneipe.lavis.unam.mx:/home/srohit/unamworkshop2021/G4PromFinder_outputs .
```

## Test Negative control

**How To Generate Randomized Sequence Based On Sequence Already Known?**

> Example: https://www.biostars.org/p/69756/

```python
import random
dna = list('ATATTCATGAGTACCGTA'); 
print('Original:',''.join(dna))
random.shuffle(dna); 
print('Random:  ',''.join(dna))
```



Instead of python, lets to use the home-made script [here](https://github.com/RJEGR/2nd-Workshop-in-Advanced-Bioinformatics/blob/main/shuffled_genomes.R) (based on https://stat.ethz.ch/pipermail/bioconductor/2013-March/051640.html)



Based on literature (), lets shuffle the genomes based on a sliding window (2-kmer) using `biasaway k -h`  (note: Quadruplexes can be uni-, bi- or tetramolecular (Julian L. Huppert et al, 2005).

```bash
genome=genome.fa

biasaway k --foreground $genome --nfold 1 --kmer 2 --seed 202102 > ${genome%.fa}.random.fn

for genome in $(ls *fa); do biasaway k --foreground $genome --nfold 1 --kmer 2 --seed 20210210 > ${genome%.fa}.random.fn; done
```

Then 

```bash
for i in $(ls kmer_2/*fn); do  python HMG4PromFinder.py $i; done
```



### Calculate GC content per genome

```bash
awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f\n", (gc*100)/(gc+at) }' Tenacibaculum_todarodis_gca_001889045.fa

# for i in $(ls *fa); do  awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f\n", (gc*100)/(gc+at) }' $i; done
```

### Calculate genome length

```bash
awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f\n", gc+at }' Tenacibaculum_todarodis_gca_001889045.fa

# for i in $(ls *fa); do  awk '!/^>/{gc+=gsub(/[gGcC]/,""); at+=gsub(/[aAtT]/,"");} END{ printf "%.2f\n", gc+at }' $i; done
```

| Index                                            | GC    |  Size   |
| ------------------------------------------------ | ----- | :-----: |
| Tenacibaculum_caenipelagi_gca_004363005          | 31.88 | 3266097 |
| Tenacibaculum_dicentrarchi_gca_001483385         | 31.48 | 2918253 |
| Tenacibaculum_dicentrarchi_gca_900239305         | 30.18 | 2663544 |
| Tenacibaculum_dicentrarchi_gca_900239345         | 30.09 | 2804033 |
| Tenacibaculum_discolor_gca_003664185             | 31.63 | 3376020 |
| Tenacibaculum_finnmarkense_gca_900239185         | 30.91 | 2923232 |
| Tenacibaculum_finnmarkense_gca_900239485         | 31.01 | 2955007 |
| Tenacibaculum_gallaicum_gca_003387615            | 31.74 | 3422893 |
| Tenacibaculum_holothuriorum_gca_002120225        | 31    | 3147654 |
| Tenacibaculum_jejuense_gca_900198195             | 30.31 | 4614879 |
| Tenacibaculum_litoreum_gca_003937815             | 31.92 | 3376591 |
| Tenacibaculum_lutimaris_gca_003610735            | 31.9  | 2868850 |
| Tenacibaculum_maritimum_ncimb_2154_gca_900119795 | 32.01 | 3435971 |
| Tenacibaculum_mesophilum_gca_900129475           | 31.61 | 3286619 |
| Tenacibaculum_skagerrakense_gca_004345825        | 31.25 | 3712612 |
| Tenacibaculum_soleae_gca_001693415               | 30.25 | 3006944 |
| Tenacibaculum_sp_4g03_gca_002700005              | 31.8  | 3812359 |
| Tenacibaculum_sp_bg11_29_gca_002836595           | 29.77 | 4537492 |
| Tenacibaculum_sp_dsm_106434_gca_003867015        | 32.01 | 3511704 |
| Tenacibaculum_sp_m341_gca_004337695              | 30.09 | 5123891 |
| Tenacibaculum_sp_mar_2009_124_gca_900105455      | 32.4  | 5524492 |
| Tenacibaculum_sp_mar_2010_89_gca_900105985       | 29.32 | 4172912 |
| Tenacibaculum_sp_sg_28_gca_002954385             | 34.37 | 2761692 |
| Tenacibaculum_sp_sz_18_gca_002813915             | 31.14 | 4019179 |
| Tenacibaculum_sp_tno020_gca_900239505            | 30.71 | 2452834 |
| Tenacibaculum_todarodis_gca_001889045            | 30.73 | 3019213 |

### Data-Viz

Preliminary results

![Figure 1](G4PromFinder_outputs/G4PromFinder_kmer2.png)





![Figure 1](G4PromFinder_outputs/S1_figure.png)

## Cites

Aziz Khan, Rafael Riudavets Puig, Paul Boddie, Anthony Mathelier, BiasAway: command-line and web server to generate nucleotide composition-matched DNA background sequences, *Bioinformatics*, 2020;, btaa928, https://doi.org/10.1093/bioinformatics/btaa928



Ding, Y., Fleming, A. M., & Burrows, C. J. (2018). Case studies on potential G-quadruplex-forming sequences from the bacterial orders Deinococcales and Thermales derived from a survey of published genomes. *Scientific reports*, *8*(1), 1-11. DOI:10.1038/s41598-018-33944-4

Di Salvo, M., Pinatel, E., Talà, A., Fondi, M., Peano, C., & Alifano, P. (2018). G4PromFinder: an algorithm for predicting transcription promoters in GC-rich bacterial genomes based on AT-rich elements and G-quadruplex motifs. *BMC bioinformatics*, *19*(1), 1-11. https://doi.org/10.1186/s12859-018-2049-x



Julian L. Huppert, Shankar Balasubramanian, Prevalence of quadruplexes in the human genome, *Nucleic Acids Research*, Volume 33, Issue 9, 1 May 2005, Pages 2908–2916, https://doi.org/10.1093/nar/gki609

