[![DOI](https://zenodo.org/badge/607588720.svg)](https://doi.org/10.5281/zenodo.14012589)

# The reference-free *Arabidopsis thaliana* pangenome

This repository collects all the commands and script used for the reference-free *Arabidopsis thaliana* pangenome. This work is part of the PhD project carried out by Lia Obinu.

If you use the code and methods explain in this repository, please cite our pre-print:

> The reference-free pangenome of Arabidopsis thaliana \
Lia Obinu, Andrea Guarracino, Urmi Trivedi, Andrea Porceddu \
bioRxiv 2025.06.25.661508; doi: https://doi.org/10.1101/2025.06.25.661508 


# Table of contents

- [The reference-free *Arabidopsis thaliana* pangenome](#the-reference-free-arabidopsis-thaliana-pangenome)
- [Table of contents](#table-of-contents)
- [Building the reference-free *Arabidopsis thaliana* pangenome graph](#building-the-reference-free-arabidopsis-thaliana-pangenome-graph)
  - [Get the assemblies](#get-the-assemblies)
  - [Quality check and trimming of the assemblies](#quality-check-and-trimming-of-the-assemblies)
  - [Renaming the assemblies](#renaming-the-assemblies)
  - [Format fasta headers and remove organelles from assemblies](#format-fasta-headers-and-remove-organelles-from-assemblies)
  - [Sequence Partitioning using partition-before-pggb](#sequence-partitioning-using-partition-before-pggb)
  - [Running pggb on each community](#running-pggb-on-each-community)
    - [community 0 - chr1](#community-0---chr1)
    - [community 1 - chr2](#community-1---chr2)
    - [community 2 - chr3](#community-2---chr3)
    - [community 3 - chr4](#community-3---chr4)
    - [community 4 - chr5](#community-4---chr5)
- [Pangenome annotation - reference based](#pangenome-annotation---reference-based)
  - [Get odgi untangle for annotation](#get-odgi-untangle-for-annotation)
  - [Prepare the annotation files for injection](#prepare-the-annotation-files-for-injection)
    - [Genes](#genes)
      - [community 0 - chr1](#community-0---chr1-1)
      - [community 1 - chr2](#community-1---chr2-1)
      - [community 2 - chr3](#community-2---chr3-1)
      - [community 3 - chr4](#community-3---chr4-1)
      - [community 4 - chr5](#community-4---chr5-1)
    - [Pseudogenes](#pseudogenes)
      - [community 0 - chr1](#community-0---chr1-2)
      - [community 1 - chr2](#community-1---chr2-2)
      - [community 2 - chr3](#community-2---chr3-2)
      - [community 3 - chr4](#community-3---chr4-2)
    - [community 4 - chr5](#community-4---chr5-2)
  - [Inject genes from the reference  into the pangenome graph](#inject-genes-from-the-reference--into-the-pangenome-graph)
    - [Genes](#genes-1)
      - [community 0 - chr1](#community-0---chr1-3)
      - [community 1 - chr2](#community-1---chr2-3)
      - [community 2 - chr3](#community-2---chr3-3)
      - [community 3 - chr4](#community-3---chr4-3)
      - [community 4 - chr5](#community-4---chr5-3)
    - [Pseudogenes](#pseudogenes-1)
      - [community 0 - chr1](#community-0---chr1-4)
      - [community 1 - chr2](#community-1---chr2-4)
      - [community 2 - chr3](#community-2---chr3-4)
      - [community 3 - chr4](#community-3---chr4-4)
      - [community 4 - chr5](#community-4---chr5-4)
  - [Untangle genes across all the assemblies of the pangenome](#untangle-genes-across-all-the-assemblies-of-the-pangenome)
    - [Genes](#genes-2)
      - [community 0 - chr1](#community-0---chr1-5)
      - [community 1 - chr2](#community-1---chr2-5)
      - [community 2 - chr3](#community-2---chr3-5)
      - [community 3 - chr4](#community-3---chr4-5)
      - [community 4 - chr5](#community-4---chr5-5)
    - [Pseudogenes](#pseudogenes-2)
      - [community 0 - chr1](#community-0---chr1-6)
      - [community 1 - chr2](#community-1---chr2-6)
      - [community 2 - chr3](#community-2---chr3-6)
      - [community 3 - chr4](#community-3---chr4-6)
      - [community 4 - chr5](#community-4---chr5-6)
  - [Untangle PAF processing](#untangle-paf-processing)
    - [collapse and merge](#collapse-and-merge)
      - [community 0 - chr1](#community-0---chr1-7)
      - [community 1 - chr2](#community-1---chr2-7)
      - [community 2 - chr3](#community-2---chr3-7)
      - [community 3 - chr4](#community-3---chr4-7)
      - [community 4 - chr5](#community-4---chr5-7)
    - [genes classification](#genes-classification)
      - [community 0 - chr1](#community-0---chr1-8)
      - [community 1 - chr2](#community-1---chr2-8)
      - [communty 2 - chr3](#communty-2---chr3)
      - [community 3 - chr4](#community-3---chr4-8)
      - [community 4 - chr5](#community-4---chr5-8)
    - [pseudogenes classification](#pseudogenes-classification)
      - [community 0 - chr1](#community-0---chr1-9)
      - [community 1 - chr2](#community-1---chr2-9)
      - [community 2 - chr3](#community-2---chr3-8)
      - [community 3 - chr4](#community-3---chr4-9)
      - [community 4 - chr5](#community-4---chr5-9)
- [Pangenome annotation - non-reference sequences](#pangenome-annotation---non-reference-sequences)
  - [Get the FASTA file for annotation](#get-the-fasta-file-for-annotation)
    - [community 0 - chr1](#community-0---chr1-10)
      - [1. Extract non-reference sequences](#1-extract-non-reference-sequences)
      - [2. Filter non\_reference.bed](#2-filter-non_referencebed)
      - [3. Apply bedtools slop and merge](#3-apply-bedtools-slop-and-merge)
      - [4. Extract the FASTA](#4-extract-the-fasta)
    - [community 1 - chr2](#community-1---chr2-10)
      - [1. Extract non-reference sequences](#1-extract-non-reference-sequences-1)
      - [2. Filter non\_reference.bed](#2-filter-non_referencebed-1)
      - [3. Apply bedtools slop and merge](#3-apply-bedtools-slop-and-merge-1)
      - [4. Extract the FASTA](#4-extract-the-fasta-1)
    - [community 2 - chr3](#community-2---chr3-9)
      - [1. Extract non-reference sequences](#1-extract-non-reference-sequences-2)
      - [2. Filter non\_reference.bed](#2-filter-non_referencebed-2)
      - [3. Apply bedtools slop and merge](#3-apply-bedtools-slop-and-merge-2)
      - [4. Extract the FASTA](#4-extract-the-fasta-2)
    - [community 3 - chr4](#community-3---chr4-10)
      - [1. Extract non-reference sequences](#1-extract-non-reference-sequences-3)
      - [2. Filter non\_reference.bed](#2-filter-non_referencebed-3)
      - [3. Apply bedtools slop and merge](#3-apply-bedtools-slop-and-merge-3)
      - [4. Extract FASTA](#4-extract-fasta)
    - [community 4 - chr5](#community-4---chr5-10)
      - [1. Extract non-reference sequences](#1-extract-non-reference-sequences-4)
      - [2. Filter non\_reference.bed](#2-filter-non_referencebed-4)
      - [3. Apply bedtools slop and merge](#3-apply-bedtools-slop-and-merge-4)
      - [4. Extract the FASTA](#4-extract-the-fasta-3)
  - [Annotation](#annotation)
    - [Get UniRef\_100 and format as diamond database](#get-uniref_100-and-format-as-diamond-database)
    - [Download conversion table between UniProt names and Araport](#download-conversion-table-between-uniprot-names-and-araport)
    - [Extract attributes from UniRef100 FASTA headers](#extract-attributes-from-uniref100-fasta-headers)
    - [community 0 - chr1](#community-0---chr1-11)
      - [1. Diamond](#1-diamond)
      - [2. Exonerate](#2-exonerate)
        - [split diamond results in chunks](#split-diamond-results-in-chunks)
        - [run exonerate](#run-exonerate)
      - [3. Parse exonerate results](#3-parse-exonerate-results)
    - [community 1 - chr2](#community-1---chr2-11)
      - [1. Diamond](#1-diamond-1)
      - [2. Exonerate](#2-exonerate-1)
        - [split diamond results in chunks](#split-diamond-results-in-chunks-1)
        - [run exonerate](#run-exonerate-1)
      - [3. Parse exonerate results](#3-parse-exonerate-results-1)
    - [community 2 - chr3](#community-2---chr3-10)
      - [1. Diamond](#1-diamond-2)
      - [2. Exonerate](#2-exonerate-2)
        - [split diamond results in chunks](#split-diamond-results-in-chunks-2)
        - [run exonerate](#run-exonerate-2)
      - [3. Parse exonerate results](#3-parse-exonerate-results-2)
    - [community 3 - chr4](#community-3---chr4-11)
      - [1. Diamond](#1-diamond-3)
      - [2. Exonerate](#2-exonerate-3)
        - [split diamond results in chunks](#split-diamond-results-in-chunks-3)
        - [run exonerate](#run-exonerate-3)
      - [3. Parse exonerate results](#3-parse-exonerate-results-3)
    - [community 4 - chr5](#community-4---chr5-11)
      - [1. Diamond](#1-diamond-4)
      - [2. Exonerate](#2-exonerate-4)
        - [split](#split)
        - [run exonerate](#run-exonerate-4)
      - [3. Parse exonerate results](#3-parse-exonerate-results-4)
  - [Annotation review](#annotation-review)
    - [Review results that have a corresponding TAIR annotation](#review-results-that-have-a-corresponding-tair-annotation)
    - [Review results with no corresponding TAIR annotation](#review-results-with-no-corresponding-tair-annotation)
      - [community 0 - chr1](#community-0---chr1-12)
      - [community 1 - chr2](#community-1---chr2-12)
      - [community 2 - chr3](#community-2---chr3-11)
      - [community 3 - chr4](#community-3---chr4-12)
      - [community 4 - chr5](#community-4---chr5-12)
- [Final annotation screening](#final-annotation-screening)
  - [genes](#genes-3)
    - [community 0 - chr1](#community-0---chr1-13)
    - [community 1 - chr2](#community-1---chr2-13)
    - [community 2 - chr3](#community-2---chr3-12)
    - [community 3 - chr4](#community-3---chr4-13)
    - [community 4 - chr5](#community-4---chr5-13)
  - [pseudogenes](#pseudogenes-3)
    - [community 0 - chr1](#community-0---chr1-14)
    - [community 1 - chr2](#community-1---chr2-14)
    - [community 2 - chr3](#community-2---chr3-13)
    - [community 3 - chr4](#community-3---chr4-14)
    - [chr5](#chr5)
- [Gene ontology enrichment analysis](#gene-ontology-enrichment-analysis)
- [Sequence-based pangenome analysis](#sequence-based-pangenome-analysis)
  - [Graphs characteristics](#graphs-characteristics)
  - [Obtain node length](#obtain-node-length)
  - [Node coverage matrices](#node-coverage-matrices)
      - [Matrices processing](#matrices-processing)
  - [Node-based similarity](#node-based-similarity)
  - [Growth curve simulation](#growth-curve-simulation)
- [Similarity analysis](#similarity-analysis)
  - [Gene based](#gene-based)
    - [make the figure](#make-the-figure)
  - [Pseudogene based](#pseudogene-based)
    - [make the figure](#make-the-figure-1)
  - [Node based](#node-based)
- [Other data visualisation](#other-data-visualisation)
  - [Assembly based gene classification - PAV](#assembly-based-gene-classification---pav)
    - [Obtain the files](#obtain-the-files)
    - [Make the figure](#make-the-figure-2)
  - [Assembly based gene classification - CNV](#assembly-based-gene-classification---cnv)
    - [Obtain the files](#obtain-the-files-1)
    - [Make the figure](#make-the-figure-3)
  - [Assembly based pseudogene classification - PAV](#assembly-based-pseudogene-classification---pav)
    - [Obtain the files](#obtain-the-files-2)
    - [Make the figure](#make-the-figure-4)
  - [Assembly based pseudogene classification - CNV](#assembly-based-pseudogene-classification---cnv)
    - [Obtain the files](#obtain-the-files-3)
  - [Assembly based node classification](#assembly-based-node-classification)
    - [Make the figures](#make-the-figures)
  - [Pie charts](#pie-charts)


# Building the reference-free *Arabidopsis thaliana* pangenome graph

## Get the assemblies

We found in total 93 suitable assemblies, which are at chromosome level or complete genome level. They are coming from the following BioProjects: PRJEB55353, PRJEB55632, PRJEB50694, PRJEB30763, PRJEB31147, PRJEB37252, PRJEB37257, PRJEB37258, PRJEB37260, PRJEB37261, PRJEB40125, PRJEB51511, PRJNA777107, PRJNA779205, PRJNA828135, PRJNA834751, PRJNA10719, PRJNA915353, PRJNA311266, PRJCA005809, PRJCA007112. 

## Quality check and trimming of the assemblies

We should exclude short contigs from the pangenome, but we need to set a threshold, and for this we must quality inspect the assemblies and check the length distribution. We can do this using QUAST:

```
quast accessions/*.fna.gz accessions/*.fasta.gz -o ./output_quast -r accessions/GCA_000001735.2_TAIR10.1_genomic.fna.gz -t 40 --eukaryote --large -k --plots-format png
```

Based on QUAST results, we can decided to trim the contigs that are smaller than 5kbp. The assemblies that will suffer a loss will be the following according to QUAST:

Assembly | # contigs (>= 0 bp) | # contigs (>= 5000 bp) | Total length (>= 0 bp) | Total length (>= 5000 bp) |
---------|---------------------|-------------------------|------------------------|----------------------------|
GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic | 109 | 78 | 119626746 | 119530365 |
GCA_902460265.3_Arabidopsis_thaliana_An-1_chrom_genomic | 111 | 83 | 120129838 | 120059751 |
GCA_902460275.1_Arabidopsis_thaliana_Cvi-0_genomic | 102 | 71 | 119749512 | 119657772 |
GCA_902460285.1_Arabidopsis_thaliana_Ler_genomic | 105 | 75 | 120338059 | 120246237 |
GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic | 94 | 74 | 120289901 | 120224609 |
GCA_902460305.1_Arabidopsis_thaliana_Kyo_genomic | 184 | 155 | 122202079 | 122112347 |
GCA_902460315.1_Arabidopsis_thaliana_Eri-1_genomic | 142 | 115 | 120795191 | 120728700 |

For the filtering we can use `bioawk`. 

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic.fna > GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic.trim.5k.fna

#check:
grep '>' GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic.fna | wc -l
---
109

grep '>' GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic.trim.5k.fna | wc -l
---
78
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460265.3_Arabidopsis_thaliana_An-1_chrom_genomic.fna > GCA_902460265.3_Arabidopsis_thaliana_An-1_chrom_genomic.trim.5k.fna

#check:
grep '>' GCA_902460265.3_Arabidopsis_thaliana_An-1_chrom_genomic.fna | wc -l
---
111

grep '>' GCA_902460265.3_Arabidopsis_thaliana_An-1_chrom_genomic.trim.5k.fna | wc -l
---
83
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460275.1_Arabidopsis_thaliana_Cvi-0_genomic.fna > GCA_902460275.1_Arabidopsis_thaliana_Cvi-0_genomic.trim.5k.fna

#check:
grep '>' GCA_902460275.1_Arabidopsis_thaliana_Cvi-0_genomic.fna | wc -l
---
102

grep '>' GCA_902460275.1_Arabidopsis_thaliana_Cvi-0_genomic.trim.5k.fna | wc -l
---
71
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460285.1_Arabidopsis_thaliana_Ler_genomic.fna > GCA_902460285.1_Arabidopsis_thaliana_Ler_genomic.trim.5k.fna

#check:
grep '>' GCA_902460285.1_Arabidopsis_thaliana_Ler_genomic.fna | wc -l
---
105

grep '>' GCA_902460285.1_Arabidopsis_thaliana_Ler_genomic.trim.5k.fna | wc -l
---
75
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna > GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.trim.5k.fna 

#check:
grep '>' GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.fna | wc -l
---
94

grep '>' GCA_902460295.1_Arabidopsis_thaliana_Sha_genomic.trim.5k.fna | wc -l
---
74
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460305.1_Arabidopsis_thaliana_Kyo_genomic.fna > GCA_902460305.1_Arabidopsis_thaliana_Kyo_genomic.trim.5k.fna

#check:
grep '>' GCA_902460305.1_Arabidopsis_thaliana_Kyo_genomic.fna | wc -l
---
184

grep '>' GCA_902460305.1_Arabidopsis_thaliana_Kyo_genomic.trim.5k.fna | wc -l
---
155
```

```
bioawk -c fastx '(length($seq)>5000){print ">" $name ORS $seq}' GCA_902460315.1_Arabidopsis_thaliana_Eri-1_genomic.fna > GCA_902460315.1_Arabidopsis_thaliana_Eri-1_genomic.trim.5k.fna

#check:
grep '>' GCA_902460315.1_Arabidopsis_thaliana_Eri-1_genomic.fna | wc -l
---
142

grep '>' GCA_902460315.1_Arabidopsis_thaliana_Eri-1_genomic.trim.5k.fna | wc -l
---
115
```

Now we should zip the .trim.5k.fna files:

```
for i in *trim.5k.fna ; do (echo "gzip $i"); done | bash
```

From now on, we are going to used the trimmed version of these assemblies.

## Renaming the assemblies

We can follow the [S. cerevisiae example](https://pggb.readthedocs.io/en/latest/rst/tutorials/sequence_partitioning.html#download-the-assemblies):

```
ls *.gz | cut -f 1 -d '.' | uniq | while read f; do
    echo $f
    zcat $f.* > $f.fa
done
```

To change the sequence names according to [PanSN-spec](https://github.com/pangenome/PanSN-spec), we can use [fastix](https://github.com/ekg/fastix):

```
ls *.fa | while read f; do
    sample_name=$(echo $f | cut -f 1 -d '.');
    echo ${sample_name}
    fastix -p "${sample_name}#1#" $f >> Arabidopsis.pg.in.fasta
done
```

With `"${sample_name}#1#"` we specify haplotype_id equals to 1 for all the assemblies, as they are all haploid.
The names of the scaffolds are now following the pattern:

`[sample_name][delim][haplotype_id][delim][contig_or_scaffold_name]`

## Format fasta headers and remove organelles from assemblies

We should remove all what comes after `\t` from the fasta headers as this can led to errors in subsequent steps:

```
cat Arabidopsis.pg.in.fasta | sed -E '/^>/s/( +|\t).*//' > Arabidopsis.adj.pg.in.fasta
```

We then should remove organelle scaffolds from our input file. Only three assemblies have Mt e Plt, and the corresponding scaffolds are:

```
>GCA_000001735#1#BK010421.1
>GCA_000001735#1#AP000423.1
>GCA_023115395#1#CP096029.1
>GCA_023115395#1#CP096030.1
>GCA_904420315#1#LR881472.1
>GCA_904420315#1#LR881471.1
```

We need to put these headers in a file that we will call `Mt_and_Plt_headers.txt`.

We can now exclude these scaffolds from our previous multifasta input file:

```
awk '(NR==FNR) { toRemove[$1]; next }
     /^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }
    p' Mt_and_Plt_headers.txt Arabidopsis.adj.pg.in.fasta > Arabidopsis.pg.input.fasta
```

We need to compress our input and to index it:

```
bgzip -@ 16 Arabidopsis.pg.input.fasta
samtools faidx Arabidopsis.pg.input.fasta.gz
```

## Sequence Partitioning using partition-before-pggb

We can now run `partition-before-pggb` on our multifasta (note: if it is not possible to estimate the divergence, you should try different settings to individuate the right value to set for `-p`):

```
partition-before-pggb -i Arabidopsis.pg.input.fasta.gz -o output_pbp_p90 -n 93 -t 100 -Y '#' -p 90 -s 10k -V 'GCA_000001735:#:1000' -m -D ./temp  
```

- -i input
- -o output directory
- -n number of haplotypes
- -p percent identity for mapping/alignment
- -s segment length for mapping [default: 5000]
- -V to produce vcf against ref
- -m generate MultiQC report of graphs' statistics and visualizations,automatically runs odgi stats
- -D temporary files directory


## Running pggb on each community

We can now run pggb on each community! :-)

To build the graph of each community we will use `pggb -p 95`.

### community 0 - chr1

```
pggb -i output_pbp_p90/Arabidopsis.pg.input.fasta.gz.c325321.community.0.fa \
     -o output_pbp_p90/community.0.out \
     -s 10000 -l 50000 -p 95 -n 93 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10000000 \
     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V GCA_000001735:#:1000 --multiqc --temp-dir ./temp --threads 100 --poa-threads 100
```

### community 1 - chr2

```
pggb -i output_pbp_p90/Arabidopsis.pg.input.fasta.gz.c325321.community.1.fa \
     -o output_pbp_p90/community.1.out \
     -s 10000 -l 50000 -p 95 -n 93 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10000000 \
     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V GCA_000001735:#:1000 --multiqc --temp-dir ./temp --threads 100 --poa-threads 100
```

### community 2 - chr3

```
pggb -i output_pbp_p90/Arabidopsis.pg.input.fasta.gz.c325321.community.2.fa \
     -o output_pbp_p90/community.2.out \
     -s 10000 -l 50000 -p 95 -n 93 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10000000 \
     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V GCA_000001735:#:1000 --multiqc --temp-dir ./temp --threads 100 --poa-threads 100
```

### community 3 - chr4

```
pggb -i output_pbp_p90/Arabidopsis.pg.input.fasta.gz.c325321.community.3.fa \
     -o output_pbp_p90/community.3.out \
     -s 10000 -l 50000 -p 95 -n 93 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10000000 \
     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V GCA_000001735:#:1000 --multiqc --temp-dir ./temp --threads 100 --poa-threads 100
```

### community 4 - chr5

```
pggb -i output_pbp_p90/Arabidopsis.pg.input.fasta.gz.c325321.community.4.fa \
     -o output_pbp_p90/community.4.out \
     -s 10000 -l 50000 -p 95 -n 93 -K 19 -F 0.001 -g 30 \
     -k 23 -f 0 -B 10000000 \
     -j 0 -e 0 -G 700,900,1100 -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
     -Y "#" -V GCA_000001735:#:1000 --multiqc --temp-dir ./temp --threads 100 --poa-threads 100
```

# Pangenome annotation - reference based

## Get odgi untangle for annotation

To use `odgi untangle` for annotating the pangenome graph, we will switch to `untangle_for_annotation` branch.

```
cd odgi
git checkout untangle_for_annotation && git pull && git submodule update --init --recursive
cmake -H. -Bbuild && cmake --build build -- -j 16
```

## Prepare the annotation files for injection

### Genes

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/735/GCA_000001735.2_TAIR10.1/GCA_000001735.2_TAIR10.1_genomic.gff.gz
gunzip GCA_000001735.2_TAIR10.1_genomic.gff.gz
```

Now, we will inject only the genes, so we need to extract only gene rows from the gff of TAIR10 and split them based on chromosomes.

Let's have a look to the gff file:

```
head GCA_000001735.2_TAIR10.1_genomic.gff
---
##gff-version 3
#!gff-spec-version 1.21
#!processor NCBI annotwriter
#!genome-build TAIR10.1
#!genome-build-accession NCBI_Assembly:GCA_000001735.2
##sequence-region CP002684.1 1 30427671
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702
CP002684.1      Genbank region  1       30427671        .       +       .       ID=CP002684.1:1..30427671;Dbxref=taxon:3702;Name=1;chromosome=1;ecotype=Columbia;gbkey=Src;mol_type=genomic DNA
CP002684.1      Genbank gene    3631    5899    .       +       .       ID=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC domain containing protein 1,T25K16.1,T25K16_1;locus_tag=AT1G01010
CP002684.1      Genbank mRNA    3631    5899    .       +       .       ID=rna-gnl|JCVI|mRNA.AT1G01010.1;Parent=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010;gbkey=mRNA;gene=NAC001;inference=similar to RNA sequence%2C mRNA:INSD:BT001115.1%2CINSD:AF439834.1%2CINSD:AK226863.1;locus_tag=AT1G01010;orig_protein_id=gnl|JCVI|AT1G01010.1;orig_transcript_id=gnl|JCVI|mRNA.AT1G01010.1;product=NAC domain containing protein 1
```

#### community 0 - chr1

In order to inject the genes in the pangenome, we need to extract the information we need from the gff, and modify them according to the PanSN-spec naming:

```
awk -F "\t|ID=gene-|;Dbxref" '$1 ~/^CP002684.1$/ && $3 ~/^gene$/ {print ("GCA_000001735#1#"$1"\t"$4"\t"$5"\t"$10":"$4"-"$5)}' ../../../GCA_000001735.2_TAIR10.1_genomic.gff > chr1.genes.adj.TAIR10.bed

head chr1.genes.adj.TAIR10.bed | column -t
---
GCA_000001735#1#CP002684.1  3631   5899   AT1G01010:3631-5899
GCA_000001735#1#CP002684.1  6788   9130   AT1G01020:6788-9130
GCA_000001735#1#CP002684.1  11101  11372  AT1G03987:11101-11372
GCA_000001735#1#CP002684.1  11649  13714  AT1G01030:11649-13714
GCA_000001735#1#CP002684.1  23121  31227  AT1G01040:23121-31227
GCA_000001735#1#CP002684.1  23312  24099  AT1G03993:23312-24099
GCA_000001735#1#CP002684.1  28500  28706  AT1G01046:28500-28706
GCA_000001735#1#CP002684.1  31170  33171  AT1G01050:31170-33171
GCA_000001735#1#CP002684.1  32727  33009  AT1G03997:32727-33009
GCA_000001735#1#CP002684.1  33365  37871  AT1G01060:33365-37871
```

Make sure that community 0 is chromosome 1 of TAIR10:

```
odgi paths -i Arabidopsis.pg.input.community.0.og -L | grep 'GCA_000001735'
---
GCA_000001735#1#CP002684.1
```

#### community 1 - chr2


In order to inject the genes in the pangenome, we need to extract the information we need from the gff, and modify them according to the PanSN-spec naming:

```
awk -F "\t|ID=gene-|;Dbxref" '$1 ~/^CP002685.1$/ && $3 ~/^gene$/ {print ("GCA_000001735#1#"$1"\t"$4"\t"$5"\t"$10":"$4"-"$5)}' ../../../GCA_000001735.2_TAIR10.1_genomic.gff > chr2.genes.adj.TAIR10.bed

head chr2.genes.adj.TAIR10.bed | column -t
---
GCA_000001735#1#CP002685.1  1025   3173   AT2G01008:1025-3173
GCA_000001735#1#CP002685.1  2805   3176   AT2G03855:2805-3176
GCA_000001735#1#CP002685.1  3706   5513   AT2G01010:3706-5513
GCA_000001735#1#CP002685.1  5782   5945   AT2G01020:5782-5945
GCA_000001735#1#CP002685.1  6571   6672   AT2G01021:6571-6672
GCA_000001735#1#CP002685.1  7090   7505   AT2G03865:7090-7505
GCA_000001735#1#CP002685.1  7669   9399   AT2G03875:7669-9399
GCA_000001735#1#CP002685.1  9648   9767   AT2G01023:9648-9767
GCA_000001735#1#CP002685.1  9849   10177  AT2G00340:9849-10177
GCA_000001735#1#CP002685.1  50288  51487  AT2G01035:50288-51487
```

Make sure that community 1 is chromosome 2 of TAIR10:

```
odgi paths -i Arabidopsis.pg.input.community.1.og -L | grep 'GCA_000001735'
---
GCA_000001735#1#CP002685.1
```

#### community 2 - chr3

In order to inject the genes in the pangenome, we need to extract the information we need from the gff, and modify them according to the PanSN-spec naming:

```
awk -F "\t|ID=gene-|;Dbxref" '$1 ~/^CP002686.1$/ && $3 ~/^gene$/ {print ("GCA_000001735#1#"$1"\t"$4"\t"$5"\t"$10":"$4"-"$5)}' ../../../GCA_000001735.2_TAIR10.1_genomic.gff > chr3.genes.adj.TAIR10.bed

head chr3.genes.adj.TAIR10.bed | column -t
---
GCA_000001735#1#CP002686.1  1609   4159   AT3G01015:1609-4159
GCA_000001735#1#CP002686.1  4342   4818   AT3G01010:4342-4818
GCA_000001735#1#CP002686.1  5104   6149   AT3G01020:5104-6149
GCA_000001735#1#CP002686.1  6657   7772   AT3G01030:6657-7772
GCA_000001735#1#CP002686.1  8723   12697  AT3G01040:8723-12697
GCA_000001735#1#CP002686.1  13046  15906  AT3G01050:13046-15906
GCA_000001735#1#CP002686.1  15934  16320  AT3G00970:15934-16320
GCA_000001735#1#CP002686.1  16631  18909  AT3G01060:16631-18909
GCA_000001735#1#CP002686.1  19409  20806  AT3G01070:19409-20806
GCA_000001735#1#CP002686.1  25355  27712  AT3G01080:25355-27712
```

Make sure that community 2 is chromosome 3 of TAIR10:

```
odgi paths -i Arabidopsis.pg.input.community.2.og -L | grep 'GCA_000001735'
---
GCA_000001735#1#CP002686.1
```

#### community 3 - chr4

In order to inject the genes in the pangenome, we need to extract the information we need from the gff, and modify them according to the PanSN-spec naming:

```
awk -F "\t|ID=gene-|;Dbxref" '$1 ~/^CP002687.1$/ && $3 ~/^gene$/ {print ("GCA_000001735#1#"$1"\t"$4"\t"$5"\t"$10":"$4"-"$5)}' ../../../GCA_000001735.2_TAIR10.1_genomic.gff > chr4.genes.adj.TAIR10.bed

head chr4.genes.adj.TAIR10.bed | column -t
---
GCA_000001735#1#CP002687.1  1180   1536   AT4G00005:1180-1536
GCA_000001735#1#CP002687.1  2895   10504  AT4G00020:2895-10504
GCA_000001735#1#CP002687.1  10815  13359  AT4G00026:10815-13359
GCA_000001735#1#CP002687.1  13527  14413  AT4G00030:13527-14413
GCA_000001735#1#CP002687.1  14627  16079  AT4G00040:14627-16079
GCA_000001735#1#CP002687.1  17639  20183  AT4G00050:17639-20183
GCA_000001735#1#CP002687.1  21351  29064  AT4G00060:21351-29064
GCA_000001735#1#CP002687.1  29072  31426  AT4G00070:29072-31426
GCA_000001735#1#CP002687.1  32748  33756  AT4G00080:32748-33756
GCA_000001735#1#CP002687.1  33800  33872  AT4G00085:33800-33872
```

Make sure that community 3 is chromosome 4 of TAIR10:

```
odgi paths -i Arabidopsis.pg.input.community.3.og -L | grep 'GCA_000001735'
---
GCA_000001735#1#CP002687.1
```

#### community 4 - chr5

In order to inject the genes in the pangenome, we need to extract the information we need from the gff, and modify them according to the PanSN-spec naming:

```
awk -F "\t|ID=gene-|;Dbxref" '$1 ~/^CP002688.1$/ && $3 ~/^gene$/ {print ("GCA_000001735#1#"$1"\t"$4"\t"$5"\t"$10":"$4"-"$5)}' ../../../GCA_000001735.2_TAIR10.1_genomic.gff > chr5.genes.adj.TAIR10.bed

head chr5.genes.adj.TAIR10.bed | column -t
---
GCA_000001735#1#CP002688.1  2      303    AT5G00730:2-303
GCA_000001735#1#CP002688.1  995    5156   AT5G01010:995-5156
GCA_000001735#1#CP002688.1  5256   5907   AT5G01015:5256-5907
GCA_000001735#1#CP002688.1  5339   5593   AT5G01017:5339-5593
GCA_000001735#1#CP002688.1  5917   8467   AT5G01020:5917-8467
GCA_000001735#1#CP002688.1  9780   13235  AT5G01030:9780-13235
GCA_000001735#1#CP002688.1  13128  16236  AT5G01040:13128-16236
GCA_000001735#1#CP002688.1  18086  20887  AT5G01050:18086-20887
GCA_000001735#1#CP002688.1  22684  24934  AT5G01060:22684-24934
GCA_000001735#1#CP002688.1  25012  25908  AT5G01070:25012-25908
```

Make sure that community 4 is chromosome 5 of TAIR10:

```
odgi paths -i Arabidopsis.pg.input.community.4.og -L | grep 'GCA_000001735'
---
GCA_000001735#1#CP002688.1
```

### Pseudogenes

The original file from [Mascagni et al., 2021](https://www.nature.com/articles/s41598-021-84778-6) was modified as explained in scripts/bed_to_gff3.ipynb. The coordinates of pseudogenes annotated in the reverse strain were then converted in forward strain to allow compatibility with the downstream analysis. This is how the final annotations file looks like:

```
head pseudogenes.txt
---
Chromosome	end	start	coverage	Genewise_score	position	size	size	Strand	Type	Pater	Species	Stop_codon	Frameshift
Chr1	6144	6269	0.09	56.39	UTR	129	plus	FRAG	AT1G01010.1	A.thaliana	0	0
Chr1	47087	47325	0.06	74.09	Intergenic	206	plus	FRAG	AT1G02730.1	A.thaliana	1	3
Chr1	89932	90153	0.02	24.75	Intergenic	51	plus	FRAG	AT5G41740.2	A.thaliana	1	3
Chr1	259011	259151	0.06	48.95	Intergenic	84	plus	FRAG	AT1G01690.1	A.thaliana	0	0
Chr1	266562	267495	0.99	68.96	UTR	979	plus	DUP	AT4G00525.1	A.thaliana	0	0
Chr1	426018	426176	0.06	17.28	Intergenic	36	plus	FRAG	AT3G04430.1	A.thaliana	1	1
Chr1	578305	578484	0.04	47.13	Intron_cds	102	plus	FRAG	AT1G05120.1	A.thaliana	0	0
Chr1	582401	582585	0.46	52.69	Intergenic	236	plus	SE	AT4G02160.1	A.thaliana	0	1
Chr1	582634	582971	0.57	117.33	Intergenic	247	plus	SE	AT5G61710.1	A.thaliana	0	1
```

#### community 0 - chr1

```
awk '$1 == "Chr1"' ../../../../pseudogenes.txt | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' | sed 's/Chr1/GCA_000001735#1#CP002684.1/g' > Pseudogenes_TAIR10_chr1.bed

head Pseudogenes_TAIR10_chr1.bed 
---
GCA_000001735#1#CP002684.1	6144	6269	GCA_000001735#1#CP002684.1:6144-6269
GCA_000001735#1#CP002684.1	47087	47325	GCA_000001735#1#CP002684.1:47087-47325
GCA_000001735#1#CP002684.1	89932	90153	GCA_000001735#1#CP002684.1:89932-90153
GCA_000001735#1#CP002684.1	259011	259151	GCA_000001735#1#CP002684.1:259011-259151
GCA_000001735#1#CP002684.1	266562	267495	GCA_000001735#1#CP002684.1:266562-267495
GCA_000001735#1#CP002684.1	426018	426176	GCA_000001735#1#CP002684.1:426018-426176
GCA_000001735#1#CP002684.1	578305	578484	GCA_000001735#1#CP002684.1:578305-578484
GCA_000001735#1#CP002684.1	582401	582585	GCA_000001735#1#CP002684.1:582401-582585
GCA_000001735#1#CP002684.1	582634	582971	GCA_000001735#1#CP002684.1:582634-582971
GCA_000001735#1#CP002684.1	893047	893493	GCA_000001735#1#CP002684.1:893047-893493
```

#### community 1 - chr2

```
awk '$1 == "Chr2"' ../../../../pseudogenes.txt | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' | sed 's/Chr2/GCA_000001735#1#CP002685.1/g' > Pseudogenes_TAIR10_chr2.bed

head Pseudogenes_TAIR10_chr2.bed 
---
GCA_000001735#1#CP002685.1	2629	2769	GCA_000001735#1#CP002685.1:2629-2769
GCA_000001735#1#CP002685.1	2939	3079	GCA_000001735#1#CP002685.1:2939-3079
GCA_000001735#1#CP002685.1	126312	126472	GCA_000001735#1#CP002685.1:126312-126472
GCA_000001735#1#CP002685.1	167479	167678	GCA_000001735#1#CP002685.1:167479-167678
GCA_000001735#1#CP002685.1	184200	184421	GCA_000001735#1#CP002685.1:184200-184421
GCA_000001735#1#CP002685.1	214567	214668	GCA_000001735#1#CP002685.1:214567-214668
GCA_000001735#1#CP002685.1	249464	249634	GCA_000001735#1#CP002685.1:249464-249634
GCA_000001735#1#CP002685.1	287831	287959	GCA_000001735#1#CP002685.1:287831-287959
GCA_000001735#1#CP002685.1	338842	339105	GCA_000001735#1#CP002685.1:338842-339105
GCA_000001735#1#CP002685.1	343391	343537	GCA_000001735#1#CP002685.1:343391-343537
```

#### community 2 - chr3

```
awk '$1 == "Chr3"' ../../../../pseudogenes.txt | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' | sed 's/Chr3/GCA_000001735#1#CP002686.1/g' > Pseudogenes_TAIR10_chr3.bed

head Pseudogenes_TAIR10_chr3.bed 
---
GCA_000001735#1#CP002686.1	1009	1200	GCA_000001735#1#CP002686.1:1009-1200
GCA_000001735#1#CP002686.1	430020	430130	GCA_000001735#1#CP002686.1:430020-430130
GCA_000001735#1#CP002686.1	484416	484737	GCA_000001735#1#CP002686.1:484416-484737
GCA_000001735#1#CP002686.1	630152	630235	GCA_000001735#1#CP002686.1:630152-630235
GCA_000001735#1#CP002686.1	792374	792475	GCA_000001735#1#CP002686.1:792374-792475
GCA_000001735#1#CP002686.1	799113	799223	GCA_000001735#1#CP002686.1:799113-799223
GCA_000001735#1#CP002686.1	832148	832279	GCA_000001735#1#CP002686.1:832148-832279
GCA_000001735#1#CP002686.1	863502	863939	GCA_000001735#1#CP002686.1:863502-863939
GCA_000001735#1#CP002686.1	864622	864843	GCA_000001735#1#CP002686.1:864622-864843
GCA_000001735#1#CP002686.1	987088	990672	GCA_000001735#1#CP002686.1:987088-990672
```

#### community 3 - chr4

```
awk '$1 == "Chr4"' ../../../../pseudogenes.txt | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' | sed 's/Chr4/GCA_000001735#1#CP002687.1/g' > Pseudogenes_TAIR10_chr4.bed

head Pseudogenes_TAIR10_chr4.bed 
---
GCA_000001735#1#CP002687.1	11253	11351	GCA_000001735#1#CP002687.1:11253-11351
GCA_000001735#1#CP002687.1	62136	62339	GCA_000001735#1#CP002687.1:62136-62339
GCA_000001735#1#CP002687.1	70243	70416	GCA_000001735#1#CP002687.1:70243-70416
GCA_000001735#1#CP002687.1	99056	99229	GCA_000001735#1#CP002687.1:99056-99229
GCA_000001735#1#CP002687.1	99355	99464	GCA_000001735#1#CP002687.1:99355-99464
GCA_000001735#1#CP002687.1	102138	102411	GCA_000001735#1#CP002687.1:102138-102411
GCA_000001735#1#CP002687.1	102614	103287	GCA_000001735#1#CP002687.1:102614-103287
GCA_000001735#1#CP002687.1	103362	103703	GCA_000001735#1#CP002687.1:103362-103703
GCA_000001735#1#CP002687.1	105289	105375	GCA_000001735#1#CP002687.1:105289-105375
GCA_000001735#1#CP002687.1	121936	122058	GCA_000001735#1#CP002687.1:121936-122058
```

### community 4 - chr5

```
awk '$1 == "Chr5"' ../../../../pseudogenes.txt | awk -v OFS='\t' '{print($1,$2,$3,$1":"$2"-"$3)}' | sed 's/Chr5/GCA_000001735#1#CP002688.1/g' > Pseudogenes_TAIR10_chr5.bed

head Pseudogenes_TAIR10_chr5.bed 
---
GCA_000001735#1#CP002688.1	17952	18035	GCA_000001735#1#CP002688.1:17952-18035
GCA_000001735#1#CP002688.1	18038	18085	GCA_000001735#1#CP002688.1:18038-18085
GCA_000001735#1#CP002688.1	453356	453511	GCA_000001735#1#CP002688.1:453356-453511
GCA_000001735#1#CP002688.1	453793	454342	GCA_000001735#1#CP002688.1:453793-454342
GCA_000001735#1#CP002688.1	498633	498719	GCA_000001735#1#CP002688.1:498633-498719
GCA_000001735#1#CP002688.1	595390	595590	GCA_000001735#1#CP002688.1:595390-595590
GCA_000001735#1#CP002688.1	597848	597991	GCA_000001735#1#CP002688.1:597848-597991
GCA_000001735#1#CP002688.1	660001	660384	GCA_000001735#1#CP002688.1:660001-660384
GCA_000001735#1#CP002688.1	672914	673057	GCA_000001735#1#CP002688.1:672914-673057
GCA_000001735#1#CP002688.1	679891	680126	GCA_000001735#1#CP002688.1:679891-680126
```

## Inject genes from the reference  into the pangenome graph

### Genes 

#### community 0 - chr1

```
odgi inject -i Arabidopsis.pg.input.community.0.og -b chr1.genes.adj.TAIR10.bed -o Arabidopsis.pg.community.0.inject.adj.og -t 100 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.0.inject.adj.og -L | grep 'AT1G' > chr1.genenames.txt

# security check
wc -l chr1.genenames.txt 
---
8771 chr1.genenames.txt

wc -l chr1.genes.adj.TAIR10.bed
---
8771 chr1.genes.adj.TAIR10.bed

```

#### community 1 - chr2

```
odgi inject -i Arabidopsis.pg.input.community.1.og -b chr2.genes.adj.TAIR10.bed -o Arabidopsis.pg.community.1.inject.adj.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.1.inject.adj.og -L | grep 'AT2G' > chr2.genenames.txt

# security check
wc -l chr2.genenames.txt 
---
5265 chr2.genenames.txt

wc -l chr2.genes.adj.TAIR10.bed
---
5265 chr2.genes.adj.TAIR10.bed

```

#### community 2 - chr3

```
odgi inject -i Arabidopsis.pg.input.community.2.og -b chr3.genes.adj.TAIR10.bed -o Arabidopsis.pg.community.2.inject.adj.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.2.inject.adj.og -L | grep 'AT3G' > chr3.genenames.txt

# security check
wc -l chr3.genenames.txt 
---
6544 chr3.genenames.txt

wc -l chr3.genes.adj.TAIR10.bed
---
6544 chr3.genes.adj.TAIR10.bed

```


#### community 3 - chr4


```
odgi inject -i Arabidopsis.pg.input.community.3.og -b chr4.genes.adj.TAIR10.bed -o Arabidopsis.pg.community.3.inject.adj.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.3.inject.adj.og -L | grep 'AT4G' > chr4.genenames.txt

# security check
wc -l chr4.genenames.txt 
---
5007 chr4.genenames.txt

wc -l chr4.genes.adj.TAIR10.bed
---
5007 chr4.genes.adj.TAIR10.bed

```

#### community 4 - chr5


```
odgi inject -i Arabidopsis.pg.input.community.4.og -b chr5.genes.adj.TAIR10.bed -o Arabidopsis.pg.community.4.inject.adj.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.4.inject.adj.og -L | grep 'AT5G' > chr5.genenames.txt

# security check
wc -l chr5.genenames.txt 
---
7469 chr5.genenames.txt

wc -l chr5.genes.adj.TAIR10.bed
---
7469 chr5.genes.adj.TAIR10.bed

```

### Pseudogenes

#### community 0 - chr1

```
odgi inject -i ../Arabidopsis.pg.input.community.0.og -b Pseudogenes_TAIR10_chr1.bed -o Arabidopsis.pg.community.0.inject.pseudogenes.og -t 50 -P

```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.0.inject.pseudogenes.og -L | grep 'GCA_000001735#1#CP002684.1:' > chr1.pseudogenenames.txt

# security check
wc -l chr1.pseudogenenames.txt 
---
1083 chr1.pseudogenenames.txt

wc -l Pseudogenes_TAIR10_chr1.bed
---
1083 Pseudogenes_TAIR10_chr1.bed

```

#### community 1 - chr2

```
odgi inject -i ../Arabidopsis.pg.input.community.1.og -b Pseudogenes_TAIR10_chr2.bed -o Arabidopsis.pg.community.1.inject.pseudogenes.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.1.inject.pseudogenes.og -L | grep 'GCA_000001735#1#CP002685.1:' > chr2.pseudogenenames.txt

# security check
wc -l chr2.pseudogenenames.txt 
---
904 chr2.pseudogenenames.txt

wc -l Pseudogenes_TAIR10_chr2.bed
---
904 Pseudogenes_TAIR10_chr2.bed

```

#### community 2 - chr3

```
odgi inject -i ../Arabidopsis.pg.input.community.2.og -b Pseudogenes_TAIR10_chr3.bed -o Arabidopsis.pg.community.2.inject.pseudogenes.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.2.inject.pseudogenes.og -L | grep 'GCA_000001735#1#CP002686.1:' > chr3.pseudogenenames.txt

# security check
wc -l chr3.pseudogenenames.txt 
---
965 chr3.pseudogenenames.txt

wc -l Pseudogenes_TAIR10_chr3.bed
---
965 Pseudogenes_TAIR10_chr3.bed

```

#### community 3 - chr4

```
odgi inject -i ../Arabidopsis.pg.input.community.3.og -b Pseudogenes_TAIR10_chr4.bed -o Arabidopsis.pg.community.3.inject.pseudogenes.og -t 50 -P
```


Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.3.inject.pseudogenes.og -L | grep 'GCA_000001735#1#CP002687.1:' > chr4.pseudogenenames.txt

# security check
wc -l chr4.pseudogenenames.txt 
---
804 chr4.pseudogenenames.txt

wc -l Pseudogenes_TAIR10_chr4.bed
---
804 Pseudogenes_TAIR10_chr4.bed

```

#### community 4 - chr5

```
odgi inject -i ../Arabidopsis.pg.input.community.4.og -b Pseudogenes_TAIR10_chr5.bed -o Arabidopsis.pg.community.4.inject.pseudogenes.og -t 50 -P
```

Extract gene names from injected pangenome:

```
odgi paths -i Arabidopsis.pg.community.4.inject.pseudogenes.og -L | grep 'GCA_000001735#1#CP002688.1:' > chr5.pseudogenenames.txt

# security check
wc -l chr5.pseudogenenames.txt 
---
919 chr5.pseudogenenames.txt

wc -l Pseudogenes_TAIR10_chr5.bed
---
919 Pseudogenes_TAIR10_chr5.bed

```

## Untangle genes across all the assemblies of the pangenome

### Genes 

#### community 0 - chr1

```
odgi stepindex -i Arabidopsis.pg.community.0.inject.adj.og -a 0 -o Arabidopsis.pg.community.0.inject.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.0.inject.adj.og -a Arabidopsis.pg.community.0.inject.index -R chr1.genenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.0.inject.paf
```

#### community 1 - chr2

```
odgi stepindex -i Arabidopsis.pg.community.1.inject.adj.og -a 0 -o Arabidopsis.pg.community.1.inject.index -t 50 -P

odgi untangle -i Arabidopsis.pg.community.1.inject.adj.og -a Arabidopsis.pg.community.1.inject.index -R chr2.genenames.txt -j 0 -t 50 -P -p > Arabidopsis.pg.community.1.inject.paf
```

#### community 2 - chr3

```
odgi stepindex -i Arabidopsis.pg.community.2.inject.adj.og -a 0 -o Arabidopsis.pg.community.2.inject.index -t 50 -P

odgi untangle -i Arabidopsis.pg.community.2.inject.adj.og -a Arabidopsis.pg.community.2.inject.index -R chr3.genenames.txt -j 0 -t 50 -P -p > Arabidopsis.pg.community.2.inject.paf
```

#### community 3 - chr4

```
odgi stepindex -i Arabidopsis.pg.community.3.inject.adj.og -a 0 -o Arabidopsis.pg.community.3.inject.index -t 50 -P

odgi untangle -i Arabidopsis.pg.community.3.inject.adj.og -a Arabidopsis.pg.community.3.inject.index -R chr4.genenames.txt -j 0 -t 50 -P -p > Arabidopsis.pg.community.3.inject.paf
```

#### community 4 - chr5

```
odgi stepindex -i Arabidopsis.pg.community.4.inject.adj.og -a 0 -o Arabidopsis.pg.community.4.inject.index -t 50 -P

odgi untangle -i Arabidopsis.pg.community.4.inject.adj.og -a Arabidopsis.pg.community.4.inject.index -R chr5.genenames.txt -j 0 -t 50 -P -p > Arabidopsis.pg.community.4.inject.paf
```

### Pseudogenes

#### community 0 - chr1

```
odgi stepindex -i Arabidopsis.pg.community.0.inject.pseudogenes.og -a 0 -o Arabidopsis.pg.community.0.inject.pseudogenes.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.0.inject.pseudogenes.og -a Arabidopsis.pg.community.0.inject.pseudogenes.index -R chr1.pseudogenenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.0.inject.paf
```

#### community 1 - chr2

```
odgi stepindex -i Arabidopsis.pg.community.1.inject.pseudogenes.og -a 0 -o Arabidopsis.pg.community.1.inject.pseudogenes.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.1.inject.pseudogenes.og -a Arabidopsis.pg.community.1.inject.pseudogenes.index -R chr2.pseudogenenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.1.inject.paf
```

#### community 2 - chr3

```
odgi stepindex -i Arabidopsis.pg.community.2.inject.pseudogenes.og -a 0 -o Arabidopsis.pg.community.2.inject.pseudogenes.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.2.inject.pseudogenes.og -a Arabidopsis.pg.community.2.inject.pseudogenes.index -R chr3.pseudogenenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.2.inject.paf
```

#### community 3 - chr4


```
odgi stepindex -i Arabidopsis.pg.community.3.inject.pseudogenes.og -a 0 -o Arabidopsis.pg.community.3.inject.pseudogenes.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.3.inject.pseudogenes.og -a Arabidopsis.pg.community.3.inject.pseudogenes.index -R chr4.pseudogenenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.3.inject.paf
```

#### community 4 - chr5


```
odgi stepindex -i Arabidopsis.pg.community.4.inject.pseudogenes.og -a 0 -o Arabidopsis.pg.community.4.inject.pseudogenes.index -t 100 -P

odgi untangle -i Arabidopsis.pg.community.4.inject.pseudogenes.og -a Arabidopsis.pg.community.4.inject.pseudogenes.index -R chr5.pseudogenenames.txt -j 0 -t 100 -P -p > Arabidopsis.pg.community.4.inject.paf
```


## Untangle PAF processing 

### collapse and merge 

#### community 0 - chr1

```
# genes

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.0.inject.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

```
# pseudogenes

# pseudogenes paf needs one more step before of the processing, due to the way we named pseudogenes. We have to exclude the alignments of pseudogenes against eachother. This is done automatically for genes, because their pattern doesn't match with the reference name

grep -v '^GCA_000001735#1#CP002684.1:' Arabidopsis.pg.community.0.inject.paf > Arabidopsis.pg.community.0.inject.filtered.paf

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.0.inject.filtered.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```


#### community 1 - chr2

```
# genes

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.1.inject.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

```
# pseudogenes

# pseudogenes paf needs one more step before of the processing, due to the way we named pseudogenes. We have to exclude the alignments of pseudogenes against eachother. This is done automatically for genes, because their pattern doesn't match with the reference name

grep -v '^GCA_000001735#1#CP002685.1:' Arabidopsis.pg.community.1.inject.paf > Arabidopsis.pg.community.1.inject.filtered.paf

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.1.inject.filtered.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```


#### community 2 - chr3

```
# genes 

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.2.inject.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

```
# pseudogenes
# pseudogenes paf needs one more step before of the processing, due to the way we named pseudogenes. We have to exclude the alignments of pseudogenes against eachother. This is done automatically for genes, because their pattern doesn't match with the reference name

grep -v '^GCA_000001735#1#CP002686.1:' Arabidopsis.pg.community.2.inject.paf > Arabidopsis.pg.community.2.inject.filtered.paf

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.2.inject.filtered.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```



#### community 3 - chr4

```
# genes

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.3.inject.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

```
# pseudogenes

# pseudogenes paf needs one more step before of the processing, due to the way we named pseudogenes. We have to exclude the alignments of pseudogenes against eachother. This is done automatically for genes, because their pattern doesn't match with the reference name

grep -v '^GCA_000001735#1#CP002687.1:' Arabidopsis.pg.community.3.inject.paf > Arabidopsis.pg.community.3.inject.filtered.paf

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.3.inject.filtered.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```


#### community 4 - chr5

```
# genes

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.4.inject.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

```
# pseudogenes

# pseudogenes paf needs one more step before of the processing, due to the way we named pseudogenes. We have to exclude the alignments of pseudogenes against eachother. This is done automatically for genes, because their pattern doesn't match with the reference name

grep -v '^GCA_000001735#1#CP002688.1:' Arabidopsis.pg.community.4.inject.paf > Arabidopsis.pg.community.4.inject.filtered.paf

~/software/scripts/run_collapse_and_merge.sh -i ../Arabidopsis.pg.community.4.inject.filtered.paf -c ~/software/scripts/collapse_paf.py -m ~/software/scripts/merge_paf.py -d 100 -j 20
```

### genes classification

#### community 0 - chr1

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_genes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff /GCA_000001735.2_TAIR10.1_genomic.gff -ft gene
```

To verify correctness of matrix:

```
Maximum value: 49.0 found for Feature: AT1G40104 and Assembly: GCA_949796415

grep 'AT1G40104' filter_passed_features.csv | grep 'GCA_949796415' | wc -l
49
```


#### community 1 - chr2 

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_genes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff /GCA_000001735.2_TAIR10.1_genomic.gff -ft gene
```

To verify correctness of matrix:

```
Maximum value: 129.0 found for Feature: AT2G01020 and Assembly: GCA_946412065

grep 'AT2G01020' filter_passed_features.csv | grep 'GCA_946412065' | wc -l
129
```

#### communty 2 - chr3  

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_genes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff /GCA_000001735.2_TAIR10.1_genomic.gff -ft gene
```

To verify correctness of matrix:

```
Maximum value: 42.0 found for Feature: AT3G41761 and Assembly: GCA_933208065

grep 'AT3G41761' filter_passed_features.csv | grep 'GCA_933208065' | wc -l
42
```


#### community 3 - chr4 

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_genes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff /GCA_000001735.2_TAIR10.1_genomic.gff -ft gene
```

To verify correctness of matrix:

```
Maximum value: 6.0 found for Feature: AT4G20680 and Assembly: GCA_949796755

grep 'AT4G20680' filter_passed_features.csv | grep 'GCA_949796755' | wc -l
6
```

#### community 4 - chr5 

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_genes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff /GCA_000001735.2_TAIR10.1_genomic.gff -ft gene
```

To verify correctness of matrix:

```
Maximum value: 8.0 found for Feature: AT5G05050 and Assembly: GCA_949796525

grep 'AT5G05050' filter_passed_features.csv | grep 'GCA_949796525' | wc -l
8
```


### pseudogenes classification

#### community 0 - chr1

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_pseudogenes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff pseudogenes.gff -ft pseudogene -chr 'CP002684.1' -ref 'GCA_000001735'
```

#### community 1 - chr2

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_pseudogenes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff pseudogenes.gff -ft pseudogene -chr 'CP002685.1' -ref 'GCA_000001735'
```

#### community 2 - chr3

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_pseudogenes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff pseudogenes.gff -ft pseudogene -chr 'CP002686.1' -ref 'GCA_000001735'
```

#### community 3 - chr4

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_pseudogenes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff pseudogenes.gff -ft pseudogene -chr 'CP002687.1' -ref 'GCA_000001735'
```

#### community 4 - chr5

```
mkdir screening_cov95_id95
cd screening_cov95_id95

~/software/scripts/core_dispensable_pseudogenes.py -i ../collapse_and_merge/collapsed_and_merged.paf -id 95 -jc 0.95 -cov 95 -na 93 -sl 75 -gff pseudogenes.gff -ft pseudogene -chr 'CP002688.1' -ref 'GCA_000001735'
```

# Pangenome annotation - non-reference sequences

## Get the FASTA file for annotation 

### community 0 - chr1

#### 1. Extract non-reference sequences

- find the reference paths using odgi paths + grep 'reference_name' and place it in `reference.txt`
- extract non reference sequences using odgi paths with the flag `--non-reference-ranges`

```
mkdir extract_non_reference_paths
cd extract_non_reference_paths/

odgi paths -i ../Arabidopsis.pg.input.community.0.og -L | grep 'GCA_000001735' > reference.txt

odgi paths -i ../Arabidopsis.pg.input.community.0.og --non-reference-ranges reference.txt > non_reference.bed
```

#### 2. Filter non_reference.bed

`filter_private_bed.py` script is in the scripts directory.
It takes the output from the previous command and separates sequences shorther than threshold from sequences longer than threshold. It also prints the distribution of the paths length to stdout and to the log file.


```
~/software/scripts/filter_private_bed.py -i non_reference.bed -o non_reference_filtered.bed -thr 500

---

size_category    count
            1 16361432
         2-50  3136982
        51-99   168284
      100-499    93759
      500-999    22737
    1000-4999    25706
    5000-9999     5370
  10000-19999     1905
  20000-99999     1253
100000-499999      280
500000-999999        7
    from_1Mbp        0

```

#### 3. Apply bedtools slop and merge

Bedtools needs a genome file for slop. We should use the fasta.fai of the multifasta used for building the pangenome with pggb.
I will extend regions smaller than 500 bp of 500bp in both ends.

- filtered_out.bed is the bed file containing all the regions smaller than 500 bp
- non_reference_filtered.bed is the bed file containing all the regions larger than 500 bp

I will enlarge only the sequences contained in filtered_out.bed

```
bedtools slop -i filtered_out.bed -g ../../Arabidopsis.pg.input.fasta.gz.fai -b 500 > expanded_regions.bed

head expanded_regions.bed 
---
GCA_001651475#1#CM004359.1      0       789
GCA_001651475#1#CM004359.1      0       882
GCA_001651475#1#CM004359.1      0       887
GCA_001651475#1#CM004359.1      0       912
GCA_001651475#1#CM004359.1      0       1013
GCA_001651475#1#CM004359.1      36      1037
GCA_001651475#1#CM004359.1      104     1105
GCA_001651475#1#CM004359.1      752     1753
GCA_001651475#1#CM004359.1      933     1934
GCA_001651475#1#CM004359.1      997     1998

head filtered_out.bed
---
GCA_001651475#1#CM004359.1      0       289
GCA_001651475#1#CM004359.1      381     382
GCA_001651475#1#CM004359.1      385     387
GCA_001651475#1#CM004359.1      411     412
GCA_001651475#1#CM004359.1      413     513
GCA_001651475#1#CM004359.1      536     537
GCA_001651475#1#CM004359.1      604     605
GCA_001651475#1#CM004359.1      1252    1253
GCA_001651475#1#CM004359.1      1433    1434
GCA_001651475#1#CM004359.1      1497    1498
```

Right, now the smaller regions have been expanded. Now I will join the expanded regions (expanded_regions.bed) with the regions that were already longer than 500 bp (non_reference_filtered.bed):

```
cat expanded_regions.bed non_reference_filtered.bed | bedtools sort -i - > slop_final.bed
```

Now we can apply bedtools merge to exclude overlappings and join close sequences:

```
bedtools merge -i slop_final.bed -d 100 > non_reference_chr1_for_annotation.bed
```

- -d 100 means that if there are gaps of 100bp between two sequences they are merged.

#### 4. Extract the FASTA

Now, with this bed file we should be able to extract the fasta that we want to annotate.

```
bedtools getfasta -fi ../../Arabidopsis.pg.input.fasta.gz -bed non_reference_chr1_for_annotation.bed -fo non_reference_chr1_for_annotation.fasta -name
```

```
head non_reference_chr1_for_annotation.fasta 
---
>::GCA_001651475#1#CM004359.1:0-1998
GTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTAGGGTTTTAGGTTTTAGGGTTCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCCTTAAACCCTAAACCCTAAACCCTAAATAAAGCGCTGTGGGATCAATCATTGCATTTTTCCATAGGATAGATAGGGCCGACAAGATCATCAGGGAAGAAGTCAAATCACATCCGAATTCAATTGTTCTTTTCCTAAACCCTAAACCCTAAACACTAAACCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTTATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTATGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATTTAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAAAAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGCTTGTTTTGCTTCTTTGAAGTAGTTTCTCTTTGCAAAATTCCTCTTTTTTTAGAGTGATTTGGATGATTCAAGACTTCTCGGTACTGCAAAGTTCTTCCGCCTGATTAATTATCCATTTTACCTTTGTCGTAGATATTAGGTAATCTGTAAGTCAACTCATATACAACTCATAATTTAAAAGAAAATTATGATCGACACACGTTTACACATAAAATCTGTAAATCAACTCATATACCCGTTATTCTCACAATCATATGCTTTCTAAAAGCAAAAGTATATGTCAACAATTGGTTATAAATTATTAGAAGTTTTCCACTTATGACTTAAGAACTTGTGAAGCAGAAAGTGGCAACACCCCCCACCTCCCCCCCCCCCCACCCCCCAAATTGAGAAGTCAATTTTATATAATTTAATCAAATAAATAAGTTTATGGTTAAGAGTTTTTTACTCTCTTTATTTTTCTTTTTCTTTTTGAGACATACTGAAAAAAGTTGTAATTATTAATGATAGTTCTGTGATTCCTCCATGAATCACATCTGCTTGATTTTTCTTTCATAAATTTATAAGTAATACATTCTTATAAAATGGTCAGAGAAACACCAAAGATCCCGAGATTTCTTCTCACTTACTTTTTTTCTATCTATCTAGATTATATAAATGAGATGTTGAATTAGAGGAACCTTTGATTCAATGATCATAGAAAAATTAGGTAAAGAGTCAGTGTCGTTATGTTATGGAAGATGTG
>::GCA_001651475#1#CM004359.1:6363-7364
TATATGTAATAATAATTAGTGCATCGTTTTGTGGTGTAGTTTATATAAATAAAGTGATATATAGTCTTGTATAAGAAAGGGATTTTACATGAGACCCAAATATGAGTAAAGGGTGTTGGCTCAAAGATTCATTTAGCAACCAAAGTTGCATTTGCAAGGAAATGAAAAGGTGTTAACAATGTCACTGCGTAACATGACTTCGATACAAAATCTCAAACAATACATCTTCAACTGTGGATTATGATGGACTTGGGGTTGCAGGTGATATGTCTATAGAAAAACGGGTGGGAATGGAAAATGGCTGAAGAAGAAAGCTCCAACTAAGCATCATAACTGGATTGTTTTAGAGGAGATGAGTCAAAGGAATTCACAGTGGAACTATCTCAAGAACATGATCATTGGCTTCTTATTGTTCATCTCCATCATTGGCTGGATCATTCTGGTTCAAGAGGTCAAATTATATATACATAACGGATCTAAGAAGTATAGTGTAGTCAATTAAAAACAAAACGAGACTTGAAAATAAGCATAAGTATTATTAAGTTAACCCAATATTCGTTTCAATGCTTTAGTTATCATCAGAATTCTCAACATTTTCAGATATTTAACTTGCCTTTGGTTGCTTTCTCGCCATGGTAGTAGCATTCTCCGGATAAGAATCAAGGGGAGCCTCAACTTCGGCTTCAACCGTCTCCTCGTCTTTCCTATGACATTCACTTGGTGTTGCAACAATGTGTTGATGCCCCACATTCATAGGAGAGTCCCACATATCTCCAACATTCATAGGAGGGTTCCACATATGCCCACCATTCCAAGGAGGGTACCACCACATATGCCCAGCATTCCAAGAAAGATGGAGCGACGATATCGTGAGGAACAAGGTTTTTGTAGGGCAAATAATTGTAGGTGGTTGAGTTATATCCGCCCGCACAGAGTAACAGACCAACAATTAACTTTTGATATTTTAGTAAGGTCTAATTCAATTTTTGGTGGCGATAATA
>::GCA_001651475#1#CM004359.1:8710-9736
AGCTAGAATCAGACAGGAACTAGCAATGCTTGAAATCAAGAACTTGAATTGAAATAGTTTTTTACCTGAATATTGACAGTTGCTGGATTAATTGCATTGTAGAGGACGTGTCTATATACCTTTGGTCTGTGAAGGATTAAATCGATGAAAATAATCTGCCAAAGAAAACAATTAAAGAACCAAAAACCAAAATTGGAAAGAAATAGGGAAACACCCAAAAAGGGAAAGAAAGTGATTAAAACAGACCATGCGTTCACACTCGATGTACTCATCTGCTACTTCCTTGCAATTTCCCTAAATATAACAATATGATCAAAGATGGAAACTTTGAAGAAATTTAATAGAGAATCTTATAAACCCTAATTGGGTCAAAGAAGATCCATTAATACAAAAATCTTACGCATTTCATGAGACGAATGTTACCCGGAGAGTATTGAATGAACAATGACTTTACCCTAAAACCACATCCCACGCATCTGTGTTCACTCGCCGCCATTGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCAAGAGAAGAAGAATACGGAGCAATTAGAGTCCGGGTCTGGGCTACTGTTTTAACCCTAAATGGGCTTATTCATGGGCCAAGTTTTTGAAGTCTTAACTTTAAATTTGTTAGGCCCACTTTTGCTCTAAGCCGGGGTATTTGTACCCCAAAATTTAAAAATCATATACACGTTGTAATTTATAAATAGTTCAATTTGGATCAAAATCTTGTCCATATGACATAGCATTTTAAAATGCGTAGGTTCATGAATGAAACATATTATAGGCCTCAGATAAAGATATACATATTAAGTCTAAATTATTTAGTCTTCAGAATTTACCACACTTACTGAAAAGTCTAGTGGTTCACTAATATTATTACTGTCGTGTTACTTTCTATATATAGTTCATGACTTGTGAGTTGTGATGGATAAGTTTATAAGAAAATAAATTATTTATTACAATTCAACAGTGAAGAAATTTATTTAGTTTGATTAAATAAGAAAGGTAAATAA
>::GCA_001651475#1#CM004359.1:10945-11946
TTTAACACTTAATTAACGCAATCTTACCATATATTCATGGACTACATGCAGAATAGTAATTCTCCCAACCTTTCTAGTTATTTACCTGAATGTGTTTATGTACATGGACCGGTAACCTCATGTATATATGCACATACTGACAATCTGACATACATATATATAGTAGATATGACAACAACAACAAAAAAAAAAAAAATTCCTTGTTCGTGAAGCATGATCTGAGAGTTCCTAGTTAGCATGTTGTGTGGGATCATACTTTTAATATGCTGCAAGTACCAGTCAATTTTAGTATGGGAAACTATAAACATGTATAATCAACCAATGAACACGTCAATAACCTATTGAACAGCTTAGGGTGAAAATTATGATCCGTAGAGACAGCATTTAAAAGTTCCTTACGTCCACGTAAAATAATATATCAATTTATACATATACATGTGTAAACTGTGTATATATAGGGTAGGTATATGTGTATATATATAGTAATTGACAAATGATTTTGGTTCTAACATATATTCTAAAAGTACTCATGAGTTTGTGAGATCTACACAAGATACCTGATTTGATAAAAATGGCTTCAACTTGCAATCCAAACCAAACCAAACAAAGTTAATAACCAAGGGTTAATAACAAAAACAAGAATCTAGAATTAGTAAAAAAATGAGAAATTAATGAACCTGTGATCATAAAAAAAGTCAAACAATGTGAAAACATATCATACCTTTTGTTCTTTTTAATATAATAATTGAATTACTAAATGGATGGATCAGTCCTTCTTCCATAGCTAGCTTCTCTTTATTTTCTCTGCCCATAACCTGCAGAAAATCTCTTTAACCGAGCAAAATTACAAGAGTAACCAAACAAACAAAATAGGACCATCAGAGAGAGAGAAAAGAGTGCCTTTTTTTGGACCTGCCCATTAGCTTAGAATGTACCATGAAACACTTTGTCAGTGTAGGGAGATAGGCACAGAGAGTACAATTCATGAAATTTATAAGCTT
>::GCA_001651475#1#CM004359.1:13387-14389
CTGCTGCTTGCTCCGATGTTGGAGGTTAATTCTTGGTCTCTGTCTTGTTCTTGGTCGGAACTTGTTGTTGTTGTCGGAGCCAGGGATAGATCCATTATGACTTATCTGATTTTCCTTGGAGAAATCTTGATTTACGAGATGGTTCTTGAATCTTTTTTTTTTTTTTTTTGGTTTTTGGCAAAACCTTGCTCCGTTTAGTTTTGAGTTGGCTCTTAGTGGGTCTTAGTGAATCTTGATGAAGTTTAAACCTTTCTCTTGTTCTTCAATTGGTGAATTGTTGCGGATGTTTTTCAAAGCCATTGACGATGATGATGACTCAACTCGAGATCTTCAACATACATTAGAGAGATAATATTAAAAAGAGGAAACGAAACAAAAATCTAGCCACAATTATAAAAAAGGGTAGATGAATTAGTGAAAAAAAATTCAAAGAATATATTAAAAGAAGGGAAGAAGGAAAAAGAAAATCAAAAGTTCTTGTGGTCTCACCTTATAATATAAGAGAGAGAGAGAGGGAGAGAGAGATGTTGATATTGAAGTCTTCTTTCTCTTGCTGGTTCGATTGTTTCTTAGATTTTCTTCTTTGACCATGAAAATTTAGAGGAAATATGAAAACCCTAGAATCGGAAGAAAACTATATATGTATATCTTTCCGTTGACTTTATATAGAATGAAATCAAGGAAAGAAAAGAGCTAAGCAAATCAGGACTTAGAACTCAAATTGGGTTCTTGCCAAACAAGAGGATCTCTCTCTTTCTCTTTATGAAGATGAAGAAAAAATGGTGGCTTGAGATTCTTGATAGAGTTTCAAGACTTAGAGAGAGAGAGCGAGATTGGTACCGAGTAATGCGTAGGGGCCCATGAATGATATGGTTCGGGTCATGTGAGCTTCTTGTCCTTTTCATTAGGCTTCGAATTAGCTCTTTGTATTGAAATTTCAAAATCTCTTTATACATATGTGTATGTATATGTTTAGTTGTCGACGGTGTTAATTTGAGAATC
```

This is the fasta for annotation.

remove :: from the fasta headers (they have been automatically added):

```
sed -i 's/^>::/>/' non_reference_chr1_for_annotation.fasta
```


### community 1 - chr2

#### 1. Extract non-reference sequences

- find the reference paths using odgi paths + grep 'reference_name' and place it in `reference.txt`
- extract non reference sequences using odgi paths with the flag `--non-reference-ranges`

```
mkdir extract_non_reference_paths
cd extract_non_reference_paths/

odgi paths -i ../Arabidopsis.pg.input.community.1.og -L | grep 'GCA_000001735' > reference.txt

odgi paths -i ../Arabidopsis.pg.input.community.1.og --non-reference-ranges reference.txt > non_reference.bed
```

#### 2. Filter non_reference.bed

`filter_private_bed.py` script is in the scripts folder.
It takes the output from the previous command and separates sequences shorther than a given threshold from longer sequences. It also prints the distribution of the paths length to stdout and to the log file.

```
~/software/scripts/filter_private_bed.py -i non_reference.bed -o non_reference_filtered.bed -thr 500

---
size_category    count
            1 11130257
         2-50  2121266
        51-99   119889
      100-499    58584
      500-999    15457
    1000-4999    19663
    5000-9999     4421
  10000-19999     1694
  20000-99999     1209
100000-499999      256
500000-999999       11
    from_1Mbp        0
```

#### 3. Apply bedtools slop and merge

Bedtools needs a genome file for slop. We can use the fasta.fai of the multifasta used for building the pangenome with pggb.
I will extend regions smaller than 500 bp of 500bp in both ends.

- filtered_out.bed is the bed file containing all the regions smaller than 500 bp
- non_reference_filtered.bed is the bed file containing all the regions larger than 500 bp

I will enlarge only the sequences contained in filtered_out.bed

```
bedtools slop -i filtered_out.bed -g ../../Arabidopsis.pg.input.fasta.gz.fai -b 500 > expanded_regions.bed

head expanded_regions.bed 
---
GCA_001651475#1#CM004360.1      0       614
GCA_001651475#1#CM004360.1      0       703
GCA_001651475#1#CM004360.1      0       780
GCA_001651475#1#CM004360.1      0       796
GCA_001651475#1#CM004360.1      0       966
GCA_001651475#1#CM004360.1      284     1285
GCA_001651475#1#CM004360.1      516     1517
GCA_001651475#1#CM004360.1      564     1565
GCA_001651475#1#CM004360.1      674     1675
GCA_001651475#1#CM004360.1      972     1973

head filtered_out.bed
---
GCA_001651475#1#CM004360.1      113     114
GCA_001651475#1#CM004360.1      202     203
GCA_001651475#1#CM004360.1      279     280
GCA_001651475#1#CM004360.1      295     296
GCA_001651475#1#CM004360.1      465     466
GCA_001651475#1#CM004360.1      784     785
GCA_001651475#1#CM004360.1      1016    1017
GCA_001651475#1#CM004360.1      1064    1065
GCA_001651475#1#CM004360.1      1174    1175
GCA_001651475#1#CM004360.1      1472    1473
```

Right, now the smaller regions have been expanded. Now I will join the expanded regions (expanded_regions.bed) with the regions that were already above 500 bp (non_reference_filtered.bed):

```
cat expanded_regions.bed non_reference_filtered.bed | bedtools sort -i - > slop_final.bed
```

Now we can apply bedtools merge to exclude overlappings and join close sequences:

```
bedtools merge -i slop_final.bed -d 100 > non_reference_chr2_for_annotation.bed
```

- -d 100 means that if there are gaps of 100bp between two sequences they are merged.

#### 4. Extract the FASTA

Now, with this bed file we should be able to extract the fasta that we want to annotate.

```
bedtools getfasta -fi ../../Arabidopsis.pg.input.fasta.gz -bed non_reference_chr2_for_annotation.bed -fo non_reference_chr2_for_annotation.fasta -name
```

remove :: from the fasta headers (they have been automatically added):

```
sed -i 's/^>::/>/' non_reference_chr2_for_annotation.fasta
```

### community 2 - chr3

#### 1. Extract non-reference sequences

- find the reference paths using odgi paths + grep 'reference_name' and place it in `reference.txt`
- extract non reference sequences using odgi paths with the flag `--non-reference-ranges`

```
mkdir extract_non_reference_paths
cd extract_non_reference_paths/

odgi paths -i ../Arabidopsis.pg.input.community.2.og -L | grep 'GCA_000001735' > reference.txt

odgi paths -i ../Arabidopsis.pg.input.community.2.og --non-reference-ranges reference.txt > non_reference.bed
```

#### 2. Filter non_reference.bed

`filter_private_bed.py` script is in the scripts folder.
It takes the output from the previous command and separates sequences shorther than a given threshold from longer sequences. It also prints the distribution of the paths length to stdout and to the log file.

```
~/software/scripts/filter_private_bed.py -i non_reference.bed -o non_reference_filtered.bed -thr 500

---
size_category    count
            1 12646306
         2-50  2437703
        51-99   107639
      100-499    59027
      500-999    17275
    1000-4999    23841
    5000-9999     5898
  10000-19999     2585
  20000-99999     1596
100000-499999      324
500000-999999       12
    from_1Mbp        1
```

#### 3. Apply bedtools slop and merge

Bedtools needs a genome file for slop. We can use the fasta.fai of the multifasta used for building the pangenome with pggb.

I will extend regions smaller than 500 bp of 500bp in both ends.

- filtered_out.bed is the bed file containing all the regions smaller than 500 bp
- non_reference_filtered.bed is the bed file containing all the regions larger than 500 bp

I will enlarge only the sequences contained in filtered_out.bed

```
bedtools slop -i filtered_out.bed -g ../../Arabidopsis.pg.input.fasta.gz.fai -b 500 > expanded_regions.bed

head expanded_regions.bed 
---
GCA_001651475#1#CM004361.1      0       594
GCA_001651475#1#CM004361.1      0       677
GCA_001651475#1#CM004361.1      0       734
GCA_001651475#1#CM004361.1      0       742
GCA_001651475#1#CM004361.1      0       800
GCA_001651475#1#CM004361.1      0       808
GCA_001651475#1#CM004361.1      0       816
GCA_001651475#1#CM004361.1      0       823
GCA_001651475#1#CM004361.1      0       837
GCA_001651475#1#CM004361.1      0       841

head filtered_out.bed
---
GCA_001651475#1#CM004361.1      92      94
GCA_001651475#1#CM004361.1      168     177
GCA_001651475#1#CM004361.1      233     234
GCA_001651475#1#CM004361.1      241     242
GCA_001651475#1#CM004361.1      292     300
GCA_001651475#1#CM004361.1      306     308
GCA_001651475#1#CM004361.1      315     316
GCA_001651475#1#CM004361.1      322     323
GCA_001651475#1#CM004361.1      336     337
GCA_001651475#1#CM004361.1      340     341
```

Right, now the smaller regions have been expanded. Now I will join the expanded regions (expanded_regions.bed) with the regions that were already above 500 bp (non_reference_filtered.bed):

```
cat expanded_regions.bed non_reference_filtered.bed | bedtools sort -i - > slop_final.bed
```

Now we can apply bedtools merge to exclude overlappings and join close sequences:

```
bedtools merge -i slop_final.bed -d 100 > non_reference_chr3_for_annotation.bed
```

- -d 100 means that if there are gaps of 100bp between two sequences they are merged.

#### 4. Extract the FASTA

Now, with this bed file we should be able to extract the fasta that we want to annotate.

```
bedtools getfasta -fi ../../Arabidopsis.pg.input.fasta.gz -bed non_reference_chr3_for_annotation.bed -fo non_reference_chr3_for_annotation.fasta -name
```

remove :: from the fasta headers (they have been automatically added):

```
sed -i 's/^>::/>/' non_reference_chr3_for_annotation.fasta
```

### community 3 - chr4

#### 1. Extract non-reference sequences

- find the reference paths using odgi paths + grep 'reference_name' and place it in `reference.txt`
- extract non reference sequences using odgi paths with the flag `--non-reference-ranges`

```
mkdir extract_non_reference_paths
cd extract_non_reference_paths/

odgi paths -i ../Arabidopsis.pg.input.community.3.og -L | grep 'GCA_000001735' > reference.txt

odgi paths -i ../Arabidopsis.pg.input.community.3.og --non-reference-ranges reference.txt > non_reference.bed
```

#### 2. Filter non_reference.bed

`filter_private_bed.py` script is in the scripts folder.
It takes the output from the previous command and separates sequences shorther than a given threshold from longer sequences. It also prints the distribution of the paths length to stdout and to the log file.

```
~/software/scripts/filter_private_bed.py -i non_reference.bed -o non_reference_filtered.bed -thr 500

---
size_category    count
            1 11148357
         2-50  2200442
        51-99   112329
      100-499    70063
      500-999    15581
    1000-4999    19389
    5000-9999     4439
  10000-19999     1948
  20000-99999     1331
100000-499999      392
500000-999999       25
    from_1Mbp        2
```

#### 3. Apply bedtools slop and merge

Bedtools needs a genome file for slop. We can use the fasta.fai of the multifasta used for building the pangenome with pggb.
I will extend regions smaller than 500 bp of 500bp in both ends.

- filtered_out.bed is the bed file containing all the regions smaller than 500 bp
- non_reference_filtered.bed is the bed file containing all the regions larger than 500 bp

I will enlarge only the sequences contained in filtered_out.bed

```
bedtools slop -i filtered_out.bed -g ../../Arabidopsis.pg.input.fasta.gz.fai -b 500 > expanded_regions.bed

head expanded_regions.bed 
---
GCA_001651475#1#CM004362.1      4567    5570
GCA_001651475#1#CM004362.1      4571    5574
GCA_001651475#1#CM004362.1      4575    5584
GCA_001651475#1#CM004362.1      4586    5612
GCA_001651475#1#CM004362.1      4614    5616
GCA_001651475#1#CM004362.1      4617    5619
GCA_001651475#1#CM004362.1      4620    5623
GCA_001651475#1#CM004362.1      4624    5628
GCA_001651475#1#CM004362.1      4630    5662
GCA_001651475#1#CM004362.1      4673    5674

head filtered_out.bed
---
GCA_001651475#1#CM004362.1      5067    5070
GCA_001651475#1#CM004362.1      5071    5074
GCA_001651475#1#CM004362.1      5075    5084
GCA_001651475#1#CM004362.1      5086    5112
GCA_001651475#1#CM004362.1      5114    5116
GCA_001651475#1#CM004362.1      5117    5119
GCA_001651475#1#CM004362.1      5120    5123
GCA_001651475#1#CM004362.1      5124    5128
GCA_001651475#1#CM004362.1      5130    5162
GCA_001651475#1#CM004362.1      5173    5174
```

Right, now the smaller regions have been expanded. Now I will join the expanded regions (expanded_regions.bed) with the regions that were already above 500 bp (non_reference_filtered.bed):

```
cat expanded_regions.bed non_reference_filtered.bed | bedtools sort -i - > slop_final.bed
```

Now we can apply bedtools merge to exclude overlappings and join close sequences:

```
bedtools merge -i slop_final.bed -d 100 > non_reference_chr4_for_annotation.bed
```

- -d 100 means that if there are gaps of 100bp between two sequences they are merged.

#### 4. Extract FASTA

Now, with this bed file we should be able to extract the fasta that we want to annotate.

```
bedtools getfasta -fi ../../Arabidopsis.pg.input.fasta.gz -bed non_reference_chr4_for_annotation.bed -fo non_reference_chr4_for_annotation.fasta -name
```

remove :: from the fasta headers (they have been automatically added):

```
sed -i 's/^>::/>/' non_reference_chr4_for_annotation.fasta
```

### community 4 - chr5

#### 1. Extract non-reference sequences

- find the reference paths using odgi paths + grep 'reference_name' and place it in `reference.txt`
- extract non reference sequences using odgi paths with the flag `--non-reference-ranges`

```
mkdir extract_non_reference_paths
cd extract_non_reference_paths/

odgi paths -i ../Arabidopsis.pg.input.community.4.og -L | grep 'GCA_000001735' > reference.txt

odgi paths -i ../Arabidopsis.pg.input.community.4.og --non-reference-ranges reference.txt > non_reference.bed
```

#### 2. Filter non_reference.bed

`filter_private_bed.py` script is in the scripts folder.
It takes the output from the previous command and separates sequences shorther than a given threshold from longer sequences. It also prints the distribution of the paths length to stdout and to the log file.

```
~/software/scripts/filter_private_bed.py -i non_reference.bed -o non_reference_filtered.bed -thr 500

---
size_category    count
            1 15088021
         2-50  2535317
        51-99    98407
      100-499   158317
      500-999    21568
    1000-4999    23741
    5000-9999     6136
  10000-19999     2140
  20000-99999     1024
100000-499999      180
500000-999999       10
    from_1Mbp        1
```

#### 3. Apply bedtools slop and merge

Bedtools needs a genome file for slop. We can use the fasta.fai of the multifasta used for building the pangenome with pggb.
I will extend regions smaller than 500 bp of 500bp in both ends.

- filtered_out.bed is the bed file containing all the regions smaller than 500 bp
- non_reference_filtered.bed is the bed file containing all the regions larger than 500 bp

I will enlarge only the sequences contained in filtered_out.bed

```
bedtools slop -i filtered_out.bed -g ../../Arabidopsis.pg.input.fasta.gz.fai -b 500 > expanded_regions.bed

head expanded_regions.bed 
---
GCA_001651475#1#CM004363.1      0       501
GCA_001651475#1#CM004363.1      0       503
GCA_001651475#1#CM004363.1      0       508
GCA_001651475#1#CM004363.1      0       510
GCA_001651475#1#CM004363.1      0       515
GCA_001651475#1#CM004363.1      0       517
GCA_001651475#1#CM004363.1      0       522
GCA_001651475#1#CM004363.1      0       524
GCA_001651475#1#CM004363.1      0       529
GCA_001651475#1#CM004363.1      0       531

head filtered_out.bed
---
GCA_001651475#1#CM004363.1      0       1
GCA_001651475#1#CM004363.1      2       3
GCA_001651475#1#CM004363.1      6       8
GCA_001651475#1#CM004363.1      9       10
GCA_001651475#1#CM004363.1      14      15
GCA_001651475#1#CM004363.1      16      17
GCA_001651475#1#CM004363.1      21      22
GCA_001651475#1#CM004363.1      23      24
GCA_001651475#1#CM004363.1      28      29
GCA_001651475#1#CM004363.1      30      31
```

Right, now the smaller regions have been expanded. Now I will join the expanded regions (expanded_regions.bed) with the regions that were already above 500 bp (non_reference_filtered.bed):

```
cat expanded_regions.bed non_reference_filtered.bed | bedtools sort -i - > slop_final.bed
```

Now we can apply bedtools merge to exclude overlappings and join close sequences:

```
bedtools merge -i slop_final.bed -d 100 > non_reference_chr5_for_annotation.bed
```

- -d 100 means that if there are gaps of 100bp between two sequences they are merged.

#### 4. Extract the FASTA

Now, with this bed file we should be able to extract the fasta that we want to annotate.

```
bedtools getfasta -fi ../../Arabidopsis.pg.input.fasta.gz -bed non_reference_chr5_for_annotation.bed -fo non_reference_chr5_for_annotation.fasta -name
```

remove :: from the fasta headers (they have been automatically added):

```
sed -i 's/^>::/>/' non_reference_chr5_for_annotation.fasta
```

## Annotation

### Get UniRef_100 and format as diamond database

```
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz

#get taxonomy files 

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
```

```
# format db

diamond makedb --in uniref100.fasta --taxonmap prot.accession2taxid.FULL.gz --taxonnodes nodes.dmp --taxonnames names.dmp -d uniref100 --threads 20

```

### Download conversion table between UniProt names and Araport

We need to download the conversion table between UniProt names and Araport from TAIR10 website: https://www.arabidopsis.org/download/file?path=Proteins/Id_conversions/TAIR2UniprotMapping.txt.We will need this table to use the `parse_exonerate_results.py` script.

### Extract attributes from UniRef100 FASTA headers

For the same script, we also need to extract informations from the uniref100 headers:

```
grep '>' uniref100.fasta > headers_uniref100
sed -i 's/^>//' headers_uniref100
awk '{split($0, a, " "); printf "%s\t", a[1]; for (i=2; i<=NF; i++) printf "%s%s", $i, (i==NF ? "" : "_"); print ""}' headers_uniref100 > info_headers.txt
rm headers_uniref100
```

### community 0 - chr1

#### 1. Diamond

```
diamond blastx -d uniref100 -q non_reference_chr1_for_annotation.fasta -o diamond_chr1_top0_results.tsv --range-culling --top 0 -F 15 -f 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames sskingdoms skingdoms sphylums -t ./ -p 50 --log
```

#### 2. Exonerate

We want to refine the annotations found by diamond blastx, and not repeat the whole analysis with exonerate.
To do so, we should prepare some input files and run exonerate using the 'PW_exonerate.pl' script.

In order to speed up the analysis, we are going to split the diamond results in chunks, and with the help of seqtk we will create the protein and sequence fasta files, in order to parallelise exonerate (which doesn't allow multithreading).

For each chunk we need to create:

- a file containing the whole sequences where diamond found proteins
- a file containing the proteins that matched with our sequences from the database.

To automate this step I created the `split_and_extract_for_exonerate.sh` script.
It splits diamond results in chunks and for each chunks seqtk extracts the fasta of the proteins and the fasta of the sequences.

##### split diamond results in chunks

Here's how we can use it:
```
~/software/scripts/split_and_extract_for_exonerate.sh -i ../diamond_chr1_top0_results.tsv -l 100000 -s ../non_reference_chr1_for_annotation.fasta -p uniref100.fasta
```

##### run exonerate

Then, we can run exonerate in parallel on the different chunks with the following command for each chunk (change suffix for each chunk):

```
perl ~/software/scripts/PW_exonerate.pl --protein_fasta diamond_chunk_aa_proteins.fasta --genome_fasta diamond_chunk_aa_seqs.fasta --diamond_results diamond_chunk_aa --suffix aa

# join all the results together

cat exonerate_results_aa.gff exonerate_results_ab.gff exonerate_results_ac.gff exonerate_results_ad.gff exonerate_results_ae.gff exonerate_results_af.gff exonerate_results_ag.gff exonerate_results_ah.gff exonerate_results_ai.gff exonerate_results_aj.gff exonerate_results_ak.gff > chr1_exonerate_results.gff

# remove chunks to save space

rm exonerate_results_*
```

#### 3. Parse exonerate results

To parse the results of exonerate we are going to use `parse_exonerate.py`.

```
# extract C4 alignments

awk '/C4 Alignment:/,/vulgar/{if (!/vulgar/) print}' chr1_exonerate_results.gff > chr1_C4_alignments.txt

# Parse the results
 
~/software/scripts/parse_exonerate.py -e chr1_exonerate_results.gff -c chr1_C4_alignments.txt -t TAIR2UniprotMapping.txt -i info_headers.txt
```

### community 1 - chr2

#### 1. Diamond

```
diamond blastx -d uniref100 -q non_reference_chr2_for_annotation.fasta -o diamond_chr2_top0_results.tsv --range-culling --top 0 -F 15 -f 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames sskingdoms skingdoms sphylums -t ./ -p 50 --log
```

#### 2. Exonerate 

##### split diamond results in chunks
```
~/software/scripts/split_and_extract_for_exonerate.sh -i ../diamond_chr2_top0_results.tsv -l 100000 -s ../non_reference_chr2_for_annotation.fasta -p uniref100.fasta
```

##### run exonerate

Then, we can run exonerate in parallel on the different chunks with the following command for each chunk (changing suffix for each chunk):
```
perl ~/software/scripts/PW_exonerate.pl --protein_fasta diamond_chunk_aa_proteins.fasta --genome_fasta diamond_chunk_aa_seqs.fasta --diamond_results diamond_chunk_aa --suffix aa

cat exonerate_results_aa.gff exonerate_results_ab.gff exonerate_results_ac.gff exonerate_results_ad.gff exonerate_results_af.gff exonerate_results_ag.gff exonerate_results_ah.gff > chr2_exonerate_results.gff

rm exonerate_results_*
```

#### 3. Parse exonerate results

```
awk '/C4 Alignment:/,/vulgar/{if (!/vulgar/) print}' chr2_exonerate_results.gff > chr2_C4_alignments.txt

~/software/scripts/parse_exonerate.py -e chr2_exonerate_results.gff -c chr2_C4_alignments.txt -t ~/Arabidopsis/pangenome/ara_pan_from_April/TAIR2UniprotMapping.txt -i info_headers.txt
```


### community 2 - chr3

#### 1. Diamond

```
diamond blastx -d uniref100 -q non_reference_chr3_for_annotation.fasta -o diamond_chr3_top0_results.tsv --range-culling --top 0 -F 15 -f 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames sskingdoms skingdoms sphylums -t ./ -p 50 --log
```

#### 2. Exonerate 

##### split diamond results in chunks

```
~/software/scripts/split_and_extract_for_exonerate.sh -i ../diamond_chr3_top0_results.tsv -l 100000 -s ../non_reference_chr3_for_annotation.fasta -p uniref100.fasta
```

##### run exonerate

Then, we can run exonerate in parallel on the different chunks with the following command for each chunk (change suffix for each chunk):
```
perl ~/software/scripts/PW_exonerate.pl --protein_fasta diamond_chunk_aa_proteins.fasta --genome_fasta diamond_chunk_aa_seqs.fasta --diamond_results diamond_chunk_aa --suffix aa

cat exonerate_results_aa.gff exonerate_results_ab.gff exonerate_results_ac.gff exonerate_results_ad.gff exonerate_results_ae.gff exonerate_results_af.gff exonerate_results_ag.gff exonerate_results_ah.gff exonerate_results_ai.gff > chr3_exonerate_results.gff

rm exonerate_results_*
```

#### 3. Parse exonerate results

```

awk '/C4 Alignment:/,/vulgar/{if (!/vulgar/) print}' chr3_exonerate_results.gff > chr3_C4_alignments.txt

~/software/scripts/parse_exonerate.py -e chr3_exonerate_results.gff -c chr3_C4_alignments.txt -t TAIR2UniprotMapping.txt -i info_headers.txt
```

### community 3 - chr4

#### 1. Diamond

```
diamond blastx -d uniref100 -q non_reference_chr4_for_annotation.fasta -o diamond_chr4_top0_results.tsv --range-culling --top 0 -F 15 -f 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames sskingdoms skingdoms sphylums -t ./ -p 50 --log
```

#### 2. Exonerate

##### split diamond results in chunks
```
~/software/scripts/split_and_extract_for_exonerate.sh -i ../diamond_chr4_top0_results.tsv -l 100000 -s ../non_reference_chr4_for_annotation.fasta -p uniref100.fasta
```

##### run exonerate

Then, we can run exonerate in parallel on the different chunks with the following command for each chunk (changing suffix for each chunk):

```
perl ~/software/scripts/PW_exonerate.pl --protein_fasta diamond_chunk_aa_proteins.fasta --genome_fasta diamond_chunk_aa_seqs.fasta --diamond_results diamond_chunk_aa --suffix aa

cat exonerate_results_aa.gff exonerate_results_ab.gff exonerate_results_ac.gff exonerate_results_ad.gff exonerate_results_ae.gff exonerate_results_af.gff exonerate_results_ag.gff > chr4_exonerate_results.gff

rm exonerate_results_*
```

#### 3. Parse exonerate results

```
awk '/C4 Alignment:/,/vulgar/{if (!/vulgar/) print}' chr4_exonerate_results.gff > chr4_C4_alignments.txt

~/software/scripts/parse_exonerate.py -e chr4_exonerate_results.gff -c chr4_C4_alignments.txt -t TAIR2UniprotMapping.txt -i info_headers.txt
```

### community 4 - chr5

#### 1. Diamond

```
diamond blastx -d uniref100 -q non_reference_chr5_for_annotation.fasta -o diamond_chr5_top0_results.tsv --range-culling --top 0 -F 15 -f 6  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen staxids sscinames sskingdoms skingdoms sphylums -t ./ -p 50 --log
```

#### 2. Exonerate

##### split
```
~/software/scripts/split_and_extract_for_exonerate.sh -i ../diamond_chr5_top0_results.tsv -l 100000 -s ../non_reference_chr5_for_annotation.fasta -p uniref100.fasta
```

##### run exonerate

Then, we can run exonerate in parallel on the different chunks with the following command for each chunk (changing suffix for each chunk):

```
perl ~/software/scripts/PW_exonerate.pl --protein_fasta diamond_chunk_aa_proteins.fasta --genome_fasta diamond_chunk_aa_seqs.fasta --diamond_results diamond_chunk_aa --suffix aa

cat exonerate_results_aa.gff exonerate_results_ab.gff exonerate_results_ac.gff exonerate_results_ad.gff exonerate_results_ae.gff exonerate_results_af.gff exonerate_results_ag.gff exonerate_results_ah.gff exonerate_results_ai.gff exonerate_results_aj.gff > chr5_exonerate_results.gff

rm exonerate_results_*
```

#### 3. Parse exonerate results

```
awk '/C4 Alignment:/,/vulgar/{if (!/vulgar/) print}' chr5_exonerate_results.gff > chr5_C4_alignments.txt

~/software/scripts/parse_exonerate.py -e chr5_exonerate_results.gff -c chr5_C4_alignments.txt -t TAIR2UniprotMapping.txt -i info_headers.txt
```

## Annotation review

These are the files of the results that we get from `parse_exonerate.py`:

- `corresponding_to_tair.tsv` --> these are the UniRef100 annotations that have surely a correspondent TAIR10 annotation. From these, we will need to exclude UniRef genes/pseudogenes which have not a unique correspondent TAIR locus. This is needed as itg may lead to uncorrect calculations of statistics and downstream analysis.
- `new_results_not_tair.tsv` --> These are the "new genes" that we found with the analysis, meaning that they don't have a correspondent in TAIR10 annotations. To have very robust results, we need to filter them to keep only the results corresponding to reviewed proteins (i.e. those coming from UniProtKB/Swiss-Prot), as UniRef also contains proteins which were not reviewed. 


### Review results that have a corresponding TAIR annotation

Since the gene names were grouped by UniProt ID, and the corresponding TAIR names were concatenated in a single field separated by a `,` when more than one, we can simply remove those lines which contain a comma in the field containing TAIR names:

```
# chr1 

awk -F'\t' '$26 !~ /,/' corresponding_to_tair.tsv > reviewed_corresponding_to_tair.tsv


# chr2

awk -F'\t' '$26 !~ /,/' corresponding_to_tair.tsv > reviewed_corresponding_to_tair.tsv


# chr3

awk -F'\t' '$26 !~ /,/' corresponding_to_tair.tsv > reviewed_corresponding_to_tair.tsv


# chr4

awk -F'\t' '$26 !~ /,/' corresponding_to_tair.tsv > reviewed_corresponding_to_tair.tsv

# chr5

awk -F'\t' '$26 !~ /,/' corresponding_to_tair.tsv > reviewed_corresponding_to_tair.tsv
```

### Review results with no corresponding TAIR annotation

We need to review `new_results_not_tair.tsv` by looking in uniprot if they have been confirmed or not. 

We will proceed in this way:

- extraction of the unique UniRef IDs from `new_results_not_tair.tsv` (for genes and pseudogenes, unclassified will be excluded)
- go to [UniProt ID mapping](https://www.uniprot.org/id-mapping)
- select from UniRef100 to UniProtKB/Swiss-Prot --> this will give us automatically only reviewed results 
- upload IDs list, submit the job
- costumise results (TSV format) to get also Gene Ontology info (we will need this later)
- Be sure to include info about proteins review (it's under miscellaneous options)
- include also taxonomic lineage 


#### community 0 - chr1

```
grep 'gene' new_results_not_tair.tsv | cut -f 1 | sort | uniq > IDs_chr1.txt

```

The script `review_exonerate_results.py` will drop all the results coming from proteins that were not reviewed.

```
~/software/scripts/review_exonerate_results.py -i chr1_reviewed_proteins.tsv -a new_results_not_tair.tsv -o reviewed_not_tair_results.tsv
```




#### community 1 - chr2

```
grep 'gene' new_results_not_tair.tsv | cut -f 1 | sort | uniq > IDs_chr2.txt


~/software/scripts/review_exonerate_results.py -i chr2_reviewed_proteins.tsv -a new_results_not_tair.tsv -o reviewed_not_tair_results.tsv
```


#### community 2 - chr3

```
grep 'gene' new_results_not_tair.tsv | cut -f 1 | sort | uniq > IDs_chr3.txt

~/software/scripts/review_exonerate_results.py -i chr3_reviewed_proteins.tsv -a new_results_not_tair.tsv -o reviewed_not_tair_results.tsv

```

#### community 3 - chr4

```
grep 'gene' new_results_not_tair.tsv | cut -f 1 | sort | uniq > IDs_chr4.txt

~/software/scripts/review_exonerate_results.py -i chr4_reviewed_proteins.tsv -a new_results_not_tair.tsv -o reviewed_not_tair_results.tsv

```

#### community 4 - chr5

```
grep 'gene' new_results_not_tair.tsv | cut -f 1 | sort | uniq > IDs_chr5.txt

~/software/scripts/review_exonerate_results.py -i chr5_reviewed_proteins.tsv -a new_results_not_tair.tsv -o reviewed_not_tair_results.tsv

```

# Final annotation screening

Now, we can screen the whole pangenome, using all the data coming from `exonerate` and `odgi untangle`.

## genes

We'll need the following files:

- `filter_passed_features.csv` -->  obtained from `core_dispensable_genes.py`
- `pangenome_screening.csv` -->  obtained from `core_dispensable_genes.py`
- `reviewed_corresponding_to_tair.tsv` --> obtained from `parse_exonerate.py` and reviewed with `review_exonerate_results.py`
- `reviewed_not_tair_results.tsv` --> obtained from `parse_exonerate.py` and reviewed with `review_exonerate_results.py`

### community 0 - chr1

```
~/software/scripts/final_genes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```

### community 1 - chr2

```
~/software/scripts/final_genes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


### community 2 - chr3

```
~/software/scripts/final_genes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


### community 3 - chr4

```
~/software/scripts/final_genes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


### community 4 - chr5

```
~/software/scripts/final_genes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


## pseudogenes

We'll need the following files:

- `filter_passed_features.csv` -->  obtained from `core_dispensable_pseudogenes.py`
- `pangenome_screening.csv` -->  obtained from `core_dispensable_pseudogenes.py`
- `reviewed_corresponding_to_tair.tsv` --> obtained from `parse_exonerate.py` and reviewed with `review_exonerate_results.py`
- `reviewed_not_tair_results.tsv` --> obtained from `parse_exonerate.py` and reviewed with `review_exonerate_results.py`


### community 0 - chr1

```
~/software/scripts/final_pseudogenes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


### community 1 - chr2 

```
~/software/scripts/final_pseudogenes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```

### community 2 - chr3

```
~/software/scripts/final_pseudogenes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```

### community 3 - chr4

```
~/software/scripts/final_pseudogenes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```


### chr5

```
~/software/scripts/final_pseudogenes_screening.py -u filter_passed_features.csv -p pangenome_screening.csv -t reviewed_corresponding_to_tair.tsv -n reviewed_not_tair_results.tsv -a 93 -s 75
```

# Gene ontology enrichment analysis

The universe for A. thaliana was obtained from: https://www.arabidopsis.org/download/list?dir=GO_and_PO_Annotations

The analysis was then performed on using TopGO as shown in `scripts/TopGO.R`.



# Sequence-based pangenome analysis

## Graphs characteristics

To obtain information about the structure of each graph, we can use the command `odgi stats`:

```
odgi stats -i Arabidopsis.pg.input.community.0.og -m > chr1_odgi_stats.yaml

odgi stats -i Arabidopsis.pg.input.community.1.og -m > chr2_odgi_stats.yaml

odgi stats -i Arabidopsis.pg.input.community.2.og -m > chr3_odgi_stats.yaml 

odgi stats -i Arabidopsis.pg.input.community.3.og -m > chr4_odgi_stats.yaml

odgi stats -i Arabidopsis.pg.input.community.4.og -m > chr5_odgi_stats.yaml
```

## Obtain node length

We can use `odgi paths --coverage-levels` as in the following example:

```
odgi paths -i Arabidopsis.pg.input.community.0.og --coverage-levels 2,75,93 > chr1_nodes_classification.txt
```

`odgi paths --coverage-levels` classifies the nodes in genomic classes. In our case, we couldn't use this script for classification because the number of paths did not correspond to the number of assemblies included in the pangenome.



## Node coverage matrices 

```
odgi paths -i Arabidopsis.pg.input.community.0.og -H -D "#" -p1 -t 100 > chr1_nodes_composition.txt

odgi paths -i Arabidopsis.pg.input.community.1.og -H -D "#" -p1 -t 100 > chr2_nodes_composition.txt

odgi paths -i Arabidopsis.pg.input.community.2.og -H -D "#" -p1 -t 100 > chr3_nodes_composition.txt

odgi paths -i Arabidopsis.pg.input.community.3.og -H -D "#" -p1 -t 100 > chr4_nodes_composition.txt

odgi paths -i Arabidopsis.pg.input.community.4.og -H -D "#" -p1 -t 100 > chr5_nodes_composition.txt
```

#### Matrices processing

```
~/software/scripts/node_matrix_processing.py -i chr1_nodes_composition.txt -p chr1 -o chr1_matrix.csv -na 93 -sl 75

~/software/scripts/node_matrix_processing.py -i chr2_nodes_composition.txt -p chr2 -o chr2_matrix.csv -na 93 -sl 75

~/software/scripts/node_matrix_processing.py -i chr3_nodes_composition.txt -p chr3 -o chr3_matrix.csv -na 93 -sl 75

~/software/scripts/node_matrix_processing.py -i chr4_nodes_composition.txt -p chr4 -o chr4_matrix.csv -na 93 -sl 75

~/software/scripts/node_matrix_processing.py -i chr5_nodes_composition.txt -p chr5 -o chr5_matrix.csv -na 93 -sl 75
```

## Node-based similarity

We will join the graph with `odgi squeeze` and subsequently we will run `odgi similarity` to obtain a similarity matrix.

`odgi squeeze` needs a file with the list of the graphs as input

```
nano graphs.txt

Arabidopsis.pg.input.community.0.og
Arabidopsis.pg.input.community.1.og
Arabidopsis.pg.input.community.2.og
Arabidopsis.pg.input.community.3.og
Arabidopsis.pg.input.community.4.og
```

Squeeze the graphs:

```
odgi squeeze -f graphs.txt -o squeezed_pangenome.og -t 100 -P
```

Obtain the matrix:

```
odgi similarity -i squeezed_pangenome.og -D '#' -d -t 100 -P > similarity_nodes.tsv
```

## Growth curve simulation

```
# get dataset from gene matrix

perl shuffle_count.pl --file concatenated_matrix.csv --iterations 10 --campionamenti 10-20-30-40-50-60-70-80-90-93 --mode NOT_REI --output out_shuffling.txt

# get the figure

python3 curve_plot.py out_shuffling.txt # for the zoom of the "All" curve
python3 curve_plot2.py out_shuffling.txt
```


# Similarity analysis

## Gene based

Join not transposed matrices

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/all

mkdir gene_cluster
cd gene_cluster

ln -s ../../chr1/matrix_chr1.txt
ln -s ../../chr2/matrix_chr2.txt
ln -s ../../chr3/matrix_chr3.txt
ln -s ../../chr4/matrix_chr4.txt
ln -s ../../chr5/matrix_chr5.txt

./combine_notransp_matrices.py

mv concatenated_matrix.csv genes_concatenated_matrix.csv

cp genes_concatenated_matrix.csv ~/Arabidopsis_pangenome/R_plots/Cluster/
```

### make the figure

```
R

library("ggplot2")
library("vegan")
library(viridis)
library(ggtree)
library(ape)

setwd("/home/lia/Arabidopsis_pangenome/R_plots")

# Load data
my_data <- read.csv("Cluster/genes_concatenated_matrix.csv", header = TRUE, sep = "\t")
head(my_data)

# Extract assembly names
assembly_names <- colnames(my_data)[2:94]

# Load locations
locations <- read.csv("genebank_country_1.csv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(locations)

colnames(locations) <- c("Assembly", "Location") 

country_to_colour <- c(
'Afghanistan' = "salmon2",
'Belgium' = "gold3",
'Cape Verde' = "darkgreen",
'Czech Republic' = "grey46",
'France' = "red", 
'Germany' = "magenta",
'Italy' = "brown",
'Japan' = "springgreen4",
'Lithuania' = "black",
'Netherlands' = "green",
'Poland' = "orange",
'Madeira' = "darkviolet",
'Romania' = "pink",
'Spain' = 'deepskyblue3',
'Sweden' = "olivedrab",
'Tanzania' = 'blue',
'United Kingdom' = 'mediumpurple',
'USA' = 'cyan')


# Create a mapping from assembly names to countries
assembly_to_country <- setNames(locations$Location, locations$Assembly)

# Map assembly names to colors through their countries
assembly_to_color <- sapply(assembly_names, function(assembly) {
  country <- assembly_to_country[assembly]
  color <- country_to_colour[country]
  return(color)
})


# Calculate Jaccard distance
dist.jaccard <- vegdist(t(my_data[, 2:94]), method = "jaccard")
write.csv(as.matrix(dist.jaccard), "genes_PAV_jaccard_distance_matrix.csv", row.names = TRUE)

# Perform hierarchical clustering
res.hc <- hclust(d = dist.jaccard, method = "ward.D2")

# Calculate the cophenetic distance matrix from the clustering result
res.coph<-cophenetic(res.hc)

# Compute the cophenetic correlation
cophenetic_corr <- cor(as.dist(dist.jaccard), as.dist(res.coph))

# Print cophenetic correlation to check the quality of the clustering
print(cophenetic_corr)

phylo_tree <- as.phylo(res.hc)

# Prepare a data frame that ggtree can use to map colors to tips
# The row names of 'my_data' are assumed to be the assembly names used in 'phylo_tree$tip.label'
# Ensure these names match those in 'assembly_to_color'
tip_colors <- data.frame(label = names(assembly_to_color), color = assembly_to_color)


# Plot the tree with colored tips
p <- ggtree(phylo_tree, layout="circular") +
  geom_tiplab(aes(color = label), size = 4) +
  scale_color_manual(values = tip_colors$color) +
  theme_tree2() +
  theme(legend.position = "none",  # Adjust legend position
        plot.background = element_blank(),  # Customize background
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text = element_blank(),  # Hide axis text
        axis.ticks = element_blank())  # Hide axis ticks

p

# Save the plot
ggsave(filename = "Cluster/phylo_jaccard_ward_PAV_genes.png", plot = p, width = 16, height = 16, dpi = 300)
```

## Pseudogene based

Join not transposed matrices

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/all

mkdir pseudogene_cluster
cd pseudogene_cluster

ln -s ../../chr1/matrix_chr1.txt
ln -s ../../chr2/matrix_chr2.txt
ln -s ../../chr3/matrix_chr3.txt
ln -s ../../chr4/matrix_chr4.txt
ln -s ../../chr5/matrix_chr5.txt

./combine_notransp_matrices.py

mv concatenated_matrix.csv pseudogenes_concatenated_matrix.csv

cp pseudogenes_concatenated_matrix.csv ~/Arabidopsis_pangenome/R_plots/Cluster/
```

### make the figure

```
mamba activate R
R

library("ggplot2")
library("vegan")
library(ggtree)
library(ape) 

setwd("/home/lia/Arabidopsis_pangenome/R_plots")

# Load data
my_data <- read.csv("Cluster/pseudogenes_concatenated_matrix.csv", header = TRUE, sep = "\t")
head(my_data)

# Extract assembly names
assembly_names <- colnames(my_data)[2:94]

# Load locations
locations <- read.csv("genebank_country_1.csv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(locations)

colnames(locations) <- c("Assembly", "Location") 

country_to_colour <- c(
'Afghanistan' = "salmon2",
'Belgium' = "gold3",
'Cape Verde' = "darkgreen",
'Czech Republic' = "grey46",
'France' = "red", 
'Germany' = "magenta",
'Italy' = "brown",
'Japan' = "springgreen4",
'Lithuania' = "black",
'Netherlands' = "green",
'Poland' = "orange",
'Madeira' = "darkviolet",
'Romania' = "pink",
'Spain' = 'deepskyblue3',
'Sweden' = "olivedrab",
'Tanzania' = 'blue',
'United Kingdom' = 'mediumpurple',
'USA' = 'cyan')


# Create a mapping from assembly names to countries
assembly_to_country <- setNames(locations$Location, locations$Assembly)

# Map assembly names to colors through their countries
assembly_to_color <- sapply(assembly_names, function(assembly) {
  country <- assembly_to_country[assembly]
  color <- country_to_colour[country]
  return(color)
})


# Calculate Jaccard distance
dist.jaccard <- vegdist(t(my_data[, 2:94]), method = "jaccard")
write.csv(as.matrix(dist.jaccard), "pseudogenes_PAV_jaccard_distance_matrix.csv", row.names = TRUE)


# Hierarchical clustering
res.hc <- hclust(d = dist.jaccard, method = "ward.D2")


# Calculate the cophenetic distance matrix from the clustering result
res.coph<-cophenetic(res.hc)

# Compute the cophenetic correlation
cophenetic_corr <- cor(as.dist(dist.jaccard), as.dist(res.coph))

# Print cophenetic correlation to check the quality of the clustering
print(cophenetic_corr)

# Convert to phylo object
phylo_tree <- as.phylo(res.hc)

# Prepare a data frame that ggtree can use to map colors to tips
# The row names of 'my_data' are assumed to be the assembly names used in 'phylo_tree$tip.label'
# Ensure these names match those in 'assembly_to_color'
tip_colors <- data.frame(label = names(assembly_to_color), color = assembly_to_color)


# Plot the tree
p <- ggtree(phylo_tree, layout="circular") +
  geom_tiplab(aes(color = label), size = 4) +
  scale_color_manual(values = tip_colors$color) +
  theme_tree2() +
  theme(legend.position = "none",  # Adjust legend position
        plot.background = element_blank(),  # Customize background
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text = element_blank(),  # Hide axis text
        axis.ticks = element_blank())  # Hide axis ticks

p

# Save the plot
ggsave(filename = "Cluster/phylo_jaccard_ward_PAV_pseudogenes.png", plot = p, width = 16, height = 16, dpi = 300)
```

## Node based 
For this figure we generated the input file using `odgi similarity` after joining all the graphs from each chromosomes with `odgi squeeze`.
To generate the phylogenetic tree we will follow these instructions: https://hackmd.io/@AndreaGuarracino/SyhbiKuE2#Primate-chromosome-6

```
mamba activate R

R

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ape)
library(ggtree)
library(ggplot2)

setwd ("/home/lia/Arabidopsis_pangenome/FINAL/nodes_similarity")

path_dist_tsv <- 'similarity_nodes.tsv'

# Read sparse matrix
sparse_matrix_df <- read_tsv(path_dist_tsv)

# Prepare distance matrix
jaccard_dist_df <- sparse_matrix_df %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")

dist.jaccard <- as.dist(jaccard_dist_df)
write.csv(as.matrix(dist.jaccard), "nodes_jaccard_distance_matrix.csv", row.names = TRUE)

# Perform hierarchical clustering
res.hc <- hclust(d = dist.jaccard)

# Calculate the cophenetic distance matrix from the clustering result
res.coph<-cophenetic(res.hc)

# Compute the cophenetic correlation
cophenetic_corr <- cor(as.dist(dist.jaccard), as.dist(res.coph))

# Print cophenetic correlation to check the quality of the clustering
print(cophenetic_corr)

# Convert hclust object to phylo
phylo_tree <- as.phylo(res.hc)

# Load locations
locations <- read.csv("genebank_country_1.csv", sep="\t", stringsAsFactors = FALSE, header = TRUE)
head(locations)

colnames(locations) <- c("Assembly", "Location") 

#extract assembly names
assembly_names <- locations$Assembly

country_to_colour <- c(
'Afghanistan' = "salmon2",
'Belgium' = "gold3",
'Cape Verde' = "darkgreen",
'Czech Republic' = "grey46",
'France' = "red", 
'Germany' = "magenta",
'Italy' = "brown",
'Japan' = "springgreen4",
'Lithuania' = "black",
'Netherlands' = "green",
'Poland' = "orange",
'Madeira' = "darkviolet",
'Romania' = "pink",
'Spain' = 'deepskyblue3',
'Sweden' = "olivedrab",
'Tanzania' = 'blue',
'United Kingdom' = 'mediumpurple',
'USA' = 'cyan')


# Create a mapping from assembly names to countries
assembly_to_country <- setNames(locations$Location, locations$Assembly)


# Map assembly names to colors through their countries
assembly_to_color <- sapply(assembly_names, function(assembly) {
  country <- assembly_to_country[assembly]
  color <- country_to_colour[country]
  return(color)
})

# Prepare a data frame that ggtree can use to map colors to tips
# The row names of 'my_data' are assumed to be the assembly names used in 'phylo_tree$tip.label'
# Ensure these names match those in 'assembly_to_color'
tip_colors <- data.frame(label = names(assembly_to_color), color = assembly_to_color)


# Remove country information from tip_colors labels
tip_colors$label <- sub("\\..*", "", tip_colors$label)

# Ensure tip_colors is a named vector
tip_colors_named <- setNames(tip_colors$color, tip_colors$label)

# Plot the tree with color mapping
p <- ggtree(phylo_tree, layout="circular") +
  geom_tiplab(aes(label = label, color = label), size = 4) +
  scale_color_manual(values = tip_colors_named) +
  theme_tree2() +
  theme(legend.position = "none",  # Adjust legend position
        plot.background = element_blank(),  # Customize background
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  # Remove grid lines
        axis.text = element_blank(),  # Hide axis text
        axis.ticks = element_blank())  # Hide axis ticks

# Display the plot
print(p)

# Save the plot
ggsave(filename = "phylo_jaccard_nodes.png", plot = p, width = 16, height = 16, dpi = 300)
```

# Other data visualisation

## Assembly based gene classification - PAV

### Obtain the files

We will need the following files:

- `presence_absence_matrix.csv` with the results from all the chromosomes
- A file called `loci_classes.txt`, where we have the `Feature/Pater` and `Label` columns from `final_screening_GO.tsv` with the results from all the chromosomes
- `matrix_transposed.csv`, to be obtained with `transpose.pl`
- The dataset rearranged for R, to be obtained with `comp_node_path.pl`

1. Final screening and Presence-absence matrix per chromosome

```
# chr1

cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/chr1
cut -f 1,4 final_genes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' loci_classes.txt > chr1_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr1.txt

# transpose matrix
perl ../transpose.pl matrix_chr1.txt > chr1_transp_matrix.txt


# chr2

cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/chr2
cut -f 1,4 final_genes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' loci_classes.txt > chr2_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr2.txt

# transpose matrix
perl ../transpose.pl matrix_chr2.txt > chr2_transp_matrix.txt



# chr3

cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/chr3
cut -f 1,4 final_genes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' loci_classes.txt > chr3_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr3.txt

# transpose matrix
perl ../transpose.pl matrix_chr3.txt > chr3_transp_matrix.txt


# chr4

cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/chr4
cut -f 1,4 final_genes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' loci_classes.txt > chr4_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr4.txt

# transpose matrix
perl ../transpose.pl matrix_chr4.txt > chr4_transp_matrix.txt


# chr5

cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/chr5
cut -f 1,4 final_genes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' loci_classes.txt > chr5_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr5.txt

# transpose matrix
perl ../transpose.pl matrix_chr5.txt > chr5_transp_matrix.txt
```

2. Join Transposed matrices

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/all

ln -s ../chr1/chr1_transp_matrix.txt
ln -s ../chr2/chr2_transp_matrix.txt
ln -s ../chr3/chr3_transp_matrix.txt
ln -s ../chr4/chr4_transp_matrix.txt
ln -s ../chr5/chr5_transp_matrix.txt

../combine_transposed_matrices.py
```

3. Join loci classes

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/genes/all

ln -s ../chr1/chr1_loci_classes.txt
ln -s ../chr2/chr2_loci_classes.txt
ln -s ../chr3/chr3_loci_classes.txt
ln -s ../chr4/chr4_loci_classes.txt
ln -s ../chr5/chr5_loci_classes.txt

for f in chr*_loci_classes.txt; do tail -n +2 "$f"; done >> combined_loci_classes.txt
```

3. R dataset

```
perl ../comp_node_path_PAV.pl combined_loci_classes.txt concatenated_matrix.csv > loci_comp_classes.txt
```

This last script produces also `data_set_R.txt` which is the file we are going to use for the figure. Let's rename the file:

```
mv data_set_R.txt R_data_genes_PAV_total.txt
```

### Make the figure

Install required packages

```
mamba activate R
mamba install r-ggplot2
mamba install r-factoextra
mamba install conda-forge::r-vegan
mamba install bioconda::bioconductor-ggtree
mamba install conda-forge::r-ape
mamba install bioconda::bioconductor-ggtreeextra
```

Open R

```
R
```

```
setwd("/home/lia/Arabidopsis_pangenome/R_plots")
```

```
data<-read.table(file="./Assembly_based_genes/R_data_genes_PAV_total.txt", sep="\t", header=FALSE)
head(data)
```

```
colnames(data)<-c("Assembly", "Class", "Count")
head(data)
```

```
library("ggplot2")
```

```
custom_colors <- c("private" = "black", "dispensable" = "beige", "softcore" = "olivedrab", "core" = "firebrick")  # Replace with your actual classes and desired colors

# Order the classes in a custom order
data$Class <- factor(data$Class, levels = c("private", "dispensable", "softcore", "core"))  # Replace with your actual class order

# Assign names to the colors based on the levels of Class
names(custom_colors) <- levels(data$Class)

# Create the plot
gene_plot_impr <- ggplot(data = data, aes(x = Assembly, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +  # Set a white background
  theme(text = element_text(size = 6),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Ensure the background is white
        plot.background = element_blank()) +
  scale_fill_manual(values = custom_colors)  # Apply the custom colors

# Print the plot
print(gene_plot_impr)

```

```
ggsave(filename = "Assembly_based_genes/gene_plot_PAV_total.png", plot = gene_plot_impr, width = 10, height = 7, dpi = 300)
```

## Assembly based gene classification - CNV

This will be based on untangle results only.

### Obtain the files

We will need the following files:

- `matrix.csv` with the results from all the chromosomes
- A file called `loci_classes.txt`, where we have the `Feature/Pater` and `Label` columns from `pangenome_screening.csv` with the results from all the chromosomes
- `matrix_transposed.csv`, to be obtained with `transpose.pl`
- The dataset rearranged for R, to be obtained with `comp_node_path.pl`

1. Pangenome screening and matrix per chromosome

```
# chr1

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/chr1
cut -f 2,8 pangenome_screening.csv > chr1_loci_classes.txt

# transpose matrix
perl ../transpose.pl matrix.csv > chr1_transp_matrix.txt


# chr2

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/chr2
cut -f 2,8 pangenome_screening.csv > chr2_loci_classes.txt

# transpose matrix
perl ../transpose.pl matrix.csv > chr2_transp_matrix.txt


# chr3

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/chr3
cut -f 2,8 pangenome_screening.csv > chr3_loci_classes.txt

# transpose matrix
perl ../transpose.pl matrix.csv > chr3_transp_matrix.txt


# chr4

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/chr4
cut -f 2,8 pangenome_screening.csv > chr4_loci_classes.txt

# transpose matrix
perl ../transpose.pl matrix.csv > chr4_transp_matrix.txt


# chr5

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/chr5
cut -f 2,8 pangenome_screening.csv > chr5_loci_classes.txt

# transpose matrix
perl ../transpose.pl matrix.csv > chr5_transp_matrix.txt
```

2. Join Transposed matrices

```
cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/
mkdir all
cd all

ln -s ../chr1/chr1_transp_matrix.txt
ln -s ../chr2/chr2_transp_matrix.txt
ln -s ../chr3/chr3_transp_matrix.txt
ln -s ../chr4/chr4_transp_matrix.txt
ln -s ../chr5/chr5_transp_matrix.txt

../combine_transposed_matrices.py
```

3. Join loci classes

```
cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/genes/all

ln -s ../chr1/chr1_loci_classes.txt
ln -s ../chr2/chr2_loci_classes.txt
ln -s ../chr3/chr3_loci_classes.txt
ln -s ../chr4/chr4_loci_classes.txt
ln -s ../chr5/chr5_loci_classes.txt

for f in chr*_loci_classes.txt; do tail -n +2 "$f"; done >> combined_loci_classes.txt
```

3. R dataset

```
perl ../comp_node_path_CNV.pl combined_loci_classes.txt concatenated_matrix.csv > loci_comp_classes.txt
```

This last script produces also `data_set_R.txt` which is the file we are going to use for the figure. Let's rename the file:

```
mv data_set_R.txt R_data_genes_CNV_untangle.txt
```

### Make the figure

```{r}
setwd("/home/lia/Arabidopsis_pangenome/R_plots")
```

```{r}
data<-read.table(file="Assembly_based_genes/R_data_genes_CNV_untangle.txt", sep="\t", header=FALSE)
head(data)
```

```{r}
colnames(data)<-c("Assembly", "Class", "Count")
head(data)
```

```{r}
library("ggplot2")
```

```{r}
custom_colors <- c("private" = "black", "dispensable" = "beige", "softcore" = "olivedrab", "core" = "firebrick")  # Replace with your actual classes and desired colors

# Order the classes in a custom order
data$Class <- factor(data$Class, levels = c("private", "dispensable", "softcore", "core"))  # Replace with your actual class order

# Assign names to the colors based on the levels of Class
names(custom_colors) <- levels(data$Class)

# Create the plot
gene_plot_impr <- ggplot(data = data, aes(x = Assembly, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +  # Set a white background
  theme(text = element_text(size = 6),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Ensure the background is white
        plot.background = element_blank()) +
  scale_fill_manual(values = custom_colors)  # Apply the custom colors

# Print the plot
print(gene_plot_impr)

```

```{r}
ggsave(filename = "Assembly_based_genes/gene_plot_CNV_untangle.png", plot = gene_plot_impr, width = 10, height = 7, dpi = 300)
```

## Assembly based pseudogene classification - PAV

### Obtain the files

We will need the following files:

- `presence_absence_matrix.csv` with the results from all the chromosomes
- A file called `loci_classes.txt`, where we have the `Feature/Pater` and `Label` columns from `final_screening_GO.tsv` with the results from all the chromosomes
- `matrix_transposed.csv`, to be obtained with `transpose.pl`
- The dataset rearranged for R, to be obtained with `comp_node_path.pl`

1. Final screening and Presence-absence matrix per chromosome

```
# chr1

cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/chr1
cut -f 1,4 final_pseudogenes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' loci_classes.txt > chr1_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr1.txt

# transpose matrix
perl ../transpose.pl matrix_chr1.txt > chr1_transp_matrix.txt


# chr2

cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/chr2
cut -f 1,4 final_pseudogenes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' loci_classes.txt > chr2_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr2.txt

# transpose matrix
perl ../transpose.pl matrix_chr2.txt > chr2_transp_matrix.txt



# chr3

cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/chr3
cut -f 1,4 final_pseudogenes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' loci_classes.txt > chr3_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr3.txt

# transpose matrix
perl ../transpose.pl matrix_chr3.txt > chr3_transp_matrix.txt


# chr4

cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/chr4
cut -f 1,4 final_pseudogenes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' loci_classes.txt > chr4_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr4.txt

# transpose matrix
perl ../transpose.pl matrix_chr4.txt > chr4_transp_matrix.txt


# chr5

cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/chr5
cut -f 1,4 final_pseudogenes_screening_GO.tsv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' loci_classes.txt > chr5_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' presence_absence_matrix.csv > matrix_chr5.txt

# transpose matrix
perl ../transpose.pl matrix_chr5.txt > chr5_transp_matrix.txt
```

2. Join Transposed matrices

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/

mkdir all
cd all

ln -s ../chr1/chr1_transp_matrix.txt
ln -s ../chr2/chr2_transp_matrix.txt
ln -s ../chr3/chr3_transp_matrix.txt
ln -s ../chr4/chr4_transp_matrix.txt
ln -s ../chr5/chr5_transp_matrix.txt

../combine_transposed_matrices.py
```

3. Join loci classes

```
cd ~/Arabidopsis_pangenome/FINAL/final_screening/pseudogenes/all

ln -s ../chr1/chr1_loci_classes.txt
ln -s ../chr2/chr2_loci_classes.txt
ln -s ../chr3/chr3_loci_classes.txt
ln -s ../chr4/chr4_loci_classes.txt
ln -s ../chr5/chr5_loci_classes.txt

for f in chr*_loci_classes.txt; do tail -n +2 "$f"; done >> combined_loci_classes.txt
```

3. R dataset

```
perl ../comp_node_path_PAV.pl combined_loci_classes.txt concatenated_matrix.csv > loci_comp_classes.txt
```

This last script produces also `data_set_R.txt` which is the file we are going to use for the figure. Let's rename the file:

```
mv data_set_R.txt R_data_pseudogenes_PAV_total.txt
```

### Make the figure

Install required packages

```
mamba activate R
mamba install r-ggplot2
mamba install r-factoextra
```

Open R

```
R
```

```
setwd("/home/lia/Arabidopsis_pangenome/R_plots")
```

```
data<-read.table(file="./Assembly_based_pseudogenes/R_data_pseudogenes_PAV_total.txt", sep="\t", header=FALSE)
head(data)
```

```
colnames(data)<-c("Assembly", "Class", "Count")
head(data)
```

```
library("ggplot2")
```

```
custom_colors <- c("private" = "black", "dispensable" = "beige", "softcore" = "olivedrab", "core" = "firebrick")  # Replace with your actual classes and desired colors

# Order the classes in a custom order
data$Class <- factor(data$Class, levels = c("private", "dispensable", "softcore", "core"))  # Replace with your actual class order

# Assign names to the colors based on the levels of Class
names(custom_colors) <- levels(data$Class)

# Create the plot
gene_plot_impr <- ggplot(data = data, aes(x = Assembly, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +  # Set a white background
  theme(text = element_text(size = 6),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Ensure the background is white
        plot.background = element_blank()) +
  scale_fill_manual(values = custom_colors)  # Apply the custom colors

# Print the plot
print(gene_plot_impr)

```

```
ggsave(filename = "Assembly_based_pseudogenes/pseudogene_plot_PAV_total.png", plot = gene_plot_impr, width = 10, height = 7, dpi = 300)
```

## Assembly based pseudogene classification - CNV

This will be based on untangle results only.

### Obtain the files

We will need the following files:

- `matrix.csv` with the results from all the chromosomes
- A file called `loci_classes.txt`, where we have the `Feature/Pater` and `Label` columns from `pangenome_screening.csv` with the results from all the chromosomes
- `matrix_transposed.csv`, to be obtained with `transpose.pl`
- The dataset rearranged for R, to be obtained with `comp_node_path.pl`

1. Pangenome screening and matrix per chromosome

```
# chr1

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/chr1
cut -f 2,8 pangenome_screening.csv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' loci_classes.txt > chr1_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr1"; print $0 }' OFS='\t' matrix.csv > matrix_chr1.txt

# transpose matrix
perl ../transpose.pl matrix_chr1.txt > chr1_transp_matrix.txt


# chr2

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/chr2
cut -f 2,8 pangenome_screening.csv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' loci_classes.txt > chr2_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr2"; print $0 }' OFS='\t' matrix.csv > matrix_chr2.txt

# transpose matrix
perl ../transpose.pl matrix_chr2.txt > chr2_transp_matrix.txt


# chr3

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/chr3
cut -f 2,8 pangenome_screening.csv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' loci_classes.txt > chr3_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr3"; print $0 }' OFS='\t' matrix.csv > matrix_chr3.txt

# transpose matrix
perl ../transpose.pl matrix_chr3.txt > chr3_transp_matrix.txt


# chr4

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/chr4
cut -f 2,8 pangenome_screening.csv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' loci_classes.txt > chr4_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr4"; print $0 }' OFS='\t' matrix.csv > matrix_chr4.txt

# transpose matrix
perl ../transpose.pl matrix_chr4.txt > chr4_transp_matrix.txt


# chr5

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/chr5
cut -f 2,8 pangenome_screening.csv > loci_classes.txt

# append _chr* to gene names
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' loci_classes.txt > chr5_loci_classes.txt
awk -F'\t' 'NR==1 {print; next} { $1 = $1 "_chr5"; print $0 }' OFS='\t' matrix.csv > matrix_chr5.txt

# transpose matrix
perl ../transpose.pl matrix_chr5.txt > chr5_transp_matrix.txt

2.  Join Transposed matrices

```

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/
mkdir all
cd all

ln -s ../chr1/chr1_transp_matrix.txt
ln -s ../chr2/chr2_transp_matrix.txt
ln -s ../chr3/chr3_transp_matrix.txt
ln -s ../chr4/chr4_transp_matrix.txt
ln -s ../chr5/chr5_transp_matrix.txt

../combine_transposed_matrices.py

```

3.  Join loci classes

```

cd ~/Arabidopsis_pangenome/FINAL/untangle_based_results/pseudogenes/all

ln -s ../chr1/chr1_loci_classes.txt
ln -s ../chr2/chr2_loci_classes.txt
ln -s ../chr3/chr3_loci_classes.txt
ln -s ../chr4/chr4_loci_classes.txt
ln -s ../chr5/chr5_loci_classes.txt

for f in chr*_loci_classes.txt; do tail -n +2 "$f"; done >> combined_loci_classes.txt

```

3.  R dataset

```

perl ../comp_node_path_CNV.pl combined_loci_classes.txt concatenated_matrix.csv > loci_comp_classes.txt

```

This last script produces also `data_set_R.txt` which is the file we are going to use for the figure. Let's rename the file:

```

mv data_set_R.txt R_data_pseudogenes_CNV_untangle.txt

```

```

mamba activate R
R

```

### Make the figure

```{r}
setwd("/home/lia/Arabidopsis_pangenome/R_plots")
```

```{r}
data<-read.table(file="Assembly_based_pseudogenes/R_data_pseudogenes_CNV_untangle.txt", sep="\t", header=FALSE)
head(data)
```

```{r}
colnames(data)<-c("Assembly", "Class", "Count")
head(data)
```

```{r}
library("ggplot2")
```

```{r}
custom_colors <- c("private" = "black", "dispensable" = "beige", "softcore" = "olivedrab", "core" = "firebrick")  # Replace with your actual classes and desired colors

# Order the classes in a custom order
data$Class <- factor(data$Class, levels = c("private", "dispensable", "softcore", "core"))  # Replace with your actual class order

# Assign names to the colors based on the levels of Class
names(custom_colors) <- levels(data$Class)

# Create the plot
gene_plot_impr <- ggplot(data = data, aes(x = Assembly, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +  # Set a white background
  theme(text = element_text(size = 6),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Ensure the background is white
        plot.background = element_blank()) +
  scale_fill_manual(values = custom_colors)  # Apply the custom colors

# Print the plot
print(gene_plot_impr)

```

```{r}
ggsave(filename = "Assembly_based_pseudogenes/pseudogene_plot_CNV_untangle.png", plot = gene_plot_impr, width = 10, height = 7, dpi = 300)
```

## Assembly based node classification

we need to obtain a huge matrix:

```
cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/all_graphs

ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr1/node_matrix/chr1_matrix.csv
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr2/node_matrix/chr2_matrix.csv
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr3/node_matrix/chr3_matrix.csv
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr4/node_matrix/chr4_matrix.csv
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr5/node_matrix/chr5_matrix.csv

screen -S combine_matrices

mamba activate pandas

combine_transposed_matrices.py


cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr1/node_matrix/
cut -f 1,4 PAV_CNV.csv > chr1_loci_classes.txt

cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr2/node_matrix/
cut -f 1,4 PAV_CNV.csv > chr2_loci_classes.txt

cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr3/node_matrix/
cut -f 1,4 PAV_CNV.csv > chr3_loci_classes.txt

cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr4/node_matrix/
cut -f 1,4 PAV_CNV.csv > chr4_loci_classes.txt

cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr5/node_matrix/
cut -f 1,4 PAV_CNV.csv > chr5_loci_classes.txt

cd /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/all_graphs

ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr1/node_matrix/chr1_loci_classes.txt
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr2/node_matrix/chr2_loci_classes.txt
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr3/node_matrix/chr3_loci_classes.txt
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr4/node_matrix/chr4_loci_classes.txt
ln -s /home/edg01/edg01/lia/Arabidopsis/pangenome/ara_pan_from_April/chr5/node_matrix/chr5_loci_classes.txt


for f in chr*_loci_classes.txt; do tail -n +2 "$f"; done >> combined_loci_classes.txt

screen -r combine_matrices

./comp_node_path_PAV.pl combined_loci_classes.txt concatenated_matrix.csv > loci_comp_classes_pav.txt

mv data_set_R.txt pav_data_set_R.txt
```


### Make the figures

```

mamba activate R
R

```

```{r}
setwd("/home/lia/Arabidopsis_pangenome/R_plots")
```

```{r}
data<-read.table(file="Assembly_based_nodes/pav_data_set_R.txt", sep="\t", header=FALSE)
head(data)
```

```{r}
colnames(data)<-c("Assembly", "Class", "Count")
head(data)
```

```{r}
library("ggplot2")
library(scale)
```

```{r}
custom_colors <- c("private" = "black", "dispensable" = "beige", "softcore" = "olivedrab", "core" = "firebrick")  # Replace with your actual classes and desired colors

# Order the classes in a custom order
data$Class <- factor(data$Class, levels = c("private", "dispensable", "softcore", "core"))  # Replace with your actual class order

# Assign names to the colors based on the levels of Class
names(custom_colors) <- levels(data$Class)

# Create the plot
node_plot_pav <- ggplot(data = data, aes(x = Assembly, y = Count, fill = Class)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +  # Set a white background
  theme(text = element_text(size = 6),
        #panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_blank(),  # Ensure the background is white
        plot.background = element_blank()) +
  scale_fill_manual(values = custom_colors) +  # Apply the custom colors
  scale_y_continuous(labels = label_number())  # Format y-axis labels to avoid scientific notation

# Print the plot
print(node_plot_pav)

```

```{r}
ggsave(filename = "Assembly_based_nodes/nodes_plot_PAV.png", plot = node_plot_pav, width = 10, height = 7, dpi = 300)
```


## Pie charts

Refer to `scripts/python_plot.ipynb`
