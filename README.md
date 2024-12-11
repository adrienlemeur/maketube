
<h1 align="center"> Maketube </h1>
<p align="center">
   <img src="/maketube_logo.png" width="200" height="100">
</p>

maketube is an artificial genome generator for Mycobacterium tuberculosis

### Table of contents

<!--ts-->
   - [Installation](#install)
   - [Quickstart](#quickstart)
   - [Pipeline](#pipeline)
   - [Usage](#usage)

### <a name="install"></a>Installation
#### Maketube :
  ##### Dependencies:

  - [R](https://www.r-project.org/) (v4.1.2)
  - [ape](https://cran.r-project.org/web/packages/ape/index.html) (v5.8)
  - [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) (v4.2-36)
  - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (v2.62.0)
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (v1.1.4)
  - [jackalope](https://cran.r-project.org/web/packages/jackalope/index.html) (v1.1.5)

   ##### Installation procedure

   Install the dependencies & >>
   ```
   git clone https://github.com/adrienlemeur/maketube.git
   cd maketube
   echo "PATH=\"$(pwd)\"/:$PATH" >> ~/.bashrc && source ~/.bashrc
   ```
#### vcf2metrics.py
  ##### Dependencies
  - [cyvcf2](https://brentp.github.io/cyvcf2/) (0.30.18)
  - [numpy](https://numpy.org/) (1.26.4) (⚠ 2.0 produces error)
  - [argparse](https://pypi.org/project/argparse/) (3.2)
  - math, re, sys, gc & os

#### Container
Don't want to install all these pesky packages and their dependencies ? There is a [container](https://hub.docker.com/r/alemeur/maketube) !
#### Singularity
```
singularity pull maketube.img docker://alemeur/maketube:latest
#run maketube.R
singularity run maketube.img
#run vcf2metrics.R
singularity exec maketube.img "/usr/local/bin/maketube/vcf2metrics.py"
```
##### Docker (sill in progress)
```
docker pull alemeur/maketube:latest
#run maketube.R
docker run maketube:latest
#run vcf2metrics.R
# in /usr/local/bin/maketube/vcf2metrics.py
```

### <a name="quickstart"></a>Quickstart
#### Building a set of genomes from H37Rv (minimal input)
```
cd maketube
gunzip REF/nonH37Rv_pool_sequence.fasta.gz

Rscript maketube.R \
   --reference REF/H37Rv.fasta \
   --transposon REF/H37Rv_transposon.bed \
   --nonhomoseq_pool REF/nonH37Rv_pool_sequence.fasta
```
#### Comparing a test VCF to the reference VCF
```
#Comparing variants found in Haplotype 1 (H1) of population 1

bcftools view --samples H1 maketube_run/SV1/SV1_pop1.vcf.gz > my_reference_vcf.vcf.gz

vcf2metrics.py -i my_sample_vcf.vcf.gz --reference my_reference_vcf.vcf.gz \
		--backtrack maketube_run/SV1/SV1_equivalence.bed
      --bed maketube_run/SV1/SV1_SV.bed
```

#### Comparing variants to the reference

vcf2metrics.py compares a sample vcf to the reference VCF built by maketube using the backtrack partition. It can also compare a test VCF to a reference VCF.

Usage :
```
vcf2metrics [-h] -i I [I ...] --reference REFERENCE [REFERENCE ...]
                   [--bed BED [BED ...]]
                   [--backtrack BACKTRACK [BACKTRACK ...]] [--trace]
                   [--add_col ADD_COL [ADD_COL ...]] [--sample SAMPLE]
                   [-o O [O ...]]
options:
  -h, --help            show this help message and exit
  -i I [I ...]          folder containing VCF to test against the reference
  --reference REFERENCE [REFERENCE ...]
                        reference vcf
  --bed BED [BED ...]   bed file with the region annotations
  --backtrack BACKTRACK [BACKTRACK ...]
                        bed backtrack file produced by maketube
  --trace
  --add_col ADD_COL [ADD_COL ...]
                        Space separated list of columns specifying different
                        parameters of the experiment if needed
  -o O [O ...]          output directory

Written by Adrien Le Meur, v??
```

#### Maketube output
```
maketube_run/
├── my_run_arborescence.tsv		#tsv with strains information (position in the arborescence, corresponding by-products)
├── SV1					# Independant set of structural variant n°1
│   ├── my_run_unmuted.fasta		# sequence with structural variants but without snps and small indels
│   ├── pop1				# population index of SV1
│   │   ├── FASTA				# Fasta file for every strain, with structural variant Set n°1
│   │   ├── FASTQ				# Corresponding FASTQ ⚠️ BROKEN ⚠️; Delete this & run independantly art_illumina from the FASTA
│   │   ├── SV1_pop1.nwk		# Population tree
│   │   ├── SV1_pop1.vcf.gz		# Multisample VCF containing snps and indels of every haplotype (H1, H2, ..., Hn)
│   ├── SV1_equivalence.bed		# Backtrack partition used by vcf2metrics.py
│   └── SV1_SV.bed			# Position of structural variants
│
[...]
│
└── SV10					# Independant set of structural variant n°10
    ├── my_run_unmuted.fasta
    ├── pop1
    │   ├── FASTA
    │   ├── FASTQ
    │   ├── SV10_pop1.nwk
    │   ├── SV10_pop1.vcf.gz
    ├── SV10_equivalence.bed
    └── SV10_SV.bed
```

### <a name="usage"></a>Full description
```
./maketube.R \
   --reference (fasta) :      fasta files of the reference sequence that will be used as a base
   --transposon (bed) :      bed file delimitating the transposon sequences that will be cut and paste across the genome
   --nonhomoseq_pool (fasta) :      a fasta file with no sequence name delimiter (">"). Maketube will randomely sample across these sequences and concatenate them to create than. As such, it is advised to give sequences that are not present in the reference. \
   --haplotype_count :      number of strain to compute for each population. Default : 10 \
   --pop_count :      number of population to compute. Default : 2 \
   --structural_variants :      number of different set of structural variant by population. Default : 2 \
   #Total number of strain : haplotype_count x pop_count x structural_variants \
   --unmuted (flag) :      write the fasta sequence without the mutation. Default : false (no flag) \
   --duplication_region_size (int) :      the size of the region to be duplicated \
   --prefix (string) :      prefix of the generated sample \
   -o/--output (string) :      name directory to put the results in. Default : maketube_results \
   --mutation_rate (float) :      in substitution by generation by site. Default : 1.23*10^-7 \
   --effective_pop (int) :      the dreadful Ne. Default : 700 \
   --TCAGI (R vector) :      parameters of the GTR model, in the same order as specified in the flag. Default : "0.172,0.329,0.172,0.329,0" \
   --ABCDEF (R vector) :      parameters of the GTR model, in the same order as specified in the flag. Default : "0.65,0.05,0.21,0.27,0.02,0.65" \
   --scaling_factor (int) :      empirical factor to get a number of indel equivalent to 1/8 the number of SNP, as found in natural strains. Only affects the number of indels. Default : 0.125 \
   --slope (int, > 0) :      Size of the distance between 2 structural variants. Default : 300. \
   --threads (int) :      Maximum number of threads to use for generating fastq. Default : 4.
```
