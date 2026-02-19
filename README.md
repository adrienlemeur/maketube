link github
<h1 align="center"> Maketube </h1>
<p align="center">
   <img src="/maketube_logo.png" width="200" height="100">
</p>

maketube evolves a _Mycobacterium tuberculosis_ genome

## Table of contents

<!--ts-->
   - [Installation](#install)
     - [Maketube dependencies](#install_maketube)
     - [vcf2metrics dependencies](#install_vcf2metrics)
     - [Installation](#installation_procedure)
     - [Container](#container)
   - [Quickstart](#quickstart)
   - [Pipeline](#pipeline)
   - [Usage](#usage)


## <a name="install"></a>Installation

### <a name="install_maketube"></a>Maketube:
  #### Dependencies:

  - [R](https://www.r-project.org/) (v4.1.2)
  - [ape](https://cran.r-project.org/web/packages/ape/index.html) (v5.8)
  - [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) (v4.2-36)
  - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (v2.62.0)
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (v1.1.4)
  - [jackalope](https://cran.r-project.org/web/packages/jackalope/index.html) (v1.1.5)

### <a name="install_vcf2metrics"></a>Mvcf2metrics.py
  #### Dependencies
  - [cyvcf2](https://brentp.github.io/cyvcf2/) (0.30.18)
  - [numpy](https://numpy.org/) (1.26.4) (⚠ 2.0 produces error)
  - [argparse](https://pypi.org/project/argparse/) (3.2)
  - math, re, sys, gc & os

#### <a name="installation_procedure"></a>Installation procedure

Install the dependencies & >>
  ```
  git clone https://github.com/adrienlemeur/maketube.git
  cd maketube
  echo "PATH=\"$(pwd)\"/:$PATH" >> ~/.bashrc && source ~/.bashrc
  ```

### <a name="container">Container
Don't want to install all these pesky packages and their dependencies ? There is a [container](https://hub.docker.com/r/alemeur/maketube) !
### Singularity
```
singularity pull maketube.img docker://alemeur/maketube:latest
#run maketube.R
singularity run maketube.img
#run vcf2metrics.R
singularity exec maketube.img "/usr/local/bin/maketube/vcf2metrics.py"
```
#### Docker
```
docker pull alemeur/maketube:latest
#run maketube.R
docker run maketube:latest
#run vcf2metrics.R
# in /usr/local/bin/maketube/vcf2metrics.py
```

## <a name="quickstart"></a>Quickstart
### Building a set of genomes from H37Rv (minimal input)
```
cd maketube
gunzip REF/nonH37Rv_pool_sequence.fasta.gz

Rscript maketube.R \
   --reference REF/H37Rv.fasta \
   --transposon REF/H37Rv_transposon.bed \
   --nonhomoseq_pool REF/nonH37Rv_pool_sequence.fasta
```

### Comparing a test VCF to the reference VCF
```
#Comparing variants found in Haplotype 1 (H1) of population 1

bcftools view --samples H1 maketube_run/SV1/SV1_pop1.vcf.gz > my_reference_vcf.vcf.gz

vcf2metrics.py -i my_sample_vcf.vcf.gz \
		--reference my_reference_vcf.vcf.gz \
		--backtrack maketube_run/SV1/SV1_equivalence.bed \
		--bed maketube_run/SV1/SV1_SV.bed
```

## <a name="usage"></a>Complete usage

## Full description
```
Usage: ./maketube.R [options]


Options:
	--reference=CHARACTER
		reference sequence to evolve (fasta)

	--transposon=CHARACTER
		Inserting Sequences (transposon-like sequences) initial positions in the reference sequence (bed). You can get these by blasting the genome sequence for the IS/transposon-like elements.

	--nonhomoseq_pool=CHARACTER
		a list of kmers (fasta). All kmers should be the same size. Maketube subsamples this list at random and concatenates them to build ancetral-like regions (regions that are present in the evolved genome but not in the reference sequence).

	--duplication_region_size=CHARACTER
		size of the duplicated region.

	--haplotype_count=HAPLOTYPE_COUNT
		number of haplotype to compute (>=1)

	--pop_count=POP_COUNT
		number of haplotypes to compute (>=1)

	--structural_variants=STRUCTURAL_VARIANTS
		number of structural variants

	--deletion_count=DELETION_COUNT
		number of deletion region to remove

	--unmuted
		Should maketube create an unmuted sequence of the evolved genome

	-p PREFIX, --prefix=PREFIX
		artificial genome name prefix

	--mutation_rate=MUTATION_RATE
		the genome wide mutation rate

	--effective_pop=EFFECTIVE_POP
		the effective population size

	--TCAGI=TCAGI
		the rate of the different nucleotide (TCAG) + invariants (I)

	--ABCDEF=ABCDEF
		the parameters of the GTR model

	--indel_scaling_factor=INDEL_SCALING_FACTOR
		indel/snp scaling factor. Changes the number of indel.

	--slope=SLOPE
		Size of the slope between two structural variants

	--threads=THREADS
		Number of threads used for fastq generation

	-o OUTPUT, --output=OUTPUT
		output folder name

	-h, --help
		Show this help message and exit
```

### Maketube output
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


### vcf2metrics

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
