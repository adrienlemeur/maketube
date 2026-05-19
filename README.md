<h1 align="center"> Maketube </h1>
<p align="center">
   <img src="/maketube_logo.png" width="200" height="100">
</p>

maketube introduces structural variants and short variants to a _Mycobacterium tuberculosis_ reference genome sequence. Variants found by aligning these genomes to the reference genome can be compared to the original ones in order to study the impact of structural variants on variant calling pipeline performances. For exemples of benchmarking pipelines, and the associated article check the sister repo adrienlemeur/maketube\_supplemental. Code at adrienlemeur/maketube. Only maketube genome adrienlemeur/maketube\_genomes.

If you have trouble running these tools / suggestions / comments, feel free to send me a mail at : alemeur at biophylo dot com

### Table of contents

<!--ts-->
   - [Installation](#install)
     - [Maketube dependencies](#install_maketube)
     - [vcf2metrics dependencies](#install_vcf2metrics)
     - [Installation](#installation_procedure)
     - [Container](#container)
   - [Quickstart](#quickstart)
   - [Pipeline](#pipeline)
   - [Usage](#usage)


### <a name="install"></a>Installation

#### <a name="install_maketube"></a>Maketube:
  ##### Dependencies:

  - [R](https://www.r-project.org/) (v4.1.2)
  - [ape](https://cran.r-project.org/web/packages/ape/index.html) (v5.8)
  - [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) (v4.2-36)
  - [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (v2.62.0)
  - [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (v1.1.4)
  - [jackalope](https://cran.r-project.org/web/packages/jackalope/index.html) (v1.1.5)

#### <a name="install_vcf2metrics"></a>vcf2metrics.py
  ##### Dependencies
  - [cyvcf2](https://brentp.github.io/cyvcf2/) (0.30.18)
  - [numpy](https://numpy.org/) (1.26.4) (⚠ 2.0 produces error)
  - [argparse](https://pypi.org/project/argparse/) (3.2)
  - math, re, sys, gc & os

##### <a name="installation_procedure"></a>Installation procedure

Install the dependencies & >>
  ```
  git clone https://github.com/adrienlemeur/maketube.git
  cd maketube
  echo "PATH=\"$(pwd)\"/:$PATH" >> ~/.bashrc && source ~/.bashrc
  ```

#### <a name="container">Container
Don't want to install all these pesky packages and their dependencies ? There is a [container](https://hub.docker.com/r/alemeur/maketube) !

#### Singularity
```
singularity pull maketube.img docker://alemeur/maketube:latest
#run maketube.R
singularity run maketube.img
#run vcf2metrics.R
singularity exec maketube.img "/usr/local/bin/maketube/vcf2metrics.py"
```
##### Docker
```
docker pull alemeur/maketube:latest
#run maketube.R
docker run maketube:latest
#run vcf2metrics.R
# in /usr/local/bin/maketube/vcf2metrics.py
```

### <a name="quickstart"></a>Quickstart
#### Benchmark pipeline

 <img src="/maketube_pipeline.png">
A. Maketube generates artificial genomes

B. The user extract variants from these genomes using genome wide alignment software (minimap2, nucmer) OR by generating an _in silico_ sequencing run and a pipeline (alignment + variant caller).

C. vcf2metrics.py compare the list of variants to the original list of variants and provide precision and recall

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
Comparing variants found in Haplotype 1 (H1) of population 1. 
```

bcftools view --samples H1 maketube_run/SV1/SV1_pop1.vcf.gz > my_reference_vcf.vcf.gz

vcf2metrics.py -i my_sample_vcf.vcf.gz \
		--reference my_reference_vcf.vcf.gz \
		--backtrack maketube_run/SV1/SV1_equivalence.bed \
		--bed maketube_run/SV1/SV1_SV.bed
```

### <a name="usage"></a>Complete usage

#### Full description
```
Usage: ./maketube.R [options]

Options:
	--reference=file
		reference genome sequence to evolve (fasta file with header)

	--transposon=file
		Inserting Sequences (transposon-like sequences) initial positions in the reference sequence (bed) : CHROM_NAME\tTRANSPOSON_START\tTRANSPOSON_STOP
		You can get these by blasting the genome sequence for the IS/transposon-like elements.

	--nonhomoseq_pool=file
		a list of kmers (fasta). All kmers should be the same size. Maketube subsamples this list at random and concatenates them to build ancetral-like regions (regions that are present in the evolved genome but not in the reference sequence). Check REF/nonH37Rv_pool_sequence.fasta.gz for reference.

	--duplication_region_size=string
		size of the duplicated region. Default : 150kbp

	--haplotype_count=numeric
		number of haplotype (genome from the same population) to compute (>=1). Final number of haplotype is haplotype_count * pop_count * structural_variants. Default = 10.

	--pop_count=numeric
		number of populations to compute (>=1).

	--structural_variants=numeric
		number of set of independant structural variants. All the populations from the same set of structural variants share the same structural variants at the same positions.. Default = 2.

	--deletion_count=numeric
		number of deletion region to remove. Default = 3. 

	--unmuted
		Maketube will only write the unmutated evolved genome (with structural variants but no snps and no indels) if the flag is specified.

	-p PREFIX, --prefix=PREFIX
		artificial genome name prefix.

	--mutation_rate=string
		substitution/site/generation. Default to 1.23e-7.

	--effective_pop=EFFECTIVE_POP
		the effective population size. Default to 2000.

	--TCAGI=numerics
		the rate of the different nucleotides (TCAG) + invariants (I). Comma separated list of numerics. Default to "0.172,0.329,0.172,0.329,0"

	--ABCDEF=numerics
		the parameters of the GTR model (check the jackalope package for reference). Comma separated list of numerics. Default to "0.65,0.05,0.21,0.27,0.02,0.65".

	--indel_scaling_factor=numerics
		indel/snp ratio. Default to 0.125, ie. 1 indel will be created for every 8 SNPs. 0 to remove indels.

	--slope=numerics
		when writing the structural variant neighbouring positions, the number of nucleotide around each structural variants for a snps or indels to be considered "near" it. You may specify comma separated value. Default = "50,300" : maketube will write 2 entries in the bed file for structural variants neighbouring region. One for SNP <= 50bp around the structural variant and one for SNP <= 300bp.

	--threads=THREADS
		Number of threads used for fastq generation

	-o OUTPUT, --output=OUTPUT
		output folder name

	-h, --help
		Show this help message and exit
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


#### vcf2metrics

vcf2metrics.py compares a sample vcf to the reference VCF built by maketube using the backtrack partition. It can also compare a test VCF to a reference VCF without backtracking. 

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
