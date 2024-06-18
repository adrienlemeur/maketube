
# maketube
<p align="center">
<img src="/maketube_logo.png" width="200" height="100">
</p>

maketube is an artificial genome generator for Mycobacterium tuberculosis

Table of contents
=================

<!--ts-->
   * [Installation](#install)
   * [Usage](#usage)

### <a name="install"></a>Installation

Dependencies & (versions used for developpement) (bioconducters packages are marked with *):
- [R](https://www.r-project.org/) (v4.1.2)
- [ape](https://cran.r-project.org/web/packages/ape/index.html) (v5.8)
- [seqinr](https://cran.r-project.org/web/packages/seqinr/index.html) (v4.2-36)
- *[Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) (v2.62.0)
- [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (v1.1.4)
- [jackalope](https://cran.r-project.org/web/packages/jackalope/index.html) (v1.1.5)

Don't want to install all these pesky packages and their dependencies ? There is a [container](https://hub.docker.com/r/alemeur/maketube) !

### <a name="usage"></a>Usage

```
./maketube.R \
  --reference (fasta) #Fasta files of the reference sequence that will be used as a base
  --transposon (bed) #bed file delimitating the transposon sequences that will be cut and paste across the genome
  --nonhomoseq_pool (fasta) #a fasta file with no sequence name delimiter (">"). Maketube will randomely sample across these sequences and concatenate them to create than. As such, it is advised to give sequences that are not present in the reference. \
  --haplotype_count #number of strain to compute for each population. Default : 10 \
  --pop_count : number of population to compute. Default : 2 \
  --structural_variants : number of different set of structural variant by population. Default : 2 \
  #Total number of strain : haplotype_count x pop_count x structural_variants \
  --unmuted (flag) : write the fasta sequence without the mutation. Default : false (no flag) \
  --duplication_region_size (int) : the size of the region to be duplicated \
  --prefix (string) : prefix of the generated sample \
  -o/--output (string) : name directory to put the results in. Default : maketube_results \
  --mutation_rate (float) : in substitution by generation by site. Default : 1.23*10^-7 \
  --effective_pop (int) : the dreadful Ne. Default : 50000*10^10 \
  --TCAGI (R vector) : parameters of the GTR model, in the same order as specified in the flag. Default : c(0.172, 0.329, 0.172, 0.329, 0.99978) \
  --ABCDEF (R vector): parameters of the GTR model, in the same order as specified in the flag. Default : c(0.65, 0.05, 0.21, 0.27, 0.02, 0.65) \
  --scaling_factor (int) : empirical factor to get a number of indel equivalent to 1/8 the number of SNP, as found in natural strains. Only affects the number of indels. Default : 5*10^-14. \
  --slope (int, > 0) : Size of the distance between 2 structural variants. Default : 300. \
  --threads (int) : Maximum number of threads to use for generating fastq. Default : 4.
```
