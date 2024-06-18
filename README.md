# maketube

maketube is an artificial genome generator for Mycobacterium tuberculosis

```
Mandatory inputs :
--reference (fasta) : Fasta files of the reference sequence that will be used as a base
--transposon (bed) : a bed file delimitating the transposon sequences that will be cut and paste across the genome
--nonhomoseq_pool (fasta) : a fasta file with no sequence name delimiter (">"). Maketube will randomely sample across these sequences and concatenate them to create than. As such, it is advised to give sequences that are not present in the reference.
```
```
Other parameters :
--haplotype_count : number of strain to compute for each population (replicate). Default : 10
--pop_count : number of population to compute (number of replicate)
--structural_variants : number of 
--duplication_region_size (int) : the size of the region to be duplicated
```
