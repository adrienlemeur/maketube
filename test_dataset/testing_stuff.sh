#!/usr/bin/env bash

sample="H1"
bcftools norm -a -m- "SV1_pop1.vcf.gz" | bcftools view -v snps -s $sample -Ov | grep -v "	0:" > H1.reference.vcf
bgzip -f H1.reference.vcf

#genome_length=4576628
#genome_length_minus_150k=4576628

bcftools norm -a -m- "test_VCF/SV1_pop1_H1.vcf.gz" | bcftools view -v snps -Oz > "test_VCF/SV1_pop1_H1.filtered.vcf.gz"

python3 ../vcf2metrics.py -i "test_VCF/SV1_pop1_H1.filtered.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
	--subtract <(echo -e "H37Rv\t4426628\t4576628\tduplicata_region") \
	--backtrack "SV1_equivalence.bed" \
	--add_col $sample
