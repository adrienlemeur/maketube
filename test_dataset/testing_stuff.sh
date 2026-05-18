#!/usr/bin/env bash

sample="H1"
bcftools norm -a -m- "SV1_pop1.vcf.gz" | bcftools view -v snps -s $sample -Ov | grep -v "	0:" > H1.reference.vcf
bgzip -f H1.reference.vcf

#genome_length=4576628
#genome_length_minus_150k=4576628

bcftools norm -a -m- "test_VCF/SV1_pop1_H1.vcf.gz" | bcftools view -v snps -Oz > "test_VCF/SV1_pop1_H1.filtered.vcf.gz"
#python3 ../foretrack.py --reference "H1.reference.vcf.gz" --backtrack "SV1_equivalence.bed"

python3 ../foretrack.py --reference "H1.reference.vcf.gz" --backtrack "SV1_equivalence.bed" | bcftools view -Oz > "H1.reference.foretracked.vcf.gz"

#python3 ../vcf2metrics.py -i "H1.reference.foretracked.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
#	--backtrack "SV1_equivalence.bed"


python3 ../vcf2metrics.py -i "test_VCF/SV1_pop1_H1.filtered.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
	--subtract <(echo -e "H37Rv\t4426628\t4576628\tduplicata_region") \
	--backtrack "SV1_equivalence.bed" \
	--add_col $sample

exit

python3 ../foretrack.py --reference "H1.reference.vcf.gz" --backtrack "SV1_equivalence.bed"
python3 ../foretrack.py --reference "H1.reference.vcf.gz" --backtrack "SV1_equivalence.bed" | bcftools view -Oz > "H1.reference.foretracked.vcf.gz"

#python3 ../vcf2metrics.py -i "H1.reference.foretracked.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
#	--backtrack "SV1_equivalence.bed" --trace

#python3 ../vcf2metrics.py -i "test_VCF/SV1_pop1_H1.filtered.vcf.gz" --sample $sample --reference "H1.reference.foretracked.vcf.gz" \
#	--subtract <(echo -e "H37Rv\t4426628\t4576628\tduplicata_region")

#python3 ../vcf2metrics.py -i "test_VCF/SV1_pop1_H1.filtered.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
#	--backtrack "SV1_equivalence.bed" \
#	--subtract <(echo -e "H37Rv	394136	544136	duplicated_region\nH37Rv\t4426628\t4576628\tduplicata_region")

exit

rm -rf test
bcftools index test_VCF/SV1_pop1_H1.filtered.vcf.gz
bcftools index H1.reference.foretracked.vcf.gz
bcftools isec test_VCF/SV1_pop1_H1.filtered.vcf.gz H1.reference.foretracked.vcf.gz -p test
#intéressant mais pourquoi ???
exit

python3 ../vcf2metrics.py -i "test_VCF/SV1_pop1_H1.filtered.vcf.gz" --sample $sample --reference "H1.reference.vcf.gz" \
	--subtract <(echo -e "H37Rv\t4426628\t4576628\tduplicata_region") \
	--backtrack "SV1_equivalence.bed" \
	--add_col $sample
