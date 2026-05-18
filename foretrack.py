#!/usr/bin/env python3

from cyvcf2 import VCF
import numpy as np
import argparse
import gzip
import re, sys, gc, os
import math

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(
			prog="vcf2metrics",
			description="""foretrack takes the maketube reference vcf and output (stdout) the reference-equivalent VCF""",
			epilog="Written by Adrien Le Meur")

parser.add_argument('--reference', type=str, nargs='+', required=True, help='reference vcf')
parser.add_argument('--backtrack', type=str, nargs='+', required=True, help='bed backtrack file produced by maketube')

parser.add_argument('-o', type=str, nargs='+', required=False, default = "./", help='output directory')


#sort a dictionary by key
def sort_dict(my_dict):
	return(dict(sorted(my_dict.items(), key=len, reverse=True)))

#import a vcf to a dict object
#key = start_referenceVariant_alternativeVariant
#one entry for each variant
def vcf2dict(vcf, correction:int):
	vcf_dict = {}

	#already bed dict key, so de facto unique
	for v in vcf:
		#should output the variant anyway
		if v.FILTERS == [] or v.FILTERS[0] == 'PASS':
			variant_start = int(v.POS) - correction
			ALT = v.ALT[0]

			variant_key = "_".join( [ str(variant_start), v.REF, ALT])

			vcf_dict[variant_key] = {}
			vcf_dict[variant_key]['START'] = variant_start
			vcf_dict[variant_key]['CLASS'] = v.var_type
			vcf_dict[variant_key]['END'] = variant_start+len(ALT)
			vcf_dict[variant_key]['REF'] = v.REF
			vcf_dict[variant_key]['ALT'] = ALT
	return(vcf_dict)

#import a backtrack bed file into a dict object.
#key = transposon|insertion|deletion_start_stop
#for transposons and insertions :
# line 1 and 2 are the structural variant border coordinates before structural variations
# line 1 and 2 are the structural variant border coordinates after structural variations
# for deletions, there are no borders after structural variations (the fragment disappeared)
def read_backtrack(file):
	bed_dict = {}

	with open(file) as bed_lines:
		for line in bed_lines:
			if line.startswith('#'):
				continue
			line = line.rstrip().split()
			key = "_".join([ line[3], line[1], line[2] ])

			if(line[3] == "transposon"):
				key_ins = "_".join([ line[3], line[1], line[2], "ins"])
				key_del = "_".join([ line[3], line[1], line[2], "del"])
#				bed_dict[key] = {"class":line[3], "start":int(line[1]), "stop":int(line[2]), "ins_start":int(line[4]), "ins_stop":int(line[5])}
				bed_dict[key_del] = {"class":"deletion", "start":int(line[1]), "stop":int(line[2]) }
				bed_dict[key_ins] = {"class":"insertion", "start":int(line[4]), "stop":int(line[5])}
			elif(line[3] == "insertion"):
				bed_dict[key] = {"class":line[3], "start":int(line[4]), "stop":int(line[5])}
			elif(line[3] == "deletion"):
				bed_dict[key] = {"class":line[3], "start":int(line[1]), "stop":int(line[2])}
	return(bed_dict)

#def foretrack_backtrack_file(bed_dict):

#	for SV_A in bed_dict:
#		for SV_B in bed_dict:
#			if(bed_dict[SV_A]["class"] != "transposon"):
#				if(bed_dict[SV_A]["start"] <):

def foretrack(vcf_dict, backtrack_dict):
	new_dict = {}

	for variant_key in vcf_dict:
		ovPOS = vcf_dict[variant_key]["START"]

		for interval in backtrack_dict:
			start = backtrack_dict[interval]["start"]
			stop = backtrack_dict[interval]["stop"]
			#stop + 1 because of technical debt
			if( (vcf_dict[variant_key]["START"] >= stop) and (backtrack_dict[interval]["class"] == "insertion") ):
#				ovPOS = ovPOS + backtrack_dict[interval]["ins_start"] - backtrack_dict[interval]["ins_stop"] + 1
				ovPOS = ovPOS + backtrack_dict[interval]["start"] - backtrack_dict[interval]["stop"] + 1

			if( (vcf_dict[variant_key]["START"] >= stop) and (backtrack_dict[interval]["class"] == "deletion") ):
#				ovPOS = ovPOS - start + backtrack_dict[interval]["stop"] + 1
				ovPOS = ovPOS - start + backtrack_dict[interval]["stop"] + 1

			if(backtrack_dict[interval]["class"] == "transposon"):
				start = backtrack_dict[interval]["start"]
				stop = backtrack_dict[interval]["stop"]
				insertion_start = backtrack_dict[interval]["ins_start"]
				insertion_stop = backtrack_dict[interval]["ins_stop"]


				if( vcf_dict[variant_key]["START"] <= insertion_stop & vcf_dict[variant_key]["START"] >= insertion_start ):
					print("this happened")
					continue
				#I : a transposon that was after the variant jumped to a position before the variant
				# remove the size of the transposon to the variant position
				elif((vcf_dict[variant_key]["START"] > start) and (vcf_dict[variant_key]["START"] < insertion_start) ):
#				elif((vcf_dict[variant_key]["START"] > start) and (vcf_dict[variant_key]["START"] < insertion_start) and (insertion_start > start) ):
#					ovPOS = ovPOS + stop - start + 1
					ovPOS = ovPOS + stop - start + 1

				#II : a transposon that was before the variant jumped to a position after the variant
				# add the size of the transposon to the variant position
				elif((vcf_dict[variant_key]["START"] < start) and (vcf_dict[variant_key]["START"] > insertion_start)):
#				elif((vcf_dict[variant_key]["START"] < start) and (vcf_dict[variant_key]["START"] > insertion_start) and (insertion_start < start)):
					ovPOS = ovPOS - insertion_stop + insertion_start
				
			#print(vcf_dict[variant_key]["START"], ovPOS, vcf_dict[variant_key]["REF"], vcf_dict[variant_key]["ALT"], "\t".join(map(str, [backtrack_dict[interval][i] for i in backtrack_dict[interval]])), sep="\t")


		new_key = "_".join([ str(ovPOS), vcf_dict[variant_key]["REF"], vcf_dict[variant_key]["ALT"] ])
		new_dict[new_key] = {}
		new_dict[new_key]["START"] = ovPOS
		new_dict[new_key]["END"] = ovPOS + len(vcf_dict[variant_key]["ALT"])
		new_dict[new_key]["CLASS"] = vcf_dict[variant_key]["CLASS"]
		new_dict[new_key]["REF"] = vcf_dict[variant_key]["REF"]
		new_dict[new_key]["ALT"] = vcf_dict[variant_key]["ALT"]
		new_dict[new_key]["PREVIOUS_START"] = vcf_dict[variant_key]["START"]

	return(new_dict)

def vcf_dict2vcf(vcf_dict, vcf):
	new_dict = {}

	for variant_key in vcf_dict:
		ovPOS = vcf_dict[variant_key]["START"]

#THE MAIN STARTS 

args = parser.parse_args()

output_directory = str(args.o[0])

if(args.o):
	output_directory = str(args.o[0])
	if not os.path.exists(output_directory):
		os.mkdir(output_directory)

backtrack_dict = read_backtrack(args.backtrack[0])

reference_vcf = VCF(args.reference[0])
reference_chromosome = reference_vcf.seqnames[0]
reference_dictionnary = sort_dict(vcf2dict(reference_vcf, 0))

foretracked_dict = foretrack(reference_dictionnary, backtrack_dict)

handle = open(args.reference[0], "rb")

with gzip.open(args.reference[0],'rt') as f:
	for line in f:
		if "#" in line:
			print(line.rstrip())

for variant in foretracked_dict:
	print(reference_chromosome, foretracked_dict[variant]["START"], ".", foretracked_dict[variant]["REF"], foretracked_dict[variant]["ALT"], 441453, "PASS", "NS=10;AC=1;AN=1", "GT:GQ", "1:441453", sep="\t")

