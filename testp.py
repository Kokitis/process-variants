

from pprint import pprint
import shlex

string = """python3 {script} \
		--tsv-in {infile} \
		--vcf-out {outfile} \
		--normal-sample-name {normalid} \
		--tumor-sample-name {tumorid} \
		--individual-mutation-tools {mutation-tools} \
		--emit-all \
		--phred-scale"""
string1 = [i.strip() for i in string.split('\t\t')]
pprint(string1)
pprint(shlex.split(string))