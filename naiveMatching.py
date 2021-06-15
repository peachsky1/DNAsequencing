#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 20:41:20 2021

@author: jasonlee
"""

import math
import wget
import matplotlib.pyplot as plt
import collections
import random


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def readFastq(filename):
    sequence = []
    qualities = []
    with open(filename) as filehandle:
        while True:
            filehandle.readline()
            seq = filehandle.readline().rstrip()
            filehandle.readline()
            qual = filehandle.readline().rstrip()
            if len(seq) == 0:
                break
            sequence.append(seq)
            qualities.append(qual)
    return sequence, qualities


# generate random genome using random 
def generateReads(genome, numReads, readLen):
	reads = []
	for _ in range(numReads):
		start = random.randint(0, len(genome)-readLen) -1
		reads.append(genome[start: start+readLen])
	return reads
		

# import numpy as np
def downloadFile(url):
    filename = wget.download(url)
    print("{} has been downloaded".format(filename))
    return filename


def readGenome(filename):
    genome = ''

    # >gi|9626243|ref|NC_001416.1| Enterobacteria phage lambda, complete genome
    # GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCG
    # TCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGGC
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                # rstrip removes trailing white space from the ends of the string. in this case it will trip off the new line at the end
                genome += line.rstrip()
    return genome

# p=read sequence
# t=genome
def naive(pattern, textToSearch):
	p = pattern
	t = textToSearch
	occurrences = []
	for i in range(len(t)-len(p)+1):
		match = True
		for j in range(len(p)):
			if not t[i+j] == p[j]:
				match = False
				break
			
		if match:
			occurrences.append(i)
	return occurrences

 
		
	
def main():
	url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa"
# 	downloadFile(url)
	filename = "phix.fa"
	genome = readGenome(filename)
	t = 'AGTTWEASGAGAGERTWETAG'
	p = 'AG'
	print(naive(p,t))
	
	reads = generateReads(genome, 100, 100)
	numMatched = 0
	for r in reads:
		matches = naive(r, genome)
		if len(matches) > 0:
			numMatched += 1
	
	url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR266411_1.first1000.fastq"
# 	downloadFile(url)
	filename = "ERR266411_1.first1000.fastq"
	
	phix_reads, _ = readFastq(filename)
	numMatched = 0
	n = 0
	for r in phix_reads:
# 		possible sequencing error might be contains. use only 30 bases per read
		r = r[:30]
		matches = naive(r, genome)
		if len(matches) > 0:
			numMatched +=1
# 	numMatched
# Out[47]: 476
# Why this is less than half? 2 possible reasons
# 1. sequencing error.
# 2. genome is double stranded. check the reverse complement.
	
	numMatched = 0
	n = 0
	for r in phix_reads:
# 		possible sequencing error might be contains. use only 30 bases per read
		r = r[:30]
		matches = naive(r, genome)
# 		add reverse matches
		matches.extend(naive(reverseComplement(r),genome))
		if len(matches) > 0:
			numMatched +=1
	
# 	numMatched
# Out[53]: 932
# Good improvement!
if __name__ == '__main__':
	main()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	