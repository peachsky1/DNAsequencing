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

def naive_2mm(pattern, textToSearch):
	p = pattern
	t = textToSearch
	occurrences = []
	for i in range(len(t)-len(p)+1):
		match = True
		tolarence = 0
		for j in range(len(p)):
			if not t[i+j] == p[j]:
				tolarence +=1
				if tolarence > 2:
					match = False
					break
		if match:
			occurrences.append(i)
	return occurrences

def findGCByPos(reads):
	gc = [0]*100
	tot = [0]*100
	for r in reads:
		for i in range(len(r)):
			if r[i] =='C' or r[i] == 'G':
				gc[i]+=1
			tot[i]+=1
	for i in range(len(gc)):
		if tot[i] !=0:
			gc[i] /= float(tot[i])
	return gc
		
	
def main():
	
	lambdaURL = "https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa"
	downloadFile(lambdaURL) 
	filename = "lambda_virus.fa"
	genome = readGenome(filename)
	
	
	numMatched = 0
	n = 0
	matches = naive("AGGT", genome)
	matchesRev = naive(reverseComplement("AGGT"),genome)
# 	print(matches)
	answer1 = len(matches) + len(matchesRev)

	numMatched = 0
	n = 0
	matches = naive("TTAA", genome)
	matchesRev = naive(reverseComplement("TTAA"),genome)
# 	print(matches)
	answer2 = len(matches) 

	
	tempSeq = "ACTAAGT"
	matches = naive(tempSeq,genome)
	matchesRev = naive(reverseComplement(tempSeq),genome)
	matches[0]
	matchesRev[0]
	answer3 = matchesRev[0]
	
	tempSeq = "TTCAAGCC"
	answer5 = len(naive_2mm(tempSeq, genome))
	
	tempSeq = "AGGAGGTT"
	answer6 = naive_2mm(tempSeq, genome)[0]
	
	q6URL = "https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq"
	downloadFile(q6URL)
	filename = "ERR037900_1.first1000.fastq"
	human_reads, _ = readFastq(filename)
	
	gc = findGCByPos(human_reads)
	plt.plot(range(len(gc)),gc)
	plt.show()
	for i in range(len(gc)):
		if 0.1 > gc[i]:
			print(i)
	human_reads[66]
	
	




	
	url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/phix.fa"
	downloadFile(url) 
	filename = "phix.fa"
	genome = readGenome(filename)

	t = 'AGTTWEASGAGAGERTWETAG'
	p = 'AGAA'
	print(naive_2mm(p,t))
	
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
		print(matches)
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	