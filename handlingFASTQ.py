import math
import wget
import matplotlib.pyplot as plt
# import numpy as np
def downloadFile(url):
    filename = wget.download(url)
    print("{} has been downloaded".format(filename))
    return filename


def QtoPhred33(Q):
    return (int(math.ceil(Q)) + 33)


def phred33toQ(C):
    return ord(C) - 33


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


# it gives frequency of quality base
def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33toQ(phred)
            hist[q] += 1
    return hist



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
    url = "http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/SRR835775_1.first1000.fastq"
    # downloadFile(url)
    filename = "SRR835775_1.first1000.fastq"
    seqs, quals = readFastq(filename)
    print(seqs[:5])
    print(quals[:5])
    # print(phred33toQ('@'))
    lst = createHist(quals)
    print(lst)
    print(range(50))
    plt.bar(range(len(lst)),lst)
    plt.show()

    gc = findGCByPos(seqs)
	plt.plot(range(len(gc)),gc)
	plt.show()



if __name__ == '__main__':
    main()
