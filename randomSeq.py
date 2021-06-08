import random
import wget


def randomSeqGen(n):
    dnaMolecule = 'ATGC'
    seq = ''
    for _ in range(n):
        # random.seed(8)
        # this causes bug since every char will be fixed by seed
        seq += random.choice(dnaMolecule)

    # print(seq)
    return seq


def longestCommonPrefix(s1, s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i += 1
    return s1[:i]


def matchVerification(s1, s2):
    if not len(s1) == len(s2):
        return False
    for i in range(len(s1)):
        if not s1[i] == s2[i]:
            return False
    return True


def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


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


def main():
    random.seed(13)
    seq1 = randomSeqGen(10)
    print(seq1)

    # random.seed(13)
    seq2 = ''.join([random.choice('ATGC') for _ in range(10)])
    print(seq2)
    #     These 2 methods(inline and method) will generate same 'pattern' of the sequence due to the same seed.
    #     CATTTAGCGT
    #     Systematically, as long as I don't give new seed value, the 'patter' will continue.
    commonPrefix = longestCommonPrefix('GGTTTTTTAT', 'GGTTTAAAAAAA')
    print(commonPrefix)

    matchBool1 = matchVerification(seq2, seq2)
    matchBool2 = matchVerification(seq1, seq2)
    print(matchBool1)
    print(matchBool2)

    # reverseComplement - use Dictionary
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    print(complement['C'])
    complementMethod = reverseComplement(seq1)
    # imagine there is a double stranded DNA,
    # seq is one of those strands from top to bottom,
    # complement reverse is the other strand from bottom to top
    print(seq1)
    print(complementMethod)

    url = 'http://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa'
    # filename = wget.download(url)
    filename = 'lambda_virus.fa'
    print(filename)
    genome = readGenome('lambda_virus.fa')
    print(genome[:100])
    print(len(genome))
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for base in genome:
        counts[base] += 1
    print(counts)
    # {'A': 12334, 'T': 11986, 'G': 12820, 'C': 11362}
    import collections
    cc = collections.Counter(genome)
    print(cc)


if __name__ == '__main__':
    main()
