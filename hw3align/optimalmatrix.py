import sys

def readFasta(filename):
	'''
	reads in a protein sequence from a fasta file.
	returns a string.
	'''
	seq = ""
	for line in open(filename, "r"):
		if line[0] == ">":
			continue
		seq = seq + line.rstrip().upper()
	return seq

seqnames = []
for line in open("hw3align/seqnamesprep.txt", "r"):
	seqnames.append(line.rstrip())

sequences = {}
for filename in seqnames:
	sequences[filename] = readFasta("hw3align/"+filename)
positives = [line.rstrip().split() for line in open("hw3align/Pospairs.txt")]
positives = [(sequences[fileA], sequences[fileB]) for fileA, fileB in positives]

negatives = [line.rstrip().split() for line in open("hw3align/Negpairs.txt")]
negatives = [(sequences[fileA], sequences[fileB]) for fileA, fileB in negatives]
