import glob
import itertools
import os
from .utils import AASeq

#This top section of functions help me read in my sequences from the fasta files.
#The final section is to get the lists of tuples of sequence strings belonging to the postiive and negative pairs.


#for all sequences
def read_aa_sequences(dir):
    """
    Read in all of the protein sequences from the given directory.

    Input: directory
    Output: list of protein sequences
    """
    files = glob.glob(dir + '/*.fa')

    aa_sequences = []
    # iterate over each .fa file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):

        aa_sequences.append(read_aa_sequence(filepath))

    print("Read in %d amino acid sequences"%len(aa_sequences))

    return aa_sequences 

#for all sequences' names
def read_aa_sequences_names(dir):
    """
    Read in all of the protein sequences from the given directory.

    Input: directory
    Output: list of protein sequences
    """
    files = glob.glob(dir + '/*.fa')

    aa_sequences = []
    # iterate over each .fa file in the given directory
    for filepath in glob.iglob(os.path.join(dir, "*.fa")):

        aa_sequences.append(read_aa_sequence_name(filepath))

    print("Read in %d amino acid sequences' names"%len(aa_sequences))

    return aa_sequences

#for one  sequence
def read_aa_sequence(filepath):
    """
    Read in a single amino acid sequence given a fasta file

    Input: fasta file path
    Output: amino acid sequence instance
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".fa":
        raise IOError("%s is not a fasta file"%filepath)

    aa_sequence = AASeq(name[0])


    # open .fa file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in itertools.islice(f, 0, 5):
            part_sequence = f.read()
  
    part_sequence = ''.join(part_sequence.split())    
    aa_sequence.partialsequence.append(part_sequence)
       

    return aa_sequence.partialsequence 

#for one  sequence name
def read_aa_sequence_name(filepath):
    """
    Read in a single amino acid sequence given a fasta file

    Input: fasta file path
    Output: amino acid sequence instance
    """
    basename = os.path.basename(filepath)
    name = os.path.splitext(basename)

    if name[1] != ".fa":
        raise IOError("%s is not a fasta file"%filepath)

    aa_sequence = AASeq(name[0])


    # open .fa file
    with open(filepath, "r") as f:
        # iterate over each line in the file
        for line in itertools.islice(f, 0, 5):
            part_sequence = f.read() 
    part_sequence = ''.join(part_sequence.split())    
    aa_sequence.partialsequence.append(part_sequence)
       

    return aa_sequence



