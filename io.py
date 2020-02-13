import glob
import os
from .utils import Atom, Residue, ActiveSite


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
        part_sequence = f.readlines()


    return active_site


def write_clustering(filename, clusters):
    """
    Write the clustered ActiveSite instances out to a file.

    Input: a filename and a clustering of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusters)):
        out.write("\nCluster %d\n--------------\n" % i)
        for j in range(len(clusters[i])):
            out.write("%s\n" % clusters[i][j])

    out.close()


def write_mult_clusterings(filename, clusterings):
    """
    Write a series of clusterings of ActiveSite instances out to a file.

    Input: a filename and a list of clusterings of ActiveSite instances
    Output: none
    """

    out = open(filename, 'w')

    for i in range(len(clusterings)):
        clusters = clusterings[i]

        for j in range(len(clusters)):
            out.write("\nCluster %d\n------------\n" % j)
            for k in range(len(clusters[j])):
                out.write("%s\n" % clusters[j][k])

    out.close()
