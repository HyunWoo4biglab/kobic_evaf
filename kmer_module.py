# A PYTHON MODULE HANDLING KMERS
# HYUN WOO KIM
# 03-11-2024
import re

def get_kmers_with_count(kmer_table, max_count):
    """
    | make a kmer dictionary as {kmer_bit : kmer_count} from a kmc dump file
    : kmer_table : input kmc dump.txt file
    : max_count : maximum kmer count threshold
    """

    kmers = dict()
    print('**** in progress of getting kmer count from dump ****')
    print(':)')
    with open(kmer_table, 'r') as infile:
        for line in infile:
            kmer = re.split('\s+', line.strip())[0]
            count = int(line.strip().split('\t')[-1])
            if count < int(max_count):
                bkmer = convert_to_bit(kmer)
                kmers[bkmer] = count
        print('**** getting kmer count from dump completed! ****')
    return kmers


def save_kmerdb_with_pickle(kmer_dump, outputdir, max_count_cutoff=10000):
    """
    | Save a kmer dictionary as a pickle object from {@get_kmers_with_count}
    """

    kmer_db = get_kmers_with_count(kmer_dump, max_count_cutoff)
    output = f'{outputdir}roi_repeatkmer_intersect.pickle'

    print("----generating kmer db pickle object")
    with open(output, 'wb') as outfile:
        pickle.dump(kmer_db, outfile)
    print('kmer database object saved at {}'.format(output))

def lexicographical_comparison(str1, str2):
    """
    | Compare the lexicographical order of two strings and return the smaller one
    """
    # Iterate over the characters of the two strings simultaneously
    for char1, char2 in zip(str1, str2):
        # Compare the characters using the built-in ord() function
        if ord(char1) < ord(char2):
            return str1
        elif ord(char1) > ord(char2):
            return str2

    # If the two strings are equal up to the length of the shorter string,
    # return the shorter string
    if len(str1) < len(str2):
        return str1
    else:
        return str2


def convert_to_bit(seq):
    """
    | Convert a 31-mer nucleotide sequence to a 62-bit integer
    """
    seq = seq.upper()
    # Define a dictionary to map nucleotides to bits
    bit_dict = {'A': 0b00, 'C': 0b01, 'G': 0b10, 'T': 0b11, 'N' : 0b00}

    # Initialize an empty 62-bit integer
    int_seq = 0b0

    # Iterate over the sequence and add the bits to the integer
    for i in range(len(seq)):
        int_seq <<= 2  # shift the bits left by 2 to make room for the next nucleotide
        int_seq |= bit_dict[seq[i]]

    return int_seq

