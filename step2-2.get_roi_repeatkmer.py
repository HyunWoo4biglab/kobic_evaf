# A script to extract ROI refkmers identical with repeatMasker-annotated repeat sequences
# Written by Hyun Woo Kim
# update-note
# updated to use kmc intersect

def make_roi_repeatKmer(repeatkmer_db, roikmer_db, outputdir):
    # run KMC intersect to extract ROI repeat kmers
    kmc_path = cf.kmc
    kmc_intersect_cmd = f'{kmc_path}_tools simple {roikmer_db} {repeatkmer_db} \
                        intersect {outputdir}roi_repeatkmer_intersect -ocleft'
    print("----KMC intersection running")
    print(kmc_intersect_cmd)
    md.execute_subprocess(kmc_intersect_cmd)

    # make roi repeat kmer dump file
    kmc_dump_cmd = f'{kmc_path}_tools transform {outputdir}roi_repeatkmer_intersect \
                    dump {outputdir}roi_repeatkmer_intersect_dump.txt'
    print("----KMC dump running")
    print(kmc_dump_cmd)
    md.execute_subprocess(kmc_dump_cmd)

    return f'{outputdir}roi_repeatkmer_intersect_dump.txt'



def save_kmerdb_with_pickle(kmer_dump, outputdir, max_count_cutoff=10000):
    kmer_db = kmer_md.get_kmers_with_count(kmer_dump, max_count_cutoff)
    output = f'{outputdir}roi_repeatkmer_intersect.pickle'

    print("----generating kmer db pickle object")
    with open(output, 'wb') as outfile:
        pickle.dump(kmer_db, outfile)
    print('kmer database object saved at {}'.format(output))



def main(repeatkmer_db, roikmer_db, outputdir):
    if not outputdir.endswith('/'):
        outputdir += '/'
    kmer_dump = make_roi_repeatKmer(repeatkmer_db, roikmer_db, outputdir)
    save_kmerdb_with_pickle(kmer_dump, outputdir)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for finding kmers having identical sequence with repeat regions")
    parser.add_argument("-i", "--roi_kmer", help = "input kmer KMC db prefix", required = True)
    parser.add_argument("-o", "--output", help = "output directory", required = True)
    parser.add_argument("-r", "--repeat_kmer", help = "repeat kmer KMC db prefix ", required = True)


    args = parser.parse_args()

    import sys, os, pysam, pickle, re, subprocess
    import config as cf
    import dnsv_module as md
    import kmer_module as kmer_md
    from Bio.Seq import reverse_complement

    main(args.repeat_kmer, args.roi_kmer, args.output)
