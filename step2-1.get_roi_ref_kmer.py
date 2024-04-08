"""A script to extract ROI(region of interest) reference kmers from the input VCF file, reference genome file, construct a reference kmer fasta to intersect with sample total kmer using kmc simple intersect"""
# written by Hyun Woo kim
# ver1.0 2023-10-11


def get_roi(vcf, window=50):
    roiD = dict() # in to form of {region1 : [(), ()], ... regionN : [(),()]

    exclude_contig = ['chrM', 'chrX', 'chrY']
    #global window
    vcf_file = pysam.VariantFile(vcf, "r")
    for variant in vcf_file:
        contig = variant.contig
        info = variant.info
        mate_contig = info.get('CHR2')# this line can only be applied to ETCHING-generated vcfs
        filt = list(variant.filter)[0]
        if filt == 'PASS' and contig not in exclude_contig and mate_contig not in exclude_contig and not contig.startswith('GL') and not mate_contig.startswith('GL'):
            sv_id = variant.id
            sv_type = info.get('SVTYPE')

            #--- Catching SV 101bp window for each breakends ---#
            breakpoint1 = variant.start
            if breakpoint1 >= window:
                #-- resizing window for copy-number affected types (DEL, DUP)
                if sv_type == 'DUP':
                    start_pos = breakpoint1 +1
                else:
                    start_pos = breakpoint1 - window
            else:
                start_pos = breakpoint1

            if sv_type == 'DEL':
                end_pos = breakpoint1
            else:
                end_pos = breakpoint1 + window + 1
            bnd1 = (contig, start_pos, end_pos)

            breakpoint2 = int(info.get('REPATH').split(':')[-1].split('(')[0]) -1 #this can be only applied to ETCHINg-generated vcfs
            if breakpoint2 > window:
                #-- resizing window for copy-number affected types (DEL, DUP)
                if sv_type == 'DEL':
                    mate_start_pos = breakpoint2 +1
                else:
                    mate_start_pos = breakpoint2 - window
            else:
                mate_start_pos = breakpoint2
            if sv_type == 'DUP':
                mate_end_pos = breakpoint2 -1
            else:
                mate_end_pos = breakpoint2 + window + 1
            bnd2 = (mate_contig, mate_start_pos, mate_end_pos)

            region = [bnd1, bnd2]
            roiD[sv_id] = region

    return roiD

def get_reference_sequence(refobj, contig, start, end):
    """simply get the reference genome sequence inside the given window"""

    ref_seq = refobj.fetch(contig, start, end)

    return ref_seq


def get_roi_reference_kmers(vcf, genome, k=31):
    #roi_ref_kmerD = dict()# in to form of {region1_bp1 : [kmer_sequences], ... regionN_bp2 : [kmer_sequence]}

    roiD = get_roi(vcf)
    ref_genome = pysam.FastaFile(genome)

    roi_ref_kmers = set()
    for var in roiD:
        for i in range(len(roiD[var])):
           # kmers = list()
            contig = roiD[var][i][0]; start_pos = roiD[var][i][1]; end_pos = roiD[var][i][2]
            ref_seq = get_reference_sequence(ref_genome, contig, start_pos, end_pos)
            for j in range(len(ref_seq)-k+1):
                kmer = ref_seq[j:j+k]
                if 'N' in kmer: continue
                rc_kmer = reverse_complement(kmer)
                kmer = kmer_md.lexicographical_comparison(kmer, rc_kmer)
                roi_ref_kmers.add(kmer)
                #kmers.append(kmer)
                #roi_ref_kmerD[f'{var}_bp{i}'] = kmers

    return roi_ref_kmers


def make_reference_fasta(output):
    # use enumerate to make kmer id
    global roi_ref_kmers

    with open(output, 'w') as outfile:
        for i, kmer in enumerate(roi_ref_kmers):
            outfile.write(f'>kmer_{i}\n{kmer}\n')

    return None

### this must be changed to a shell coded process
def kmc_intersect_and_dump(roi_ref_fasta, total_kmer, outputdir):
    # use python subprocess Popen to perform "kmc_tools simple intersect"

    kmc_path = '/home/hyunwoo/programs/kmc/bin/kmc'
    # make kmc db with roi reference
    kmc_cmd = f'{kmc_path} -v -fa -k31 -ci1 -t10 -m8 {roi_ref_fasta} {outputdir}/roi_refkmer {outputdir}'
    md.execute_subprocess(kmc_cmd)

    # intersect with total kmc db
    kmc_intersect_cmd = f'{kmc_path}_tools simple {total_kmer} {outputdir}/roi_refkmer intersect {outputdir}/roi_refkmer_intersect -ocleft'
    md.execute_subprocess(kmc_intersect_cmd)

    # dump into text file
    kmc_dump_cmd = f'{kmc_path}_tools transform {outputdir}/roi_refkmer_intersect dump {outputdir}/roi_refkmer_intersect_dump.txt'
    md.execute_subprocess(kmc_dump_cmd)

    return None

def save_roi_reference_kmerdb_with_pickle(kmer_dump, output, max_count):

    kmer_db = get_kmers_with_count(kmer_dump, max_count)

    with open(output, 'wb') as outfile:
        pickle.dump(kmer_db, outfile)

    print('kmer databse pickle object saved at {}'.format(output))

    return None


def save_kmerdb_with_pickle(kmer_dump, kmer_db_output, max_count_cutoff):
    #kmer_db = get_kmers(kmer_db)
    kmer_db = get_kmers_with_count(kmer_dump, max_count_cutoff)

    with open(kmer_db_output, 'wb') as outfile:
        pickle.dump(kmer_db, outfile)
    print('kmer database object save at {}'.format(kmer_db_output))

    #with open(kmer_sample_output, 'wb') as outfile:
    #    pickle.dump(sample_kmer, outfile)
    #print('sample kmer database object save at {}'.format(kmer_sample_output))


def get_kmers_with_count(kmer_table, max_count):
    kmers = dict()
    print('**** in progress of getting kmer count from dump ****')
    print(':)')
    with open(kmer_table, 'r') as infile:
        for line in infile:
            kmer = line.split('\t')[0]
            count = int(line.strip().split('\t')[-1])
            if count < int(max_count):
                #kmer =  lexicographical_comparison(kmer, reverse_complement(kmer))
                bkmer = kmer_md.convert_to_bit(kmer)
                kmers[bkmer] = count
        print('**** getting kmer count from dump completed! ****')
    return kmers



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for intersecting kmc db with ROI reference kmers ")
    parser.add_argument("-c", "--max_count", help = " maximum kmer count cutoff ", default = 10000,type = int, required = False)
    parser.add_argument("-v", "--vcf", help = "input vcf file", required = True)
    parser.add_argument("-g", "--genome", help = "input reference genome file", required = True)
    parser.add_argument("-f", "--fasta_output", help = "outputput ROI reference fasta file prefix", required = True)
    parser.add_argument("-k", "--total_kmerdb", help = " sample's total kmer kmc file prefix ", required = True)
    parser.add_argument("-o", "--outputdir", help = " output directory ", required = True)
    parser.add_argument("-p", "--output_pickle", help = " output pickle file prefix ", required = True)


    args = parser.parse_args()

    import sys, os, pysam, pickle, subprocess
    from Bio.Seq import reverse_complement
    import dnsv_module as md
    import kmer_module as kmer_md

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    output_fasta = args.outputdir + args.fasta_output
    output_pickle = args.outputdir + args.output_pickle

    roi_ref_kmers = get_roi_reference_kmers(args.vcf, args.genome, 31)
    make_reference_fasta(output_fasta)

    kmc_intersect_and_dump(output_fasta, args.total_kmerdb, args.outputdir)
    
    dump_file = args.outputdir + '/roi_refkmer_intersect_dump.txt'
    save_roi_reference_kmerdb_with_pickle(dump_file, output_pickle, args.max_count)

