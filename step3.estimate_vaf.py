
def estimate_vaf(vcf, bam, ref, output, model, mode, depth=None, kmer=None, repeat_kmer=None):
    """
    | A FUNCTION TO ESTIMATE/CALCULATE VARIANT ALLELE FREQUENCY FOR A GIVEN VARAINT BP WINDOW
    : vcf : input vcf file
    : bam : input bam file (ETCHING-filtered bam or total bam)
    : ref : reference genome fasta file
    : mode : VAF calculation mode (read / kmer)
    """

    lowmap_bed = '/home/hyunwoo/data_set/project/svProject/genomic_context/giab_genome_stratification/GRCh38@all/Mappability/GRCh38_nonunique_l100_m2_e1.bed'
    repeat_bed = '/home/hyunwoo/data_set/project/svProject/genomic_context/repeatMasker/hg38_nonHeader.fa.out'
    mappability = '/home/hyunwoo/data_set/project/svProject/mappability_score/genmap/grch38.p14/grch38.p14_genmap_score.txt'
    #mappability = '/home/hyunwoo/data_set/project/svProject/mappability_score/genmap/hg38_gatk_bundle/hg38_gatk_bundle_genmapscore.txt'
    lcr = ''
    #outfile = open(output, 'w')
    
    outvcf = vcf.split('.vcf')[0] + '.addVAF.vcf'
    dnoutvcf = outvcf.split('.vcf')[0] + '.dnsv.vcf'
    outvcffile = open(outvcf, 'w')
    dnoutvcffile = open(dnoutvcf, 'w')

    if mode == 'kmer':
        print('STARTING KMER DATABASE LOADING ---')
        with open(kmer, 'rb') as infile:
            kmercount_d = pickle.load(infile)
        print('KMER DATABASE LOADED FROM {}'.format(kmer))
        
        print('STARTING MAPPABILITY SCORE LOADING ---')
        mapscore_d = md.get_mappability_score(mappability)
        print('MAPPABILITY SCORE LOADED FROM {}'.format(mappability))

        print('STARTING REPEAT KMER LOADING ---')
        with open(repeat_kmer, 'rb') as repeat_kmers:
            repeatkmer_l = pickle.load(repeat_kmers)

        lowmap_d = md.get_difficult_genomic_context(lowmap_bed)
        repeat_d = md.get_repeat_region(repeat_bed)

        print('ALL REQUIRED INPUT DATA LOADED! --- STARTING VAF ESTIMATION.')

    if mode == 'read':
        depth_d = md.get_depth(depth)

    refobj = pysam.FastaFile(ref)
    bamobj = pysam.AlignmentFile(bam, 'rb')
    parser = md.VariantParser(vcf)#, ref, mode, 30)
    parser.parse()#ref, kmercount_d, mapscore_d, repeatkmer_l, repeat_d, lowmap_d, 31)
    v_l = parser.get_variant_list()
    header_l = parser.get_header()
    for h in header_l:
        outvcffile.write(h)
        dnoutvcffile.write(h)

    if mode == 'kmer':
        df = pd.DataFrame(columns=['sv_id', 'pair_type', 'average_read_coverage', 'split_read', 'sv_length', 'vaf_prior'])
    elif mode == 'read':
        df = pd.DataFrame(columns=['sv_id', 'pair_type', 'average_read_coverage', 'split_read', 'sv_length', 'vaf_prior', 'discordant'])

    for v in v_l: #v is each variant class object
        attr_d = vars(v)
        vid = attr_d['vid']; vtype = attr_d['vtype']; 

        vwindow_l = md.get_variant_window(v, win=30)
        v.set_variant_window(vwindow_l)

        vaf_l = list(); vcov_l = list()
        split_l = list(); discordant_l = list(); vread_l = list(); vread_count_l = list()
        call_flag_l = list()

        for i, w in enumerate(vwindow_l):
            if i == 0:  pair_type = 'BP1'
            else:   pair_Type = 'BP2'

            if mode == 'read':
                vcov = md.get_variant_read_coverage(w, depth_d)
                split, discordant = md.extract_vreads(bamobj, w, v, True)
                vread = split.union(discordant)
                call_flag = 'callable'

            elif mode == 'kmer':
                kmer_l, blackkmer_l, count_l, mapscore_l = md.get_varaint_window_kmers(w, refobj, kmercount_d, mapscore_d, repeatkmer_l, repeat_d, lowmap_d, 31)
                vcov = md.estimate_coverage(count_l, mapscore_l)
                #print(v.vid, v.vtype, vcov)
                if len(count_l) == 0:
                    print(v.vid, v.vtype, w)
                split = md.extract_vreads(bamobj, w, v, False)
                vread = split

                #print(v.vid, w, len(kmer_l))
                try:
                    blackkmerprop = len(blackkmer_l) / len(kmer_l)
                except ZeroDivisionError:
                    blackkmerprop = 0.0

                if (blackkmerprop) > 0.8:               
                    call_flag = 'uncallable'
                else:
                    call_flag = 'callable'

                v.call_flag_l.append(call_flag)

            vread_count = len(vread)
#            vread_l.append(vread)
#            vread_count = len(vread_l)
#            vread_count_l.append(vread_count)
#            vcov_l.append(vcov)
#            split_l.append(split); discordant_l.append(discordant)

            try:
                vaf = round(vread_count / vcov + vread_count, 4)
                #vaf = round(vread_count / vcov, 4)
            except ZeroDivisionError:
                vaf = 0.0

            v.vaf_l.append(vaf)
            v.vcov_l.append(vcov)
            #v.vread_l.append(vread)
            #v.vread_count_l.append(vread_count)
            v.vsplit_l.append(len(split))#


            if mode == 'read':
                v.vdiscordant_l.append(discordant)#

            ##--  write_df for RF regression
            #df = df.append({'coverage': vcov, 'split_read' : vread_count, 'sv_length' : v.vlen, 'vaf_prior' : vaf})
            if mode == 'kmer':
                data = {'sv_id' : vid, 'pair_type' : pair_type, 'average_read_coverage': vcov, 'split_read' : len(split), 'sv_length' : v.vlen, 'vaf_prior' : vaf}
            elif mode == 'read':
                data = {'sv_id' : vid, 'pair_type' : pair_type, 'average_read_coverage': vcov, 'split_read' : len(split), 'sv_length' : v.vlen, 'vaf_prior' : vaf, 'discordant_read' : discordant}
            df.loc[len(df)] = data

    ##-- PREDICT DISCORDANT READ WITH RANDOM FOREST
    print(df.head(10))
    print(df.isnull().any())
    if mode == 'kmer':
        regressor = load(model)
        discordant_d = md.predict_discordant(regressor, df)# {sv_id : [bp1_discordant_count, bp2_discordant_count], ...}
        for v in v_l:
            discordant_l = discordant_d[v.vid] # ['bp1_discordant_count', 'bp2_discordant_count']
            #v.vdiscordant_l.append(discordant_l)
            v.vdiscordant_l = discordant_l
            for i, d in enumerate(discordant_l):
                vread_count = v.vsplit_l[i]
                #print('predicted discordant', d)
                #print('coverage', v.vcov_l[i])
                #print('split read', vread_count)
                try:
                    vaf_update = round((vread_count + d) / (v.vcov_l[i] + vread_count + d), 4)
                except ZeroDivisionError:
                    vaf_update = 0.0
                v.vaf_l[i] = vaf_update

        write_table(v_l, output)
        write_vcf(v_l, outvcffile, dnoutvcffile)


    #df.to_csv(output)
    #outfile.close()
    outvcffile.close()
    refobj.close()
    bamobj.close()

#    return vaf_df


def write_vcf(v_l, outfile, dnoutfile):
    for v in v_l:
        v.add_info('eVAF', v.vaf_l)
        v.add_info('SR', v.vsplit_l)
        v.add_info('DR', v.vdiscordant_l)
        v.add_info('eCOV', v.vcov_l)
        print(v.vinfo)
        line = v.make_vcf_line()
        outfile.write(f'{line}\n')
        vaf_filt_l = list(filter(lambda x: x >= 0.2 and x <= 0.8, v.vaf_l))
        #vaf_filt_l = [af for af in v.vaf_l if (af >= 0.2) and (af <= 0.8)]
        if len(vaf_filt_l) == 2:
            dnoutfile.write(f'{line}\n')


def write_table(v_l, output):
    with open(output, 'w') as outfile:
        outfile.write('SVID\tSVTYPE\tSVLEN\tBP1_SR\tBP1_COV\tBP1_DR\tBP1_VAF\tBP2_SR\tBP2_COV\tBP2_DR\tBP2_VAF\n')
        for v in v_l:
            vid = v.vid; vtype = v.vtype; vlen = v.vlen
            split_l = v.vsplit_l; cov_l = v.vcov_l; discordant_l = v.vdiscordant_l
            vaf_l = v.vaf_l
            outfile.write(f'{vid}\t{vtype}\t{vlen}\t{split_l[0]}\t{cov_l[0]}\t{discordant_l[0]}\t{vaf_l[0]}\t{split_l[1]}\t{cov_l[1]}\t{discordant_l[1]}\t{vaf_l[1]}\n')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for processing VCF and calculating performance ")
    parser.add_argument("-i", "--bam", help = " input bam file ", required = True)
    parser.add_argument("-v", "--vcf", help = " input vcf file ", required = True)
    parser.add_argument("-o", "--output", help = " output table ", required = True)
    parser.add_argument("-d", "--depth", help = " input depth file ", required = False)
    #parser.add_argument("-O", "--outtable", help = " output table file ", required = False)
    parser.add_argument("-g", "--genome", help = " input reference genome fasta file ", required = True)
    parser.add_argument("-m", "--mode", help = " VAF calculation mode (read / kmer) ", default = 'kmer', required = False)
    parser.add_argument("-r", "--model", help = " pre-trained RF regression model to predict discordant reads", required = False)
    parser.add_argument("-k", "--kmer", help = " a pre-selected ROI kmer database pickle object ", required = False)
    parser.add_argument("-R", "--repeat_kmer", help = " a pre-selected ROI repeat kmer database pickle object ", required = False)
    #subparsers = parser.add_subparsers(title="mode")
    #parser_read = subparsers.add_parser('read')
    #parser_read.set_defaults(func=estimate_vaf)

    args = parser.parse_args()

    import sys, os, pysam, pickle, re, math
    from Bio.Seq import reverse_complement
    import numpy as np
    import dnsv_module as md
    import kmer_module as kmer_md
    from joblib import dump, load

    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_score
    import pandas as pd
    from sklearn.model_selection import GridSearchCV
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error

    if args.mode == 'read':
        estimate_vaf(args.vcf, args.bam, args.genome, args.depth, args.output, args.model, args.mode)
    elif args.mode == 'kmer':
        estimate_vaf(args.vcf, args.bam, args.genome, args.output, args.model, args.mode, None, args.kmer, args.repeat_kmer)
