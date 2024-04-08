# A SCRIPT TO FILTER VCF VARIANTS HAVING THEIR BREAKPOINT WINDOW OVERLAPPING WITH LOW-COMPLEXITY REGIONS
# HYUN WOO KIM
# 03-25-2024

def filter_variant(vcf, bed, output, minlen=50):
    """
    | A FUNCTION TO FILTER VARAINTS OVERLAPPING WITH LCR
    : vcf : input VCF file
    : bed : LCR bed file
    : output : output VCF file
    """

    lcr_d = md.get_bed(bed)
    parser = md.VariantParser(vcf)
    parser.parse()
    v_l = parser.get_variant_list()
    header_l = parser.get_header()

    with open(output, 'w') as outfile:
        for h in header_l:
            outfile.write(h)
        for v in v_l:
            #-- FILTER UNPALCED/UNLOCALIZED/MITOCHONDRIAL ETC. CONTIGS and SVLEN < 50
            if re.match('(chr[1-9][0-9]?|chrX|chrY)\\b', v.contig) and re.match('(chr[1-9][0-9]?|chrX|chrY)\\b', v.mate_contig) and (v.vlen >= minlen) and (v.vfilter == 'PASS'):

                #-- FILTER INTRA-CHROMOSOMAL BND
                if v.vtype == 'BND':
                    if v.contig == v.mate_contig:
                        continue
                #v.set_variant_window(vwindow_l)
                lcr_flag = False
                vwindow_l = md.get_variant_window(v, win=30)
                for w in vwindow_l:
                    contig, mate_contig, start, end, wstart, wend = w
                    if md.check_overlap(lcr_d, contig, wstart, wend):
                        lcr_flag = True
#                if not lcr_flag:
                if True:
                    line = v.make_vcf_line()
                    outfile.write(f"{line}\n")
            
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = " A script for comparing an input VCF with low complexity region bed file ")
    parser.add_argument("-b", "--bed", help = " benchmark vcf to compare ", required = True)
    parser.add_argument("-v", "--vcf", help = " input vcf ", required = True)
    parser.add_argument("-o", "--output", help = " output vcf ", required = True)
    args = parser.parse_args()

    import re
    import dnsv_module as md

    filter_variant(args.vcf, args.bed, args.output) 
