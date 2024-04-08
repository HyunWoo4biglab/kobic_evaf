import os, sys, re, subprocess, pickle
import pysam
from Bio.Seq import reverse_complement
import numpy as np
import pandas as pd
#from sklearn.lienar_model import LinearRegression
import statsmodels.api as sm
import statsmodels.formula.api as smf
import kmer_module as kmer_md

class Variant:
    
    def __init__(self, contig, vid, mate_contig, mateid, vtype, vstart, vend, vinfo, vfilter, valt, vline):
        self.contig = contig
        self.vid = vid
        self.mate_contig = mate_contig
        self.mateid = mateid
        self.vtype = vtype
        self.vstart = vstart
        self.vend = vend
        self.vinfo = vinfo
        self.vfilter = vfilter
        self.valt = valt
        if self.vtype != 'BND':
            self.vlen = (self.vend - self.vstart)
        else:
            self.vlen = 50
        self.vaf_l = list()
        self.vcov_l = list()
        self.vsplit_l = list()
        self.vdiscordant_l = list()
        self.vline = vline
        self.call_flag_l = list()
        #self.vread_l = list()
        #self.vread_count_l = list()


    def set_contig(self, contig):
        self.contig = contig
    def set_vid(self, id):
        self.vid = id
    def set_mateid(self, mate):
        self.mateid = mate
    def set_vtype(self, svtype):
        self.vtype = svtype
    def set_vstart(self, position):
        self.vstart = position
    def set_vend(self, position):
        self.vend = position
    def set_vlen(self, length):
        self.vlen = length
    def set_vinfo(self, info):
        self.vinfo = info
    def set_filter(self, filt):
        self.vfilter = filt
    def set_alt(self, alt):
        self.valt = alt

    @staticmethod
    def description():
        print('A class for each variant')



    #@classmethod
    def set_variant_window(self, vwindow_l):
        self.vwindow_l = vwindow_l
    #def set_variant_window(self, start, end, mate_start, mate_end):
    #    self.vwindow = [(start, end), (mate_start, mate_end)]

    def set_variant_kmer(self, kmer_l, blacklistkmer_l, count_l, mapscore_l):
        self.kmer_l = kmer_l
        self.blacklistkmer_l = blacklistkmer_l
        self.count_l = count_l
        self.mapscore_l = count_l

    #@classmethod
    def set_variant_cov(self, vcov_l):
        self.vcov_l = vcov_l

    def set_variant_adjcov(self, adjcov):
        self.adjcov = adjcov

    def set_variant_supporting_read(self, vread_l):
        self.vread_l = vread_l

    def set_varant_supporting_readcount(self, vread_count_l):
        self.vread_count_l = vread_count_l

    def set_variant_af(self, vaf_l):
        self.vaf_l = vaf_l

    def set_vaf(self, vaf):
        self.vaf = vaf


    def get_variant_window(self):
        return self.vwindow

    def get_variant_cov(self):
        return self.cov

    def get_variant_kmer(self):
        return self.kmer_l

    def get_variant_kmer_count(self):
        return self.count_l

    def get_variant_kmer_mapscore(self):
        return self.mapscore_l

    def get_variant_blacklistkmer(self):
        return self.blacklistkmer_l

    def add_info(self, key, value):

        key = key.upper()
        if type(value) == list:
            for i in range(len(value)):
                self.vinfo += f';{key}{i+1}={value[i]}'
        else:
            self.vinfo += f';{key}={value}'

    def make_vcf_line(self):
        col = re.split('\s+', self.vline)
        col[7] = self.vinfo
        self.vline = '\t'.join(col)
        return self.vline


#    def set_variant_read_cov(self, depth):
#        windowcov = list()
#        for window in self.vwindow:
#            depth_l = list()
#            win_start, win_end = window
#            for i in range(win_start, win_end):
#                d = depth[self.contig][i]
#                depth_l.append(d)
#            average_cov = round(np.mean(depth_l), 4)
#            windowcov.append(average_cov)

#        self.vcov = widowcov


class VariantParser:
    #def __init__(self, vcf, ref, mode, win):
    def __init__(self, vcf):
        self.vcf = vcf
        #self.mode = mode # vaf calculation mode : readVAF | kmerVAF
        #self.win = win # variant window size
        #self.ref = ref
        self.variant_l = list()
        self.header_l = list()


    #def parse(self, ref, count_d, mappscore_d, repeatkmer_l, repeat_d, lowmap_d, k=31):
    def parse(self):
        """read the file and return """
        #refobj = pysam.FastaFile(ref)
        with open(self.vcf, 'r') as infile:
            for line in infile:
                if line.startswith('#'):
                    self.header_l.append(line)
                    continue
                col = re.split('\s+', line)
                contig = col[0]; start = int(col[1])-1; vid = col[2]; ref = col[3]; alt = col[4]; filt = col[6]; info = col[7]

                mateid = self.get_mateid(vid, info)
                vtype = self.get_vtype(info)
                mate_contig, end = self.get_vend(contig, info)

                variant = Variant(contig, vid, mate_contig, mateid, vtype, start, end, info, filt, alt, line.strip())

                vaf = self.get_vaf(info)
                if vaf:
                    variant.set_vaf(vaf)
                self.variant_l.append(variant)

#                variant_window = get_variant_window(vtype, start, end)
#                variant.set_variant_window(variant_window)
#
#                if self.mode == 'read':
#                    variant_cov = self.get_variant_read_cov(variant_window, count_d) #count_d is depth_d for read mode
#                elif self.mode == 'kmer':
#                    kmer_l, blacklistkmer_l, count_l, mappscore_l = self.get_varaint_window_kmers(self, window_l, refobj, count_d, mappscore_d, repeatkmer_l, repeat_d, lowmap_d, k)
#                    
#                    #kmer_count_l, kmer_adjcount_l = get_kmer_counts(kmer_l, count_d) # count_d is kmc_dump_d for kmer mode
#                    
#                    variant_cov = self.estimate_cov(count_l, mappscore_l)
#                    variant_adjcov = np.mean(kmer_ajdcount_l)
#                    variant.set_variant_adjcov(variant_adjcov)
#                    variant.set_variant_kmer(kmer_l, blacklistkmer_l, count_l, mappscore_l)
#
#                variant.set_variant_cov(variant_cov)
#                self.variant_l.append(variant)


    #@classmethod
    def get_variant_list(self):
        return self.variant_l

    def get_header(self):
        return self.header_l

    def get_vend(self, contig, info):
        match = re.search('=([a-zA-Z]+[1-9][0-9]?):(\d+)\([+-]\)-([a-zA-Z]+[1-9][0-9]?):(\d+)\([+-]\)', info)
        if match:
            mate_contig = match.group(3)
            #start = int(match.group(2)) -1
            end = int(match.group(4)) -1
        else:
            mate_contig = contig
            match = re.search('(;)?END=(\d+)(;)?', info)
            end = int(match.group(2)) -1
        return mate_contig, end

    def get_mateid(self, vid, info):
        match = re.search('(;)?MATEID=(\S+)(;)?', info)
        if match:
            mateid = match.group(2)
        else:
            mateid = f'{vid}_mate'
        return mateid

    def get_vtype(self, info):
        vtype = re.search('(;)?SVTYPE=([A-Z]+)(;)?', info).group(2)
        return vtype

    def get_vaf(self, info):
        vaf_match = re.search(';VAF=([+-]?([0-9]*[.])?[0-9]+)', info)
        if vaf_match:
            vaf = vaf_match.group(1)
        else:
            vaf = None
        return vaf


    def get_varaint_window(self, vtype, contig, start, mate_contig, end):#
        #variant.set_variant_window(150)
        """
        | get variant's breakpoint window
        : variant : Class object Variant
        : vtype : variant type
        : start : varaint start position
        : end : variant end position
        : win : variant window 1-side size
        """
        windows = list()
        ##-- BND1 START POS
        if start >= self.win:
            if vtype == 'DUP':
                bp1_start = start
            else:
                bp1_start = start - self.win
        else:
            bp1_start = 0

        ##--BND1 END POS
        if vtype == 'DEL':
            bp1_end = start
        else:
            bp1_end = start + self.win +1
        windows.append((contig, mate_contig, bp1_start, bp1_end))

        ##-- BND2 START POS
        if end >= self.win:
            if vtype == 'DEL':
                bp2_start = end +1
            else:
                bp2_start = end - self.win
        else:
            bp2_start = 0

         ##-- BND2 END POS
        if vtype == 'DUP':
            bp2_end = end +1
        else:
            bp2_end = end + self.win +1
        windows.append((mate_contig, contig, bp2_start, bp2_end))

        return windows

    def get_variant_read_coverage(self, window_l, depth):#
        """
        | GET VARIANT WINDOW READ COVERAGE FROM SAMTOOLS DEPTH DICTIONARY
        : window_l : a list of variant window
        : depth : a dictionary of ROI read depth calculated by samtools depth
        """
        windowcov = list()
        for w in window_l:
            depth_l = list()
            contig, win_start, win_end = w
            for i in range(win_start, win_end):
                d = depth[self.contig][i]
                depth_l.append(d)
            average_cov = round(np.mean(depth_l), 4)
            windowcov.append(average_cov)
 
        return windowcov

    def get_varaint_window_kmers(self, window_l, refobj, count_d, mapscore_d, repeatkmer_l, repeat_d, lowmap_d, k=31):#
        """
        | GET REFERENCE KMERS, KMER COUNTS, MAPPABILITY SCORE AND BLACKLIST KMERS FROM A GIVEN BP WINDOW
        : contig : contig of the start(BP1)
        : mate_contig : contig of the end(BP2)
        : window_l : a list of window returned from self.get_varaint_window()
        : refobj : reference genome sequence object from pysam.FastaFile()
        : k : kmer size; default=31mer
        """        
        kmer_l = list() # in the form of [[bp1_kmers, ...], [bp2_kmers, ...]]
        count_l = list()
        mapscore_l = list()
        blacklistkmer_l = list()
        for w in window_l:
            contig, start, end = w
            kmers = list(); counts = list(); map_scores = list(); blacklistkmers = list()
            seq = get_ref_seq(refobj, contig, start, end)
            seq.upper()
            for i in range(len(seq) - k+1):
                kmer = seq[i:i+k]
                if 'N' in kmer:
                    continue
                rc_kmer = reverse_complement(kmer)
                kmer = kmer_md.lexicographical_comparison(kmer, rc_kmer)
                kmers.append(kmer)
                if (kmer in repeatkmer_l) or (genomic_overlap(repeat_d, contig, start +i, start +i +k)) or (genomic_overlap(lowmap_d, contig, start +i, start +i +k)):
                    score = mapscore_d[contig][start +i]
                else:
                    score = 1.0

                bit_kmer = kmer_md.convert_to_bit(kmer)
                try:
                    count = count_d[bit_kmer]
                    #adjcount = count * map_score
                except KeyError:
                    count = 0.0

                if (kmer in repeatkmer_l) or (count >= 200):
                    blacklistkmers.append(kmer)

                counts.append(count)
                map_scores.append(score)

            kmer_l.append(kmers)
            blacklistkmer_l.append(blacklistkmers)
            count_l.append(counts)
            mapscore_l.append(map_scores)

        return kmer_l, blacklistkmer_l, count_l, mapscore_l



    def estimate_coverage(self, count_l, mapscore_l):
        """
        | ESTIMATE VARIANT WINDOW COVERAGE USING REFERENCE KMER AND ITS CORRESPONDING MAPPABILITY SCORE
        : count_l : kmer count list(containing counts of kmer not in the kmerdb ;value set as 0) from self.get_variant_window_kmer()
        : mappscore_l : mappability score list from self.get_variant_window_kmer()
        """

        nonzero_index_l = [i for i, e in enumerate(count_l) if e != 0]
        adjcount_l = [count_l[i]*mapscore_l[i] for i in nonzero_index_l]
        estimated_cov = np.mean(adjcount_l)
        #adjcount_l = list(np.multiply(count_l, mapscore_l))

        return estimated_cov

def get_variant_window(v, win):
    #variant.set_variant_window(150)
    """
    | GET VARIANT'S BREAKPOINT WINDOW
    : v : vairant class object
    : win : variant window 1-sided size
    """
    window_l = list()#in the form of [(bp1_window_info), (bp2_window_info)]

    ##-- BND1 START POS
    if v.vstart >= win:
        if v.vtype == 'DUP':
            bp1_start = v.vstart

        elif v.vtype == 'DEL':
            bp1_start = v.vstart - win*2
        else:
            bp1_start = v.vstart - win
    else:
        bp1_start = 0

    ##--BND1 END POS
    if v.vtype == 'DEL':
        bp1_end = v.vstart
    elif v.vtype == 'DUP':
        bp1_end = v.vstart + win*2 +1
    else:
        bp1_end = v.vstart + win +1
    window_l.append((v.contig, v.mate_contig, v.vstart, v.vend, bp1_start, bp1_end))

    ##-- BND2 START POS
    if v.vend >= win:
        if v.vtype == 'DEL':
            bp2_start = v.vend +1
        elif v.vtype == 'DUP':
            bp2_start = v.vend - win*2
        else:
            bp2_start = v.vend - win
    else:
        bp2_start = 0

     ##-- BND2 END POS
    if v.vtype == 'DUP':
        bp2_end = v.vend +1
    elif v.vtype == 'DEL':
        bp2_end = v.vend + win*2 +1
    else:
        bp2_end = v.vend + win +1
    window_l.append((v.mate_contig, v.contig, v.vend, v.vstart, bp2_start, bp2_end))

    return window_l


def get_variant_read_coverage(w, depth_d):
    """
    | GET VARIANT WINDOW READ COVERAGE FROM SAMTOOLS DEPTH DICTIONARY
    : window_l : a list of variant window
    : depth : a dictionary of ROI read depth calculated by samtools depth
    """
    contig, mate_contig, wstart, wend = w
    for i in range(wstart, wend):
        d = depth_d[contig][i]
        depth_l.append(d)
    average_cov = round(np.mean(depth_l), 4)
    #cov_l.append(average_cov)

    return average_cov

def estimate_coverage(count_l, mapscore_l):
    """
    | ESTIMATE VARIANT WINDOW COVERAGE USING REFERENCE KMER AND ITS CORRESPONDING MAPPABILITY SCORE
    : count_l : kmer count list(containing counts of kmer not in the kmerdb ;value set as 0) from get_variant_window_kmer()
    : mappscore_l : mappability score list from get_variant_window_kmer()
    """
    print('count_list', count_l)
    nonzero_index_l = [i for i, e in enumerate(count_l) if e != 0]
    adjcount_l = [count_l[i]*mapscore_l[i] for i in nonzero_index_l]
    print('adjusted_count_list', adjcount_l)
    if len(adjcount_l) != 0:
        estimated_cov = round(np.mean(adjcount_l), 4)
    else:
        estimated_cov = 0.0
    #adjcount_l = list(np.multiply(count_l, mappscore_l))

    return estimated_cov


def get_varaint_window_kmers(w, refobj, count_d, mapscore_d, repeatkmer_l, repeat_d, lowmap_d, k=31):
    """
    | GET REFERENCE KMERS, KMER COUNTS, MAPPABILITY SCORE AND BLACKLIST KMERS FROM A GIVEN BP WINDOW
    : contig : contig of the start(BP1)
    : mate_contig : contig of the end(BP2)
    : window_l : a list of window returned from v.get_varaint_window()
    : refobj : reference genome sequence object from pysam.FastaFile()
    : k : kmer size; default=31mer
    """
    kmer_l = list() # in the form of [[bp_kmers, ...]]
    count_l = list()
    mapscore_l = list()
    blacklistkmer_l = list()

    contig, mate_contig, start, end, wstart, wend = w
    #kmers = list(); counts = list(); mapp_scores = list(); blacklistkmers = list()

    seq = get_ref_seq(refobj, contig, wstart, wend)
    seq = seq.upper()

    count_check_nz = 0; count_check_z = 0
    for i in range(len(seq) - k+1):
        kmer = seq[i:i+k]
        if 'N' in kmer:
            continue
        rc_kmer = reverse_complement(kmer)
        kmer = kmer_md.lexicographical_comparison(kmer, rc_kmer)
        kmer_l.append(kmer)
        if (kmer in repeatkmer_l) \
        or (check_overlap(repeat_d, contig, wstart +i, wstart +i +k)) \
        or (check_overlap(lowmap_d, contig, wstart +i, wstart +i +k)):
#        if True:
            score = mapscore_d[contig][wstart +i]
        else:
            score = 1.0

        bit_kmer = kmer_md.convert_to_bit(kmer)
        try:
            count_check_nz += 1
            count = count_d[bit_kmer]
            #adjcount = count * mapp_score
        except KeyError:
            count_check_z += 1
            count = 0.0


        if (kmer in repeatkmer_l) or (count >= 200):
            blacklistkmer_l.append(kmer)

        count_l.append(count)
        mapscore_l.append(score)

    print('window kmer check', count_check_nz, count_check_z)

    return kmer_l, blacklistkmer_l, count_l, mapscore_l


 
def variant_match(var, true_var, offset=10):
    """
    | A function to check whether a given variant overlapps with true variant within a given offset
    : var : class object of input candidate variant
    " true_var : class object of true benchmark varaint
    """
    
    if var.contig != true_var.contig:
        return False
    if var.vtype != true_var.vtype:
        return False
    if (var.vstart < true_var.vstart - offset) or (var.vstart > true_var.vstart + offset):
        return False
    if (var.vend < true_var.vend - offset) or (var.vend > true_var.vend + offset):
        return False
    return True
    

def fit_discordant_read_pair_lm(dr, sr, cov):
    """
    | fit a linear model to predict discordant read pair counts 
    using split read, estimated coverage and their interaction term as explanatory variables
    : dr : discordant read count
    : sr : split read count
    : cov : estimated/average read depth of coverage
    """
    
    #lm = smf.ols(formula="dr ~ sr + cov + sr:cov", data=df).fit()
    lm = smf.glm(formula="dr ~ sr + cov + sr:cov", data=df).fit()
    lm = smf.rlm(formula="dr ~ sr + cov + sr:cov", data=df).fit()
    lm = LinearRegression().fit(x, y)
    #r_sq = lm.score(x, y)

    print(lm.summary())
    print(f'coefficient of determination : {lm.rsquared}')
    print(f'adjusted coefficient of determination : {lm.rsquared_adj}')
    #print(f'coefficients: {model.coef_}')
    #print(f'intercept : {model.intercept_}')
    return lm


def predict_discordant(regressor, df):
    """
    | A FUNCITON TO PREDICT DISCORDANT READ COUNTS WITH A GIVEN REGRESSION MODEL
    : regressor : pre-trained linear / random forest regression model
    : x : input pandas dataframe
    """

    print('INPUT DATAFRAME COLUMN NAMES', list(df.columns))
    discordant_d = dict()
    labels = df.iloc[:, 0]
    #pairs = df.iloc[:, 1]
    #X = df.iloc[:, 2:].values #features 2: SV_LENGTH, 3:SPLIT_READ, 4:COV, 5:VAF_PRIOR
    X = df.iloc[:, 2:] #features 2: SV_LENGTH, 3:SPLIT_READ, 4:COV, 5:VAF_PRIOR
    predicted = regressor.predict(X)

    for i, d in enumerate(predicted):
        vid = labels[i]
        #pair = pairs[i]
        if not vid in discordant_d:
            discordant_d[vid] = list()
        discordant_d[vid].append(int(d))

    return discordant_d

def train_regressor(n_tree, n_depth, seed, df):
    """
    | A FUNCTION TO TRAIN A RANDOM FOREST REGRESSION MODEL
    : df : input dataframe having {sv_id, pair_id, cov, split_read, sv_length, vaf_prior, discordant_readcount} as columns
    """

    X = df.iloc[:, 2:6].values
    y = df.iloc[:, 6].values
    regressor = RandomForestRegressor(n_estimators=n_tree, criterion ='squared_error', max_depth =n_depth, random_state=seed, oob_score=True)
    regressor.fit(X, y)
    print('OOB score : ', regressor.oob_score_)
    print('R_squred :', regressor.score(X, y))

    return regressor


def execute_subprocess(cmd):
    print(f'process : {cmd}')
    p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE)
    p.wait()
    print(f'process : {cmd} DONE!')



def get_ref_seq(refobj, contig, start, end):
    """
    | get reference genome seqeunce with a given window
    : refobj : pysam object with input reference genome
    : contig : target contig
    : start, end : starting & ending position of the target window
    """

    ref_seq = refobj.fetch(contig, start, end)

    return ref_seq

def check_overlap(region_dict, contig, start, end):
    """
    | check wether a region(SV window) overlaps with a given window
    : region_dict : a dictionary containing input regions
    : contig : target contig
    : start, end : starting & ending position of the window
    """
    overlap = False
    for region in region_dict[contig]:
        if start <= region[-1] and end >= region[0]:
            overlap = True
            break

    return overlap

def extract_vreads(bam, w, v, discordant_bool):#, cigar):
    """
    | EXTRACT BREAKPOINT-SPANNING SPLIT READS / SV-SPANNING DISCORDANT READ PAIRS FROM A GIVEN VARIANT WINDOW
    : bam : input bam file pysam object 
    : w : breakpoint window (contig, mate_contig, wstart, wend)
    : v : variant object containing the variant information
    : discordant_bool : boolean (True or False) for whether to extract discordant reads from the given window
    """
    #sv_id, sv_type, contig, mate_contig, breakpoint, matepoint, win_start, win_end, pair, sv_length = sv_info
    contig, mate_contig, breakpoint, matepoint, wstart, wend = w

    split_reads = set(); split_counts = dict()
    if discordant_bool:
        discordant_reads = dict()
        discordant_reads['orient'] = set(); discordant_reads['size'] = set(); discordant_reads['total'] = set()
        dis_offset = 250

    for read in bam.fetch(contig, wstart, wend):
        if (not read.is_unmapped) and (read.mapping_quality >= 20) and (read.mapping_quality != 255):
            ##-- GET SPLIT READS
            if read.has_tag('SA'):
                lclip_pos = read.reference_start
                rclip_pos = read.reference_end -1

                if ((lclip_pos >= breakpoint -1) and (lclip_pos <= breakpoint +1)) \
                    or ((rclip_pos >= breakpoint -1) and (rclip_pos <= breakpoint +1)):
                #FOR NOT ALLOWING 1BP OFFSET
                #if lclip_pos == breakpoint or rcplip_pos == breakpoint:

                    split_aligns = read.get_tag('SA').split(';')[:-1]
                    if len(split_aligns) <= 2:
                        #splits = list()
                        for sa in split_aligns:
                            sa_contig, sa_pos, sa_strand, sa_cigar, sa_mq, sa_nm = sa.split(',')
                            sa_pos = int(sa_pos) -1
                            if 'S' in sa_cigar:
                                clip_cigar = 'S'
                            elif 'H' in sa_cigar:
                                clip_cigar = 'H'
                            else:
                                clip_cigar = 'S'
                            matches = re.findall(r'(\d+)([A-Z])', sa_cigar)
                            n_match = [int(match[0]) for match in matches if match[1] != 'S' and match[1] != 'H' and match[1] != 'I']
                            len_match = sum(n_match)

                            if v.vtype in ['DEL', 'DUP']:
                                if ((not read.is_reverse) and (sa_strand == '+')) or\
                                (read.is_reverse and (sa_strand == '-')):
                                    if max(sa_cigar.rfind(i) for i in 'SH') > sa_cigar.find('M'):
                                        if sa_cigar.find(clip_cigar) < sa_cigar.find('M'):
                                            #splits.append(sa_pos)
                                            pass
                                        else:
                                            sa_pos += len_match -1
                                            #splits.append(sa_pos + len_match -1)
                                    else:
                                        pass
                                        #splits.append(sa_pos)
                                else:
                                    sa_pos = -1#make is as NA
                                    #splits.apend(sa_pos)
                            elif v.vtype == 'INV':
                                if ((not read.is_reverse) and (sa_strand == '-')) or \
                                ((read.is_reverse) and (sa_strand == '+')):
                                    if max(sa_cigar.rfind(i) for i in 'SH') > sa_cigar.find('M'):
                                        if sa_cigar.find(clip_cigar) < sa_cigar.find('M'):
                                            pass
                                            #splits.append(sa_pos)
                                        else:
                                            sa_pos += len_match -1
                                            #splits.append(sa_pos + len_match -1)
                                    else:
                                        pass
                                        #splits.append(sa_pos)
                                else:
                                    sa_pos = -1
                            elif v.vtype == 'BND':
                                if max(sa_cigar.rfind(i) for i in 'SH') > sa_cigar.find('M'):
                                    if sa_cigar.find(clip_cigar) < sa_cigar.find('M'):
                                        if sa_pos < matepoint:
                                            sa_pos += len_match -1
                                            #splits.append(sa_pos + len_match -1)
                                        else:
                                            pass
                                            #splits.append(sa_pos)
                                    else:
                                        sa_pos += len_match -1
                                        #splits.append(sa_pos + len_match -1)
                                else:
                                    pass
                                    #splits.append(sa_pos)

                    #for pos in splits:
                    if not sa_pos in split_counts:
                        split_counts[sa_pos] = set()
                    if sa_pos == matepoint:
                    #if ((pos >= matepoint -1) and (pos <= matepoint +1)):
                        split_counts[sa_pos].add(read.query_name)

            #---GET DISCORDANT READ PAIRS        
            if discordant_bool:

                discordant = False; mate_match = False
                if read.next_reference_name == mate_contig:
                    if (read.flag & 0x1) and (not read.flag & 0x4) and (not read.flag & 0x8) and (not read.flag & 0x100) and (not read.flag & 0x400) and \
                    ((not read.flag & 0x2) or \
                    ((read.flag & 0x10 and read.flag & 0x20) or ((not read.flag & 0x10) and (not read.flag & 0x20))) or \
                    (read.next_reference_start < read.reference_start)):
                    #if ((read.is_paired) and (not read.mate_is_unmapped) and (not read.is_proper_pair)) or \
                    #((read.is_paired) and (not read.mate_is_unmapped) and (read.is_reverse == read.mate_is_reverse)):
                        if sv_type == 'INV':
                            #if read.is_reverse == read.mate_is_reverse:
                            if (read.flag & 0x10 and read.flag & 0x20) or \
                            ((not read.flag & 0x10) and (not read.flag & 0x20)):
                                discordant_bool = True
                                discordant_type = 'orient'

                        elif sv_type == 'DEL':
                            if (not read.flag & 0x2) or \
                            (abs(read.template_length) > 850):
                                discordant_bool = True
                                discordant_type = 'size'

                        elif sv_type == 'DUP':
                            rl_pair = False

                            try:
                                mate_read = bamobj.mate(read)
                            except Exception as e:
                            #except ValueError:
                                print('an exception occured, ', e)
                                if read.flag & 0x40: # first in pair
                                    if (read.next_reference_start < read.reference_start):
                                        rl_pair = True
                                elif read.flag & 0x80: # second in pair
                                    if (read.reference_start < read.next_reference_start):
                                        rl_pair = True
                            else:
                                if read.flag & 0x40: # first in pair
                                    if (mate_read.reference_start < read.reference_start):
                                        rl_pair = True
                                elif read.flag & 0x80: # second in pair
                                    if (read.reference_start < mate_read.reference_start):
                                        rl_pair = True
                                else:
                                    print(read.query_name, read.flag)
                            finally:
                                if rl_pair or \
                                (not read.flag & 0x2) or \
                                (abs(read.template_length) > 850) or (abs(read.template_length) < 250):
                                    discordant_bool = True
                                    discordant_type = 'size'
                        try:
                            mate_read = bamobj.mate(read)
                        except Exception as e:
                            if (read.next_reference_start >= breakpoint_mate - breakend_extend_size) and \
                            (read.next_reference_end <= breakpoint_mate + breakend_extend_size):
                                bp_match = True
                        else:
                            if (mate_read.reference_start >= breakpoint_mate - breakend_extend_size) and \
                            (mate_read.reference_end <= breakpoint_mate + breakend_extend_size):
                                bp_match = True
                        finally:
                            if discordant and bp_match:
                                if dis_type == 'orient':
                                    discordants['orient'].add(read.query_name)
                                elif disc_type == 'size':
                                    discordants['size'].add(read.query_name)
                                    read_type_dict[sv_identifier]['discordant_s'].add(read.query_name)
                                    read_dict[sv_id]['discordant_s'].add(read.query_name)
                                    read_dict[mate_id]['discordant_s'].add(read.query_name)


    split_counts_filt = [split_counts[sa] for sa in split_counts if len(split_counts[sa]) >= 2]

    for spr in split_counts_filt:
        split_reads.update(spr)

    if discordant_bool:
        return split_reads, discordant_reads
    else:
        return split_reads


def get_depth(inputfile):
    """
    | get depth for each postion in a given variant window from samtools depth result
    : inputfile : samtools depth result file
    """
    depth_d = dict()
    with open(inputfile, 'r') as infile:
        for line in infile:
            col = re.split('\s+', line.strip())
            contig, pos, depth = col
            if not contig in depth_d:
                depth_d[contig] = dict()
            depth_d[contig][int(pos) -1] = int(depth)

    return depth_d

def get_mappability_score(mapscore_file):
    contig_list = set()
    map_score_dict = dict()
    with open(mapscore_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                contig = line.strip()[1:]
                contig_list.add(contig)
            else:
#                if not contig in ['chrM', 'chrX', 'chrY'] and not contig.startswith('GL'):
                if re.match('(chr[1-9][0-9]?|chrX|chrY)\\b', contig):
                    map_scores = [float(s) for s in line.strip().split()]
                    map_score_dict[contig] = map_scores
    print(contig_list, 'contigs in you mappability score file')
    return map_score_dict

def get_difficult_genomic_context(genomic_context_bed):
    """Load GIAB genome stratification - difficult region bed annotation"""
    diff_region_dict = dict()
    with open(genomic_context_bed, 'r') as infile:
        regions = [line.strip().split('\t') for line in infile if not line.startswith('#')]
        for region in regions:
            contig, start, end = region
            if not contig.startswith('chr'):
                contig = 'chr' + contig
            if not contig in diff_region_dict:
                diff_region_dict[contig] = list()
            diff_region_dict[contig].append((int(start), int(end)))

    return diff_region_dict

def get_repeat_region(repeat_annotation_file):
    repeat_dict = dict()
    with open(repeat_annotation_file, 'r') as infile:
        lines = [line for line in infile if line.startswith(' ')]
        for line in lines:
            col = line.split(' ')
            col = list(filter(lambda x: x != '', col))
            contig = col[4]; start = col[5]; end = col[6];
            if not contig in repeat_dict:
                repeat_dict[contig] = list()
            repeat_dict[contig].append((int(start) -1, int(end) -1))

    return repeat_dict

def get_bed(genomic_bed):
    """Load genomic bed file - region bed annotation"""
    out_d = dict()
    with open(genomic_bed, 'r') as infile:
        regions = [line.strip().split('\t') for line in infile if not line.startswith('#')]
        for region in regions:
            contig, start, end = region
            if not contig.startswith('chr'):
                contig = 'chr' + contig
            if not contig in out_d:
                out_d[contig] = list()
            out_d[contig].append((int(start), int(end)))

    return out_d

