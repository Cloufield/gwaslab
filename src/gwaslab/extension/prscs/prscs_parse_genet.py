#!/usr/bin/env python

"""
Parse the reference panel, summary statistics, and validation set.

"""

import os
import numpy as np
from scipy.stats import norm
from scipy import linalg
import h5py
import pandas as pd


def parse_ref(ref_file, chrom, log):
    # ref_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[]}
    ref_dict = pd.read_csv(ref_file ,sep="\t")
    ref_dict = ref_dict.loc[ref_dict["CHR"]==chrom,:]
    return ref_dict


def parse_bim(bim_file, chrom):

    vld_dict = pd.read_csv(bim_file + '.bim' ,sep="\t", usecols=[1,3,4])
    vld_dict.columns=["SNP","A1","A2"]
    return vld_dict


def parse_sumstats(ref_dict, vld_dict, sst_file, n_subj, log):
    
    n_sqrt = np.sqrt(n_subj)

    sst_file.dropna()
    sst_file["CHR"] = sst_file["CHR"].astype("int64")
    sst_file["BP"] = sst_file["BP"].astype("int64")
    sst_file["EA"] = sst_file["EA"].astype("string")
    sst_file["NEA"] = sst_file["NEA"].astype("string")

    sst_file = pd.merge(sst_file, ref_dict, on=["SNP","CHR","BP"],how="inner")
    
    is_flipped = ((sst_file["NEA"] == sst_file["A1"]) &(sst_file["EA"] == sst_file["A2"]))
    is_valid = ((sst_file["EA"] == sst_file["A1"]) & (sst_file["NEA"] == sst_file["A2"]))| is_flipped
    
    sst_file = sst_file.loc[is_valid,:]
    
    
    sst_file.loc[is_flipped, "MAF"] = 1 - sst_file.loc[is_flipped, "MAF"]
    sst_file["BETA"] = sst_file["BETA"] / sst_file["SE"] / n_sqrt

    sst_file.loc[~is_flipped, "BETA"] = 1 * sst_file.loc[~is_flipped, "BETA"]
    sst_file.loc[is_flipped, "BETA"] = -1 * sst_file.loc[is_flipped, "BETA"]

    sst_file["FLP"] = 1
    sst_file.loc[is_flipped, "FLP"] = -1
    log.write(" -Number of common SNPs:{}".format(len(sst_file)))
    sst_dict= sst_file[['CHR', 'SNP', 'BP', 'A1', 'A2', 'MAF', 'BETA', 'FLP']].to_dict("list")
    
    return sst_dict


def parse_ldblk(ldblk_dir, sst_dict, chrom, log):
    log.write('... parse reference LD on chromosome %d ...' % chrom)

    if '1kg' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_1kg_chr' + str(chrom) + '.hdf5'
    elif 'ukbb' in os.path.basename(ldblk_dir):
        chr_name = ldblk_dir + '/ldblk_ukbb_chr' + str(chrom) + '.hdf5'

    hdf_chr = h5py.File(chr_name, 'r')
    n_blk = len(hdf_chr)
    ld_blk = [np.array(hdf_chr['blk_'+str(blk)]['ldblk']) for blk in range(1,n_blk+1)]

    snp_blk = []
    for blk in range(1,n_blk+1):
        snp_blk.append([bb.decode("UTF-8") for bb in list(hdf_chr['blk_'+str(blk)]['snplist'])])

    blk_size = []
    mm = 0
    for blk in range(n_blk):
        idx = [ii for (ii, snp) in enumerate(snp_blk[blk]) if snp in sst_dict['SNP']]
        blk_size.append(len(idx))
        if idx != []:
            idx_blk = range(mm,mm+len(idx))
            flip = [sst_dict['FLP'][jj] for jj in idx_blk]
            ld_blk[blk] = ld_blk[blk][np.ix_(idx,idx)]*np.outer(flip,flip)

            _, s, v = linalg.svd(ld_blk[blk])
            h = np.dot(v.T, np.dot(np.diag(s), v))
            ld_blk[blk] = (ld_blk[blk]+h)/2            

            mm += len(idx)
        else:
            ld_blk[blk] = np.array([])

    return ld_blk, blk_size


