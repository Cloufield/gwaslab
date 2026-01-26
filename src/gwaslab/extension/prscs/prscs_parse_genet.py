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
    # Read bim file and filter by chromosome (matching original PRScs logic)
    vld_dict = pd.read_csv(bim_file + '.bim', sep="\t", header=None, usecols=[0, 1, 4, 5])
    vld_dict.columns = ["CHR", "SNP", "A1", "A2"]
    # Filter by chromosome
    vld_dict = vld_dict.loc[vld_dict["CHR"] == chrom, ["SNP", "A1", "A2"]]
    return vld_dict


def parse_sumstats(ref_dict, vld_dict, sst_file, n_subj, log):
    """
    Parse sumstats using original PRScs logic with set-based matching.
    This matches the original PRScs behavior exactly.
    """
    import tempfile
    
    # Convert DataFrame to file-like format if needed
    # Original PRScs expects a file path, but GWASLab passes a DataFrame
    # We'll handle both cases
    if isinstance(sst_file, pd.DataFrame):
        # Save to temporary file to match original logic
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp_file:
            # Write header - original PRScs expects: SNP A1 A2 BETA P (or SE)
            # Map EA/NEA to A1/A2
            header_cols = ['SNP', 'A1', 'A2']
            if 'BETA' in sst_file.columns:
                header_cols.append('BETA')
            if 'P' in sst_file.columns:
                header_cols.append('P')
            elif 'SE' in sst_file.columns:
                header_cols.append('SE')
            tmp_file.write('\t'.join(header_cols) + '\n')
            # Write data - map EA/NEA to A1/A2
            for _, row in sst_file.iterrows():
                # Get A1 and A2 from EA/NEA or A1/A2 columns
                a1 = str(row.get('EA', row.get('A1', '')))
                a2 = str(row.get('NEA', row.get('A2', '')))
                line = [str(row.get('SNP', '')), a1, a2]
                if 'BETA' in sst_file.columns:
                    line.append(str(row['BETA']))
                if 'P' in sst_file.columns:
                    line.append(str(row['P']))
                elif 'SE' in sst_file.columns:
                    line.append(str(row['SE']))
                tmp_file.write('\t'.join(line) + '\n')
            sst_file_path = tmp_file.name
    else:
        # Already a file path
        sst_file_path = sst_file
    
    # Original PRScs logic
    ATGC = ['A', 'T', 'G', 'C']
    mapping = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    
    # Read sumstats file to extract SNP, A1, A2
    sst_dict_temp = {'SNP':[], 'A1':[], 'A2':[]}
    with open(sst_file_path) as ff:
        header = next(ff)
        for line in ff:
            ll = (line.strip()).split()
            if len(ll) >= 3 and ll[1] in ATGC and ll[2] in ATGC:
                sst_dict_temp['SNP'].append(ll[0])
                sst_dict_temp['A1'].append(ll[1])
                sst_dict_temp['A2'].append(ll[2])
    
    log.write('... %d SNPs read from sumstats ...' % len(sst_dict_temp['SNP']))
    
    # Build sets for matching (original PRScs logic)
    # Convert ref_dict and vld_dict to lists (they're DataFrames)
    ref_snp_list = ref_dict['SNP'].tolist()
    ref_a1_list = ref_dict['A1'].tolist()
    ref_a2_list = ref_dict['A2'].tolist()
    
    vld_snp_list = vld_dict['SNP'].tolist()
    vld_a1_list = vld_dict['A1'].tolist()
    vld_a2_list = vld_dict['A2'].tolist()
    
    vld_snp = set(zip(vld_snp_list, vld_a1_list, vld_a2_list))
    
    ref_snp = set(zip(ref_snp_list, ref_a1_list, ref_a2_list)) | \
              set(zip(ref_snp_list, ref_a2_list, ref_a1_list)) | \
              set(zip(ref_snp_list, [mapping[aa] for aa in ref_a1_list], [mapping[aa] for aa in ref_a2_list])) | \
              set(zip(ref_snp_list, [mapping[aa] for aa in ref_a2_list], [mapping[aa] for aa in ref_a1_list]))
    
    sst_snp = set(zip(sst_dict_temp['SNP'], sst_dict_temp['A1'], sst_dict_temp['A2'])) | \
              set(zip(sst_dict_temp['SNP'], sst_dict_temp['A2'], sst_dict_temp['A1'])) | \
              set(zip(sst_dict_temp['SNP'], [mapping[aa] for aa in sst_dict_temp['A1']], [mapping[aa] for aa in sst_dict_temp['A2']])) | \
              set(zip(sst_dict_temp['SNP'], [mapping[aa] for aa in sst_dict_temp['A2']], [mapping[aa] for aa in sst_dict_temp['A1']]))
    
    comm_snp = vld_snp & ref_snp & sst_snp
    log.write('... %d common SNPs in the reference, sumstats, and validation set ...' % len(comm_snp))
    
    # Compute beta_std using original PRScs logic
    n_sqrt = np.sqrt(n_subj)
    sst_eff = {}
    with open(sst_file_path) as ff:
        header = (next(ff).strip()).split()
        header = [col.upper() for col in header]
        for line in ff:
            ll = (line.strip()).split()
            if len(ll) < 4:
                continue
            snp = ll[0]; a1 = ll[1]; a2 = ll[2]
            if a1 not in ATGC or a2 not in ATGC:
                continue
            
            # Original PRScs logic: determine flip during conversion
            if (snp, a1, a2) in comm_snp or (snp, mapping[a1], mapping[a2]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = np.log(float(ll[3]))
                
                if 'SE' in header and len(ll) > 4:
                    se = float(ll[4])
                    beta_std = beta/se/n_sqrt
                elif 'P' in header and len(ll) > 4:
                    p = max(float(ll[4]), 1e-323)
                    beta_std = np.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                else:
                    continue
                
                sst_eff.update({snp: beta_std})
            
            elif (snp, a2, a1) in comm_snp or (snp, mapping[a2], mapping[a1]) in comm_snp:
                if 'BETA' in header:
                    beta = float(ll[3])
                elif 'OR' in header:
                    beta = np.log(float(ll[3]))
                
                if 'SE' in header and len(ll) > 4:
                    se = float(ll[4])
                    beta_std = -1*beta/se/n_sqrt
                elif 'P' in header and len(ll) > 4:
                    p = max(float(ll[4]), 1e-323)
                    beta_std = -1*np.sign(beta)*abs(norm.ppf(p/2.0))/n_sqrt
                else:
                    continue
                
                sst_eff.update({snp: beta_std})
    
    # Build sst_dict in reference panel order (original PRScs logic)
    sst_dict = {'CHR':[], 'SNP':[], 'BP':[], 'A1':[], 'A2':[], 'MAF':[], 'BETA':[], 'FLP':[]}
    
    # Get reference data as lists (ref_dict is a DataFrame)
    ref_chr_list = ref_dict['CHR'].tolist()
    ref_bp_list = ref_dict['BP'].tolist()
    ref_maf_list = ref_dict['MAF'].tolist()
    
    for (ii, snp) in enumerate(ref_snp_list):
        if snp in sst_eff:
            sst_dict['SNP'].append(snp)
            sst_dict['CHR'].append(ref_chr_list[ii])
            sst_dict['BP'].append(ref_bp_list[ii])
            sst_dict['BETA'].append(sst_eff[snp])
            
            a1 = ref_a1_list[ii]; a2 = ref_a2_list[ii]
            # Determine A1, A2, MAF, FLP based on which tuple matches (original logic)
            if (snp, a1, a2) in comm_snp:
                sst_dict['A1'].append(a1)
                sst_dict['A2'].append(a2)
                sst_dict['MAF'].append(ref_maf_list[ii])
                sst_dict['FLP'].append(1)
            elif (snp, a2, a1) in comm_snp:
                sst_dict['A1'].append(a2)
                sst_dict['A2'].append(a1)
                sst_dict['MAF'].append(1-ref_maf_list[ii])
                sst_dict['FLP'].append(-1)
            elif (snp, mapping[a1], mapping[a2]) in comm_snp:
                sst_dict['A1'].append(mapping[a1])
                sst_dict['A2'].append(mapping[a2])
                sst_dict['MAF'].append(ref_maf_list[ii])
                sst_dict['FLP'].append(1)
            elif (snp, mapping[a2], mapping[a1]) in comm_snp:
                sst_dict['A1'].append(mapping[a2])
                sst_dict['A2'].append(mapping[a1])
                sst_dict['MAF'].append(1-ref_maf_list[ii])
                sst_dict['FLP'].append(-1)
    
    # Clean up temporary file if we created one
    if isinstance(sst_file, pd.DataFrame) and 'sst_file_path' in locals():
        try:
            os.unlink(sst_file_path)
        except:
            pass
    
    log.write(" -Number of common SNPs:{}".format(len(sst_dict['SNP'])))
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


