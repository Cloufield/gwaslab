
## extract lead variants

## for each lead variant 
    ## extract snp list from sumstats

    ## check available snps with reference file

    ## matching alleles

    ## extract available SNPs and calculate ld matrix using PLINK

plink \
  --bfile ${plinkFile} \
  --keep-allele-order \
  --r square \
  --extract sig_locus.snplist \
  --out sig_locus_mt

plink \
  --bfile ${plinkFile} \
  --keep-allele-order \
  --r2 square \
  --extract sig_locus.snplist \
  --out sig_locus_mt_r2
