import os
import subprocess


def clump(insumstats, vcf, bfile=None, chrom=None, clump_p1=5e-8, clump_p2=5e-8, clump_r2=0.2, clump_kb=250):
    ## process reference
    sumstats = insumstats.loc[insumstats["P"]<clump_p1,:].copy()
    
    if bfile is None:
        if vcf is not None:
            print("Processing VCF...")
            for i in sumstats["CHR"].unique():
                bfile_gwaslab = vcf.replace(".vcf.gz","") + ".{}".format(i)
                print("Processing VCF for CHR {}...".format(i))
                if not os.path.exists(bfile_gwaslab+".bed"):
                    script_vcf_to_bfile = """
                    plink \
                        --vcf {} \
                        --chr {} \
                        --make-bed \
                        --keep-allele-order \
                        --list-duplicate-vars suppress-first \
                        --out {}
                    """.format(vcf, i, bfile_gwaslab)

                    subprocess.run(script_vcf_to_bfile,capture_output=True, shell=True)

                    
                    exclude = bfile_gwaslab + ".dupvar"
                    script_vcf_to_bfile = """
                    plink \
                        --bfile {} \
                        --exclude {}\
                        --make-bed \
                        --keep-allele-order \
                        --out {}
                    """.format(bfile_gwaslab, exclude, bfile_gwaslab)

                    subprocess.run(script_vcf_to_bfile,capture_output=True,shell=True)

                    
            
    ## process sumstats by CHR
    for i in sumstats["CHR"].unique():
        print("Processing sumstats for CHR {}...".format(i))
        sumstats.loc[sumstats["CHR"]==i,["SNPID","P"]].to_csv("_gwaslab_tmp.{}.SNPIDP".format(i),index=None,sep="\t")
    
    
    results = pd.DataFrame()
    ## clump using plink
    for i in sumstats["CHR"].unique():
        if bfile is None:
            bfile_gwaslab = vcf.replace(".vcf.gz","") + ".{}".format(i)

        chrom = i
        clump = "_gwaslab_tmp.{}.SNPIDP".format(chrom)
        clump_p1 = 0.00001 
        clump_p2 = 0.001
        clump_r2 = 0.1
        clump_kb = 250
        out="_gwaslab_tmp.{}".format(chrom)
        
        if bfile is None:
            bfile_to_use = bfile_gwaslab
        else:
            bfile_to_use = bfile
        print("Clumping for CHR {}...".format(i))
        script = """
        plink \
            --bfile {}\
            --chr {} \
            --clump {} \
            --clump-field P \
            --clump-snp-field SNPID \
            --clump-p1 {} \
            --clump-p2 {} \
            --clump-r2 {} \
            --clump-kb {} \
            --out {}
        """.format(bfile_to_use, chrom, clump, clump_p1, clump_p2, clump_r2, clump_kb, out)

        os.system(script)
    
        clumped = pd.read_csv("_gwaslab_tmp.{}.clumped".format(chrom),usecols=[2,0,3,4],sep="\s+")
        results = pd.concat([results,clumped],ignore_index=True)
    return results