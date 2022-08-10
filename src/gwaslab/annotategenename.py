
def closest_gene(x,chrom="CHR",pos="POS",maxiter=20000,step=50):
        gene = data.gene_names_at_locus(contig=x[chrom], position=x[pos])
        if len(gene)==0:
            i=0
            while i<maxiter:
                distance = i*step
                gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance)
                gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance)
                if len(gene_u)>0 and len(gene_d)>0:
                    for j in range(0,step,2):
                        distance = (i-1)*step
                        gene_u = data.gene_names_at_locus(contig=x[chrom], position=x[pos]-distance-j)
                        gene_d = data.gene_names_at_locus(contig=x[chrom], position=x[pos]+distance+j)
                        if len(gene_u)>0:
                            return -1,",".join(gene_u)
                        else:
                            return 1,",".join(gene_d)
                elif len(gene_u)>0:
                    return -1,",".join(gene_u)
                elif len(gene_d)>0:
                    return 1,",".join(gene_d)
                else:
                    i+=1
            return 9,pd.NA
        else:
            return 0,",".join(gene)