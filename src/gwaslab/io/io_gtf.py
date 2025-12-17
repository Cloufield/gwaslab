import pandas as pd
from os import path
from gwaslab.g_Log import Log
from gtfparse import read_gtf
from gwaslab.bd.bd_download import check_and_download

def read_gtf_file(gtf_path):
    return read_gtf(
        gtf_path,
        usecols=[
            "seqname",
            "start",
            "end",
            "strand",
            "feature",
            "gene_biotype",
            "gene_id",
            "gene_name",
        ],
    )

def get_gtf(chrom, build="19", source="ensembl"):
    gtf = None
    if source == "ensembl":
        if build == "19":
            data_path = check_and_download("ensembl_hg19_gtf")
            gtf = read_gtf(
                data_path,
                usecols=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    "feature",
                    "gene_biotype",
                    "gene_id",
                    "gene_name",
                ],
            )
            gtf = gtf.loc[gtf["seqname"] == chrom, :]
        if build == "38":
            data_path = check_and_download("ensembl_hg38_gtf")
            gtf = read_gtf(
                data_path,
                usecols=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    "feature",
                    "gene_biotype",
                    "gene_id",
                    "gene_name",
                ],
            )
            gtf = gtf.loc[gtf["seqname"] == chrom, :]
    if source == "refseq":
        from gwaslab.bd.bd_common_data import get_chr_to_NC
        chrom_NC = get_chr_to_NC(build=build)[str(chrom)]
        if build == "19":
            data_path = check_and_download("refseq_hg19_gtf")
            gtf = read_gtf(
                data_path,
                usecols=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    "feature",
                    "gene_biotype",
                    "gene_id",
                    "gene_name",
                ],
            )
            gtf = gtf.loc[gtf["seqname"] == chrom_NC, :]
        if build == "38":
            data_path = check_and_download("refseq_hg38_gtf")
            gtf = read_gtf(
                data_path,
                usecols=[
                    "seqname",
                    "start",
                    "end",
                    "strand",
                    "feature",
                    "gene_biotype",
                    "gene_id",
                    "gene_name",
                ],
            )
            gtf = gtf.loc[gtf["seqname"] == chrom_NC, :]
        gtf["seqname"] = str(chrom)
    if gtf is None:
        gtf = pd.DataFrame(
            columns=[
                "seqname",
                "start",
                "end",
                "strand",
                "feature",
                "gene_biotype",
                "gene_id",
                "gene_name",
            ]
        )
    return gtf

def gtf_to_protein_coding(gtfpath, log=Log(), verbose=True):
    protein_coding_path = gtfpath[:-6] + "protein_coding.gtf.gz"
    if not path.isfile(protein_coding_path):
        log.write(
            " - Extracting protein_coding genes from {}".format(gtfpath),
            verbose=verbose,
        )
        gtf = read_gtf(
            gtfpath,
            usecols=["feature", "gene_biotype", "gene_id", "gene_name"],
        )
        gene_list = (
            gtf.loc[
                (gtf["feature"] == "gene") & (gtf["gene_biotype"] == "protein_coding"),
                "gene_id",
            ]
            .values
        )
        log.write(
            " - Loaded {} protein_coding genes.".format(len(gene_list)),
            verbose=verbose,
        )
        gtf_raw = pd.read_csv(gtfpath, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[gtf_raw["_gene_id"].isin(gene_list), :]
        gtf_raw = gtf_raw.drop("_gene_id", axis=1)
        log.write(
            " - Extracted records are saved to : {} ".format(protein_coding_path),
            verbose=verbose,
        )
        gtf_raw.to_csv(protein_coding_path, header=None, index=None, sep="\t")
    return protein_coding_path

def gtf_to_all_gene(gtfpath, log=Log(), verbose=True):
    all_gene_path = gtfpath[:-6] + "all_genes.gtf.gz"
    if not path.isfile(all_gene_path):
        log.write(" - Extracting genes from {}".format(gtfpath), verbose=verbose)
        gtf = read_gtf(
            gtfpath,
            usecols=["feature", "gene_biotype", "gene_id", "gene_name"],
        )
        gene_list = gtf.loc[gtf["feature"] == "gene", "gene_id"].values
        log.write(" - Loaded {} genes.".format(len(gene_list)), verbose=verbose)
        gtf_raw = pd.read_csv(gtfpath, sep="\t", header=None, comment="#", dtype="string")
        gtf_raw["_gene_id"] = gtf_raw[8].str.extract(r'gene_id "([\w\.-]+)"')
        gtf_raw = gtf_raw.loc[gtf_raw["_gene_id"].isin(gene_list), :]
        gtf_raw = gtf_raw.drop("_gene_id", axis=1)
        log.write(
            " - Extracted records are saved to : {} ".format(all_gene_path),
            verbose=verbose,
        )
        gtf_raw.to_csv(all_gene_path, header=None, index=None, sep="\t")
    return all_gene_path
