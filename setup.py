import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gwaslab",
    version="3.3.4",
    author="Yunye",
    author_email="yunye@gwaslab.com",
    description="A collection of handy tools for GWAS",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Cloufield/gwaslab",
    project_urls={
        "gwaslab": "https://github.com/Cloufield/gwaslab",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    package_data={'gwaslab':['data/formatbook.json','data/high_ld/high_ld_hla_hg*.bed.gz','data/hapmap3_SNPs/*','data/RefGene/hg19/hg19.refGene.gtf.gz','data/RefGene/hg19/GRCh37/GRCh37_latest_genomic.protein_coding.gtf.gz','data/RefSeq/GRCh38/protein_coding.gtf.gz','data/Ensembl/release107/Homo_sapiens.GRCh38.107.protein_coding.chr.gtf.gz','data/Ensembl/release75/Homo_sapiens.GRCh37.75.protein_coding.gtf.gz']},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
