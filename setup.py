import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gwaslab",
    version="3.0.2",
    author="Yunye",
    author_email="610935659@qq.com",
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
    package_data={"gwaslab":["data/hapmap3_rs_chr_pos_a1_a2/*"]},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
