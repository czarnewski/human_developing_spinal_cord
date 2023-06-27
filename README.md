# human_developing_spinal_cord.  

***
### Citation
This is a public repository containing scripts used in the publication:

Li X et al. 
Decoding spatiotemporal gene expression of the developing human spinal cord and implications for ependymoma origin

**Nature Neuroscience** 2023 (https://www.nature.com/articles/s41593-023-01312-9)

***
### Running the analysis

The analysis scripts are available in the `code` folder  and is empty by default.
The `data` folder is used to host the raw counts and other files used for the analysis.
Some datasets require manual download and should be placed in the corresponding folder. The `results` folder,
also empty by default, will store the output from analysis performed.
The analysis done herein can be reproduced by installing CONDA and running:

1. Clone this repository\
```
git clone https://github.com/czarnewski/human_developing_spinal_cord.git
```

2. Create and activate the conda environment\
```
cd human_developing_spinal_cord

conda activate base
conda install -c conda-forge mamba

mamba env create -n dev_sc -f env_dev_sc.yml
conda activate dev_sc
```
***
### Datasets

The list of all datasets used in the manuscript are depicted below:

| Technology | Dataset | source publication | Accession no |
|------------|---------|--------------------|--------------|
| 10X scRNAseq | Human developing spinal cord | this manuscript | [GSE219122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219122) |](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219122) |
| 10X Visium | Human developing spinal cord | this manuscript | [GSE219122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219122) |](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219122) |
| 10X scRNAseq | Mouse adult spinal cord | [Sathyamurthy et al. 2018 Cell Rep] (https://pubmed.ncbi.nlm.nih.gov/29466745/)| [GSE103892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103892) |
| 10X scRNAseq | Mouse juvenile spinal cord | [Zeisel et al. 2018 Cell] (https://pubmed.ncbi.nlm.nih.gov/30096314/)| [SRP135960](https://www.ncbi.nlm.nih.gov/sra/SRP135960) |
| SplitSeq | Mouse postnatal | [Rosenberg et al. 2018 Science] (https://pubmed.ncbi.nlm.nih.gov/29545511/)| [GSE110823](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110823) |
| 10X scRNAseq | Mouse adult spinal cord | [Blum et al. 2021  Nat Neurosci](https://pubmed.ncbi.nlm.nih.gov/33589834/) | [GSE161621](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161621) |
| 10X scRNAseq | Mouse adult spinal cord | [Alkaslasi et al. 2021 Nat Comm](https://pubmed.ncbi.nlm.nih.gov/33931636/)| [GSE167597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167597) |
| 10X scRNAseq | Mouse developing spinal cord | [Delile et al 2019 Development](https://pubmed.ncbi.nlm.nih.gov/30846445/) | [E-MTAB-7320](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7320/files) |
| 10X scRNAseq | Human developing spinal cord | [Rayon et al 2021 Development](https://pubmed.ncbi.nlm.nih.gov/34351410/) | [GSE171892](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171892) |
| 10X scRNAseq | Mouse juvenile spinal cord | [Milich et al 2021 J Exp Med](https://pubmed.ncbi.nlm.nih.gov/34132743/) | [GSE162610](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162610) |
| 10X scRNAseq | Human developing spinal cord | [Zhang et al 2021 EMBO Rep](https://pubmed.ncbi.nlm.nih.gov/34605607/) | [GSE136719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136719) |
| 10X scRNAseq | Human ependymomas | [Gojo et al 2020 Cancer Cell](https://pubmed.ncbi.nlm.nih.gov/32663469/) | [GSE136719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141460) |
