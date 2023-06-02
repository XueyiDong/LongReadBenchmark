# Benchmarking long-read RNA-sequencing analysis tools using in *in silico* mixtures

This repository contains the code used to perform the analysis and generate the figures in this preprint:

**Benchmarking long-read RNA-sequencing analysis tools using in silico mixtures**
Xueyi Dong, Mei R. M. Du, Quentin Gouil, Luyi Tian, Jafar S. Jabbari, Rory Bowden, Pedro L. Baldoni, Yunshun Chen, Gordon K. Smyth, Shanika L. Amarasinghe, Charity W. Law, Matthew E. Ritchie
bioRxiv 2022.07.22.501076; doi: [https://doi.org/10.1101/2022.07.22.501076](https://doi.org/10.1101/2022.07.22.501076)

![Experimental design](ExpeDesign.png)

## Data Availability

Our RNA-seq data are available from Gene Expression Omnibus (GEO) under accession number [GSE172421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172421) (main benchmarking dataset) and [GSE227000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227000) (lab-based mixture of replicate 1). 

Please cite [our preprint](https://doi.org/10.1101/2022.07.22.501076) if you use our data and/or scripts in your studies.

## Index

### Pilot study

All scripts are available at [pilot](pilot)

### *In silico* mixture

ONT: [ONT/mix_prepare](ONT/mix_prepare)

Illumina: [illumina/mix_prepare](illumina/mix_prepare)

### Mapping and quantification

ONT: [ONT/preprocess](ONT/preprocess)

Illumina: [illumina/salmon_map.sh](illumina/salmon_map.sh)

### Quality control

ONT-specific: [ONT/QC](ONT/QC)

General: [longvsshort/overdisp.R](longvsshort/overdisp.R), [longvsshort/qc.R](longvsshort/qc.R) and [longvsshort/sequinCPMvsAbundance.R](longvsshort/sequinCPMvsAbundance.R)

### Isoform detection

Scripts to run softwares: [ONT/isoform_detection/methods](ONT/isoform_detection/methods)

Results comparison: [ONT/isoform_detection/analysis](ONT/isoform_detection/analysis)

### Differential transcript expression

ONT: [ONT/DE_mix.Rmd](ONT/DE_mix.Rmd)

Illumina: [illumina/DE_mix.Rmd](illumina/DE_mix.Rmd)

Results comparison: [longvsshort/DEmixres.R](longvsshort/DEmixres.R)

### Differential transcript usage

ONT: [ONT/DTU_mix.Rmd](ONT/DTU_mix.Rmd)

Illumina: [illumina/DTU_mix.Rmd](illumina/DTU_mix.Rmd)

Results comparison: [longvsshort/DTUmixres.R](longvsshort/DTUmixres.R) and [longvsshort/DTUmixrestx.R](longvsshort/DTUmixrestx.R) 
