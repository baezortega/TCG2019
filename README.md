## Analysis scripts accompanying publication Baez-Ortega _et al._, 2019

_Adrian Baez-Ortega  
Transmissible Cancer Group, University of Cambridge  
2018_

This repository contains custom R and Python scripts which can be used to replicate the analyses performed in Baez-Ortega _et al._, 2019 [**UNPUBLISHED**]. The set-up and usage of the data and scripts are explained below.

---

### Set-up

In order to run the scripts, the repository must first be cloned into a directory in your local machine (replace the destination path `~/Desktop/TCG2019` below with your own one). This can be done from the terminal using the command below.

```
git clone https://github.com/baezortega/TCG2019.git ~/Desktop/TCG2019
```

Then, we need to download the necessary supporting data from the University of Cambridge Repository, which can be found at https://doi.org/10.17863/CAM.24962 [**DATA NOT YET AVAILABLE**]. All the supporting data files in this repository must be downloaded and **placed in the `TCG2019/data/original` directory**.

The scripts are written to run on **[Python](https://www.python.org/)** version 2.7.10 or later (not Python 3) and **[R](https://www.r-project.org/)** version 3.3.3 or later. You should be able to check your current version of R by running one of the commands below (depending on your installation):

```
R --version
/Library/Frameworks/R.framework/Resources/bin/R --version
```

The current R version is also shown when opening RStudio or the R Console. Similarly, you should be able to check your current version of Python by running the command `python --version`.

The scripts require the following R packages: [`ape`](https://cran.r-project.org/web/packages/ape/index.html), [`RColorBrewer`](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [`stringr`](https://cran.r-project.org/web/packages/stringr/index.html), [`sigfit`](https://github.com/kgori/sigfit).

Although care has been taken to make the code distribution-independent, it is possible that some of the scripts work only on Unix (MacOS/Linux) systems, and need to be modified in order to run on Windows systems.

Please note that some of the scripts may take considerable time to run, and some of the intermediate files generated will occupy up to a few gigabytes.

Finally, before running the steps indicated below, we first need to open the terminal and navigate to the cloned directory, `TCG2019`, using the command below (replace the path below with your own one). 

```
cd ~/Desktop/TCG2019
```

---

### Step 1. Extraction of variant information from VCF files

Before importing the variant calls from Somatypus into R, we will extract the information we are interested in, namely the variant metadata (fields `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`), the numbers of total reads and supporting reads at the variant locus (subfields `NR` and `NV` in the `INFO` field, respectively). Extracting the information of interest in this way will avoid the need for reading massive VCF files directly from R.

This step is carried out by the Python script `1_ExtractVcfData.py`, which is located in the [`scripts`](scripts) directory and can be run as follows.

```
scripts/1_ExtractVcfData.py <input-variants.vcf.gz> <output-dir>
```

We run this script for each of the two VCF files produced by Somatypus, which should be located in the `data/original` directory, storing the output files in `data/processed`.

```
scripts/1_ExtractVcfData.py data/original/Somatypus_CTVT_SNVs_1052.vcf.gz data/processed
scripts/1_ExtractVcfData.py data/original/Somatypus_CTVT_Indels_1052.vcf.gz data/processed
```

The script produces three tab-delimited text files for each VCF, labelled with the suffixes `Metadata`, `NR` and `NV`. These are the files that will be read into R in the next step.

---

### Step 2. Variant import and classification

This step is carried out by the R script `2_ImportVariants.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. (Depending on your installation, you may have to replace `Rscript` with **`/Library/Frameworks/R.framework/Resources/bin/Rscript`**.)

```
Rscript scripts/2_ImportVariants.R
```

This script produces an output RData file (`data/processed/VariantTables.RData`) containing the data tables for SNVs and indels, together with indices to distinguish somatic and germline variants.

---

### Step 3. Variant annotation processing

This step is carried out by the R script `3_ProcessAnnotation.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/3_ProcessAnnotation.R
```

This script produces an output RData file (`data/processed/VariantAnnotation.RData`) containing the processed annotation tables for somatic and germline SNVs and indels.

---

### Step 4. Exploratory variant analyses

This step is carried out by the R script `4_ExploratoryAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows.

```
Rscript scripts/4_ExploratoryAnalyses.R
```

This script produces an output RData file (`data/processed/Coverage_Purity_QCindex.RData`) containing coverage, purity and 'QC index' values for every tumour sample. It also produces histograms of the coverage and VAF in tumours and hosts (`output/Coverage_Histograms_Tumours.pdf`,`output/Coverage_Histograms_Hosts.pdf`, `output/VAF_Histograms_Tumours.pdf`, `output/VAF_Histograms_Hosts.pdf`) and variant summary tables for tumours and hosts (`output/Variant_Summary_Tumours.tsv`, `output/Variant_Summary_Hosts.tsv`).

---

### Step 5. Definition of phylogenetic tumour groups and group-specific variants

This step is carried out by the R script `5_PhyloGroups.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/5_PhyloGroups.R
```

This script produces an output RData file (`data/processed/Phylo_Tree_Groups.RData`) containing the CTVT phylogenetic tree and matrix indices of the phylogenetic tumour groups and group-unique variant sets. It also produces plots of the phylogenetic tree (`output/CTVT_Tree_RAxML.pdf`) and of the mutational spectra of group-unique somatic SNVs in each phylogenetic group (`output/Phylo_Groups_Spectra.pdf`).

---

### Step 6. Somatic mutation prevalence estimation

This step is carried out by the R script `6_PrevalenceAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/6_PrevalenceAnalyses.R
```

This script produces plots of the somatic mutation prevalence in CTVT and a set of human cancer types (`output/Somatic_Mutation_Prevalence.pdf`).

---

### Step 7. Mutational spectrum and mutational signature analyses

This step is carried out by the R script `7_SignatureAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/7_SignatureAnalyses.R
```

This script produces ...

---

### Step 8. Mutational spectrum analysis of dinucleotide substitutions

This step is carried out by the R script `8_DinucAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/8_DinucAnalyses.R
```

This script produces ...

---

### Step 9. Inference of associations between signature 7 exposure, CC>TT mutations and latitude

This step is carried out by the R script `9_UV-LatAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/9_UV-LatAnalyses.R
```

This script produces ...

---

### Step 10. Processing and analysis of recurrent mutations

This step is carried out by the R script `10_RecMutAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/10_RecMutAnalyses.R
```

This script produces ...

---

### Step 11. Gene expression analyses

This step is carried out by the R script `11_ExpressionAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/11_ExpressionAnalyses.R
```

This script produces ...

---

### Step 12. Somatic mutation burden analyses incorporating gene expression and transcriptional strand information

This step is carried out by the R script `12_BurdenAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/12_BurdenAnalyses.R
```

This script produces ...

---

### Step 13. Selection analyses

This step is carried out by the R script `13_SelectionAnalyses.R`, which is located in the [`scripts`](scripts) directory and can be run either from RStudio (following the instructions at the beginning of the script), or from the terminal using the `Rscript` command as follows. 

```
Rscript scripts/13_SelectionAnalyses.R
```

This script produces ...
