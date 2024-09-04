# ROHbin 
![image140](https://github.com/user-attachments/assets/ec4048d2-5f7d-4f47-8503-188ddfa8ce1d)

A probabilistic framework to test for runs of homozygosity (ROHs) between groups of populations. 
Developed specifically for the whooping crane paper. If you use it, please cite:

**Fontsere et al.,** "Genomics of the population crash, successful recovery and captive breeding of the whooping crane." bioRxiv (2024). DOI:

## Usage
The main script 000_ROHbin.R reads the metadata for a set of individuals and the heterozygosity for a given window size. Then, define the comparisons to be applied (ALL WILD vs. ALL CAPTIVE, ALL WILD vs. LATE CAPTIVE, etc...) in a design matrix and runs empirical Bayes (see this nice tutorial for an explanation on how eBayes works http://varianceexplained.org/r/empirical_bayes_baseball/). This will give the heterozygostiy windows that are significantly different between groups.
The script will also generate plots for each variable bin in heterozygosity, centered in the significant bin or bins. DBSCAN approach is used to cluster together significant, continuous bins.
