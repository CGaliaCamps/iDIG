# iDlG: individual Detection of linkage by Genotyping
Our new approach of plotting individual genotypes along the genome greatly facilitates the detection of chromosome inversions, allowing the identification of individual-level presence of inversions, their area of influence, and chromosome-inversion specific karyotype, all at once.

# Usage
To use this software, you just need to gather individuals mapped against a chromosome-level reference genome and their genotypes in 012 format. 
**You can obtain the 012 file** from an standard vcf file by **using the function --012 of the software vcftools** (https://vcftools.sourceforge.net/man_latest.html).

Thus, you will need 3 different files with different extensions:
 

_"*.012" -> This file contains your genotypes for each individual inlcuded into the analysis_

_"*.012.indv" -> This file contains the labels for your individuals_

_"*.012.pos" -> This file contains the position of all your SNPs_
 

Copy the script **"Detect_Inversions.R"** into your R terminal, and modify the working directory and filenames to be uploaded. After that, just run it.

# Citation
Please, if you found this script useful cite:

Galià-Camps C., Schell T., Pegueroles C., Baranski, D., Ben Hamadou A., Turon X., Pascual M., Greve C., Carreras C. (2013). Genomic richness enables worldwide invasive success. _Research Square_ 10.21203/rs.3.rs-3902873/v1

