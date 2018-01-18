# GWAS_with_R
R code performing Genome-Wide-Association-Study and following a preset segregation model. The fit of the segregation model is tested on the dataset using a Chi2 Test.

Input format for the data is FASTA, with the IUPAC nucleotide characters.

## Example of input for three individuals at three SNP loci:

 >ind1 \
 ATY \
 >ind2 \
 AAY \
 >ind3 \
 TAN \

## Exemple of output:

site	di-/multiall	hyp	Chi2			P	\
1	diallelic	g	0.00698597000483793	0.99651310837294	\
2	diallelic	g	0.091387599469495	0.955334441069135	\
3	diallelic	t	0.091387599469495	0.955334441069135	\
4	diallelic	a	2.97826886015311	0.225567815781874	\
5	diallelic	c	3.39279321355224	0.183342992073151	\
6	diallelic	t	0.816850208892629	0.664696253421223 \


## Extra information:

The code create a function called association(). *This code only works for bi-allelic loci*

Within the code is the model you want to test. Here the model is fitted to a scenario where one population would be entirely heterozygous at a given locus while the other population would be homozygous at the locus.


To run the function: copy and paste the code you will have fitted to your data in R.
Then type association()
This will prompt you to select your input file.

The header of the code file contains other information for you to implement correctly the script.
