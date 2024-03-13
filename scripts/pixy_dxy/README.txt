1) compare_means_dxy.R visualizes the results of the boostrap comparison of means between dXY of introgression tracts and species tree tracts, performed in pixy_dxy.py. 

2) pixy_dxy.py calculates dXY for all introgression tracts and a random sample of intervals of the same size and length as that introgression tracts that were confidently called for the species tree. It then performs a bootstrap comparison of means.

3) pixy_dxy.R visualizes the distribution of dXY values for introgressed and non-introgressed intervals. 

4) popfile.txt is a required input file for running pixy to compute dXY. 

5) scaffolds.tsv is a tab-delimited text file that contains the name and length of each of the 21 largest scaffolds. 