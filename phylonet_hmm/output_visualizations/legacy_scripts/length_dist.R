setwd("/Users/matt/desktop/hmm_R_analysis")
x = read.csv("tract_dist.csv",header=F)
y=unlist(x)
h <- hist(y, main = "", xlab = "length (bases)", col = "lightblue",xlim=c(0,200000),ylim=c(0,250))
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

# main = "Introgression Tract Length Distribution \nFor Tracts Greater Than 10kb",