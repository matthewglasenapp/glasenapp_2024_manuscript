# cat out.min6.used_sites.tsv | grep -v POS | awk '{ print $2 }' > coordinates

setwd("/Users/matt/desktop/hmm_R_analysis/")
library("rjson")
dt=fromJSON(file="rawOutput.json")
summary(dt)
y=dt$posteriorProbabilityOfSpeciesTrees[[2]]

# Need to re download coordinate file
ds=read.table("coordinates",header=F)
x=as.vector(ds[,1])

plot(y[1011120:1011280]~x[1011120:1011280],col='red',type='l',lty=1, xlab="position on scaffold NW_022145612.1", ylab = 
       "introgressed probability")


data=data.frame(x,y)
length(data$y[data$y>0.9])
length(data$y)

length(data$y[data$y>0.99]) /length(data$y)
