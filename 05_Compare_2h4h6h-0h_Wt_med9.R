library(tidyverse)
library(magrittr)

#Set Date
Date <- gsub("-","",Sys.Date())

#Read files and Date dependent files
AllIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_AllIntronGenes.csv",sep=""))
OtherIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_OtherIntronGenes.csv",sep=""))
MediatorGenes <- read_csv(file=paste("./InOutFiles/", Date,"_MediatorGenes.csv",sep=""))

Annotation <- read_csv(file=paste("./InOutFiles/", Date,"_Annotation.csv",sep=""))

#Set saved image	(8.27x11.69
png(file=paste("./Outpng/", Date,"_med9-vs-Wt_2h4h6h-0h.png",sep=""), width=11.69, height=7.27, units="in", res=100)

#Set to 3x3 combiplot
par(mfrow=c(3,3))

source(file = "051_Compare_2h-0h_med9-vs-Wt.R")
source(file = "052_Compare_4h-0h_med9-vs-Wt.R")
source(file = "053_Compare_6h-0h_med9-vs-Wt.R")

dev.off()
