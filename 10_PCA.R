library(tidyverse)
library(magrittr)
library(broom)
library(cowplot)
#Set Date
Date <- gsub("-","",Sys.Date())

#Read files and Date dependent files
AllIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_AllIntronGenes.csv",sep=""))
OtherIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_OtherIntronGenes.csv",sep=""))
MediatorGenes <- read_csv(file=paste("./InOutFiles/", Date,"_MediatorGenes.csv",sep=""))

#Read debatch files and Date dependent files
AllIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_AllIntronGenesDebatch.csv",sep=""))
OtherIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_OtherIntronGenesDebatch.csv",sep=""))
MediatorGenes <- read_csv(file=paste("./InOutFiles/", Date,"_MediatorGenesDebatch.csv",sep=""))

Annotation <- read_csv(file=paste("./InOutFiles/20221001_Annotation.csv",sep=""))

#Set saved image	(8.27x11.69
#png(file=paste("./Outpng/", Date,"_med9-vs-Wt_2h4h6h-0h.png",sep=""), width=11.69, height=7.27, units="in", res=100)
#Set to 3x3 combiplot
#par(mfrow=c(3,3))



#Do the plotting
#dev.off()

Test <- matrix(c(c(1,2,3),
                 c(4,6,7),
                 c(8,3,2)),
byrow = TRUE,ncol=3)
row.names(Test)
