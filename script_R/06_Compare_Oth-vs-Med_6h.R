library(tidyverse)
library(magrittr)

#Set Date
Date <- gsub("-","",Sys.Date())

#Read files and Date dependent files
AllIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_AllIntronGenes.csv",sep=""))
OtherIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_OtherIntronGenes.csv",sep=""))
MediatorGenes <- read_csv(file=paste("./InOutFiles/", Date,"_MediatorGenes.csv",sep=""))

Annotation <- read_csv(file=paste("./InOutFiles/20221001_Annotation.csv",sep=""))

#Set saved image to 1 in smaller than A4  	(8.27x11.69
#png(file=paste("./Outpng/", Date,"_mut-vs-Wt_6h.png",sep=""), width=10.69, height=7.27, units="in", res=100)

#Set to 2x2 combiplot
par(mfrow=c(2,3))

Layout <- rbind(c(1,2,3),
                c(4,5,0))
layout(Layout)

#####Wt
###Check Mediator intron genes.
WtMedRelExon <- AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(contains("Wt"))%>% select(contains("30C")) %>% select(contains("6h"))
WtMedRelExon %<>% mutate(Mean=rowMeans(.))
WtMedRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtMedRelExon)
###Check Other intron genes.
WtOthRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("Wt"))%>% select(contains("30C")) %>% select(contains("6h"))
WtOthRelExon %<>% mutate(Mean=rowMeans(.))
WtOthRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtOthRelExon)

MeanMediator <- signif(mean(WtMedRelExon[["Mean"]],),2)
nMediator <- dim(WtMedRelExon)[1]
MeanOther <- signif(mean(WtOthRelExon[["Mean"]]),2)
nOther<- dim(WtOthRelExon)[1]
MediatorOther <- signif(t.test(WtOthRelExon[["Mean"]],WtMedRelExon[["Mean"]], alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(WtOthRelExon[["Mean"]],WtMedRelExon[["Mean"]], main=paste("SpliceRatio Wt 30C 6h. t.test pval", MediatorOther), names=c(paste(nOther, "Other u", MeanOther),paste(nMediator, " Mediator u", MeanMediator)))
pValues <- c(pValues,MediatorOther)


#####med9
###Check Mediator intron genes.
med9MedRelExon <- AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(contains("med9"))%>% select(contains("30C")) %>% select(contains("6h"))
med9MedRelExon %<>% mutate(Mean=rowMeans(.))
med9MedRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med9MedRelExon)
###Check Other intron genes.
med9OthRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med9"))%>% select(contains("30C")) %>% select(contains("6h"))
med9OthRelExon %<>% mutate(Mean=rowMeans(.))
med9OthRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med9OthRelExon)

MeanMediator <- signif(mean(med9MedRelExon[["Mean"]],),2)
nMediator <- dim(WtMedRelExon)[1]
MeanOther <- signif(mean(med9OthRelExon[["Mean"]]),2)
nOther<- dim(WtOthRelExon)[1]
MediatorOther <- signif(t.test(med9OthRelExon[["Mean"]],med9MedRelExon[["Mean"]], alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(med9OthRelExon[["Mean"]],med9MedRelExon[["Mean"]], main=paste("SpliceRatio med9 30C 6h. t.test pval", MediatorOther), names=c(paste(nOther, "Other u", MeanOther),paste(nMediator,"Mediator u", MeanMediator)))
pValues <- c(pValues,MediatorOther)


#####cdk8
###Check Mediator intron genes.
cdk8MedRelExon <- AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(contains("cdk8"))%>% select(contains("30C")) %>% select(contains("6h"))
cdk8MedRelExon %<>% mutate(Mean=rowMeans(.))
cdk8MedRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),cdk8MedRelExon)
###Check Other intron genes.
cdk8OthRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("cdk8"))%>% select(contains("30C")) %>% select(contains("6h"))
cdk8OthRelExon %<>% mutate(Mean=rowMeans(.))
cdk8OthRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),cdk8OthRelExon)

MeanMediator <- signif(mean(cdk8MedRelExon[["Mean"]],),2)
nMediator <- dim(WtMedRelExon)[1]
MeanOther <- signif(mean(cdk8OthRelExon[["Mean"]]),2)
nOther<- dim(WtOthRelExon)[1]
MediatorOther <- signif(t.test(cdk8OthRelExon[["Mean"]],cdk8MedRelExon[["Mean"]], alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(cdk8OthRelExon[["Mean"]],cdk8MedRelExon[["Mean"]], main=paste("SpliceRatio cdk8 30C 6h. t.test pval", MediatorOther), names=c(paste(nOther, "Other u", MeanOther),paste(nMediator,"Mediator u", MeanMediator)))
pValues <- c(pValues,MediatorOther)


#####med16
###Check Mediator intron genes.
med16MedRelExon <- AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(contains("med16"))%>% select(contains("30C")) %>% select(contains("6h"))
med16MedRelExon %<>% mutate(Mean=rowMeans(.))
med16MedRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med16MedRelExon)
###Check Other intron genes.
med16OthRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med16"))%>% select(contains("30C")) %>% select(contains("6h"))
med16OthRelExon %<>% mutate(Mean=rowMeans(.))
med16OthRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med16OthRelExon)

MeanMediator <- signif(mean(med16MedRelExon[["Mean"]],),2)
nMediator <- dim(WtMedRelExon)[1]
MeanOther <- signif(mean(med16OthRelExon[["Mean"]]),2)
nOther<- dim(WtOthRelExon)[1]
MediatorOther <- signif(t.test(med16OthRelExon[["Mean"]],med16MedRelExon[["Mean"]], alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(med16OthRelExon[["Mean"]],med16MedRelExon[["Mean"]], main=paste("SpliceRatio med16 30C 6h. t.test pval", MediatorOther), names=c(paste(nOther, "Other u", MeanOther),paste(nMediator,"Mediator u", MeanMediator)))
pValues <- c(pValues,MediatorOther)


#####med18
###Check Mediator intron genes.
med18MedRelExon <- AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(contains("med18"))%>% select(contains("30C")) %>% select(contains("6h"))
med18MedRelExon %<>% mutate(Mean=rowMeans(.))
med18MedRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% MediatorGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med18MedRelExon)
###Check Other intron genes.
med18OthRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med18"))%>% select(contains("30C")) %>% select(contains("6h"))
med18OthRelExon %<>% mutate(Mean=rowMeans(.))
med18OthRelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med18OthRelExon)

MeanMediator <- signif(mean(med18MedRelExon[["Mean"]],),2)
nMediator <- dim(WtMedRelExon)[1]
MeanOther <- signif(mean(med18OthRelExon[["Mean"]]),2)
nOther<- dim(WtOthRelExon)[1]
MediatorOther <- signif(t.test(med18OthRelExon[["Mean"]],med18MedRelExon[["Mean"]], alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(med18OthRelExon[["Mean"]],med18MedRelExon[["Mean"]], main=paste("SpliceRatio med18 30C 6h. t.test pval", MediatorOther), names=c(paste(nOther, "Other u", MeanOther),paste(nMediator,"Mediator u", MeanMediator)))
pValues <- c(pValues,MediatorOther)


