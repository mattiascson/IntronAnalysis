library(tidyverse)
library(magrittr)

#Set Date
Date <- gsub("-","",Sys.Date())

#Read files and Date dependent files
AllIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_AllIntronGenesDebatch.csv",sep=""))
OtherIntronGenes <- read_csv(file=paste("./InOutFiles/", Date,"_OtherIntronGenesDebatch.csv",sep=""))
MediatorGenes <- read_csv(file=paste("./InOutFiles/", Date,"_MediatorGenesDebatch.csv",sep=""))

Annotation <- read_csv(file=paste("./InOutFiles/20221001_Annotation.csv",sep=""))

#Set to 2x2 combiplot
par(mfrow=c(2,2))

Layout <- rbind(c(1,2),
                c(3,4))
layout(Layout)

#####Wt vs med9
###Check Wt intron genes.
WtRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("Wt"))%>% select(contains("30C")) %>% select(contains("2h"))
WtRelExon %<>% mutate(Mean=rowMeans(.))
WtRelExon <- tibble(AllIntronGenes  %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtRelExon)
###Check Mediator intron genes.
med9RelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med9"))%>% select(contains("30C")) %>% select(contains("2h"))
med9RelExon %<>% mutate(Mean=rowMeans(.))
med9RelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med9RelExon)
###Make Calculations and Plot
MeanWt <- signif(mean(WtRelExon[["Mean"]],),2)
nWt <- dim(WtRelExon)[1]
Meanmed9 <- signif(mean(med9RelExon[["Mean"]]),2)
nmed9<- dim(med9RelExon)[1]
Wtmed <- signif(t.test(WtRelExon[["Mean"]],med9RelExon[["Mean"]], alternative = "two.sided", paired = TRUE)$p.value, 2)
boxplot(WtRelExon[["Mean"]],med9RelExon[["Mean"]], main=paste("SpliceRatio Oth Wt vs med9 30C 2h. t.test pair pval", Wtmed), names=c(paste(nWt, "Wt u", MeanWt),paste(nmed9, " med9 u", Meanmed9)))
pValues <- c(pValues,Wtmed)

#####Wt vs cdk8
###Check Wt intron genes.
WtRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("Wt"))%>% select(contains("30C")) %>% select(contains("2h"))
WtRelExon %<>% mutate(Mean=rowMeans(.))
WtRelExon <- tibble(AllIntronGenes  %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtRelExon)
###Check Mediator intron genes.
cdk8RelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("cdk8"))%>% select(contains("30C")) %>% select(contains("2h"))
cdk8RelExon %<>% mutate(Mean=rowMeans(.))
cdk8RelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),cdk8RelExon)
###Make Calculations and Plot
MeanWt <- signif(mean(WtRelExon[["Mean"]],),2)
nWt <- dim(WtRelExon)[1]
Meancdk8 <- signif(mean(cdk8RelExon[["Mean"]]),2)
ncdk8<- dim(cdk8RelExon)[1]
Wtmed <- signif(t.test(WtRelExon[["Mean"]],cdk8RelExon[["Mean"]], alternative = "two.sided", paired = TRUE)$p.value, 2)
boxplot(WtRelExon[["Mean"]],cdk8RelExon[["Mean"]], main=paste("SpliceRatio Oth Wt vs cdk8 30C 2h. t.test pair pval", Wtmed), names=c(paste(nWt, "Wt u", MeanWt),paste(ncdk8, " cdk8 u", Meancdk8)))
pValues <- c(pValues,Wtmed)

#####Wt vs med16
###Check Wt intron genes.
WtRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("Wt")) %>% select(contains("30C")) %>% select(contains("2h"))
WtRelExon %<>% mutate(Mean=rowMeans(.))
WtRelExon <- tibble(AllIntronGenes  %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtRelExon)
###Check Mediator intron genes.
med16RelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med16"))%>% select(contains("30C")) %>% select(contains("2h"))
med16RelExon %<>% mutate(Mean=rowMeans(.))
med16RelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med16RelExon)
###Make Calculations and Plot
MeanWt <- signif(mean(WtRelExon[["Mean"]],),2)
nWt <- dim(WtRelExon)[1]
Meanmed16 <- signif(mean(med16RelExon[["Mean"]]),2)
nmed16<- dim(med16RelExon)[1]
Wtmed <- signif(t.test(WtRelExon[["Mean"]],med16RelExon[["Mean"]], alternative = "two.sided", paired = TRUE)$p.value, 2)
boxplot(WtRelExon[["Mean"]],med16RelExon[["Mean"]], main=paste("SpliceRatio Oth Wt vs med16 30C 2h. t.test pair pval", Wtmed), names=c(paste(nWt, "Wt u", MeanWt),paste(nmed16, " med16 u", Meanmed16)))
pValues <- c(pValues,Wtmed)

#####Wt vs med18
###Check Wt intron genes.
WtRelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("Wt")) %>% select(contains("30C")) %>% select(contains("2h"))
WtRelExon %<>% mutate(Mean=rowMeans(.))
WtRelExon <- tibble(AllIntronGenes  %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),WtRelExon)
###Check Mediator intron genes.
med18RelExon <- AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(contains("med18"))%>% select(contains("30C")) %>% select(contains("2h"))
med18RelExon %<>% mutate(Mean=rowMeans(.))
med18RelExon <- tibble(AllIntronGenes %>% filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>% select(GeneID,GeneName,LengthDiff),med18RelExon)
###Make Calculations and Plot
MeanWt <- signif(mean(WtRelExon[["Mean"]],),2)
nWt <- dim(WtRelExon)[1]
Meanmed18 <- signif(mean(med18RelExon[["Mean"]]),2)
nmed18<- dim(med18RelExon)[1]
Wtmed <- signif(t.test(WtRelExon[["Mean"]],med18RelExon[["Mean"]], alternative = "two.sided", paired = TRUE)$p.value, 2)
boxplot(WtRelExon[["Mean"]],med18RelExon[["Mean"]], main=paste("SpliceRatio Oth Wt vs med18 30C 2h. t.test pair pval", Wtmed), names=c(paste(nWt, "Wt u", MeanWt),paste(nmed18, " med18 u", Meanmed18)))
pValues <- c(pValues,Wtmed)