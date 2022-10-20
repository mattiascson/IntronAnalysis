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
png(file=paste("./Outpng/", Date,"_mut-vs-Wt_0h.png",sep=""), width=10.69, height=7.27, units="in", res=100)

#Set to 2x2 combiplot
par(mfrow=c(2,2))

#####med9 vs. Wt
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_med9_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
  med9RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med9")) %>% select(contains("0h")) %>% select(contains("30C"))
  WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
  # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
  Test <- try(t.test(as.numeric(med9RelExon),as.numeric(WtRelExon)))
  if (!inherits(Test, "try-error"))
  {
    Result <- tibble(GeneID=Gene, Mean_med9_RelExon = mean(as.numeric(med9RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med9RelExon),as.numeric(WtRelExon))$p.value)
  } else
  {
    Result <- tibble(GeneID=Gene, Mean_med9_RelExon = mean(as.numeric(med9RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
  }
  Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med9-Wt" = Mean_med9_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID %in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_med9-Wt"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_med9-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairMediator <- signif(t.test(Results[["Mean_med9_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_med9-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_med9vsWt_0h.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_med9_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    med9RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med9")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(med9RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_med9_RelExon = mean(as.numeric(med9RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med9RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_med9_RelExon = mean(as.numeric(med9RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med9-Wt" = Mean_med9_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_med9-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairOther <- signif(t.test(Results[["Mean_med9_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_med9-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_med9vsWt_0h.csv",sep=""))

MediatorOther <- signif(t.test(Results[["Difference_med9-Wt"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_med9-Wt"]],MediatorGenesDifferance, main=paste("SpliceRatio med9-Wt 0h. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u", MeanOther, "Pair pval", PairOther),paste("u", MeanMediator, "Pair pval", PairMediator)))

#####cdk8 vs. Wt
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_cdk8_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    cdk8RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("cdk8")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(cdk8RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_cdk8_RelExon = mean(as.numeric(cdk8RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(cdk8RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_cdk8_RelExon = mean(as.numeric(cdk8RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_cdk8-Wt" = Mean_cdk8_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID %in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_cdk8-Wt"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_cdk8-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairMediator <- signif(t.test(Results[["Mean_cdk8_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_cdk8-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_cdk8vsWt_0h.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_cdk8_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    cdk8RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("cdk8")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(cdk8RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_cdk8_RelExon = mean(as.numeric(cdk8RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(cdk8RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_cdk8_RelExon = mean(as.numeric(cdk8RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_cdk8-Wt" = Mean_cdk8_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_cdk8-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairOther <- signif(t.test(Results[["Mean_cdk8_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_cdk8-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_cdk8vsWt_0h.csv",sep=""))

MediatorOther <- signif(t.test(Results[["Difference_cdk8-Wt"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_cdk8-Wt"]],MediatorGenesDifferance, main=paste("SpliceRatio cdk8-Wt 0h. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u", MeanOther, "Pair pval", PairOther),paste("u", MeanMediator, "Pair pval", PairMediator)))


#####med16 vs. Wt
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_med16_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    med16RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med16")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(med16RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_med16_RelExon = mean(as.numeric(med16RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med16RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_med16_RelExon = mean(as.numeric(med16RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med16-Wt" = Mean_med16_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID %in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_med16-Wt"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_med16-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairMediator <- signif(t.test(Results[["Mean_med16_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_med16-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_med16vsWt_0h.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_med16_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    med16RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med16")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(med16RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_med16_RelExon = mean(as.numeric(med16RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med16RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_med16_RelExon = mean(as.numeric(med16RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med16-Wt" = Mean_med16_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_med16-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairOther <- signif(t.test(Results[["Mean_med16_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_med16-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_med16vsWt_0h.csv",sep=""))

MediatorOther <- signif(t.test(Results[["Difference_med16-Wt"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_med16-Wt"]],MediatorGenesDifferance, main=paste("SpliceRatio med16-Wt 0h. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u", MeanOther, "Pair pval", PairOther),paste("u", MeanMediator, "Pair pval", PairMediator)))


#####med18 vs. Wt
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_med18_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    med18RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med18")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(med18RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_med18_RelExon = mean(as.numeric(med18RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med18RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_med18_RelExon = mean(as.numeric(med18RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med18-Wt" = Mean_med18_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID %in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_med18-Wt"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_med18-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairMediator <- signif(t.test(Results[["Mean_med18_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_med18-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_med18vsWt_0h.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_med18_RelExon = numeric(),Mean_Wt_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    med18RelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("med18")) %>% select(contains("0h")) %>% select(contains("30C"))
    WtRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("Wt")) %>% select(contains("0h")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(med18RelExon),as.numeric(WtRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_med18_RelExon = mean(as.numeric(med18RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=t.test(as.numeric(med18RelExon),as.numeric(WtRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_med18_RelExon = mean(as.numeric(med18RelExon)),Mean_Wt_RelExon=mean(as.numeric(WtRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_med18-Wt" = Mean_med18_RelExon - Mean_Wt_RelExon, .after = Mean_Wt_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_med18-Wt"]], mu = 0, alternative = "two.sided")$p.value,2)
PairOther <- signif(t.test(Results[["Mean_med18_RelExon"]],Results[["Mean_Wt_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_med18-Wt"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_med18vsWt_0h.csv",sep=""))

MediatorOther <- signif(t.test(Results[["Difference_med18-Wt"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_med18-Wt"]],MediatorGenesDifferance, main=paste("SpliceRatio med18-Wt 0h. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u", MeanOther, "Pair pval", PairOther),paste("u", MeanMediator, "Pair pval", PairMediator)))

dev.off()
