####Wt
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_T2h_RelExon = numeric(),Mean_T0h_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
  T2hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("2h")) %>% select(contains("Wt")) %>% select(contains("30C"))
  T0hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("0h")) %>% select(contains("Wt")) %>% select(contains("30C"))
  # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
  Test <- try(t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon)))
  if (!inherits(Test, "try-error"))
  {
    Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon))$p.value)
  } else
  {
    Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=1)
  }
  Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_T2h-T0h" = Mean_T2h_RelExon - Mean_T0h_RelExon, .after = Mean_T0h_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_T2h-T0h"]]
MediatorGenesID <- Results[["GeneID"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_T2h-T0h"]], mu = 0, alternative = "two.sided")$p.value,2)
PairedMediator <- signif(t.test(Results[["Mean_T2h_RelExon"]],Results[["Mean_T0h_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_T2h-T0h"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_T2hvsT0h_Wt.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_T2h_RelExon = numeric(),Mean_T0h_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    T2hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("2h")) %>% select(contains("Wt")) %>% select(contains("30C"))
    T0hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("0h")) %>% select(contains("Wt")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_T2h-T0h" = Mean_T2h_RelExon - Mean_T0h_RelExon, .after = Mean_T0h_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_T2h-T0h"]], mu = 0, alternative = "two.sided")$p.value,2)
PairedOther <- signif(t.test(Results[["Mean_T2h_RelExon"]],Results[["Mean_T0h_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_T2h-T0h"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_T2hvsT0h_Wt.csv",sep=""))
MediatorOther <- signif(t.test(Results[["Difference_T2h-T0h"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_T2h-T0h"]],MediatorGenesDifferance, main=paste("SR T2h-T0h Wt. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u",MeanOther , "Pair pval", PairedOther),paste("u", MeanMediator, "Pair pval", PairedMediator)))
###To compare between Wt and mutants. Run only for Wt
GeneIDWt <- MediatorGenesID
MediatorGenesDifferanceWt <- MediatorGenesDifferance
MeanMediatorWt <- MeanMediator
MeanOtherWt <- MeanOther

####med9
###Check Mediator intron genes.
Results <- tibble(GeneID = character(),Mean_T2h_RelExon = numeric(),Mean_T0h_RelExon = numeric(),ttest_pval = numeric())
for (Gene in MediatorGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    T2hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("2h")) %>% select(contains("med9")) %>% select(contains("30C"))
    T0hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("0h")) %>% select(contains("med9")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_T2h-T0h" = Mean_T2h_RelExon - Mean_T0h_RelExon, .after = Mean_T0h_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
MediatorGenesDifferance <- Results[["Difference_T2h-T0h"]]
MediatorGenesID <- Results[["GeneID"]]
DiffZeroMediator <- signif(t.test(Results[["Difference_T2h-T0h"]], mu = 0, alternative = "two.sided")$p.value,2)
PairedMediator <- signif(t.test(Results[["Mean_T2h_RelExon"]],Results[["Mean_T0h_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanMediator <- signif(mean(Results[["Difference_T2h-T0h"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_MediatorGenes_T2hvsT0h_med9.csv",sep=""))

###Check Other intron genes.
Results <- tibble(GeneID = character(),Mean_T2h_RelExon = numeric(),Mean_T0h_RelExon = numeric(),ttest_pval = numeric())
for (Gene in OtherIntronGenes[["GeneID"]])
{
  if(any(Gene==AllIntronGenes[["GeneID"]]))
  {
    T2hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("2h")) %>% select(contains("med9")) %>% select(contains("30C"))
    T0hRelExon <- AllIntronGenes %>% filter(GeneID == Gene) %>% select(contains("0h")) %>% select(contains("med9")) %>% select(contains("30C"))
    # set pval to 1 if comparing e.g. 1,1,1,1 to 1,1,1,1
    Test <- try(t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon)))
    if (!inherits(Test, "try-error"))
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=t.test(as.numeric(T2hRelExon),as.numeric(T0hRelExon))$p.value)
    } else
    {
      Result <- tibble(GeneID=Gene, Mean_T2h_RelExon = mean(as.numeric(T2hRelExon)),Mean_T0h_RelExon=mean(as.numeric(T0hRelExon)),ttest_pval=1)
    }
    Results %<>% add_row(Result)
  }
}
Results %<>% arrange(ttest_pval)
pvalBHadjust <- p.adjust(Results[["ttest_pval"]], method = "BH")
Results %<>% add_column(pvalBHadjust)
Results %<>% mutate("Difference_T2h-T0h" = Mean_T2h_RelExon - Mean_T0h_RelExon, .after = Mean_T0h_RelExon )
ResultsGeneNames <- Annotation %>% filter(GeneID%in% Results[["GeneID"]])
Results %<>% inner_join(ResultsGeneNames,by = "GeneID") %>% relocate(GeneName, .after = GeneID)
DiffZeroOther <- signif(t.test(Results[["Difference_T2h-T0h"]], mu = 0, alternative = "two.sided")$p.value,2)
PairedOther <- signif(t.test(Results[["Mean_T2h_RelExon"]],Results[["Mean_T0h_RelExon"]], alternative = "two.sided", paired =TRUE)$p.value,2)
MeanOther <- signif(mean(Results[["Difference_T2h-T0h"]]),2)
write_csv(Results, file=paste("./Outcsv/", Date,"_OtherIntronGenes_T2hvsT0h_med9.csv",sep=""))
MediatorOther <- signif(t.test(Results[["Difference_T2h-T0h"]],MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE)$p.value, 2)
boxplot(Results[["Difference_T2h-T0h"]],MediatorGenesDifferance, main=paste("SR T2h-T0h med9. Other vs. Mediator t.test pval", MediatorOther), names=c(paste("u",MeanOther , "Pair pval", PairedOther),paste("u", MeanMediator, "Pair pval", PairedMediator)))

#Test if genes are the same then paired t.test can be run
if(all(GeneIDWt %in% MediatorGenesID))
{
  #Reorder
  MediatorGenesDifferanceWtOrdered<-MediatorGenesDifferanceWt[match(MediatorGenesID,GeneIDWt)]
  #Paired
  Paired <- signif(t.test(MediatorGenesDifferanceWtOrdered,MediatorGenesDifferance, alternative = "two.sided", paired = TRUE)$p.value,2)
} else {
  Paired <- NaN
}
#Non paired
NonPaired <- signif(t.test(MediatorGenesDifferanceWt,MediatorGenesDifferance, alternative = "two.sided", var.equal = FALSE, paired = FALSE)$p.value,2)
#boxplot(MediatorGenesDifferanceWt,MediatorGenesDifferance, main=paste("SR T2h-T0h. Wt vs. med9 pval", NonPaired, "Pair pval", Paired), names=c("Mediator Wt","Mediator med9"))
boxplot(MediatorGenesDifferanceWt,MediatorGenesDifferance, main=paste("SR T2h-T0h. Wt vs. med9 pval", NonPaired, "Pair pval", Paired), names=c(paste("Mediator Wt u", MeanMediatorWt),paste("Mediator med9 u", MeanMediator)))
