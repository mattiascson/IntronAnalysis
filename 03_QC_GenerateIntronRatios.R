library(tidyverse)
library(magrittr)

#Set save QC files 
SaveQC <- "Yes" 

#Set save out files 
SaveOutFiles <- "Yes"

#Set Date
Date <- gsub("-","",Sys.Date())

#Read infiles
#Read annotation
Annotation <- read_csv(file=paste("./InOutFiles/20221001_Annotation.csv",sep=""))
#Read batch file. This file was not 100% correct manually corrected it.
PrepClust <- read_csv("./InFiles/SF-2245_clusters (003).csv") %>% select(Sample,Cluster)
#Intronfiles
SGDIntronGenes <- read_tsv("./InFiles/20220914_SGD_Introns.tsv", col_names = FALSE)
IntronGenesStd <- read_tsv("./InFiles/20220900_IntronGenes.csv")

#Set saved image to 1 in smaller than A4  	(8.27x11.69
png(file=paste("./Outpng/", Date,"_Spliceratio-ACT1.png",sep=""), width=10.69, height=7.27, units="in", res=100)

#Set to 2x2 combiplot
par(mfrow=c(2,2))

####No debatch
ExonFiles <- list.files("./InFiles/CountFiles","exon") %>% grep("summary",.,value = TRUE,invert = TRUE) %>% paste("./InFiles/CountFiles/",.,sep="")
GeneFiles <- list.files("./InFiles/CountFiles","gene") %>% grep("summary",.,value = TRUE,invert = TRUE) %>% paste("./InFiles/CountFiles/",.,sep="")
#Set no of samples
NOSamples <- length(GeneFiles)

#Number of genes 7127
ExonTibble <- tibble(.rows = 7127)
Samplenames <- vector(mode = "character")
loopstart <- 1
for (ExonFile in ExonFiles)
{
  if(loopstart == 1)
  {
    loopstart <- 0
    ExonTable<-read_tsv(ExonFile,skip=1)
    ExonTibble %<>% add_column(GeneID = ExonTable[["Geneid"]], Length = ExonTable[["Length"]])
  }
  ExonTable<-read_tsv(ExonFile,skip=1)
  Samplename<-substr(ExonFile,22,nchar(ExonFile)-15)
  Samplenames <- c(Samplenames,Samplename)
  ExonTable %<>% select(Geneid, contains("bam")) %>% rename(GeneID = Geneid,  !!Samplename := contains("bam"))
  #Extra security
  if(all(ExonTibble[["GeneID"]]==ExonTable[["GeneID"]]))
  {
    ExonTibble %<>% left_join(ExonTable)
  } else
  {
    Print("Break")
    break()
  }
}

GeneTibble <- tibble(.rows = 7127)
loopstart <- 1
for (GeneFile in GeneFiles)
{
  if(loopstart == 1)
  {
    loopstart <- 0
    GeneTable<-read_tsv(GeneFile,skip=1)
    GeneTibble %<>% add_column(GeneID = GeneTable[["Geneid"]], Length = GeneTable[["Length"]])
  }
  GeneTable<-read_tsv(GeneFile,skip=1)
  Samplename<-substr(GeneFile,22,nchar(GeneFile)-15)
  GeneTable %<>% select(Geneid, contains("bam")) %>% rename(GeneID = Geneid,  !!Samplename := contains("bam"))
  #Extra security
  if(all(GeneTibble[["GeneID"]]==GeneTable[["GeneID"]]))
  {
    GeneTibble %<>% left_join(GeneTable)
  } else
  {
    Print("Break")
    break()
  }
}

#Assign samples in rRNA batch. 1 low rRNA 197. 2 high rRNA 283
PrepClust %<>% slice(match(names(GeneTibble[,3:(NOSamples+2)]),PrepClust[["Sample"]]))

#Remove genes with any less than 50 counts
#seq 3,14 error?
ReadsP50 <- apply(as.matrix(ExonTibble[,seq(3,(NOSamples+2))]>50),1,all)&apply(as.matrix(GeneTibble[,seq(3,(NOSamples+2))]>50),1,all)
length(which(ReadsP50))
Plus50CountGenes <- ExonTibble %>% filter(ReadsP50) %>% select("GeneID") %>% left_join(Annotation)
if(SaveQC == "Yes")
{
write_csv(Plus50CountGenes, file=paste("./Outcsv/", Date,"_Plus50CountGenes.csv",sep=""))
}

#Ratio by intron counts
RatioTibble <- GeneTibble %>% select(GeneID)
IntronCountsTibble <- GeneTibble %>% select(GeneID)
for (Samplename in Samplenames)
{
  IntronPerGene <- GeneTibble[[Samplename]]-ExonTibble[[Samplename]]
  #Negative intron counts are not allowed set to 0
  IntronCountsTibble %<>% add_column(!!Samplename := IntronPerGene)
  IntronPerGene[IntronPerGene<0] <- 0
  IntronLengthPerGene <- GeneTibble[["Length"]]-ExonTibble[["Length"]]
  #Check if ever Gene is shorter than exons
  if(any(IntronLengthPerGene<0)) print("Less zero")  
  #Division fails if IntronLengthPerGene = 0
  IntronPerBase <- IntronPerGene/(IntronLengthPerGene*sum(GeneTibble[[Samplename]]))
  #Set to 0 if division fails
  IntronPerBase[is.nan(IntronPerBase)] <- 0
  ExonPerBase <- ExonTibble[[Samplename]]/(ExonTibble[["Length"]]*sum(ExonTibble[[Samplename]]))
  Ratio <- (ExonPerBase - IntronPerBase) / ExonPerBase
  #If ratio < 0 it means that by sampling introns got more counts. Set to 0
  Ratio[Ratio<0] <- 0
  RatioTibble %<>% add_column(!!Samplename := Ratio)
}
#Which genes has negative counts
LogiocIntronCounts <- IntronCountsTibble %>% select(!GeneID) %>% .[,]<(-2)
length(which(apply(LogiocIntronCounts,1,any)))
LogiocIntronCountsVec <- apply(LogiocIntronCounts,1,any)
NegIntronCountsTibble <- IntronCountsTibble  %>% left_join(Annotation) %>% relocate(GeneName, .after = GeneID) %>% filter(LogiocIntronCountsVec)
if(SaveQC == "Yes")
{
write_csv(NegIntronCountsTibble, file=paste("./Outcsv/", Date,"_NegCountsGenes.csv",sep=""))
}

#Does not keep order but not necessary
IntronGenes <- Annotation %>% filter(GeneName %in% IntronGenesStd[["GeneID"]]) %>% select(GeneID)
AllIntronGenes <- RatioTibble %>% add_column(LengthDiff = IntronLengthPerGene, .after = "GeneID") %>%  filter(IntronLengthPerGene  > 30)
#Filter away genbes not in SGD
AllIntronGenes %<>% filter(GeneID %in% SGDIntronGenes[["X1"]])
#Filter away genes less 50 counts
AllIntronGenes %<>% filter(GeneID %in% Plus50CountGenes[["GeneID"]])
#Filter away genes genes with less than -2 intron count.
AllIntronGenes %<>% filter(! GeneID %in% NegIntronCountsTibble[["GeneID"]])
#Filter away genes with 100% splicing in all samples
AllIntronGenes %<>% filter(! apply(AllIntronGenes[, 3:(NOSamples+2)] == 1,1,all))

NotInAllSet <- IntronGenes[["GeneID"]][!IntronGenes[["GeneID"]] %in% AllIntronGenes[["GeneID"]]]
#Check length missing genes
#Remove not in set
IntronGenes %<>% filter(! IntronGenes[["GeneID"]] %in% NotInAllSet)
OtherIntronGenes <- AllIntronGenes %>% filter(! AllIntronGenes[["GeneID"]] %in% IntronGenes[["GeneID"]])

AllIntronGenes %<>% left_join(Annotation) %>% relocate(GeneName, .after = "GeneID")
any(is.na(AllIntronGenes))
#any(is.nan(AllIntronGenes))
#any(is.infinite(AllIntronGenes))
OtherIntronGenes %<>% left_join(Annotation) %>% relocate(GeneName, .after = "GeneID")
MediatorGenes <- AllIntronGenes %>% filter(GeneID %in% IntronGenes[["GeneID"]] )

#Check Actin YFL039C
ACT1SpliceRatio <- RatioTibble %>% filter(GeneID == "YFL039C") %>% .[,2:(NOSamples+1)] %>% unlist(.)
ACT1Counts <- GeneTibble %>% filter(GeneID == "YFL039C") %>% .[,3:(NOSamples+2)]  %>% unlist(.)
hist(ACT1SpliceRatio, main = "YFL039C ACT1 ", xlab = "Spliceratio distribution 480 samples")
#Check Actin with rRNAbatch
plot(ACT1Counts,ACT1SpliceRatio, col = PrepClust[["Cluster"]], main = "ACT1 spliceratio and counts 480 samples", xlab = "Samples ACT1 Genecounts")
legend("bottomright", legend = c("Clust 1 low <3%", "Clust 2 high >19%"), text.col = c("black","red"), title = "By batch rRNA content")

if(SaveOutFiles == "Yes")
{
  write_csv(AllIntronGenes, file=paste("./InOutFiles/", Date,"_AllIntronGenes.csv",sep=""))
  write_csv(OtherIntronGenes, file=paste("./InOutFiles/", Date,"_OtherIntronGenes.csv",sep=""))
  write_csv(MediatorGenes, file=paste("./InOutFiles/", Date,"_MediatorGenes.csv",sep=""))
}

#### Debatch
ExonFiles <- list.files("./InFiles/CountFiles","exon") %>% grep("summary",.,value = TRUE,invert = TRUE) %>% paste("./InFiles/CountFiles/",.,sep="")
GeneFiles <- list.files("./InFiles/CountFiles","gene") %>% grep("summary",.,value = TRUE,invert = TRUE) %>% paste("./InFiles/CountFiles/",.,sep="")

#Number of genes 7127
ExonTibble <- tibble(.rows = 7127)
Samplenames <- vector(mode = "character")
loopstart <- 1
for (ExonFile in ExonFiles)
{
  if(loopstart == 1)
  {
    loopstart <- 0
    ExonTable<-read_tsv(ExonFile,skip=1)
    ExonTibble %<>% add_column(GeneID = ExonTable[["Geneid"]], Length = ExonTable[["Length"]])
  }
  ExonTable<-read_tsv(ExonFile,skip=1)
  Samplename<-substr(ExonFile,22,nchar(ExonFile)-15)
  Samplenames <- c(Samplenames,Samplename)
  ExonTable %<>% select(Geneid, contains("bam")) %>% rename(GeneID = Geneid,  !!Samplename := contains("bam"))
  #Extra security
  if(all(ExonTibble[["GeneID"]]==ExonTable[["GeneID"]]))
  {
    ExonTibble %<>% left_join(ExonTable)
  } else
  {
    Print("Break")
    break()
  }
}

GeneTibble <- tibble(.rows = 7127)
loopstart <- 1
for (GeneFile in GeneFiles)
{
  if(loopstart == 1)
  {
    loopstart <- 0
    GeneTable<-read_tsv(GeneFile,skip=1)
    GeneTibble %<>% add_column(GeneID = GeneTable[["Geneid"]], Length = GeneTable[["Length"]])
  }
  GeneTable<-read_tsv(GeneFile,skip=1)
  Samplename<-substr(GeneFile,22,nchar(GeneFile)-15)
  GeneTable %<>% select(Geneid, contains("bam")) %>% rename(GeneID = Geneid,  !!Samplename := contains("bam"))
  #Extra security
  if(all(GeneTibble[["GeneID"]]==GeneTable[["GeneID"]]))
  {
    GeneTibble %<>% left_join(GeneTable)
  } else
  {
    Print("Break")
    break()
  }
}

#Assign samples in rRNA batch. 1 low rRNA 197. 2 high rRNA 283
PrepClust %<>% slice(match(names(GeneTibble[,3:(NOSamples+2)]),PrepClust[["Sample"]]))

#Run Combat debatch without filtering any genes on ExonTiibble and GeneTibble
library("edgeR")
library("sva")
#Debatch ExonTibble
Temp <- as.data.frame(ExonTibble[,3:(NOSamples+2)])
rownames(Temp) <- ExonTibble[["GeneID"]]
Temp2<-ComBat_seq(as.matrix(Temp),as.vector(unlist(PrepClust[,"Cluster"])))
Temp2<-matrix(unlist(Temp2),nrow=7127,ncol=NOSamples)
rownames(Temp2) <- rownames(Temp)
colnames(Temp2) <- colnames(Temp)
ExonTibble <- tibble(ExonTibble[,1:2],as_tibble(Temp2))

#Debatch GeneTibble
Temp <- as.data.frame(GeneTibble[,3:(NOSamples+2)])
rownames(Temp) <- GeneTibble[["GeneID"]]
Temp2<-ComBat_seq(as.matrix(Temp),as.vector(unlist(PrepClust[,"Cluster"])))
Temp2<-matrix(unlist(Temp2),nrow=7127,ncol=NOSamples)
rownames(Temp2) <- rownames(Temp)
colnames(Temp2) <- colnames(Temp)
GeneTibble <- tibble(GeneTibble[,1:2],as_tibble(Temp2))


#Remove genes with any less than 50 counts
ReadsP50 <- apply(as.matrix(ExonTibble[,seq(3,(NOSamples+2))]>50),1,all)&apply(as.matrix(GeneTibble[,seq(3,(NOSamples+2))]>50),1,all)
length(which(ReadsP50))
Plus50CountGenes <- ExonTibble %>% filter(ReadsP50) %>% select("GeneID") %>% left_join(Annotation)
if(SaveQC == "Yes")
{
write_csv(Plus50CountGenes, file=paste("./Outcsv/", Date,"_Plus50CountGenesDebatch.csv",sep=""))
}

#Ratio by intron counts
RatioTibble <- GeneTibble %>% select(GeneID)
IntronCountsTibble <- GeneTibble %>% select(GeneID)
for (Samplename in Samplenames)
{
  IntronPerGene <- GeneTibble[[Samplename]]-ExonTibble[[Samplename]]
  #Negative intron counts are not allowed set to 0
  IntronCountsTibble %<>% add_column(!!Samplename := IntronPerGene)
  IntronPerGene[IntronPerGene<0] <- 0
  IntronLengthPerGene <- GeneTibble[["Length"]]-ExonTibble[["Length"]]
  #Check if ever Gene is shorter than exons
  if(any(IntronLengthPerGene<0)) print("Less zero")  
  #Division fails if IntronLengthPerGene = 0
  IntronPerBase <- IntronPerGene/(IntronLengthPerGene*sum(GeneTibble[[Samplename]]))
  #Set to 0 if division fails
  IntronPerBase[is.nan(IntronPerBase)] <- 0
  ExonPerBase <- ExonTibble[[Samplename]]/(ExonTibble[["Length"]]*sum(ExonTibble[[Samplename]]))
  Ratio <- (ExonPerBase - IntronPerBase) / ExonPerBase
  #If ratio < 0 it means that by sampling introns got more counts. Set to 0
  Ratio[Ratio<0] <- 0
  RatioTibble %<>% add_column(!!Samplename := Ratio)
}
#Which genes has negative counts
LogiocIntronCounts <- IntronCountsTibble %>% select(!GeneID) %>% .[,]<(-2)
length(which(apply(LogiocIntronCounts,1,any)))
LogiocIntronCountsVec <- apply(LogiocIntronCounts,1,any)
NegIntronCountsTibble <- IntronCountsTibble  %>% left_join(Annotation) %>% relocate(GeneName, .after = GeneID) %>% filter(LogiocIntronCountsVec)
if(SaveQC == "Yes")
{
write_csv(NegIntronCountsTibble, file=paste("./Outcsv/", Date,"_NegCountsGenesDebatch.csv",sep=""))
}

#Does not keep order but not necessary
IntronGenes <- Annotation %>% filter(GeneName %in% IntronGenesStd[["GeneID"]]) %>% select(GeneID)
AllIntronGenes <- RatioTibble %>% add_column(LengthDiff = IntronLengthPerGene, .after = "GeneID") %>%  filter(IntronLengthPerGene  > 30)
#Filter away genes not in SGD
AllIntronGenes %<>% filter(GeneID %in% SGDIntronGenes[["X1"]])
#Filter away genes less 50 counts
AllIntronGenes %<>% filter(GeneID %in% Plus50CountGenes[["GeneID"]])
#Filter away genes genes with less than -2 intron count.
AllIntronGenes %<>% filter(! GeneID %in% NegIntronCountsTibble[["GeneID"]])
#Filter away genes with 100% splicing in all samples
AllIntronGenes %<>% filter(! apply(AllIntronGenes[, 3:(NOSamples+2)] == 1,1,all))

NotInAllSet <- IntronGenes[["GeneID"]][!IntronGenes[["GeneID"]] %in% AllIntronGenes[["GeneID"]]]
#Check length missing genes
#Remove not in set
IntronGenes %<>% filter(! IntronGenes[["GeneID"]] %in% NotInAllSet)
OtherIntronGenes <- AllIntronGenes %>% filter(! AllIntronGenes[["GeneID"]] %in% IntronGenes[["GeneID"]])

AllIntronGenes %<>% left_join(Annotation) %>% relocate(GeneName, .after = "GeneID")
any(is.na(AllIntronGenes))
#any(is.nan(AllIntronGenes))
#any(is.infinite(AllIntronGenes))
OtherIntronGenes %<>% left_join(Annotation) %>% relocate(GeneName, .after = "GeneID")
MediatorGenes <- AllIntronGenes %>% filter(GeneID %in% IntronGenes[["GeneID"]] )

#Check Actin YFL039C
ACT1SpliceRatio <- RatioTibble %>% filter(GeneID == "YFL039C") %>% .[,2:(NOSamples+1)] %>% unlist(.)
ACT1Counts <- GeneTibble %>% filter(GeneID == "YFL039C") %>% .[,3:(NOSamples+2)]  %>% unlist(.)
hist(ACT1SpliceRatio, main = "YFL039C ACT1 ", xlab = "Spliceratio distribution 480 samples")
#Check Actin with rRNAbatch
plot(ACT1Counts,ACT1SpliceRatio, col = PrepClust[["Cluster"]], main = "ACT1 spliceratio and counts 480 samples", xlab = "Samples ACT1 Genecounts")
legend("bottomright", legend = c("Clust 1 low <3%", "Clust 2 high >19%"), text.col = c("black","red"), title = "By batch rRNA content")

if(SaveOutFiles == "Yes")
{
write_csv(AllIntronGenes, file=paste("./InOutFiles/", Date,"_AllIntronGenesDebatch.csv",sep=""))
write_csv(OtherIntronGenes, file=paste("./InOutFiles/", Date,"_OtherIntronGenesDebatch.csv",sep=""))
write_csv(MediatorGenes, file=paste("./InOutFiles/", Date,"_MediatorGenesDebatch.csv",sep=""))
}
dev.off()
