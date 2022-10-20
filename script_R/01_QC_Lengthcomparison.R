library(tidyverse)
library(magrittr)
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
    ExonTibble %<>% add_column(GeneID = ExonTable[["Geneid"]])
  }
  ExonTable<-read_tsv(ExonFile,skip=1)
  Samplename<-substr(ExonFile,14,nchar(ExonFile)-15)
  Samplenames <- c(Samplenames,Samplename)
  ExonTable %<>% select(Geneid, Length) %>% rename(GeneID = Geneid, !!Samplename := Length)
  ExonTibble %<>% left_join(ExonTable, by = "GeneID")
}
#Check if all lengths are equal
ExonTibble %>% select(!GeneID) %>% apply(.,1,function(x){length(unique(x)) == 1}) %>% all(.)


GeneTibble <- tibble(.rows = 7127)
loopstart <- 1
for (GeneFile in GeneFiles)
{
  if(loopstart == 1)
  {
    loopstart <- 0
    GeneTable<-read_tsv(GeneFile,skip=1)
    GeneTibble %<>% add_column(GeneID = GeneTable[["Geneid"]])
  }
  GeneTable<-read_tsv(GeneFile,skip=1)
  Samplename<-substr(GeneFile,14,nchar(GeneFile)-15)
  GeneTable %<>% select(Geneid, Length) %>% rename(GeneID = Geneid,  !!Samplename := Length)
  GeneTibble %<>% left_join(GeneTable, by = "GeneID")
}
#Check if all lengths are equal
GeneTibble %>% select(!GeneID) %>% apply(.,1,function(x){length(unique(x)) == 1}) %>% all(.)

#Check same order
all(ExonTibble[["GeneID"]]==GeneTibble[["GeneID"]])

# Is Gene ever shorter than Exon? No
which((GeneTibble[,2]<ExonTibble[,2])[,1])
