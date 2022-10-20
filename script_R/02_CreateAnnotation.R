library('org.Sc.sgd.db')

#Set Date
Date <- gsub("-","",Sys.Date())

#Read GeneIDs
GeneIDs <- read.csv(file="../InFiles/20220900_GeneIDs.csv")

Annotation <- AnnotationDbi::select(org.Sc.sgd.db, GeneIDs[["GeneID"]], c("ORF","GENENAME"))
Annotation <- Annotation[,c(1,3)]
Annotation[is.na(Annotation[,"GENENAME"]),"GENENAME"]<-Annotation[is.na(Annotation[,"GENENAME"]),"ORF"]
colnames(Annotation) <-c("GeneID", "GeneName")
write.csv(Annotation, file=paste("../InOutFiles/", Date, "_Annotation.csv",sep=""), row.names = FALSE)
