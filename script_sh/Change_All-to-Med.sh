InFiles=$(ls | grep "07")
echo $InFiles

for File in $InFiles
do
OutFile=$(echo $File | sed 's/07/08/' | sed 's/All/Med/')
cat $File | sed -E 's/(All)([[:space:]])/Med\2/g'\
| sed -E 's/AllIntronGenes[[:space:]]%>%/& filter(GeneID %in% MediatorGenes[["GeneID"]]) %>%/g'\
| sed -E 's/AllIntronGenes[[:space:]][[:space:]]%>%/& filter(GeneID %in% MediatorGenes[["GeneID"]]) %>%/g'\
> $OutFile
done
