InFiles=$(ls | grep "07")
echo $InFiles

for File in $InFiles
do
OutFile=$(echo $File | sed 's/07/09/' | sed 's/All/Oth/')
cat $File | sed -E 's/(All)([[:space:]])/Oth\2/g'\
| sed -E 's/AllIntronGenes[[:space:]]%>%/& filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>%/g'\
| sed -E 's/AllIntronGenes[[:space:]][[:space:]]%>%/& filter(GeneID %in% OtherIntronGenes[["GeneID"]]) %>%/g'\
> $OutFile
done
