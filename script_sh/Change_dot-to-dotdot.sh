InFiles=$(ls | grep "07\|08\|09")
echo $InFiles

for File in $InFiles
do
sed -i 's/\.\/InOut/\.\.\/InOut/g' $File
done
