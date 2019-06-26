#!/bin/sh

#settings
workdir=`pwd`
javapath="/usr/bin"
jar="$workdir/../../jar"
prefix="Fragments"
src="$workdir/../parallelize/src"

#check command line arguments
if test $# -lt 2
then
    echo "\nNo input file found!"
    echo "Provide the names of the libraries to be merged as command line arguments.\n"
    exit
fi

#make filenames compatible to SpliAndMerge
# creating a list of files with names like Fragments_input_0_1376546271084.sdf
i=0
for inplib in "$@"
do
    i=$(( $i + 1 ))
    echo "Getting file $inplib ($i)"
    newname=$prefix"_merge_"$i"_1000000000000.sdf"
    cp $inplib $prefix"_merge_"$i"_1000000000000.sdf"
done

#Splitting libraries of fragments according tomolecular weight
echo "Preparing files for splitting..."
$javapath/java -jar $jar/SplitAndMerge.jar split -p$workdir -jp$javapath/java -GM3DFp$jar/GM3DFragmenter.jar -queue187 -sw$prefix

#Run tasks
echo "Splitting fragment's libraries on MW..."
./RunAllSPLITTING_*
echo "Done!"
echo " "

#Merge MW-bins
echo "Preparing for mergning MW-bins..."
$javapath/java -jar $jar/SplitAndMerge.jar merge -p$workdir -jp$javapath/java -GM3DFp$jar/GM3DFragmenter.jar -queue187 -swmMWBin_

#Run tasks
echo "Merging MW-bins..."
./RunAllMERGING_*
echo "Done!"
echo " "

#Make one library
$src/makeSingleLib.sh
#sort -gk 3 all_BINs_Counting.dat > all_BINs_Counting_ordered.dat 
#gnuplot -persist $src/plotBins.gnuplot


