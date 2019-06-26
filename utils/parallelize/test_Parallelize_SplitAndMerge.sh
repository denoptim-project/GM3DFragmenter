#!/bin/sh

#settings
locdir=`pwd`
javapath="/usr/bin"
src="$locdir/src"
jar="$locdir/../../jar"
data="$locdir/data"
testdir="$locdir/test"

cd $testdir

#clean up
rm -f LibFrgJob_* logFrgJob_* Run* InpFrgJob_* MolFrag-ratio_* Fragments_* CPMap_InpFrgJob_* MWBin_* IsomerCountsJob_* InpSpltJob_* logSpltJob_* MWsorted-AllMWBin_* logMrgBin_* InpMrgBin_* IsomerCounts-AllMWBin_* AllMWBin_* all_BINs_* lib_* all_MW* 

#Split the input library into smaller pieces and prepare the related tasks
$javapath/java -jar $jar/ParallelizeGM3DF.jar $data/test_fragmentation.sd -n4 -jp$javapath/java -GM3DFp$jar/GM3DFragmenter.jar -i$data/test_fragmentation.par -queue187 -p$testdir

#Run tasks
echo "Running GM3DFragmenter tasks..."
./RunAllGM3DFragmentation_*
echo "Done!"
echo " "

#Splitting libraries of fragments according tomolecular weight
echo "Preparing files for splitting..."
$javapath/java -jar $jar/SplitAndMerge.jar split -p$testdir -jp$javapath/java -GM3DFp$jar/GM3DFragmenter.jar -queue187 -swFragments_

#Run tasks
echo "Splitting fragment's libraries on MW..."
./RunAllSPLITTING_*
echo "Done!"
echo " "

#Merge MW-bins
echo "Preparing for mergning MW-bins..."
$javapath/java -jar $jar/SplitAndMerge.jar merge -p$testdir -jp$javapath/java -GM3DFp$jar/GM3DFragmenter.jar -queue187 -swmMWBin_

#Run tasks
echo "Merging MW-bins..."
./RunAllMERGING_*
echo "Done!"
echo " "

#Make one library
$src/makeSingleLib.sh
sort -gk 3 all_BINs_Counting.dat > all_BINs_Counting_ordered.dat 
gnuplot -persist $src/plotBins.gnuplot


