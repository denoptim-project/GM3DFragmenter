#!/bin/sh

#settings
locdir=`pwd`
javapath="/usr/bin"
src="$locdir/src"
jar="$locdir/../../jar"
lib="$locdir/../../lib"
CDKpath=$lib
CDKversion="1.4.17"

cd $jar

if [ -f "ParallelizeGM3DF.jar" ]; then
    rm ParallelizeGM3DF.jar
fi

# compile
$javapath/javac -cp $CDKpath/cdk-$CDKversion.jar -d . $src/ParallelizeGM3DF.java

# create manifest file
echo "Main-Class: ParallelizeGM3DF" >> manifest.mf
echo "Class-Path: $CDKpath/cdk-$CDKversion.jar" >> manifest.mf
echo >> manifest.mf

$javapath/jar cvfm ParallelizeGM3DF.jar manifest.mf *.class  $CDKpath/cdk-$CDKversion.jar

if [ ! -f "ParallelizeGM3DF.jar" ]; then
    echo "###########################################"
    echo "Cannot make ParallelizeGM3DF.jar."
    echo "###########################################"
fi

rm *.class manifest.mf

###########################################################

if [ -f "SplitAndMerge.jar" ]; then
    rm SplitAndMerge.jar
fi

# compile
$javapath/javac -cp $CDKpath/cdk-$CDKversion.jar -d . $src/SplitAndMerge.java

# create manifest file
echo "Main-Class: SplitAndMerge" >> manifest.mf
echo "Class-Path: $CDKpath/cdk-$CDKversion.jar" >> manifest.mf
echo >> manifest.mf

$javapath/jar cvfm SplitAndMerge.jar manifest.mf *.class  $CDKpath/cdk-$CDKversion.jar

if [ ! -f "SplitAndMerge.jar" ]; then
    echo "###########################################"
    echo "Cannot make SplitAndMerge.jar."
    echo "###########################################"
fi

rm *.class manifest.mf

###########################################################


