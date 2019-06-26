#!/bin/sh

#settings
locdir=`pwd`
javapath="/usr/bin"
src="$locdir/src"
jar="$locdir/jar"
lib="$locdir/../../lib"
CDKpath=$lib
CDKversion="1.4.17"


#clean previous
cd $jar
if [ -f "ExtractMostCommonConformer.jar" ]; then
    rm ExtractMostCommonConformer.jar
fi


# compile
$javapath/javac -cp $CDKpath/cdk-$CDKversion.jar -d . $src/*.java

# create manifest file
echo "Main-Class: ExtractMostCommonConformer" >> manifest.mf
echo "Class-Path: $CDKpath/cdk-$CDKversion.jar" >> manifest.mf
echo >> manifest.mf

$javapath/jar cvfm ExtractMostCommonConformer.jar manifest.mf *.class  $CDKpath/cdk-$CDKversion.jar

if [ ! -f "ExtractMostCommonConformer.jar" ]; then
    echo "###########################################"
    echo "Cannot make ExtractMostCommonConformer.jar."
    echo "###########################################"
fi


rm *.class manifest.mf


cd $old
