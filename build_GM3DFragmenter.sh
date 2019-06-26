#!/bin/sh

#settings
locdir=`pwd`
javapath="/usr/bin"
src="$locdir/src"
jar="$locdir/jar"
lib="$locdir/lib"
CDKpath=$lib
CDKversion="1.4.19"

#welcome
echo "\n####################################################################"
echo "\n          Welcome to the script BUILDING GM3DFragmenter!"
echo "\n####################################################################\n"
echo " Please, check the default parameters before building:"
echo " -> path to JAVA tools: $javapath"
echo " -> path to source code: $src"
echo " -> path to jar storage: $jar"
echo " -> path to libraries: $lib"
echo " -> CDK version ID: $CDKversion"

#need to change JAVA
while true; do
    echo "\n Do you wish to change something? [Y/N]"
    read yn
    case $yn in
        [Yy]* )
		echo "\n Please, make a copy of this script and modify the settings as needed.\n"
		exit;;
	[Nn]* )
		break;;
        * ) echo " Please answer yes/y or no/n.";;
    esac
done

#clean previous
cd $jar
if [ -f "GM3DFragmenter.jar" ]; then
    rm GM3DFragmenter.jar
fi

# compile
$javapath/javac -cp $CDKpath/cdk-$CDKversion.jar -d . $src/*.java

# create manifest file
echo "Main-Class: GM3DMain" >> manifest.mf
echo "Class-Path: $CDKpath/cdk-$CDKversion.jar" >> manifest.mf
echo >> manifest.mf

    
$javapath/jar cvfm GM3DFragmenter.jar manifest.mf *.class  $CDKpath/cdk-$CDKversion.jar

if [ ! -f "GM3DFragmenter.jar" ]; then
    echo "####################################################################"
    echo "Cannot make GM3DFragmenter.jar."
    echo "####################################################################"
else 
    echo " "
    echo "####################################################################"
    echo " "
    echo "                Building of GM3DFragmenter completed!"
    echo " "
    echo "Start working with GM3DFragmenter using:"
    echo " $javapath/java -jar $jar/GM3DFragmenter.jar <input.par>"
    echo " "
    echo "               Thanks for using Foscato's tools! :-)"
    echo " "
    echo "####################################################################"
    echo " "
fi

rm *.class manifest.mf
