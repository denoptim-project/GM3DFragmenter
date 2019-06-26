#!/bin/sh

#settings
locdir=`pwd`
testdir="$locdir"
wrkdir="/tmp/GM3DFragmenter_test2"
javapath="/usr/bin"
jarpath="$locdir/../jar"
datadir="$locdir/../data"
repOnScreen=false

#welome
echo "\n####################################################################"
echo "             Welcome to the testing suite of GM3DFragmenter!\n"
echo "           This script tries to perform a series of tasks with "
echo "             GM3DFragmenter.  There is no validation of the"
echo "             results, so you have to evaluate them yourself.\n"
echo "\n####################################################################"

echo " Check the default paths:"
echo " -> testPath: $testdir (where tests files are taken from)"
echo " -> wrkdir: $wrkdir (where tests are actually run)"
echo " -> javaPath: $javapath"
echo " -> jarPath: $jarpath"

#ask for changes
while true; do
    echo "\n Do you wish to change one of these paths? [Y/N]"
    read yn
    case $yn in
        [Yy]* ) 
                echo " Specify a new test directory or press RETURN to skip"
                read nt
                if [ ! -z $nt ]
                then 
                        testdir=$nt
                        echo " The test directory is now: $testdir"
                fi
                echo " Specify a new work directory or press RETURN to skip"
                read nt
                if [ ! -z $nt ]
                then
                        wrkdir=$nt
                        echo " The work directory is now: $wrkdir"
                fi
                echo " Insert the new path to the java executable or press RETURN to skip"
                read nt
                if [ ! -z $nt ]
                then 
                        javapath=$nt
                        echo " The path to java is now: $javapath"
                fi
                echo " Insert the new path to the jar files or press RETURN to skip"
                read nt
                if [ ! -z $nt ]
                then
                        jar=$nt
                        echo " The path to the jar files is now: $jarpath"
                fi
                break;;
        [Nn]* ) 
                break;;
        * ) echo " Please answer yes/y or no/n.";;
    esac
done

#check existence
if [ ! -d "$testdir" ] 
then
	echo "ERROR! "$testdir" not found! Are you sure this is the right path?"
	exit
fi
if [ ! -d "$wrkdir" ]
then
	mkdir -p "$wrkdir"
	if [ $? -ne 0 ]
	then
		echo "ERROR! Cannot make '$wrkdir'!"
		exit
	fi
fi
if [ ! -d "$javapath" ] 
then
        echo "ERROR! $javapath not found! Are you sure this is the right path?"
        exit
fi
if [ ! -d "$jarpath" ] 
then
        echo "ERROR! $jarpath not found! Are you sure this is the right path?"
        exit
fi

#reporting mode
while true; do
    echo "\n Do you wish to run test in low-reporting mode? (write output in log files) [Y/N]"
    read yn
    case $yn in
        [Yy]* ) 
                echo " Low-reporting mode active. Messages from GM3DFragmenter will be redirected to log files (i.e. log_test*)."
                repOnScreen=false
                break;;
        [Nn]* )
                echo " High-reporting mode. All message from GM3DFragmenter will be printen on the screen, and log files are also created (i.e. log_test*)."
                repOnScreen=true
                break;;
        * ) echo " Please answer yes/y or no/n.";;
    esac
done


# Strat doing something
echo "\n####################################################################"
echo " Removing left overs from previous run..."

#run all available tests
cp -r "$testdir" "$datadir" "$wrkdir/"
cd "$wrkdir/test"
./cleanup.sh

echo "\n####################################################################"
echo "\n    RUNNING TEST: 1. Fragmentation \n ...Please wait...\n"
if ! "$repOnScreen" 
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation.par > log_test_fragmentation
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation.par 2>&1 | tee -a log_test_fragmentation
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 2. Fix Aromaticity \n ...Please wait...\n"
if ! "$repOnScreen" 
then    
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fix_aromaticity.par > log_test_fix_aromaticity
else
       "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fix_aromaticity.par 2>&1 | tee -a log_test_fix_aromaticity
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 3. Filters \n ...Please wait...\n"
if ! "$repOnScreen" 
then    
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_filter.par > log_test_filter
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_filter.par 2>&1 | tee -a log_test_filter
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 4. Extract by AP Class \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_extractClass.par > log_test_extractClass
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_extractClass.par 2>&1 | tee -a log_test_extractClass
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 5. Extract by SMARTS \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_extractSMARTS.par > log_test_extractSMARTS
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_extractSMARTS.par 2>&1 | tee -a log_test_extractSMARTS
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 6. Fragmentation collecting selected fragments \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_collectTargets.par > log_test_collectTargets
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_collectTargets.par 2>&1 | tee -a log_test_collectTargets
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 7. Fragmentation with unique labelling of APs \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation_2.par > log_test_fragmentation_2
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation_2.par 2>&1 | tee -a log_test_fragmentation_2
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 8. Fragmentation with library of frags to ignore \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation_3.par > log_test_fragmentation_3
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_fragmentation_3.par 2>&1 | tee -a log_test_fragmentation_3
fi

echo "\n####################################################################"
echo "\n    RUNNING TEST: 9. Fragmentation with RING and OMRING options \n ...Please wait...\n"
if ! "$repOnScreen"
then
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_omring.par > log_test_omring
else
        "$javapath"/java -jar "$jarpath"/GM3DFragmenter.jar test_omring.par 2>&1 | tee -a log_test_omring
fi

#goodbye
echo "\n####################################################################"
echo "\n TESTING DONE!"
echo "To delete all the files created by the test use ./cleanup.sh (from test folder)\n"
echo "Thanks for using Foscato's sripts and softwares. "
echo "Mandi! :) \n"


