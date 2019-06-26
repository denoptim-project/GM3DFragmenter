#!/bin/sh
./build_GM3DFragmenter.sh
if [ $? -ne 0 ]; then
   echo "Compilation failed!"
   exit -1
fi
cd test
./runAllTests.sh
