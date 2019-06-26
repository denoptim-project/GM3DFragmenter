#!/bin/sh
#
#  Builds a draft of compatibility matrix from a library of fragments
#

#check command line arguments
if test $# -lt 1
then
    echo "\nNo input file found!"
    echo "Provide the name of the library of fragments from which the Compatibility Matrix should be created.\n"
    exit
fi

libName=$1
name=`basename $libName`
name=${libName%.*}
outFile="CPMap_$name.par"
cap="hyd"
scap=":1"
capb="me"
scapb=":1"

#extract classes
classes=($(grep -A 1 "CLASS" $libName | grep -v "CLASS" | grep -v "\-\-" | sed 's/[0-9]*\#/@/g' | sed 's/[ ,\,]/@/g' | tr '@' '\n' | sed 's/\:/ /g' | awk '{print $1}' | grep -v '^ *$' | grep -v '^$' | sort -u))
#write header
echo "# " > $outFile
echo "# CompatibilityMatrix for Class Based Builders" >> $outFile
echo "# " >> $outFile
# write compatibility matrix (with only complementary pairs)
for cl in ${classes[@]}; do echo "RCN $cl:0 $cl:1\nRCN $cl:1 $cl:0"  >> $outFile ; done
# write class-to-bond order section
echo "# CLASS-to-BondOrder conversion" >> $outFile
for cl in ${classes[@]}; do echo "RBO $cl 1" >> $outFile ; done
# write class-to-bond order for capping groups
#echo "RBO $cap 1" >> $outFile
#echo "RBO $capb 1" >> $outFile
# write gapping groups section
echo "# Capping" >> $outFile
for cl in ${classes[@]}; do echo "CAP $cl:0 $cap$scap\nCAP $cl:1 $cap$scap" >> $outFile  >> $outFile ; done

#rm tmp_make_CPMap*
echo "DONE! Check output: $outFile"


