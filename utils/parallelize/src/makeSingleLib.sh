#! /bin/bash
# Make a single file for all fragments
ls MWBin_* | grep "AllFrg" > all_MWBin_AllFrg
awk -F'_' '{print $2" " $0}' all_MWBin_AllFrg | sort -n | awk -F' ' '{ print $2}' | xargs cat > lib_allFrags.sdf

# Make a single file with all UNIQUE fragments
ls MWsorted-AllMWBin_InpMrgBin_* > all_MWsorterd-AllMWBin
awk -F'_' '{print $3" " $0}' all_MWsorterd-AllMWBin | sort -n | awk -F' ' '{ print $2}' | xargs cat > lib_unqFrags.sdf

# Get some statistics
echo "Total number of fragments in lib_allFrags.sdf:"
grep "\$\$\$\$" lib_allFrags.sdf | wc -l
echo " "
echo "Total number of unique fragments in lib_unqFrags.sdf:"
grep "\$\$\$\$" lib_unqFrags.sdf | wc -l





