call setbow
raxml -T %1 -s data/all_isolates_withRef.fa -m GTRCAT -p 12345 -f a -x 12345 -N 1000 -n global -w //wpia-hpc-hn/Fisher/global/output
