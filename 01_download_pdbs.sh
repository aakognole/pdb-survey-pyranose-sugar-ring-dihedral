#!/bin/bash

if [ $1 ]; then
    resname=`echo ${1} | awk '{print toupper($0)}'`
else
    echo "Usage: $ bash download_pdbs.sh <RESNAME>"
    exit
fi

# convert comma separated PDBs to seprated lines
sed -e "s~,~\n~g" list-of-pdbs-with-${resname}.txt > list.txt

# download PDBs from website
mkdir -p pdbs-with-${resname}
i=0
t=`wc -l list.txt | awk '{print $1}'`
cd pdbs-with-${resname}
echo "Downloading PDB files..."
for pdbid in `cat ../list.txt`; do
    if [ ! -e ${pdbid}.pdb ]; then
	wget https://files.rcsb.org/download/${pdbid}.pdb 2> /dev/null
    fi
    i=$((i+1))
    if [ $((i%10)) == 0 ]; then echo "$i / $t"; fi
done
cd ..
rm list.txt
echo "Finished!"

