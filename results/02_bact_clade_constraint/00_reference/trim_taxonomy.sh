#!/bin/bash

cp tax_assign.txt tax_assign_trimmed.txt

echo "input:"
#grep -oP "\t[^;]+[^;]+" tax_assign.txt | sort | uniq -c
grep -oP "\t(Bacteria;Actinobacteria|Bacteria;Cyanobacteria|Bacteria;Proteobacteria|Bacteria;Firmicutes|Bacteria;Bacteroidetes)" tax_assign.txt | sort | uniq -c

for clade in "Archaea" "Eukaryota" "Bacteria;Actinobacteria" "Bacteria;Cyanobacteria" "Bacteria;Proteobacteria" "Bacteria;Firmicutes" "Bacteria;Bacteroidetes" ; do
	#echo ${clade}
	
	sed -i.bak "s/\t${clade};.*/\t${clade}/g" tax_assign_trimmed.txt
done

# put all other bacteria also in one clade, to separate them from our subgroups.
# first sed part skipts our clades, second one does the replaement.
# see https://stackoverflow.com/questions/9053100/sed-regex-and-substring-negation

sed -i.bak2 "/\t\(Bacteria;Actinobacteria\|Bacteria;Cyanobacteria\|Bacteria;Proteobacteria\|Bacteria;Firmicutes\|Bacteria;Bacteroidetes\)/b; s/\tBacteria;.*/\tBacteria;Others/g" tax_assign_trimmed.txt

echo

# count how often each clade is there
echo "results:"
grep -oP "\t.*" tax_assign_trimmed.txt | sort | uniq -c
