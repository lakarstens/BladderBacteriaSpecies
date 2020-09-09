#! /bin/bash

# checks if the species from the KTW dataset are present in the local NCBI database

# this was a end of line nightmare!
# the file had Windows CR/LF endings, 
# and grep didn't find any matches until 
# I changed it to Unix LF endings.
species=$( find . -name "ktw_type_strains_2020-2-21.txt" )
#species=$( find . -name "ktw_type_strains_for_gg_2020-3-1.txt")
#species=$( find . -name "ktw_genus_2020-3-3.txt")
#species=$( find . -name "ktw_genus_for_gg_2020-3-3.txt")



file_base="test_genomic"

#database=$( find . -name "16SMicrobial_as_fasta.fasta" )
#database=$( find . -name "gg_13_5_taxonomy.txt" )
#database=$( find . -name "dump_ncbi_genomic.txt" )
#database=$( find . -name "SILVA_132_SSURef_Nr99_tax_silva.fasta" )
#database=$( find . -name "ref_prok_rep_genomes_ACC_taxonomy.txt" )
#database=$( find . -name "${file_base}.txt")
#database=$( find . -name "dump_custom_genomic_2020-3-5.txt" )
database=$( find . -name "complete_tkw_2020-3-5.taxonomy")

output="${file_base}_attendance.csv"

if [ $species ] && [ $database ]
    then
        echo -ne "found resource files\nsearching $database\n"
        touch $output
        while read -r line; do
            echo "looking for $line"
            # running grep to search the $database var for the $line var
            # using the exit status of this line as conditional
            chunk=$(grep -wFc "$line" $database)
            # if exit status of chunk is 1 (false), then grep didn't find
            # a match of $line in $database
            if [ $? -eq 1 ]
                then
                    echo -e "\t$line is not in the taxonomy file"
                    echo "$line,0" >> $output
                else 
                    # in a genomic database, this is counting the number of contigs 
                    # in the assembly of the genome, not the number of records
                    echo "$line,$chunk" >> $output
            fi 
        done < $species
else
    # -e enables interpretation of backslash escapes
    echo -e "\nI can't find all the needed resource files. Make sure the files\n\n\tktw_type_strains_2020-2-21.txt\n\t$file_base\n\nare in the same directory.\n"
fi
