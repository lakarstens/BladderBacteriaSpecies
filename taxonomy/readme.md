





# Assigning taxonomy

blah

## primer abbreviations

There are 78 sequences generated for each predicted amplicon.
Abbreviations for the predicted amplicons were used in the filenames for
convienience, and are listed in the table below.

| region | primer                   | abbreviation |
| ------ | ------------------------ | ------------ |
| V1-V3  | A17F-515R                | k17          |
| V4-V6  | 515F-1114R               | k515         |
| V2-V3  | 16S\_BV2f-16S\_BV3r      | bbv          |
| V3-V5  | MiCSQ\_343FL-MiCSQ\_806R | b646         |
| V4     | F515-R806                | cap          |
| V3-V4  | V3F-V4R                  | gras         |
| V6     | v6\_1183F-v6\_1410R      | v6           |
| V3     | v3\_579F-v3\_779R        | v3           |

``` r
amp_plot / swa_plot + plot_annotation(tag_levels = 'A')
```

![](taxonomy_readme_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Classification

### BLCA

All predicted amplicon data was classified in batch by using a bash
script and run in a Linux environment. In general, BLCA was called in
this manner when using the default NCBI 16S
    database:

    python [path to BLCA script] -i [path to input FASTA file] -o [new name to save the output file under]

For the remaining databases, the path to the database and taxonomy file
had to be
    included:

    python [path to BLCA script] -i [path to input FASTA file] --tax [path to taxonomy file] --db [path to the database] -o [new name to save the output file under]

Just as an example, part of the `bash` file is reproduced below. You’ll
have to change your bash file to reflect where your resource files are
located.

``` bash
#!/usr/bin/env bash

# batch script for KTW type strains

# the BLCA default 16S database
python 2.blca_main.py -i extracted_16s_regions_2019-11-26/extracted_cap_tpstr_2019-11-26_1550.fna -o cap_tpstr_default_2019-11-26.outfile
python 2.blca_main.py -i extracted_16s_regions_2019-11-26/extracted_gras_tpstr_2019-11-26_1551.fna -o gras_tpstr_default_2019-11-26.outfile
#   remaining 6 predicted amplicon files would be continued here 

# the greengenes database
python 2.blca_main.py -i extracted_16s_regions_2019-11-26/extracted_cap_tpstr_2019-11-26_1550.fna --tax gg/gg_13_5_taxonomy.taxonomy --db gg/gg_13_5 -o cap_tpstr_gg_2019-11-26.outfile
python 2.blca_main.py -i extracted_16s_regions_2019-11-26/extracted_gras_tpstr_2019-11-26_1551.fna --tax gg/gg_13_5_taxonomy.taxonomy --db gg/gg_13_5 -o gras_tpstr_gg_2019-11-26.outfile
#   remaining 6 predicted amplicon files would be continued here 
```

#### reformat results

The evaluation steps need the classification results to be in a
different format than the standard BLCA output. We used the script
`blcaout_to_csv.py` to do this reformatting, where
`../raw_data/blca_tpstr_2019-11-26` is the path to the folder of
classified
    results.

    (base) C:\Users\Carter\Documents\thesis\thesis_git\taxonomy\src_files>python blcaout_to_csv.py -i ../raw_data/blca_tpstr_2019-11-26

### Qiime

Qiime needed some extra formatting and importing steps before the
predicted amplicons could be classified. The classification was also
done in a Linux environment, and all work was done in a folder that
included the data, compiled classifier, and bash scripts.

#### data prep

All the FASTA files had to be imported as a Qiime object before the
classifier would work. We used this bash script to import all the files
in a directory. One important note is that the sequences must be DNA. If
you’re using RNA sequences, replace the Uracil (U) nucleotides with
Thymine (T) (we just used the search and replace function in a plain
text program, like Window’s Notebook).

We stuck to a filename convention when dealing with all the predicted
amplicons, for example `extracted_b646_tpstr_2019-11-26_1552_dna.fna`.
The bash script will read in all file in a directory expecting this
convention, and name the output file based on the second and third
positions in the snake-case filename. So the output would be
`b646_tpstr_sequence.qza`.

``` bash
#!/bin/bash

# get the files from command line arg
seq_files=($1*.fna)

for x in ${seq_files[@]}
do
    # backticks return the commands as a variable
    # cut splits a string, returns the 2&3 element
    get_basename=`basename "$x" | cut -d'_' -f2,3`
    
    qiime tools import --type 'FeatureData[Sequence]' --input-path $x --output-path $1$get_basename"_sequence"
done
```

The bash script was called like
this:

``` shellsession
(qiime2-2019.10):~/qiime_classification_2019-12-18$ ./batch_qiime_import.bash extract_vr_dna_2019-12-18_50/
```

#### classify

The Naive Bayes classifier requires a trained database. “Trained” in
this context means a FASTA multirecord file of sequences that have been
broken into a bazillion 8-mers, and the number of specific 8-mers found
in each taxon enumerated and a probabiity calculated. “Training” the
Silva database took 4 days on the OHSU Advanced Computing Cluster. The
example bash script that performs the classification of all predicted
amplicons using the NCBI 16S database is below. When using any other
database, change the `--i-classifier` option to reflect the path to the
trained database.

``` bash
#!/bin/bash

# get the files from command line arg
seq_files=($1*.qza)

for x in ${seq_files[@]}
do
    # backticks return the commands as a variable
    # cut splits a string, returns the 2&3 element
    get_basename=`basename "$x" | cut -d'_' -f1,2`
    qiime feature-classifier classify-sklearn --i-classifier ncbi16s_classifier_2019-12-17.qza --i-reads $x --o-classification $1$get_basename"_classified" --p-confidence 0
done
```

#### extract qza files

Qiime’s .qza files are simply .zip files. They have several extra files
for metadata, but the buisiness end of it is the file `taxonomy.tsv`.
Again, using a bash script is slightly more convienient than renaming
the .qza to .zip and then going through the unzipping process by hand.

``` bash
#!/bin/bash

# get the files from command line arg
seq_files=($1*_classified.qza)

for x in ${seq_files[@]}
do
    qiime tools extract --input-path $x --output-path $2
done
```

The bash script was called like
this:

``` shellsession
(qiime2-2019.10):~/qiime_classification_2019-12-18$ ./batch_qiime_extract.bash extract_vr_dna_2019-12-18_50/ final_data/
```

#### reformat results

Like BLCA, the evaluation step needs the classification output to be in
a specific format. The `blcaout_to_csv.py` script was modified to do
this with the Qiime results, and is called `qiimeout_to_csv.py`.
`../raw_data/qiime_outfiles_2019-12-18/` is the path to the output Qiime
results.

    (base) >qiimeout_to_csv.py -i ../raw_data/qiime_outfiles_2019-12-18/

## Synonyms

Finally, the reformatted results had to be checked for bacterial species
synonyms.

### Prokaryotic Nomenclature Up-to-Date

Prokaryotic synonyms were downloaded from the [Prokaryotic Nomenclature
Up-to-Date
website](https://www.dsmz.de/services/online-tools/prokaryotic-nomenclature-up-to-date/downloads)
run by the DSMZ-German Collection of Microorganisms and Cell Cultures,
GmbH. Strains are out of scope, so entries like *Enterobacter cloacae
cloacae* and *Enterobacter cloacae dissolvens* are considered synonyms
of *Enterobacter cloacae*. Misspellings that don’t map to known synonyms
are assigned as a mismatch, but unless you’re typing your own taxonomy
into FASTA records, the chance of misspellings are low.

### in general

Any genus or species taxon labled as “unidentified” or “uncultured” is
designated as “unavailable”. “Unclassified” is reserved for when the
classification algo thinks there’s no match for the record.
“Unavailable” just means there’s no information for that record, and
are assigned as true non-matches.

Input is the formatted *csv* file as written by the python scripts
`blcaout_to_csv.py` and `qiimeout_to_csv.py`, for example:

    id,taxon,rank,confidence
    CP011325.22226.23780,Bacteria,domain,100.0

the `--map_file` option (`-m`) takes a path to a file that can be used
to map accession numbers to taxonomy by
`newick_tools.inhouse_to_species.py`, found in the `resources`
directory.

  - master\_wellcome\_list\_2019-05-07.csv

  - type\_strains\_map\_2019-11-26.csv

<!-- end list -->

    >python validate_match.py -i path/to/formatted/assigned_taxonomy_file.csv -m path/to/map/file

The screen output is a little verbose, but unless the script finds a
discrepency it will only print the output filename. An example of a
successfully resolving a synonym
is:

``` shellsession
        the result name 'Streptococcus thermophilus' was found to be a synonym linked to the current name 'Streptococcus salivarius'
        comparing the mapped name 'Streptococcus salivarius' to the current name 'Streptococcus salivarius'
        and found to be the same, writing [AY188352.1.1546,Streptococcus salivarius subsp. salivarius,Streptococcus salivarius,33.2333333333,1] to file

writing to ../processed_files/validated_outfiles_2020-01-05_47/validated_formatted_k515_tpstr_default_2019_2020-01-05_47.csv
```

and one that did
not:

``` shellsession
        the result name 'Streptococcus thermophilus' was found to be a synonym linked to the current name 'Streptococcus salivarius'
        comparing the mapped name 'Streptococcus anginosus' to the current name 'Streptococcus salivarius'
        and found to be different, writing [AFIM01000033.206.1753,Streptococcus anginosus SK52 = DSM 20563,Streptococcus salivarius,53.5833333333,0] to file
```
