
# Species level resolution of female bladder microbiota from 16S rRNA amplicon sequencing


## In general

This study was done in two parts. The first part computationally
compares the phylogenetic resolution that is achieved by combining 3
currently available databases, 2 taxonomic classifiers, and 7 subsequences
from the 16S rRNA gene sequence into 42 classification schemes.

The second part of the study was to validate one of the classification
schemes by comparing the computational outcome with data generated from
targeted amplicon sequencing of bacterial DNA obtained from urine
samples.

## File description

    * `evaluation/` - code for calculating the recall, precision and F-measure of the classification schemes' results
    * `handle_newick/` - helper functions
    * `resources/` - resource files for the repository
    * `sliding_window/` - code for performing the sliding window analysis of the 16S rRNA gene sequence
    * `split_to_vr/` - code for extracting simulated amplicons of the 16S rRNA gene sequence
    * `taxonomy/` - code for handling the assigned taxonomy from the classification schemes
    * `testing/` - python and R tests

