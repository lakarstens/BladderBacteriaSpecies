# Background

This package is intended to be used after taxonomy has been assigned by the classification schemes and synonyms have been replaced by the python script `validate_match_batch.py`. The output of that script is a CSV file formatted in the following manner:

| id | query | blca | confidence | match |
|---------|---------|---------|---------|---------|
|AB618790.1.1522 | Actinomyces naeslundii | Actinomyces naeslundii | 0.50332 | 1|
|AB920570.1.1436 | Gordonia terrae | Gordonia terrae | 0.20453 | 1|
|ACFH01000038.241.1767 | Actinomyces urogenitalis DSM 15434 | Actinomyces urogenitalis | 0.81324 | 1|
|AF003929.1.1520 | Streptococcus mitis | Streptococcus BS35a | 0.11467 | 0|
|AF133538.1.1486 | Oligella urethralis | Oligella urethralis | 0.86878 | 1|

The data dictionary for this type of file is:

|column name | description |
|-------|-----------------------|
|id | the acession number |
|query | the genus and species of the sample as identified by whole genome sequencing |
|blca | the taxonomy assignment of the classification scheme |
|confidence | the confidence score of the classification scheme |
|match | binary indication of whether the "query" and "blca" names match, disregarding strain information. 1 for match, 0 otherwise |

One item of importance is that the column name "blca" is the result of an entrenched piece of code. During development, one column was needed to hold the known genus and species name of the sample that was fed into the classification scheme (the "query" column), and one column to hold the result after the classification scheme had assigned taxonomy to the sample (the "blca" column). "blca" was used because it happened to be the first classifier that was used in the study. The pathname of the file output by `validate_match_batch.py` holds all the information about the classification scheme used, but this suffered from early development choices as well. For example, the file "validated_qiime_silva_outfiles_2020-01-05_04/validated_formatted_cap_tpstr_full_silva_classified_2020_2020-01-05_04.csv" is the very wordy way that the output of the classification scheme composed of the Silva database, 16S V4 identifier, and Naive Bayes classifier implemented by Qiime has been reformatted for the study pipeline and checked for synonyms.

The `centralEvaluation` package uses the file pathname to extract information about the classification scheme, and the information held in the columns to build the final evaluated dataframe. This dataframe was then used to make the graphs and results for the manuscript. 

# Order of functions

## list_of_df()

Once the classification scheme output files have been checked for synonyms and collected into a directory, the `list_of_df()` function reads those files and converts them into a list of dataframes. 

```{r, eval=FALSE}
blca_gg_results <- list_of_df("path/to/directory/of/input/files/", "blca")
qiime_gg_results <- list_of_df("path/to/directory/of/input/files/", "qiime")
```

Each of these lists were collected into another list by hand, and was then ready for evaluation.

```{r, eval=FALSE}
all_species_16s_results <- c(blca_gg_results, blca_silva_results, blca_ncbi16_results, qiime_gg_results, qiime_silva_results, qiime_ncbi16_results)
```

## total_results()

This function is a nested `for` loop that calls `f1_records()`.

```{r, eval=FALSE}
all_species_16s_evaluation <- total_results(all_species_16s_results)
```

### f1_records()

This function was not intended to be used alone by itself, because the output dataframe was intended to hold the evaluations of all the classification schemes used in the study. It is also a bit too long, and in the future should be split into smaller functions that do more specific tasks. 

The Recall, Precision and F-measure were calculated based on these definitions.

  * True match - All record pairs assigned as a match that have identical genus and species labels.
  * False match - All record pairs assigned as a match that did not have identical genus and species labels.
  * False non-match - If a record representing a species in the Thomas-White dataset was present in the database, but was not assigned as a match, the record was evaluated as a false non-match. 
  * True non-match - All records in the reference database that were not in the Thomas-White dataset. While records assigned to this category were not used in evaluating the classification schemes in this manuscript, the definition is still included for completeness.


Output dataframe dictionary:

| colname | description |
|---------|-------------------------|
| true_match | (numeric) Number of true matches (TM) in the dataset |
| false_match | (numeric) Number of false matches (FM) in the dataset |
| cellC | (numeric) The number of missed matches (MM) |
| recall | (numeric) Recall is calculated by TM/(MM + TM) |
| precision | (numeric) Precision is calculated by TM/(FM + TM) |
| fmeasure | (numeric) F1 is calculated by (weight * R) + ((1-weight) * P) |
| weight | (numeric) F1 can be weighted toward either Precision or Recall, but is equally weighted in this study |
| weight_r | (numeric) Intended to use, but wound up ignoring |
| weight_p | (numeric) Intended to use, but wound up ignoring |
| database | (character) Database used |
| classifier | (character) Classifier used |
| region | (character) Identifier used |
| confidence | (numeric) The number of times the same identification occured under random sampling (bootstrapping) |

### %notin%

Inside the `f1_records()` function is a section that deals with long lists of species names. It was easier to make a helper function that would exclude a few names, instead of using the `%in%` magrittr-style function and having to type a slightly-less-long list of names.
