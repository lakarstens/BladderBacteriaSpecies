#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
NULL

#---------------------------------------#
# see central_evaluation_plain.RMD
#
# and
#
#  confusion_matrix_functions_2020-3-7
#  new_taxonomy_results_2020-3-14
#---------------------------------------#

#' Helper function
#'
#' @description This little helper function is from
#' stackoverflow https://stackoverflow.com/questions/38351820/negation-of-in-in-r.
#' It was easier to define a not-in function than redo all those lists
#'
#' Function definition is based on this stackoverflow answer:
#' https://stackoverflow.com/questions/5831794/opposite-of-in-exclude-rows-with-values-specified-in-a-vector
#'
#' @param lhs The left hand side argument
#' @param rhs The right hand side argument
#'
#' @returns Elements of the `lhs` argument that are not in the `rhs` argument
#'
#' @export
`%notin%`<- function(lhs,rhs) {
  Negate(`%in%`)(lhs,rhs)
}

#' Make a list of dataframes
#'
#' @description I got tired of loading all the csv files individually,
#' so this function does it in batch. Reads all filenames in a directory,
#' then creates a dataframe from the `csv` files, then adds names to the growing list of dataframes.
#'
#' @param path_to_dir The path to the directory containing the output from the classification schemes.
#' @param prefix Which classifier was used, e.g. BLCA, Naive Bayes
#'
#' @returns A list of dataframes
#'
#' @export
list_of_df <- function(path_to_dir, prefix) {

  get_filenames <- list.files(path_to_dir, full.names = TRUE)
  build_csv_df <- lapply(get_filenames, readr::read_csv)

  grow_names <- vector()
  for (g in get_filenames) {
    new_name <- paste(prefix, strsplit(basename(g), "_")[[1]][3], strsplit(basename(g), "_")[[1]][5], sep="_")
    grow_names <- c(grow_names, new_name)
  }
  names(build_csv_df) <- grow_names
  return(build_csv_df)
}

#' Evaluate the performance of classification schemes
#'
#' @description A very long and convoluted function that I'll most
#' likely split into smaller chunks at a later time. Meant to be called from the function `total_results()`
#'
#' During development, I got a little lazy with the column names after the python script
#' is done reformatting classification results. The classifer column
#' is labled 'blca', whether or not the actual classifier used was BLCA
#' or Naive Bayes. I probably will not change that.
#'
#' This function also includes lists of bacteria species that are missing
#' from the various databases used in the study. Until I set up these
#' lists as imported data, they live here.
#'
#' @param results The results of the classification scheme
#' @param db The database used
#' @param vr The identifier used
#' @param clsfr The classifier used
#'
#' @returns A very long list of dataframes.
#'
#' Dataframe column dictionary
#' * true_match (numeric): Number of true matches (TM) in the dataset
#' * false_match (numeric): Number of false matches (FM) in the dataset
#' * cellC (numeric): The number of missed matches (MM)
#' * recall (numeric): Recall is calculated by TM/(MM + TM)
#' * precision (numeric): Precision is calculated by TM/(FM + TM)
#' * fmeasure (numeric): F1 is calculated by (weight * R) + ((1-weight) * P)
#' * weight (numeric): F1 can be weighted toward either Precision or Recall, but is equally weighted in this study
#' * weight_r (numeric): Intended to use, but wound up ignoring
#' * weight_p (numeric): Intended to use, but wound up ignoring
#' * database (character): Database used
#' * classifier (character): Classifier used
#' * region (character): Identifier used
#' * confidence (numeric): The number of times the same identification occured under random sampling (bootstrapping)

#'
#' @export
f1_records <- function(results, db, vr, clsfr) {

  # these are lists of the species missing from each reference database
  # in order to calculate the number of missed matches correctly
  gg_missing <- c('Actinotignum schaalii', 'Actinomyces naeslundii', 'Actinomyces neuii subsp. anitratus', 'Actinomyces odontolyticus', 'Actinomyces turicensis', 'Actinomyces urogenitalis DSM 15434', 'Aerococcus christensenii', 'Aerococcus sanguinicola', 'Aerococcus urinae', 'Alloscardovia omnicolens', 'Anaerococcus octavius', 'Bacillus idriensis', 'Bacillus infantis', 'Brevibacterium ravenspurgense', 'Corynebacterium amycolatum', 'Corynebacterium coyleae', 'Corynebacterium matruchotii ATCC 14266', 'Corynebacterium pyruviciproducens ATCC BAA-1742', 'Corynebacterium riegelii', 'Corynebacterium tuscaniense', 'Dermabacter hominis', 'Enterobacter asburiae', 'Enterobacter cloacae subsp. cloacae ATCC 13047', 'Enterococcus faecalis', 'Escherichia coli', 'Facklamia hominis CCUG 36813', 'Facklamia ignava', 'Gardnerella vaginalis ATCC 14018 = JCM 11026', 'Globicatella sanguinis', 'Gordonia terrae', 'Klebsiella pneumoniae', 'Kytococcus schroeteri', 'Lactobacillus crispatus', 'Lactobacillus fermentum', 'Lactobacillus gasseri', 'Lactobacillus jensenii', 'Lactobacillus johnsonii', 'Lactobacillus rhamnosus', 'Moraxella osloensis', 'Neisseria macacae', 'Neisseria perflava', 'Oligella urethralis', 'Propionibacterium avidum ATCC 25577', 'Proteus mirabilis', 'Pseudomonas aeruginosa', 'Staphylococcus hominis subsp. hominis', 'Staphylococcus saprophyticus subsp. saprophyticus ATCC 15305', 'Staphylococcus simulans', 'Staphylococcus warneri', 'Streptococcus equinus', 'Streptococcus gordonii', 'Streptococcus mitis', 'Streptococcus oralis ATCC 35037', 'Streptococcus parasanguinis', 'Streptococcus salivarius subsp. salivarius', 'Streptococcus sanguinis', 'Trueperella bernardiae', 'Varibaculum cambriense')

  custom_missing <- c('Bacillus idriensis', 'Corynebacterium species')

  ncbi_gen_missing <- c('Actinomyces naeslundii', 'Anaerococcus octavius', 'Bacillus idriensis', 'Corynebacterium amycolatum', 'Dermabacter hominis', 'Enterobacter asburiae', 'Globicatella sanguinis', 'Kytococcus schroeteri', 'Lactobacillus johnsonii', 'Neisseria macacae', 'Neisseria perflava', 'Neisseria subflava', 'Streptococcus oralis ATCC 35037', 'Corynebacterium species')


  # set up some variables
  ranger <- data.frame()

  # need to multiply the qiime confidence levels by 100
  if (clsfr=="blca") {
    results$confidence <- results$confidence / 100
  }

  for (x in seq(0,1,.001)) {
    #print(sprintf("working on %s", x))
    # seperate the confusion matrix cells first

    # True Matches:
    #   if the query and reference are the same name,
    #   and the confidence value is greater than the threshold
    #   if species isn't in DB, it never gets assigned as match
    cell_d <- ifelse((results$match==1 & results$confidence >= x), 1, 0)

    # False Matches:
    #   if the query and reference are not the same name,
    #   and the confidence value is greater than the threshold
    #   doesn't matter if the species is in the DB or not
    cell_b <- ifelse((results$match==0 & results$confidence >= x), 1, 0)

    # updating this
    # instead of creating a logical vector, loading the missing species
    # list into a variable that gets used when calculating cellC
    if (db=="gg") {
      # greengenes database
      # make a column showing if the query is in the database or not
      #results$attendence <- results$query %notin% gg_missing
      check_for_attendence <- gg_missing

    } else if (db=="genomic") {
      # ncbi genomic database

      # make a column showing if the query is in the database or not
      #results$attendence <- results$query %notin% ncbi_gen_missing
      check_for_attendence <- ncbi_gen_missing

    } else if (db=="custom") {
      # the custom genomic database

      # make a column showing if the query is in the database or not
      #results$attendence <- results$query %notin% custom_missing
      check_for_attendence <- custom_missing

    } else {
      #results$attendence <- TRUE
      check_for_attendence <- c("allgood")
    }

    # now for these two conditionals
    # checking if an element in the logical vector was failing and I didn't see it
    # the %notin% will evaluate to true or false

    # There are only two ways for a missed match
    #    1) a known true match is missed and another record is assigned as a match
    #    2) a true match that is below the confidence score, AKA a false match when the species is in the database

    # missed match below threshold:
    #   if the species is in the DB
    #   and match==1
    #   and below the confidence score
    missed_match_below <- ifelse((results$confidence < x & results$match==1 & (results$query %notin% check_for_attendence)), 1,0)

    # missed match above threshold:
    #   if the query is in the database
    #   and match==0
    #   and above the confidence score,
    #   then there needs to be a corresponding missed match
    #   but if it's below the conf score, it's only a MM
    missed_match_above <- ifelse((results$match==0 & (results$query %notin% check_for_attendence)), 1, 0)

    cell_c <- missed_match_below + missed_match_above

    # add cell counts to results dataframe
    results_confusion <- cbind(results, cell_b,cell_c,cell_d)

    w <- (sum(results_confusion$cell_d) + sum(results_confusion$cell_c)) / (sum(results_confusion$cell_b) + sum(results_confusion$cell_c) + 2*sum(results_confusion$cell_d))
    get_recall <- sum(results_confusion$cell_d)/(sum(results_confusion$cell_c) + sum(results_confusion$cell_d))
    get_precision <- sum(results_confusion$cell_d)/(sum(results_confusion$cell_b) + sum(results_confusion$cell_d))
    wr <- w*get_recall
    wp <- (1-w)*get_precision
    f1 <- (w*get_recall) + ((1-w)*get_precision)

    # replace 'region' names with more informative ones
    if (vr=="v3") {
      variable_region <- "V3"
    } else if (vr=="v6") {
      variable_region <- "V6"
    } else if (vr=="k17") {
      variable_region <- "V1-V3"
    } else if (vr=="bbv") {
      variable_region <- "V2-V3"
    } else if (vr=="cap") {
      variable_region <- "V4"
    } else if (vr=="k515") {
      variable_region <- "V4-V6"
    } else if (vr=="b646") {
      variable_region <- "V3-V5"
    } else if (vr=="ffh") {
      variable_region <- "Ffh V1-V2"
    } else if (vr=="rpob") {
      variable_region <- "RpoB V1"
    } else {
      variable_region <- "V3-V4"
    }

    # replace database
    if (db=="gg") {
      the_db <- "Greengenes"
    } else if (db=="silva") {
      the_db <- "Silva"
    } else if (db=="ncbi16") {
      the_db <- "NCBI 16S"
    } else if (db=="genomic") {
      the_db <- "NCBI Genomic"
    } else {
      the_db <- "Custom Genomic"
    }

    # replace classifier
    if (clsfr=="blca") {
      the_clsfr <- "BLCA"
    } else {
      the_clsfr <- "Naive Bayes"
    }


    small_range <- data.frame(weight=w, recall=get_recall, precision=get_precision, fmeasure=f1, weight_r=wr, weight_p=wp, true_match=sum(results_confusion$cell_d) , false_match=sum(results_confusion$cell_b), database=the_db, classifier=the_clsfr, region=variable_region, confidence=x, cellC=sum(results_confusion$cell_c))
    small_range$confidence <- x
    ranger <- rbind(ranger, small_range)
  }

  return(ranger)
}

#' total results
#'
#' @description Iteriavly feeds a list of dataframes into `f1_records()`
#'
#' @param each_result The output from `list_of_df()`
#'
#' @returns A dataframe of evaluated classification scheme results.
#'
#' Dataframe column dictionary
#' * true_match (numeric): Number of true matches (TM) in the dataset
#' * false_match (numeric): Number of false matches (FM) in the dataset
#' * cellC (numeric): The number of missed matches (MM)
#' * recall (numeric): Recall is calculated by TM/(MM + TM)
#' * precision (numeric): Precision is calculated by TM/(FM + TM)
#' * fmeasure (numeric): F1 is calculated by (weight * R) + ((1-weight) * P)
#' * weight (numeric): F1 can be weighted toward either Precision or Recall, but is equally weighted in this study
#' * weight_r (numeric): Intended to use, but wound up ignoring
#' * weight_p (numeric): Intended to use, but wound up ignoring
#' * database (character): Database used
#' * classifier (character): Classifier used
#' * region (character): Identifier used
#' * confidence (numeric): The number of times the same identification occured under random sampling (bootstrapping)
#'
#' @export
total_results <- function(each_result) {

  each_name <- names(each_result)

  final_db <- data.frame()
  for (x in each_name) {
    the_classifier <- strsplit(x, split="_")[[1]][1]
    the_var_region <- strsplit(x, split="_")[[1]][2]
    the_db <- strsplit(x, split="_")[[1]][3]
    #message(sprintf("sending to f1: name=%s, database=%s, vr=%s, classifier=%s", x, the_db, the_var_region, the_classifier))
    holder <- f1_records(each_result[x][[1]], the_db, the_var_region, the_classifier)
    final_db <- rbind(final_db, holder)
  }
  return(final_db)
}



