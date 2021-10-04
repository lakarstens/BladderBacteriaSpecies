context("Total results")

blca_silva <- list_of_df("blca_silva/", "blca")
nb_silva <- list_of_df("nb_silva/", "qiime")

all_results_list <- c(blca_silva, nb_silva)

all_results_evaluation <- total_results(all_results_list)

##########
# blca
##########

blca_results <- all_results_evaluation %>% dplyr::filter(confidence==0 & classifier=="BLCA" & region=="V4")

test_that("MM is equal to FM", {
    expect_equal(blca_results$false_match, blca_results$cellC)
})

test_that("recall is right", {
    expect_equal(blca_results$recall, (blca_results$true_match/(blca_results$cellC + blca_results$true_match)))
})

test_that("precision is right", {
    expect_equal(blca_results$precision, (blca_results$true_match/(blca_results$false_match + blca_results$true_match)))
})

##########
# nb
##########

nb_results <- all_results_evaluation %>% dplyr::filter(confidence==0 & classifier=="Naive Bayes" & region=="V4")

test_that("MM is equal to FM", {
    expect_equal(nb_results$false_match, nb_results$cellC)
})

test_that("recall is right", {
    expect_equal(nb_results$recall, (nb_results$true_match/(nb_results$cellC + nb_results$true_match)))
})

test_that("precision is right", {
    expect_equal(nb_results$precision, (nb_results$true_match/(nb_results$false_match + nb_results$true_match)))
})
