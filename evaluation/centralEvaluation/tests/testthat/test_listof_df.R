context("List of dataframes")

###########
# blca
###########

blca_silva <- list_of_df("blca_silva/", "blca")

test_that("output is a list", {
    expect_type(blca_silva, "list")
})

test_that("list is 3 elements long", {
    expect_equal(length(blca_silva), 3)
})


# sorta tricky, df is an S3 class
test_that("element 1 is a dataframe", {
    expect_s3_class(blca_silva[[1]], "data.frame")
})

test_that("element 2 is a dataframe", {
    expect_s3_class(blca_silva[[2]], "data.frame")
})

test_that("element 3 is a dataframe", {
    expect_s3_class(blca_silva[[3]], "data.frame")
})

get_names <- names(blca_silva)
refer_names <- c("blca_cap_silva", "blca_v3_silva",  "blca_v6_silva" )

test_that("names are blca+identifier+silva", {
    expect_equal(get_names, refer_names)
})

get_df_names <- colnames(blca_silva[[1]])
df_refer <- c("id", "query", "blca", "confidence", "match")

test_that("colnames are the right ones", {
    expect_equal(get_df_names, df_refer)
})

###########
# nb
###########

nb_silva <- list_of_df("nb_silva/", "qiime")

test_that("output is a list", {
    expect_type(nb_silva, "list")
})

test_that("list is 3 elements long", {
    expect_equal(length(nb_silva), 3)
})

test_that("element 1 is a dataframe", {
    expect_s3_class(nb_silva[[1]], "data.frame")
})

test_that("element 2 is a dataframe", {
    expect_s3_class(nb_silva[[2]], "data.frame")
})

test_that("element 3 is a dataframe", {
    expect_s3_class(nb_silva[[3]], "data.frame")
})

get_names <- names(nb_silva)
refer_names <- c("qiime_cap_silva", "qiime_v3_silva",  "qiime_v6_silva" )

test_that("names are qiime+identifier+silva", {
    expect_equal(get_names, refer_names)
})

nb_get_df_names <- colnames(nb_silva[[1]])
nb_df_refer <- c("id", "query", "blca", "confidence", "match")

test_that("nb colnames are the right ones", {
    expect_equal(nb_get_df_names, nb_df_refer)
})
