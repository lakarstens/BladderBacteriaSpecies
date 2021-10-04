context("testing F1")

#####################################
# mock classification scheme results
#####################################

results_one <- data.frame(match=c(1,1,1,0,0,0,0,0,0), confidence=c(.9,.6,.3,.9,.6,.3,.9,.6,.3), query=c("a", "b", "c", "d", "e", "f", "g", "h", "i"), blca=c("a", "b", "c", "e", "d", "e", "d", "e", "f"))

# results_one
#       match    confidence  query  blca
#---------------------------------------------
# 1     1        0.9         a      a
# 2     1        0.6         b      b
# 3     1        0.3         c      c
# 4     0        0.9         d      e
# 5     0        0.6         e      d
# 6     0        0.3         f      e
# 7     0        0.9         g      d
# 8     0        0.6         h      e
# 9     0        0.3         i      f

results_two <- data.frame(match=c(1,1,1,0,0,0,1,1,1), confidence=c(.9,.6,.3,.9,.6,.3,.9,.6,.3), query=c("a", "b", "c", "d", "e", "f", "g", "h", "i"), blca=c("a", "b", "c", "e", "d", "e", "g", "h", "i"))


# results_two
#       match    confidence  query  blca
#----------------------------------------------
# 1     1        0.9         a      a
# 2     1        0.6         b      b
# 3     1        0.3         c      c
# 4     0        0.9         d      e
# 5     0        0.6         e      d
# 6     0        0.3         f      e
# 7     1        0.9         g      g
# 8     1        0.6         h      h
# 9     1        0.3         i      i

# using Silva because it has a complete KTW species set
test_one_evaluation <- f1_records(results_one, "silva", "v3", "qiime")

# test_one_evaluation
# colname       value
#-------------------------
# weight        0.50
# recall        0.33
# precision     0.33
# fmeasure      0.33
# weight_r      0.17
# weight_p      0.17
# true_match    3
# false_match   6
# database      Silva
# classifier    Naive Bayes
# region        V3
# confidence    0.000
# cellC         6

##############################
# ignoring confidence score
##############################

test_conf_0 <- test_one_evaluation %>% dplyr::filter(confidence==0)

test_that("database is silva", {
    expect_equal(test_conf_0$database, "Silva")
})

test_that("classifier is NB", {
    expect_equal(test_conf_0$classifier, "Naive Bayes")
})

test_that("identifier is V3", {
    expect_equal(test_conf_0$region, "V3")
})

test_that("recall is .33", {
    expect_equal(test_conf_0$recall, 1/3)
})

test_that("precision is .33", {
    expect_equal(test_conf_0$precision, 1/3)
})

test_that("fmeasure is .33", {
    expect_equal(test_conf_0$fmeasure, 1/3)
})

test_that("NB confidence levels are < 1", {
    expect_lt(test_conf_0$confidence, 1)
})

##############################
# confidence == .3
##############################

test_two_evaluation <- f1_records(results_two, "silva", "v3", "qiime")

test_conf_33 <- test_two_evaluation %>% dplyr::filter(confidence == .3)

# test_conf_33
# colname       value
#-------------------------
# weight        0.5
# recall        0.67
# precision     0.67
# fmeasure      0.67
# weight_r      0.33
# weight_p      0.33
# true_match    6
# false_match   3
# database      Silva
# classifier    Naive Bayes
# region        V3
# confidence    .3
# cellC         3

test_that("database is silva", {
    expect_equal(test_conf_0$database, "Silva")
})

test_that("classifier is NB", {
    expect_equal(test_conf_0$classifier, "Naive Bayes")
})

test_that("identifier is V3", {
    expect_equal(test_conf_0$region, "V3")
})

test_that("recall is .66", {
    expect_equal(test_conf_33$recall, 2/3)
})

test_that("precision is .66", {
    expect_equal(test_conf_33$precision, 2/3)
})

test_that("fmeasure is .66", {
    expect_equal(test_conf_33$fmeasure, 2/3)
})

test_that("NB confidence levels are < 1", {
    expect_lt(test_conf_33$confidence, 1)
})
