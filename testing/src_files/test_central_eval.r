# I also need to write some tests so I dont' lose my mind.

# Want two databases, one with all query records and one 
# that is missing 3. 3 levels of confidence. Correct and incorrect matches.

# Setup one has a database that's missing 3 query records, 
# setup two has all query records available.

db_one <- c("g", "h", "i")
db_two <- c("z")

test_results_one <- data.frame(match=c(1,1,1,0,0,0,0,0,0), confidence=c(.9,.6,.3,.9,.6,.3,.9,.6,.3), query=c("a", "b", "c", "d", "e", "f", "g", "h", "i"), blca=c("a", "b", "c", "e", "d", "e", "d", "e", "f"))

test_results_two <- data.frame(match=c(1,1,1,0,0,0,1,1,1), confidence=c(.9,.6,.3,.9,.6,.3,.9,.6,.3), query=c("a", "b", "c", "d", "e", "f", "g", "h", "i"), blca=c("a", "b", "c", "e", "d", "e", "g", "h", "i"))

test_that("test f1_records with results one",{

    test_one_evaluation <- f1_records(test_results_one, "silva", "v3", "qiime")
  
    test_conf_0 <- test_one_evaluation %>% filter(confidence==0)
    
    expect_equal(test_conf_0$recall, 1/3)
    expect_equal(test_conf_0$precision, 1/3)
    expect_equal(test_conf_0$fmeasure, 1/3)
})

test_that("test f1_records with results two",{

    test_two_evaluation <- f1_records(test_results_two, "silva", "v3", "qiime")
  
    test_conf_0 <- test_two_evaluation %>% filter(confidence == .3)
    
    expect_equal(test_conf_0$recall, 2/3)
    expect_equal(test_conf_0$precision, 2/3)
    expect_equal(test_conf_0$fmeasure, 2/3)
})