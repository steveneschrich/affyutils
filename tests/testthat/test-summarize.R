test_matrix<-matrix(1:12, nrow = 3)
test_exprs <- readRDS("../dilution.rds")
test_that(".summarize_mean", {
  expect_equal(.summarize_mean(1:2, test_matrix), c(1.5, 4.5, 7.5, 10.50))
  expect_equal(.summarize_mean(2, test_matrix), test_matrix[2,])
})

test_that(".summarize_median", {
  expect_equal(.summarize_median(1:3, test_matrix), c(2,5,8,11) )
  expect_equal(.summarize_median(1:3, test_matrix), .summarize_mean(1:3,test_matrix))
  expect_equal(.summarize_median(2, test_matrix), test_matrix[2,])
})

test_that(".summarize_maxmedian", {
  expect_equal(.summarize_maxmedian(1:3, test_matrix), test_matrix[3,])
  expect_equal(.summarize_maxmedian(2, test_matrix), test_matrix[2,])
})

test_that(".summarize_minmedian", {
  expect_equal(.summarize_minmedian(1:3, test_matrix), test_matrix[1,])
  expect_equal(.summarize_minmedian(2, test_matrix), test_matrix[2,])
})

test_that(".summarize_tukey_biweight", {
  expect_equal(.summarize_tukey_biweight(1:3, test_matrix), test_matrix[2,])
  expect_equal(.summarize_tukey_biweight(1:4, rbind(test_matrix, test_matrix[1,])),
                                         c(1.5484027,4.5484027,7.5484027,10.5484027) )
})

s <- annotate_with_symbol(test_exprs)
egfr <- s[which(s %in% "EGFR")]

test_that("summarize_by", {
  expect_equal(
    Biobase::exprs(summarize_by(test_exprs, by = "SYMBOL", method="median"))["EGFR",],
    apply(Biobase::exprs(test_exprs)[names(egfr),],2,median)
  )
  expect_equal(
    Biobase::exprs(summarize_by_symbol(test_exprs, method="median"))["EGFR",],
    apply(Biobase::exprs(test_exprs)[names(egfr),],2,median)
  )
  expect_equal(
    Biobase::exprs(summarize_by(test_exprs, by = "SYMBOL", method="mean"))["EGFR",],
    apply(Biobase::exprs(test_exprs)[names(egfr),],2,mean)
  )
  expect_equal(
    Biobase::exprs(summarize_by(test_exprs, by = "SYMBOL", method="maxmedian"))["EGFR",],
    Biobase::exprs(test_exprs)["1537_at",]
  )
  expect_equal(
    Biobase::exprs(summarize_by(test_exprs, by = "SYMBOL", method="minmedian"))["EGFR",],
    Biobase::exprs(test_exprs)["37327_at",]
  )
})
