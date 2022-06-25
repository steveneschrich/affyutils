test_exprs <- readRDS("../dilution.rds")
s <- annotate_with_symbol(test_exprs)
egfr<-s[which(s %in% "EGFR")]

test_that("multiplication works", {
  expect_equal(
    sort(Biobase::featureNames(subset_by(test_exprs, subset = names(egfr)))),
    sort(names(egfr))
  )
  expect_equal(
    sort(Biobase::featureNames(subset_by(test_exprs, subset = egfr, by = "SYMBOL"))),
    sort(names(egfr))
  )
  expect_equal(
    Biobase::exprs(subset_by(test_exprs, subset = egfr, by = "SYMBOL")),
    Biobase::exprs(test_exprs)[names(egfr),]
  )
  expect_equal(
    Biobase::exprs(subset_by_symbol(test_exprs, subset = egfr)),
    Biobase::exprs(test_exprs)[names(egfr),]
  )
})

