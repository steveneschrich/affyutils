dilution <- readRDS("../dilution.rds")

test_that("is.ExpressionSet works", {
  expect_true(is.ExpressionSet(dilution))
  expect_false(is.ExpressionSet(iris))
})
