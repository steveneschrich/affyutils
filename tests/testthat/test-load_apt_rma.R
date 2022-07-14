test_that("load_apt_works", {
  x<-load_apt_rma("../rma")
  expect_equal(Biobase::annotation(x), "Clariom_D_Human")
  expect_equal(dim(Biobase::exprs(x)),  c(19,1))
})
