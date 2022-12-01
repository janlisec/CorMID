testthat::test_that(
  desc = "poss_local returns expected result",
  code = {
    out <- CorMID:::poss_local(vec=c(0.5,0.25,0.25), d=0.25)
    testthat::expect_true(is.matrix(out))
    testthat::expect_true(identical(dim(out), c(1L,3L)))

    out <- CorMID:::poss_local(vec=c(0.5,0.25,0.25), d=0.25, length.out=3)
    testthat::expect_true(identical(dim(out), c(7L,3L)))
    testthat::expect_true(all(apply(out,1,sum)==1))

    out <- CorMID:::poss_local(vec=c(0.5,0.25,0.25), d=0.05, by=0.01)
    testthat::expect_equal(nrow(out), 91L)

    limits <- matrix(c(0.5,0.51,0,1,0,1), nrow=2)
    out <- CorMID:::poss_local(vec=c(0.5,0.25,0.25), d=0.05, limits=limits, by=0.01)
    testthat::expect_equal(nrow(out), 21L)
  }
)
