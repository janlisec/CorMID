testthat::test_that(
  desc = "weight_errors returns expected result",
  code = {
    rMpH <- seq(0,1,0.1)
    errs <- rep(1,11)
    out <- CorMID:::weight_errors(rMpH=rMpH, errs=errs, penalize=6)
    testthat::expect_equal(length(out), 11L)
    testthat::expect_true(is.numeric(out))
    testthat::expect_equal(out[1], 10L)
    testthat::expect_equal(out[11], 1L)
    #'plot(1,1,xlim=c(0,1), ylim=c(0,11), type="n", xlab="rMpH", ylab="Weighted Error")
    #'for (i in 1:10) {
    #'lines(x=rMpH, y=CorMID:::weight_errors(rMpH=rMpH, errs=errs, penalize=i), col=i)
    #'}
    #'legend(x="topright", legend=1:10, fill=1:10, title="penalize")
  }
)
