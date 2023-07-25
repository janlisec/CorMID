# simple calculations
fml <- "C9H20O3Si2"
mid <- c(0.9,0,0,0.1)
r <- list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1)
rMID <- CorMID::recMID(mid=mid, r=r, fml=fml)
out <- CorMID::CorMID(int=rMID, fml=fml)
# avoid creating a Rplots.pdf in testthat folder
pdf(NULL)

testthat::test_that(
  desc = "CorMID returns expected result",
  code = {
    testthat::expect_true(inherits(out, "CorMID"))
    testthat::expect_true(is.numeric(out))
    testthat::expect_false(is.null(names(out)))
    testthat::expect_equal(unname(out[1:3]), c(89.06250,0.00000, 0.78125))
  }
)

testthat::test_that(
  desc = "CorMID works with fixed MID",
  code = {
    testthat::expect_equal(attr(CorMID::CorMID(int=rMID, fml=fml, mid_fix = mid),"mid_status"), "fixed")
  }
)

testthat::test_that(
  desc = "CorMID works with character vector specifying 'r'",
  code = {
    testthat::expect_equal(attr(CorMID::CorMID(int=rMID, fml=fml, r = names(r)),"ratio_status"), "estimated")
  }
)

testthat::test_that(
  desc = "CorMID works with numeric vector specifying 'r'",
  code = {
    testthat::expect_equal(attr(CorMID::CorMID(int=rMID, fml=fml, r = unlist(r)),"ratio_status"), "fixed")
  }
)

testthat::test_that(
  desc = "CorMID removes infinite intensities",
  code = {
    rMID2 <- rMID
    rMID2["M+3"] <- NA
    # removing the information on M+3 will lead to M3 to be estimated as 0
    testthat::expect_equal(unname(CorMID::CorMID(int=rMID2, fml=fml)[4]), 0)
  }
)

testthat::test_that(
  desc = "CorMID plot returns expected result",
  code = {
    vdiffr::expect_doppelganger(
      title = "CorMID_Plot",
      fig = function() plot(out)
    )
  }
)

testthat::test_that(
  desc = "CorMID plot returns expected result with alternative options",
  code = {
    vdiffr::expect_doppelganger(
      title = "CorMID_Plot_alt",
      fig = function() plot(unname(out), xlab="test", ylab="test", lwd=9, lend=2, xlim=c(0,10), ylim=c(0,2), las=1, col=2)
    )
  }
)

testthat::test_that(
  desc = "CorMID print returns expected result",
  code = {
    #browser()
    #testthat::expect_output(print(out), "\\033\[0;34mMID \[%\] \(estimated\)\\033\[0m\\n    M0    M1    M2    M3\\n 89\.06 00\.00 00\.78 10\.16\\n\\033\[0;34m\[attr\] 'r' \(estimated\)\\033\[0m\\n  M\+H   M\+   M-H   M\+H2O-CH4\\n 0\.81 0\.00  0\.10        0\.09 \\n\\033\[0;34m\[attr\] 'err'\\033\[0m\\n0\.003158\\033\[0;34m\[class\] 'CorMID'\\033\[0m")
    suppressMessages({
      # check if print returns 4 MI's
      testthat::expect_output(print(out), "    M0    M1    M2    M3")
      # check if print returns the correct err value
      testthat::expect_output(print(out), "003158")
    })
  }
)

# clean up of global objects
rm(fml, mid, r, rMID, out)
