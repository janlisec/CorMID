testthat::test_that(
  desc = "recMID returns expected result",
  code = {
    # simple calculations
    #rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2", algo = "Rdisop")
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2")
    testthat::expect_true(inherits(rMID, "recMID"))
    testthat::expect_true(is.numeric(rMID))
    testthat::expect_false(is.null(names(rMID)))
    testthat::expect_equal(unname(rMID[1:3]), c(0.06888, 0.01414, 0.55735), tolerance = 0.001)
  }
)

# avoid creating a Rplots.pdf in testthat folder
pdf(NULL)

testthat::test_that(
  desc = "recMID plot returns expected result",
  code = {
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2")
    vdiffr::expect_doppelganger(
      title = "recMID_Plot",
      fig = function() plot(rMID)
    )
  }
)

testthat::test_that(
  desc = "recMID plot returns expected result with alternative options",
  code = {
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2")
    vdiffr::expect_doppelganger(
      title = "recMID_Plot_alt",
      fig = function() plot(unname(rMID), xlab="test", ylab="test", lwd=9, lend=2, xlim=c(0,10), ylim=c(0,2), las=1)
    )
  }
)

testthat::test_that(
  desc = "recMID plot returns expected result with more alternative options",
  code = {
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2")
    names(rMID) <- 1:length(rMID)
    vdiffr::expect_doppelganger(
      title = "recMID_Plot_alt2",
      fig = function() plot(rMID, col=1)
    )
  }
)
