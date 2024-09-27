testthat::test_that(
  desc = "recMID returns expected result",
  code = {
    # simple calculations
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2", algo = "Rdisop")
    testthat::expect_true(inherits(rMID, "recMID"))
    testthat::expect_true(is.numeric(rMID))
    testthat::expect_false(is.null(names(rMID)))
    testthat::expect_equal(unname(rMID[1:3]), c(0.06881282, 0.01414192, 0.55678069))
  }
)

testthat::test_that(
  desc = "recMID plot returns expected result",
  code = {
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2", algo = "Rdisop")
    vdiffr::expect_doppelganger(
      title = "recMID_Plot",
      fig = function() plot(rMID)
    )
  }
)

testthat::test_that(
  desc = "recMID plot returns expected result with alternative options",
  code = {
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2", algo = "Rdisop")
    vdiffr::expect_doppelganger(
      title = "recMID_Plot_alt",
      fig = function() plot(unname(rMID), xlab="test", ylab="test", lwd=9, lend=2, xlim=c(0,10), ylim=c(0,2), las=1)
    )
  }
)

testthat::test_that(
  desc = "recMID plot returns expected result with more alternative options",
  code = {
    # avoid creating a Rplots.pdf in testthat folder
    pdf(NULL)
    rMID <- CorMID::recMID(mid=c(0.9,0,0,0.1), r=list("M+H"=0.8, "M-H"=0.1, "M+H2O-CH4"=0.1), fml="C9H20O3Si2", algo = "Rdisop")
    names(rMID) <- 1:length(rMID)
    vdiffr::expect_doppelganger(
      title = "recMID_Plot_alt2",
      fig = function() plot(rMID, col=1)
    )
  }
)
