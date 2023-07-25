testthat::test_that(
  desc = "FitMID returns expected result",
  code = {
    fml <- "C9H20O3Si2"; attr(fml,"nbio") <- 3
    mid <- c(0.9,0,0,0.1)
    r <- unlist(list("M+H"=0.8, "M+"=0.1, "M+H2O-CH4"=0.1))
    int <- CorMID::recMID(mid=mid, r=r, fml=fml)
    td <- CorMID::CalcTheoreticalMDV(fml=fml, nbio = attr(fml,"nbio"), nmz = attr(fml,"nbio")+3)
    r <- matrix(rep(c(0,1),3),nrow=2,dimnames=list(NULL,names(r)))
    out <- CorMID:::FitMID(md=int, td=td, r=r)
    testthat::expect_equal(length(out), 4L)
    testthat::expect_equal(sum(out), 100L)
    testthat::expect_true(all(c("err", "ratio", "ratio_status", "mid_status") %in% names(attributes(out))))

    testthat::expect_true(is.na(CorMID:::FitMID(md=0)))
    testthat::expect_equal(attr(CorMID:::FitMID(md=int, td=td, r=r, mid_fix = mid), "ratio"), unlist(list("M+H"=0.8, "M+"=0.1, "M+H2O-CH4"=0.1)))

    testthat::expect_output(CorMID:::FitMID(md=int, td=td, r=r, trace_steps = TRUE), "Testing 10 MID solutions")
  }
)
