testthat::test_that(
  desc = "CalcTheoreticalMDV returns expected result",
  code = {
    #testthat::expect_equal(unname(CorMID::CalcTheoreticalMDV(fml = "C5H6Si1", algo = "Rdisop")[1,]), c(0.87343364, 0.09358505, 0.03298131))
    testthat::expect_equal(unname(CorMID::CalcTheoreticalMDV(fml = "C5H6Si1")[1,]), c(0.87482, 0.09235, 0.03282), tolerance = 0.0001)
  }
)
