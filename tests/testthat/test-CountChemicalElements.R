testthat::test_that(
  desc = "CountChemicalElements returns expected result",
  code = {
    testthat::expect_equal(unname(CorMID::CountChemicalElements("CHC")), c(2,1))
  }
)
