test_that("hzr_argument_mapping returns required schema", {
  map <- hzr_argument_mapping()

  expect_true(is.data.frame(map))
  expect_true(nrow(map) >= 8)

  required_cols <- c(
    "sas_statement", "legacy_input", "c_concept", "r_parameter",
    "required", "expected_type", "transform_rule",
    "implementation_status", "notes"
  )
  expect_true(all(required_cols %in% names(map)))
})

test_that("mapping includes core hazard() arguments", {
  map <- hzr_argument_mapping(include_planned = FALSE)
  params <- unique(map$r_parameter)

  expect_true(all(c("time", "status", "x", "theta", "dist", "control", "...") %in% params))
})

test_that("mapping table has no duplicate legacy-to-parameter rows", {
  map <- hzr_argument_mapping()
  key <- paste(map$legacy_input, map$r_parameter, sep = "::")
  expect_equal(anyDuplicated(key), 0)
})
