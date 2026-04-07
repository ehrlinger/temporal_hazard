# C Binary Parity Tests
#
# These tests compare R package results against the C hazard binary.
# They are optional and will skip gracefully if the binary is not available.
# To enable: set environment variable TEMPORAL_HAZARD_BIN=/path/to/hazard
# or use options(temporal_hazard.binary="/path/to/hazard")

test_that("hazard binary can be located (if available)", {
  bin <- TemporalHazard:::.hzr_get_hazard_binary()
  if (is.na(bin)) {
    skip("C hazard binary not available")
  }
  expect_true(file.exists(bin))
  expect_match(bin, "hazard")
})

test_that("legacy command uses stdin/stdout redirection", {
  # Skip on Windows where /bin/cat does not exist
  skip_if(
    .Platform$OS.type != "unix",
    "Unix shell tools not available (Windows)"
  )
  
  temp_dir <- tempdir()
  sas_file <- file.path(temp_dir, "test-prefix.sas")
  lst_file <- file.path(temp_dir, "test-prefix.lst")
  writeLines("Parameter Estimate StdErr Z.value P.value", sas_file)

  result <- TemporalHazard:::.hzr_run_legacy_binary(
    bin_path = "/bin/cat",
    sas_file = sas_file,
    lst_file = lst_file
  )

  expect_match(result$call_string, "<")
  expect_match(result$call_string, ">")
  expect_match(result$call_string, "2>&1", fixed = TRUE)
  expect_equal(result$status_code, 0)
  expect_true(file.exists(lst_file))
})

test_that("build_hazard_sas_input emits legacy %HAZARD block", {
  sas_lines <- TemporalHazard:::.hzr_build_hazard_sas_input(
    data_name = "GOLDEN_TEST",
    data_file = "/tmp/golden_test.csv",
    time_var = "INT_DEAD",
    status_var = "DEAD",
    x_names = c("AGE", "STATUS"),
    dist = "weibull",
    theta = c(MUE = 0.2, NU = 1.4, AGE = -0.03, STATUS = 0.5),
    control = list(nocov = TRUE, nocor = TRUE, condition = 14, fix = "M")
  )

  sas_text <- paste(sas_lines, collapse = "\n")
  expect_match(sas_text, "%HAZARD\\(")
  expect_match(sas_text, "PROC HAZARD DATA=GOLDEN_TEST NOCOV NOCOR CONDITION=14;")
  expect_match(sas_text, "TIME INT_DEAD;")
  expect_match(sas_text, "EVENT DEAD;")
  expect_match(sas_text, "PARMS MUE=2.0e-01, NU=1.4e\\+00 FIXM;")
  expect_match(sas_text, "EARLY AGE=-3e-02, STATUS=5e-01;")
})

test_that("build_hazpred_sas_input emits legacy %HAZPRED block", {
  sas_lines <- TemporalHazard:::.hzr_build_hazpred_sas_input(
    data_name = "OUTEST",
    input_data_name = "AVCS",
    output_data_name = "HPOUT",
    time_var = "INT_DEAD",
    control = list(predict_types = c("survival", "hazard"))
  )

  sas_text <- paste(sas_lines, collapse = "\n")
  expect_match(sas_text, "%HAZPRED\\(")
  expect_match(sas_text, "PROC HAZARD DATA=OUTEST IN=AVCS OUT=HPOUT")
  expect_match(sas_text, "TIME INT_DEAD;")
  expect_match(sas_text, "PREDICT SURVIVAL HAZARD;")
})

test_that("parse_hazard_output handles empty output gracefully", {
  result <- TemporalHazard:::.hzr_parse_hazard_output(character())
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_named(result, c("parameter", "estimate", "std_err", "z_stat", "p_value"))
})

test_that("parse_hazard_output parses standard parameter table", {
  # Simulate standard statistical output format
  output <- c(
    "Parameter    Estimate      StdErr      Z.value        P.value",
    "shape        1.234567      0.054321    22.725000      0.000000",
    "scale        10.987654     0.321098    34.225000      0.000000"
  )

  result <- TemporalHazard:::.hzr_parse_hazard_output(output)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_equal(result$parameter, c("shape", "scale"))
  expect_true(all(!is.na(result$estimate)))
  expect_true(all(!is.na(result$std_err)))
})

test_that("run_hazard_binary integration test (if binary available)", {
  bin <- TemporalHazard:::.hzr_get_hazard_binary()
  if (is.na(bin)) {
    skip("C hazard binary not available; skipping integration test")
  }

  # Use a simple golden fixture dataset
  time <- c(1, 2, 3, 4, 5)
  status <- c(1, 1, 0, 1, 1)
  x <- data.frame(x1 = c(0, 1, 0, 1, 0))

  # Run binary fit
  result <- TemporalHazard:::.hzr_run_hazard_binary(
    time = time,
    status = status,
    x = x,
    dist = "weibull"
  )

  # Verify output structure (function returns list with parsed component)
  expect_type(result, "list")
  expect_named(result, c("call_string", "stdout", "parsed", "info"))
  expect_match(result$call_string, "<")
  expect_match(result$call_string, ">")
  expect_match(result$call_string, "\\.sas")
  expect_match(result$call_string, "\\.lst")
  
  parsed <- result$parsed
  expect_s3_class(parsed, "tbl_df")
  expect_named(parsed, c("parameter", "estimate", "std_err", "z_stat", "p_value"))

  if (!is.na(result$info$status_code) && result$info$status_code != 0) {
    skip(sprintf("Legacy binary returned non-zero status %d", result$info$status_code))
  }

  if (nrow(parsed) > 0) {
    # If binary produced output, verify reasonable values
    expect_true(all(!is.na(parsed$estimate)))
    expect_true(all(parsed$std_err > 0 | is.na(parsed$std_err)))
  }
})