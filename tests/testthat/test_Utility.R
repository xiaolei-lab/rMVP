test_that("MVP.Version() works fine.", {
    expect_output(MVP.Version())
    expect_output(MVP.Version(start = FALSE))
})
