# Generate a function to interpolate
func <- function(x) sin(x)

# Generate x values to interpolate
x <- seq(0, 2 * pi, length.out = 100)

# Generate interpolated y values
y <- func(x)

# Define test cases
test_that("interpolateFunction works as expected", {
    expect_equal(interpolateFunction(func, 0, 2 * pi, 100, x), y, tolerance = 1e-4)
})

test_that("interpolateFunction works as expected when values are out of range", {
    expect_equal(interpolateFunction(func, 0, 2 * pi, 100, c(0, 2 * pi)), c(y[1], y[100]), tolerance = 1e-4)
})
