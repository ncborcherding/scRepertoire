# test script for clonalSizeDistribution.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                      samples = c("P17B", "P17L", "P18B", "P18L",
                                  "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Output formats (distance matrix vs. plot) are correct", {
  # Test that exportTable = TRUE returns a matrix
  dist_matrix <- clonalSizeDistribution(combined[1:2], exportTable = TRUE)
  expect_true(is.matrix(dist_matrix))
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalSizeDistribution(combined[1:2])
  expect_s3_class(plot_output, "ggplot")
})

test_that("Distance matrix has correct properties", {
  dist_matrix <- clonalSizeDistribution(combined[1:3], exportTable = TRUE)
  num_samples <- length(combined[1:3])
  
  # Matrix should be square and have correct dimensions
  expect_equal(nrow(dist_matrix), num_samples)
  expect_equal(ncol(dist_matrix), num_samples)
  
  # Should have correct row and column names
  expect_equal(rownames(dist_matrix), names(combined)[1:3])
  expect_equal(colnames(dist_matrix), names(combined)[1:3])
  
  # Diagonal should be zero
  expect_true(all(diag(dist_matrix) == 0))
  
  # Matrix should be symmetric
  expect_equal(dist_matrix, t(dist_matrix))
  
  # All values should be non-negative
  expect_true(all(dist_matrix >= 0))
})

test_that("`group.by` parameter correctly aggregates data", {
  dist_matrix_grouped <- clonalSizeDistribution(combined, group.by = "Type", exportTable = TRUE)
  num_groups <- 2
  expect_equal(nrow(dist_matrix_grouped), num_groups)
  expect_equal(ncol(dist_matrix_grouped), num_groups)
  expect_true(all(rownames(dist_matrix_grouped) %in% c("B", "L")))
})

test_that("`method` parameter is passed to hclust correctly", {
  expect_silent(clonalSizeDistribution(combined[1:2], method = "ward.D2"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "complete"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "single"))
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  expect_silent(
    res_aa <- clonalSizeDistribution(combined[1:2], cloneCall = "aa", exportTable = TRUE)
  )
  expect_true(is.matrix(res_aa))
  
  expect_silent(
    res_tra <- clonalSizeDistribution(combined[1:2], chain = "TRA", exportTable = TRUE)
  )
  expect_true(is.matrix(res_tra))
})

test_that("Plotting output has the correct dendrogram components", {
  plot_output <- clonalSizeDistribution(combined[5:6])
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomSegment")
  expect_s3_class(plot_output$layers[[2]]$geom, "GeomText")
  expect_s3_class(plot_output$layers[[3]]$geom, "GeomPoint")
  expect_s3_class(plot_output$coordinates, "CoordFlip")
  expect_true(plot_output$scales$scales[[1]]$trans$name == "reverse")
})

#################################################################
# Section 2: Testing Core Modeling Functions
#################################################################

test_that(".fdiscgammagpd returns a correctly structured list", {
  # Generate sample data that looks like clonal counts
  set.seed(42)
  x <- rpois(100, lambda = 2) + 1 
  
  fit <- .fdiscgammagpd(x, useq = 5)
  
  expect_type(fit, "list")
  expect_named(fit, c("x", "shift", "init", "useq", "nllhuseq", "optim", "nllh", "mle", "fisherInformation"))
  expect_type(fit$mle, "double")
  expect_equal(length(fit$mle), 6)
  expect_named(fit$mle, c("phi", "shape", "rate", "thresh", "sigma", "xi"))
})

test_that("NLL functions return a single positive number", {
  # Mock parameters and data
  params <- c(1.5, 2.5) # Log-scale parameters for shape, rate
  data <- c(1, 2, 3, 4)
  
  nll_gamma <- .discgammanll(log(params), dat = data, thresh = 5, phiu = 0.1, shift = 0)
  expect_type(nll_gamma, "double")
  expect_length(nll_gamma, 1)
  expect_true(nll_gamma > 0)
  
  # Mock parameters for GPD
  gpd_params <- c(log(1), 0.1) # log(sigma), xi
  gpd_data <- c(6, 7, 8)
  
  nll_gpd <- .discgpdnll(gpd_params, dat = gpd_data, thresh = 5, phiu = 0.1)
  expect_type(nll_gpd, "double")
  expect_length(nll_gpd, 1)
  expect_true(nll_gpd > 0)
})

test_that("Discrete distribution functions work as expected", {
  # Test that the probability mass functions are non-negative
  ddg <- .ddiscgamma(x = 1:5, shape = 1, rate = 1, thresh = 6, phiu = 0.5, shift = 0)
  expect_true(all(ddg >= 0))
  
  ddgpd <- .ddiscgpd(x = 6:10, thresh = 5, sigma = 1, xi = 0.1, phiu = 0.5)
  expect_true(all(ddgpd >= 0))
  
  # Test the random number generators
  set.seed(42)
  rands <- .rdiscgamma(100, shape = 2, rate = 1, thresh = 10, shift = 0)
  expect_length(rands, 100)
  expect_true(all(rands < 10)) # All values should be below the threshold
  expect_true(all(rands == floor(rands))) # All values should be integers
})

#################################################################
# Section 3: Testing Distance Calculation Functions
#################################################################


test_that(".JS_dist returns a single non-negative number", {
  set.seed(1)
  fit1 <- .fdiscgammagpd(rpois(100, 4) + 1, useq = 10)
  fit2 <- .fdiscgammagpd(rpois(100, 5) + 1, useq = 10)
  
  jsd <- .JS_dist(fit1, fit2, grid = 0:1000)
  
  expect_type(jsd, "double")
  expect_length(jsd, 1)
  expect_gte(jsd, 0) 
})
