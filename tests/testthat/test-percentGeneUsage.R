# test script for percentGeneUsage.R - testcases are NOT comprehensive!

# Generate some data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                    "P19B","P19L", "P20B", "P20L"))

combined <- addVariable(combined, 
            variable.name = "Type", 
            variables = rep(c("B", "L"), 4))

test_that("percentGeneUsage returns a ggplot object for single gene heatmap", {
  p <- percentGeneUsage(combined, 
                        genes = "TRBV", 
                        group.by = "sample", 
                        plot.type = "heatmap")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
})

test_that("percentGeneUsage returns a ggplot object for single gene barplot", {
  p <- percentGeneUsage(combined, 
                        genes = "TRBV", 
                        group.by = "sample", 
                        plot.type = "barplot")
  expect_s3_class(p, "ggplot")
  expect_true("GeomBar" %in% class(p$layers[[1]]$geom))
})

test_that("percentGeneUsage returns a ggplot object for paired gene heatmap", {
  p <- percentGeneUsage(combined, 
                        genes = c("TRBV", "TRBJ"), 
                        group.by = "sample", 
                        plot.type = "heatmap")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
})

test_that("percentGeneUsage returns a matrix when exportTable = TRUE (single gene)", {
  mat <- percentGeneUsage(combined, 
                          genes = "TRBV", 
                          group.by = "sample", 
                          exportTable = TRUE)
  expect_type(mat, "double")
  expect_true(is.matrix(mat))
  expect_equal(colnames(mat), c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
  expect_true(any(rownames(mat) %in% c("TRBV7-2", "TRBV6-5")))
})

test_that("percentGeneUsage returns a matrix when exportTable = TRUE (paired genes)", {
  mat <- percentGeneUsage(combined, 
                          genes = c("TRBV", "TRBJ"), 
                          group.by = "sample", 
                          exportTable = TRUE)
  expect_type(mat, "double")
  expect_true(is.matrix(mat))
  expect_equal(colnames(mat), c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
})

test_that("percentGeneUsage handles different summary.fun (percent, proportion, count)", {
  mat_percent <- percentGeneUsage(combined, 
                                  genes = "TRBV",
                                  group.by = "sample",
                                  summary.fun = "percent", 
                                  exportTable = TRUE)
  mat_proportion <- percentGeneUsage(combined, 
                                     genes = "TRBV", 
                                     group.by = "sample",
                                     summary.fun = "proportion", 
                                     exportTable = TRUE)
  mat_count <- percentGeneUsage(combined, 
                                genes = "TRBV", 
                                group.by = "sample",
                                summary.fun = "count", 
                                exportTable = TRUE)
  
  # Check sums for percentages and proportions
  expect_equal(colSums(mat_percent)[1], c(P17B = 100), tolerance = 1e-6)
  expect_equal(colSums(mat_proportion)[1], c(P17B = 1), tolerance = 1e-6)
  
  expect_equal(mat_proportion * 100, mat_percent, tolerance = 1e-6)
  expect_equal(mat_count["TRBV7-2", "P17B"], 44)
  expect_equal(mat_count["TRBV6-5", "P17B"], 98)
})



test_that("percentGeneUsage handles order.by = custom vector", {
  custom_order <- c("P17L", "P17B")
  p <- percentGeneUsage(combined[1:2], 
                        genes = "TRBV", 
                        group.by = "sample", 
                        order.by = custom_order)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$Group), custom_order)
})

test_that("percentGeneUsage throws error for invalid genes parameter", {
  expect_error(percentGeneUsage(combined, genes = "InvalidGene"), 
               "Unsupported gene segment format")
  expect_error(percentGeneUsage(combined, genes = c("TRBV", "TRBJ", "TRBC")), 
               "Parameter 'genes' must be a character vector of length 1 or 2")
})

test_that("percentGeneUsage throws error for invalid group.by column", {
  expect_error(percentGeneUsage(combined, genes = "TRBV", 
                                group.by = "non_existent_col"), 
               regexp = "column 'non_existent_col' not found")
})

# --- Alias Function Tests ---

test_that("vizGenes calls percentGeneUsage correctly for single gene heatmap", {
  p <- vizGenes(combined, 
                x.axis = "TRBV", 
                group.by = "sample", 
                plot = "heatmap", 
                summary.fun = "proportion")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
  expect_equal(p$labels$fill, "Proportion")
})

test_that("vizGenes calls percentGeneUsage correctly for single gene barplot", {
  p <- vizGenes(combined, 
                x.axis = "TRBV", 
                group.by = "sample", 
                plot = "barplot")
  expect_s3_class(p, "ggplot")
  expect_true("GeomBar" %in% class(p$layers[[1]]$geom))
})

test_that("vizGenes handles y.axis as categorical variable (maps to group.by)", {
  p <- vizGenes(combined, 
                x.axis = "TRBV", 
                y.axis = "Type", 
                plot = "heatmap")
  expect_s3_class(p, "ggplot")
})

test_that("vizGenes handles y.axis as a gene (maps to paired genes)", {
  p <- vizGenes(combined, 
                x.axis = "TRBV", 
                y.axis = "TRAV", 
                group.by = "sample",
                plot = "heatmap")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
  # Check if it's a paired gene plot (Var1 and Var2 in data)
  expect_true("Var1" %in% colnames(p$data))
  expect_true("Var2" %in% colnames(p$data))
})

test_that("percentGenes calls percentGeneUsage correctly", {
  p <- percentGenes(combined, 
                    chain = "TRB", 
                    gene = "Vgene", 
                    group.by = "sample")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
  # Check if it's a single gene plot (only Var1 in data, and Group)
  expect_true("Var1" %in% colnames(p$data))
  expect_false("Var2" %in% colnames(p$data))
  expect_equal(p$labels$fill, "Percent")
})

test_that("percentGenes handles different summary.fun", {
  p_count <- percentGenes(combined, 
                          chain = "TRB", 
                          gene = "Vgene", 
                          group.by = "sample",
                          summary.fun = "count")
  expect_equal(p_count$labels$fill, "Count")
})

test_that("percentVJ calls percentGeneUsage correctly", {
  p <- percentVJ(combined, chain = "TRB", group.by = "sample")
  expect_s3_class(p, "ggplot")
  expect_true("GeomTile" %in% class(p$layers[[1]]$geom))
  # Check if it's a paired gene plot (Var1 and Var2 in data)
  expect_true("Var1" %in% colnames(p$data))
  expect_true("Var2" %in% colnames(p$data))
  expect_equal(p$labels$fill, "Percent")
})

test_that("percentVJ handles different chains correctly", {
  # For TRA, it should use TRAV and TRAJ
  p_tra <- percentVJ(combined, chain = "TRA", group.by = "sample")
  expect_s3_class(p_tra, "ggplot")
  expect_true(all(grepl("TRAV", levels(p_tra$data$Var1))))
  expect_true(all(grepl("TRAJ", levels(p_tra$data$Var2))))
})







