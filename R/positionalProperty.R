#' Examining the mean property of amino acids by position
#'
#' This function calculates the mean selected property for 
#' amino acids along the residues of the CDR3 amino acid sequence. 
#' The ribbon surrounding the individual line represents the 95% 
#' confidence interval.
#' 
#' @details
#' More information for the individual methods can be found at the following citations:
#' 
#' \strong{Atchley:} \href{https://pubmed.ncbi.nlm.nih.gov/15851683/}{citation}
#' 
#' \strong{Kidera:} \href{https://link.springer.com/article/10.1007/BF01025492}{citation}
#' 
#' \strong{stScales:} \href{https://pubmed.ncbi.nlm.nih.gov/19373543/}{citation}
#' 
#' \strong{tScales:} \href{https://www.sciencedirect.com/science/article/pii/S0022286006006314?casa_token=uDj97DwXDDEAAAAA:VZfahldPRwU1WObySJlohudtMSDwF7nJSUzcEGwPhvkY13ALLKhs08Cf0_FyyfYZjxJlj-fVf0SM}{citation}
#' 
#' \strong{VHSE:} \href{https://pubmed.ncbi.nlm.nih.gov/15895431/}{citation}
#' 
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' positionalProperty(combined, 
#'                    chain = "TRB",
#'                    method = "Atchley", 
#'                    aa.length = 20)

#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method The method to calculate the property - "Atchley", "Kidera",
#' "stScales", "tScales", or "VHSE"
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom stats qt
#' @importFrom dplyr %>% summarise n group_by 
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of line graph of diversity by position

positionalProperty <- function(input.data, 
                               chain = "TRB", 
                               group.by = NULL, 
                               order.by = NULL,
                               aa.length = 20,
                               method = "Atchley",
                               exportTable = FALSE, 
                               palette = "inferno")  {
  options( dplyr.summarise.inform = FALSE )
  if(method %!in% c("Atchley", "Kidera", "stScales", "tScales", "VHSE")) {
    stop("Please select a compatible method.")
  }
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, "CTaa", check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, "CTaa")
  
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  #Selecting Property Function
  propertyFunc <- switch(method,
                          "Atchley" = .af.ref,
                          "Kidera" = .kf.ref,
                          "stScales" = .stscales.ref,
                          "tScales" = .tscales.ref,
                          "VHSE" = .vhse.ref,
                          stop("Invalid method provided"))
  
  #Getting AA Counts
  aa.count.list <- .aa.counter(input.data, cloneCall, aa.length)
  
  #Calculating properties and melting data
  lapply(seq_along(aa.count.list), function(x) {
      lapply(seq_len(nrow(aa.count.list[[x]]))[-1], function(y) {
          pos <- aa.count.list[[x]][!is.na(aa.count.list[[x]]$AA),y]
          names(pos) <- aa.count.list[[x]][!is.na(aa.count.list[[x]]$AA),1]
          pos <- pos[pos > 0]
          lapply(seq_len(length(pos)), function(t) {
            char <- names(pos[t])
            results <- rep(propertyFunc[char,], pos[t])
            results
          }) -> output.values
          output.values <- unlist(output.values)
          df <- data.frame(group = names(output.values), 
                           value = unlist(output.values))
          
          if(nrow(df) == 0) { #For only NA positions
            summary <- data.frame(mean = rep(0, dim(propertyFunc)[2]), 
                                  ci_lower = rep(0, dim(propertyFunc)[2]),
                                  ci_upper = rep(0, dim(propertyFunc)[2]))
            
          } else {
            summary <- df %>% 
                      group_by(group) %>% 
                      summarise(mean = mean(value),
                                sd = sd(value),  # Standard deviation
                                n = n(),         # Number of observations per group
                                se = ifelse(n > 1, sd / sqrt(n), 0), # Standard error of the mean
                                ci_lower = ifelse(n > 1, mean - qt(0.975, n-1) * se, mean),
                                ci_upper = ifelse(n > 1, mean + qt(0.975, n-1) * se, mean)) %>%
                      as.data.frame()
          
           summary <- summary[,c("mean", "ci_lower", "ci_upper")]
         }
         summary$property <- colnames(propertyFunc)
         summary
      })-> output.values
      int.mat <- bind_rows(output.values, .id = "position")
      int.mat
  }) -> property.calculations
  names(property.calculations) <- names(input.data)
  mat <- bind_rows(property.calculations, .id = "group")
  mat$position <- paste0("pos", mat$position)
  mat$position <- factor(mat$position, levels = str_sort(unique(mat$position), numeric = TRUE))
  
  if(!is.null(order.by)) {
    mat <- .ordering.function(vector = order.by,
                              group.by = "group", 
                              mat)
  }
  
    plot <- ggplot(mat, aes(x = position, 
                            y = mean, 
                            group = group, 
                            color = group)) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = group), alpha = 0.5, lwd = 0) +
      geom_line(stat = "identity", alpha = 0.5) +
      geom_point() + 
      scale_color_manual(name = "Group", 
                         values = .colorizer(palette,length(unique(mat$group)))) +
      scale_fill_manual(values = .colorizer(palette,length(unique(mat$group)))) +
      xlab("Amino Acid Residues") +
      ylab("Mean Values") +
      facet_grid(property~.) + 
      guides(fill = "none", 
             color = guide_legend(override.aes=list(fill=NA))) + 
      theme_classic() + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if (exportTable == TRUE) { 
      return(mat_melt) 
    }
    return(plot)
  
}

###############################
#Amino Acid Property Matrices
###############################

.af.ref <- matrix(c(
  -0.591, -1.302, -0.733,  1.570, -0.146, # A (Alanine)
  1.538, -0.055,  1.502,  0.440,  2.897, # R (Arginine)
  0.945,  0.828,  1.299, -0.169,  0.933, # N (Asparagine)
  1.050,  0.302, -1.768,  0.276,  1.068, # D (Aspartic Acid)
  -1.343,  0.465, -0.862, -1.020, -0.255, # C (Cysteine)
  0.931, -0.179, -3.005,  0.941,  0.360, # Q (Glutamine)
  1.357, -1.453,  1.477,  0.114, -0.384, # E (Glutamic Acid)
  -0.384,  1.652,  1.330,  1.045,  2.065, # G (Glycine)
  -0.510,  0.292, -0.203, -1.378, -0.276, # H (Histidine)
  -1.006, -0.590,  1.891, -0.397,  0.412, # I (Isoleucine)
  -1.006, -0.590,  1.891, -0.397,  0.412, # L (Leucine)
  0.960, -0.181,  1.932, -0.041,  1.697, # K (Lysine)
  -1.343,  0.465, -0.862, -1.020, -0.255, # M (Methionine)
  -1.006, -0.590,  1.891, -0.397,  0.412, # F (Phenylalanine)
  0.189,  0.001, -1.756,  0.767, -0.542, # P (Proline)
  -0.228,  1.399, -4.760,  0.670, -2.647, # S (Serine)
  -0.228,  1.399, -4.760,  0.670, -2.647, # T (Threonine)
  -1.006, -0.590,  1.891, -0.397,  0.412, # W (Tryptophan)
  -1.343,  0.465, -0.862, -1.020, -0.255, # Y (Tyrosine)
  -1.006, -0.590,  1.891, -0.397,  0.412  # V (Valine)
), nrow = 20, ncol = 5, byrow = TRUE)
rownames(.af.ref) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
colnames(.af.ref) <- paste0("AF",1:5)


.kf.ref <- matrix(c(
  -1.56, -1.67, -0.97, -0.27, -0.93, -0.78, -0.20, -2.65,  0.19,  0.82, # A (Alanine)
  -2.52,  0.27, -0.85, -0.71,  1.11,  1.15,  1.07,  0.86, -0.49,  0.95, # R (Arginine)
  -1.01,  0.53,  0.06,  0.13, -1.75,  0.16,  0.42, -0.73,  0.53,  0.81, # N (Asparagine)
  -0.51,  0.07,  1.50,  0.13,  0.32,  0.42, -1.04, -0.26,  0.76,  0.22, # D (Aspartic Acid)
  0.24, -2.32,  0.60, -0.14,  1.27,  1.15,  0.21,  1.03, -0.84, -0.45, # C (Cysteine)
  -0.22, -0.04,  0.53, -0.11, -1.18,  0.45,  0.84, -0.71,  0.82, -0.62, # Q (Glutamine)
  -0.76,  0.18,  0.06, -0.42,  1.25,  0.83,  0.51,  0.44,  0.65, -0.20, # E (Glutamic Acid)
  0.81, -0.67,  2.80,  0.60,  0.88,  0.69, -0.29,  0.81, -0.92,  0.32, # G (Glycine)
  -0.47,  1.94,  0.45, -1.61,  1.49, -0.38,  1.14,  0.74, -0.72,  1.59, # H (Histidine)
  1.97, -1.73, -0.16,  0.70, -0.83,  0.69, -0.07, -0.61,  0.81, -0.03, # I (Isoleucine)
  1.41, -0.23, -0.25,  1.53, -1.50,  0.27,  0.70, -0.10,  0.45, -0.02, # L (Leucine)
  -1.27,  1.32,  0.57,  1.62, -0.30,  0.49,  0.84, -0.60,  0.71,  0.46, # K (Lysine)
  1.23, -0.08, -1.11,  0.16, -1.10,  0.84,  0.83,  1.23, -0.32,  0.51, # M (Methionine)
  1.46, -1.96, -0.23,  0.53, -0.60,  0.75,  0.20,  0.26, -0.16,  0.84, # F (Phenylalanine)
  0.59,  0.91, -0.01, -0.28,  0.22, -0.26,  0.19,  0.36, -0.33,  0.38, # P (Proline)
  -0.34,  1.15,  0.25, -1.39,  0.67, -0.76, -0.40,  0.51, -0.27,  1.56, # S (Serine)
  0.71,  0.24, -0.96,  0.67,  1.27, -1.14,  0.54,  0.43, -0.37,  0.19, # T (Threonine)
  0.85, -0.71,  2.25,  0.65, -1.09,  0.53, -0.53,  0.16,  0.93,  0.60, # W (Tryptophan)
  0.02,  0.73, -0.16,  1.08,  1.59,  0.56, -0.82,  0.29,  0.25,  0.06, # Y (Tyrosine)
  1.08, -1.81,  0.08,  0.47, -1.04,  0.06, -0.46, -0.40,  0.50, -0.05  # V (Valine)
), nrow = 20, ncol = 10, byrow = TRUE)
rownames(.kf.ref) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
colnames(.kf.ref) <- paste0("KF",1:10)

.vhse.ref <- matrix(c(
  0.15, -1.11, -1.35, -0.92,  0.02, -0.91,  0.36, -0.48, # A (Alanine)
  -1.47, 1.45,  1.24,  1.27,  1.55,  1.47,  1.30,  0.83, # C (Cysteine)
  -0.99,  0.00, -0.37,  0.69, -0.55,  0.85,  0.73, -0.80, # N (Asparagine)
  -1.15,  0.67, -0.41, -0.01, -2.68,  1.31,  0.03,  0.56, # D (Aspartic Acid)
  0.18, -1.67, -0.46, -0.21,  0.00,  1.20, -1.61, -0.19, # C (Cysteine)
  -0.96,  0.12,  0.18,  0.16,  0.09,  0.42, -0.20, -0.41, # Q (Glutamine)
  -1.18,  0.40,  0.10,  0.36, -2.16, -0.17,  0.91,  0.02, # E (Glutamic Acid)
  -0.20, -1.53, -2.63,  2.28, -0.53, -1.18,  2.01, -1.34, # G (Glycine)
  -0.43, -0.25,  0.37,  0.19,  0.51,  1.28,  0.93,  0.65, # H (Histidine)
  1.27, -0.14,  0.30, -1.80,  0.30, -1.61, -0.16, -0.13, # I (Isoleucine)
  1.36,  0.07,  0.36, -0.80,  0.22, -1.37,  0.08, -0.62, # L (Leucine)
  -1.17,  0.70,  0.70,  0.80,  1.64,  0.67,  1.63 , 0.13, # K (Lysine)
  1.01, -0.53,  0.43,  0.00,  0.23,  0.10, -0.86, -0.68, #M (Methionine)
  1.52,  0.61,  0.96, -0.16,  0.25,  0.28, -1.33, -0.20, # F (Phenylalanine)
  0.22, -0.17, -0.50,  0.05, -0.01, -1.34, -0.19, 3.56, # P (Proline)
  -0.67, -0.86, -1.07, -0.41, -0.32,  0.27, -0.64, 0.11,  # S (Serine)
  -0.34, -0.51, -0.55, -1.06, -0.06, -0.01, -0.79, 0.39, # T (Threonine)
  1.50,  2.06,  1.79,  0.75,  0.75, -0.13, -1.01, -0.85, # W (Tryptophan)
  0.61,  1.60,  1.17,  0.73,  0.53,  0.25, -0.96, -0.52, #Y (Tyrosine)
  0.76, -0.92, -0.17, -1.91,  0.22, -1.40, -0.24, -0.03 #V (Valine)
), nrow = 20, ncol = 8, byrow = TRUE)
rownames(.vhse.ref) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
colnames(.vhse.ref) <- paste0("VHSE",1:8)

.stscales.ref <- matrix(c(
  -1.552, -0.791, -0.627,  0.237, -0.461, -2.229,  0.283,  1.221,
  -0.059,  0.731, -0.013, -0.096, -0.253,  0.300,  1.256,  0.854,
  -0.888, -0.057, -0.651, -0.214,  0.917,  0.164, -0.140, -0.166,
  -0.907, -0.054, -0.781, -0.248,  1.120,  0.101, -0.245, -0.075,
  -1.276, -0.401,  0.134,  0.859, -0.196, -0.720,  0.639, -0.857,
  -0.662,  0.228, -0.193, -0.105,  0.418,  0.474,  0.172,  0.408,
  -0.629, -0.390, -0.380, -0.366,  0.635,  0.514,  0.175,  0.367,
  -1.844, -0.018, -0.184,  0.573, -0.728, -3.317,  0.166,  2.522,
  -0.225,  0.361,  0.079, -1.037,  0.568,  0.273,  1.208, -0.001,
  -0.785, -1.010, -0.349, -0.097, -0.402,  1.091, -0.139, -0.764,
  -0.826, -0.379,  0.038, -0.059, -0.625,  1.025, -0.229, -0.129,
  -0.504,  0.245,  0.297, -0.065, -0.387,  1.011,  0.525,  0.553,
  -0.693,  0.498,  0.658,  0.457, -0.231,  1.064,  0.248, -0.778,
  -0.019,  0.024,  1.080, -0.220, -0.937,  0.570, -0.357,  0.278,
  -1.049, -0.407, -0.067, -0.066, -0.813, -0.890,  0.021, -0.894,
  -1.343, -0.311, -0.917, -0.049,  0.549, -1.533,  0.166,  0.280,
  -1.061, -0.928, -0.911, -0.063,  0.538, -0.775, -0.147, -0.717,
  0.853,  0.039,  0.260,-1.163,  0.160, -0.202,  1.010,  0.195,
  0.308,  0.569,  1.100, -0.464, -0.144, -0.354, -1.099,  0.162,
  -1.133, -0.893, -0.325,  0.303, -0.561, -0.175, -0.020, -0.311
), nrow = 20, ncol = 8, byrow = TRUE)
rownames(.stscales.ref) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
colnames(.stscales.ref) <- paste0("stScales",1:8)

.tscales.ref <- matrix(c(
  -9.11, -1.63,  0.63,  1.04,  2.26,
  0.23,  3.89, -1.16, -0.39, -0.06,
  -4.62,  0.66,  1.16, -0.22,  0.93,
  -4.65,  0.75,  1.39, -0.40,  1.05,
  -7.35, -0.86, -0.33,  0.80,  0.98,
  -3.00,  1.72,  0.28, -0.39,  0.33,
  -3.03,  1.82,  0.51, -0.58,  0.43,
  -10.61, -1.21, -0.12,  0.75,  3.25,
  -1.01, -1.31,  0.01, -1.81, -0.21,
  -4.25, -0.28, -0.15,  1.40, -0.21,
  -4.38,  0.28, -0.49,  1.45,  0.02,
  -2.59,  2.34, -1.69,  0.41, -0.21,
  -4.08,  0.98, -2.34,  1.64, -0.79,
  0.49, -0.94, -0.63, -1.27, -0.44,
  -5.11, -3.54, -0.53, -0.36, -0.29,
  -7.44, -0.65,  0.68, -0.17,  1.58,
  -5.97, -0.62,  1.11,  0.31,  0.95,
  5.73, -2.67, -0.07, -1.96, -0.54,
  2.08, -0.47,  0.07, -1.67, -0.35,
  -5.87, -0.94,  0.28,  1.10,  0.48
), nrow = 20, ncol = 5, byrow = TRUE)
rownames(.tscales.ref) <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
colnames(.tscales.ref) <- paste0("tScales",1:5)

