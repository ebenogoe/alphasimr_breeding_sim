# Load required libraries
library(AlphaSimR)
library(ggplot2)
library(pegas)
library(adegenet)
library(dartR)
library(future)
library(future.apply)
library(progressr)

#----------------
# SIMULATION PARAMETERS
#----------------
Ne <- 25
n.founders <- Ne
n.chr <- 10
n.sites <- 300
n.gens <- 10
n.crosses <- 50
n.progeny <- 10
n.reps <- 10
prop.markers <- 0.5 # proportion of markers to use. Set to 1 to use all markers

#----------------
# Functions to estimate Ne
#----------------
estimate_Ne_from_LD <- function(geno_matrix, n_pairs = 1000) {
  n_ind <- nrow(geno_matrix)
  snps <- ncol(geno_matrix)
  r2_vals <- numeric()
  
  for (i in 1:n_pairs) {
    idx <- sample(1:snps, 2)
    snp1 <- geno_matrix[, idx[1]]
    snp2 <- geno_matrix[, idx[2]]
    if (var(snp1, na.rm = TRUE) > 0 && var(snp2, na.rm = TRUE) > 0) {
      r2 <- cor(snp1, snp2)^2
      r2_vals <- c(r2_vals, r2)
    }
  }
  
  mean_r2 <- mean(r2_vals, na.rm = TRUE)
  Ne_val <- 1 / (3 * (mean_r2 - (1 / n_ind)))
  list(Ne = Ne_val, r2 = mean_r2)
}

estimate_Ne_from_pegas <- function(geno_matrix, n_pairs = 1000) {
  n_ind <- nrow(geno_matrix)
  snps <- ncol(geno_matrix)
  r2_vals <- numeric()
  
  for (i in 1:n_pairs) {
    idx <- sample(1:snps, 2)
    mat <- geno_matrix[, idx]
    if (var(mat[,1], na.rm = TRUE) == 0 || var(mat[,2], na.rm = TRUE) == 0) next
    
    df <- as.loci(mat)
    tryCatch({
      res <- LD(df)
      r <- res$`Correlations among alleles`[1, 2]
      r2_vals <- c(r2_vals, r^2)
    }, error = function(e) {})
  }
  
  mean_r2 <- mean(r2_vals, na.rm = TRUE)
  Ne_val <- if (!is.na(mean_r2) && (mean_r2 - 1/n_ind) > 0) {
    1 / (3 * (mean_r2 - 1/n_ind))
  } else {
    NA
  }
  
  list(Ne = Ne_val, r2 = mean_r2)
}

estimate_Ne_from_dartR <- function(geno) {
  # Convert geno to genlight
  gl <- new("genlight", geno)
  
  #gl@ploidy <- as.integer(2)
  markers <- ncol(geno)
  gl@chromosome <- factor(rep(paste0("Chr", 1:n.chr), length.out = markers))
  
  # markers <- floor(n.sites * prop.markers)
  #gl@chromosome <- factor(rep(paste0("Chr", 1:n.chr), each = markers))
  
  indNames(gl) <- paste0("Ind", 1:nrow(geno))
  locNames(gl) <- paste0("SNP", 1:ncol(geno))
  
  # Run gl.LDNe
  result <- dartR.popgen::gl.LDNe(x = gl, plot.out = FALSE, singleton.rm = TRUE,
                                  mating = "random", verbose = 0, pairing = "separate")
  
  return(result)
}

#----------------
# parallel exec
#----------------
plan(multisession, workers = parallel::detectCores() - 2)  # or multicore if you get an error
handlers(global = TRUE)
# handlers("txtprogressbar") # DOES NOT WORK ON LOCAL MAC. 
# Causes console to freeze until 2 errors have been intentionally thrown
handlers("rstudio")

#----------------
# Parallel Simulations
#----------------
message("Running simulations with future_lapply...\n")

results.list <- data.frame()

single_rep <- function(rep) {
  #cat(sprintf("Rep %d...", rep), "\n")
  
  founders <- runMacs(nInd = n.founders,
                      nChr = n.chr,
                      segSites = n.sites,
                      inbred = FALSE,
                      #manualGenLen = 1,
                      nThreads = 1)
  
  SP <<- SimParam$new(founders)
  SP$addTraitA(nQtlPerChr = 10, mean = 0, var = 1)
  SP$addSnpChip(n.sites * 0.5, minSnpFreq = 0.05)
  SP$nThreads <- 3
  
  pop <- newPop(founders, simParam = SP)
  pop <- randCross(pop, nCrosses = n.crosses, nProgeny = n.progeny, simParam = SP)
  
  geno <- pullSegSiteGeno(pop)
  
  # Subset x% of markers. Random selecction
  n.markers <- ncol(geno)
  selected.cols <- sample(seq_len(n.markers), size = floor(n.markers * prop.markers), replace = FALSE)
  geno.sub <- geno[, selected.cols]
  
  # Estimate Ne using subsetted geno
  #est1 <- estimate_Ne_from_LD(geno.sub)
  #est2 <- estimate_Ne_from_pegas(geno.sub)
  est3 <- estimate_Ne_from_dartR(geno.sub)
  
  data.frame(
    rep = rep,
    #Ne_LD = est1$Ne,
    #r2_LD = est1$r2,
    #Ne_pegas = est2$Ne,
    #r2_pegas = est2$r2,
    Ne_glLDNe = as.numeric(est3$Pop1$`Frequency 1`[6])
  )
}

# Parallelized call to single_rep using future.apply with time tracking
start.time <- Sys.time()
with_progress({
  p <- progressor(steps = n.reps)
  results.list <- future_lapply(1:n.reps, function(rep) {
    #p(sprintf("Running rep %d", rep))
    single_rep(rep)
  }, future.seed = TRUE)
})
end.time <- Sys.time()
cat("\nStarted at", format(start.time, "%X"), "\n")
cat("Ended at", format(end.time, "%X"), "\n")

time.taken <- end.time - start.time
cat("It took:", time.taken, "\n\n")


# Non-parallelized call to single_rep to compare time taken
starting <- Sys.time()
for (i in 1:30) single_rep(i)
ending <- Sys.time()
cat("\nStarted at", format(starting, "%X"), "\n")
cat("Ended at", format(ending, "%X"), "\n")

time.spent <- ending - starting
cat("It took:", time.spent, "\n\n")

results <- do.call(rbind, results.list)


# WIP - Parallelized call to gl.LDNe using future.apply
with_progress({
  p <- progressor(steps = 1)
  results.list <- future_lapply(1, function(rep) {
    p(sprintf("Running rep %d", rep))
    
    estimate_Ne_from_dartR(geno)
  }, future.seed = TRUE)
})

#----------------
# Summary
#----------------
mean_LD <- mean(results$Ne_LD, na.rm = TRUE)
mean_pegas <- mean(results$Ne_pegas, na.rm = TRUE)
mean_glLDNe <- mean(results$Ne_glLDNe, na.rm = TRUE)

start.time <- Sys.time()
estimate_Ne_from_dartR(geno)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
cat("Started gl.LDNe analysis at", format(start.time, "%X"), "\n")
cat("Ended at", format(end.time, "%X"), "\n")
cat("Time taken to run gl.LDNe is: ", time.taken, "\n")

cat("\nAverage Ne (LD method):", mean_LD, "\n")
cat("Average Ne (pegas method):", mean_pegas, "\n")
cat("\nAverage Ne (gl.LDNe method):", mean_glLDNe, "\n")

cat("\nAll Replicate Results:\n")
print(results)

#----------------
# Plot histograms
#----------------
ggplot(results, aes(x = Ne_LD)) + 
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Ne (Manual LD Calculation)", x = "Ne (LD)", y = "Frequency") +
  theme_minimal()

ggplot(results, aes(x = Ne_pegas)) + 
  geom_histogram(binwidth = 0.5, fill = "green", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Ne (Pegas)", x = "Ne (Pegas)", y = "Frequency") +
  theme_minimal()

ggplot(results, aes(x = Ne_glLDNe)) + 
  geom_histogram(binwidth = 0.5, fill = "pink", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Ne (gl.LDNe)", x = "Ne (gl.LDNe)", y = "Frequency") +
  theme_minimal()


## Estimating QTL
# SP$traits[[1]]@lociPerChr
qtl_genos <- pullQtlGeno(pop, simParam = SP)
allele_freqs <- colMeans(qtl_genos) / 2  # For diploid 0/1/2 genos

# QTLs that are not fixed
segregating_qtls <- sum(allele_freqs > 0 & allele_freqs < 1)