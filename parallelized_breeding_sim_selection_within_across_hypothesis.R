# Load needed libraries
library(AlphaSimR)
library(ggplot2)
library(parallel)
library(dplyr)

# Set Simulation Parameters
i <- 10 # Number of founders
nChr <- 1
segSites <- 100 
species <- "MAIZE"
heritability <- 0.3
qtlPerChr <- 10
selectProp <- 0.1
numOfGens <- 5
numOfRuns <- 1000 # Number of parallel runs

#wrapped simulation into function for parallelization with mclapply
run_simulation <- function(seed) {
  set.seed(seed)
  
  # Create founders
  founderPop <- runMacs(nInd = i, nChr = nChr, segSites = segSites, species = species)
  SP <- SimParam$new(founderPop)
  assign("SP", SP, envir = .GlobalEnv)
  
  # Create single trait
  SP$addTraitA(nQtlPerChr = qtlPerChr, mean = 0, var = 1)
  founders <- newPop(founderPop)
  
  # Set heritability
  SP$setVarE(h2 = heritability)
  
  # Create all pairwise crosses
  pairs <- t(combn(1:i, 2))
  
  # Initialize F1 population
  F1 <- NULL
  for (k in 1:nrow(pairs)) {
    cross <- makeCross(pop = founders, crossPlan = matrix(pairs[k, ], ncol = 2), nProgeny = 100)
    F1 <- if (is.null(F1)) cross else c(F1, cross)
  }
  
  # variables for old global selection method
  current_pop_pheno_global <- F1
  current_pop_geno_global <- F1
  
  # variables for new cluster-based selection
  # will be initialized with real values after selection within clusters
  selected_cluster_pheno <- NULL
  selected_cluster_geno <- NULL
  
  # Perform selection within each cross
  for (k in 1:nrow(pairs)) {
    cross <- makeCross(pop = founders, crossPlan = matrix(pairs[k, ], ncol = 2), nProgeny = 100)
    
    # select top 10% from each cross
    nSelect_per_cluster <- ceiling(0.1 * cross@nInd)
    
    # phenotypic selection within cluster
    selected_pheno <- selectInd(cross, nInd = nSelect_per_cluster, trait = 1, use = "pheno")
    selected_cluster_pheno <- if(is.null(selected_cluster_pheno)) selected_pheno else c(selected_cluster_pheno, selected_pheno)
    
    # genetic selection within cluster
    selected_geno <- selectInd(cross, nInd = nSelect_per_cluster, trait = 1, use = "gv")
    selected_cluster_geno <- if(is.null(selected_cluster_geno)) selected_geno else c(selected_cluster_geno, selected_geno)
  }
  
  #Initialize results tracking with F1 gen
  results_comparison <- rbind(
    data.frame(
      Generation = 1,
      Method = "Global_Pheno",
      Num_of_Ind = current_pop_pheno_global@nInd,
      Genetic_Gain = meanG(current_pop_pheno_global),
      Genetic_Variance = varG(current_pop_pheno_global)[1]
    ),
    data.frame(
      Generation = 1,
      Method = "Global_Geno",
      Num_of_Ind = current_pop_geno_global@nInd,
      Genetic_Gain = meanG(current_pop_geno_global),
      Genetic_Variance = varG(current_pop_geno_global)[1]
    ),
    data.frame(
      Generation = 1,
      Method = "Cluster_Pheno",
      Num_of_Ind = selected_cluster_pheno@nInd,
      Genetic_Gain = meanG(selected_cluster_pheno),
      Genetic_Variance = varG(selected_cluster_pheno)[1]
    ),
    data.frame(
      Generation = 1,
      Method = "Cluster_Geno",
      Num_of_Ind = selected_cluster_geno@nInd,
      Genetic_Gain = meanG(selected_cluster_geno),
      Genetic_Variance = varG(selected_cluster_geno)[1]
    )
  )
  
  #continue selection for subsequent generations from F2
  for (gen in 2:numOfGens) {
    #Global phenotypic selection
    current_pop_pheno_global <- self(current_pop_pheno_global)
    nSelect_global <- ceiling(selectProp * current_pop_pheno_global@nInd)
    current_pop_pheno_global <- selectInd(current_pop_pheno_global, nInd = nSelect_global, trait = 1, use = "pheno")
    
    #Global genetic selection
    current_pop_geno_global <- self(current_pop_geno_global)
    nSelect_global <- ceiling(selectProp * current_pop_geno_global@nInd)
    current_pop_geno_global <- selectInd(current_pop_geno_global, nInd = nSelect_global, trait = 1, use = "gv")
    
    #Cluster phenotypic selection
    selected_cluster_pheno <- self(selected_cluster_pheno)
    nSelect_cluster <- ceiling(selectProp * selected_cluster_pheno@nInd)
    selected_cluster_pheno <- selectInd(selected_cluster_pheno, nInd = nSelect_cluster, trait = 1, use = "pheno")
    
    #Cluster genetic selection
    selected_cluster_geno <- self(selected_cluster_geno)
    nSelect_cluster <- ceiling(selectProp * selected_cluster_geno@nInd)
    selected_cluster_geno <- selectInd(selected_cluster_geno, nInd = nSelect_cluster, trait = 1, use = "gv")
    
    #track values, adding new results as rows to previously recorded observations
    results_comparison <- rbind(
      results_comparison,
      data.frame(
        Generation = gen,
        Method = "Global_Pheno",
        Num_of_Ind = current_pop_pheno_global@nInd,
        Genetic_Gain = meanG(current_pop_pheno_global),
        Genetic_Variance = varG(current_pop_pheno_global)[1]
      ),
      data.frame(
        Generation = gen,
        Method = "Global_Geno",
        Num_of_Ind = current_pop_geno_global@nInd,
        Genetic_Gain = meanG(current_pop_geno_global),
        Genetic_Variance = varG(current_pop_geno_global)[1]
      ),
      data.frame(
        Generation = gen,
        Method = "Cluster_Pheno",
        Num_of_Ind = selected_cluster_pheno@nInd,
        Genetic_Gain = meanG(selected_cluster_pheno),
        Genetic_Variance = varG(selected_cluster_pheno)[1]
      ),
      data.frame(
        Generation = gen,
        Method = "Cluster_Geno",
        Num_of_Ind = selected_cluster_geno@nInd,
        Genetic_Gain = meanG(selected_cluster_geno),
        Genetic_Variance = varG(selected_cluster_geno)[1]
      )
    )
  }
  
  return(results_comparison)
}

seeds <- 1:numOfRuns
simulation_results <- mclapply(seeds, run_simulation, mc.cores = detectCores() - 1)

# When mclapply(...) is run from within system.time with nRuns = 10000, it seems to 
# take forever to complete. Use the red stop button or Esc to abort
##### Benchmarking Start ######
#Run simulations with mclapply
#benchmark_results <- system.time({
#  seeds <- 1:numOfRuns
#  simulation_results <- mclapply(seeds, run_simulation, mc.cores = detectCores() - 1)
#})
#### End ####

#Simulation run time in minutes
#of major importance is the elapsed i.e wall clock time
#print(benchmark_results/60)

#Combine all results from the simulation
all_results <- do.call(rbind, simulation_results)

# Calculate means and standard deviations for boundaries
summary_results <- aggregate(
  cbind(Genetic_Gain, Genetic_Variance) ~ Generation + Method,
  data = all_results,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)

# some elements of summary_results are lists of lists
# doing a simple as.data.frame would not work
summary_results <- as.data.frame(sapply(summary_results, unlist))
# summary_results <- do.call(data.frame, summary_results) # more elegant
names(summary_results)[3:6] <- c("Genetic_Gain", "Genetic_Gain_SD", "Genetic_Variance", "Genetic_Variance_SD")

# Plot with grey boundaries for genetic gain
ggplot(summary_results, aes(x = Generation, y = Genetic_Gain, color = Method, fill = Method)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = Genetic_Gain - Genetic_Gain_SD, ymax = Genetic_Gain + Genetic_Gain_SD), alpha = 0.2, color = NA) +
  theme_classic() +
  labs(title = "Genetic Gain - Avg Across All Runs with Boundaries",
       y = "Genetic Gain",
       x = "Generation")

# average genetic gain comparison across all simulations -- no grey shade
# ggplot(summary_results, aes(x = Generation, y = Genetic_Gain, color = Method)) +
#  geom_line() +
#  geom_point() +
#  theme_classic() +
#  labs(title = "Genetic Gain - Avg Across All Runs",
#       y = "Genetic Gain",
#       x = "Generation")

# average genetic variance comparison across all simulations -- no grey shade
#ggplot(average_results, aes(x = Generation, y = Genetic_Variance, color = Method)) +
#  geom_line() +
#  geom_point() +
#  theme_classic() +
#  labs(title = "Comparison of Genetic Variance - Average Across All Runs",
#       y = "Genetic Variance",
#       x = "Generation")

# Print stats
print(aggregate(Genetic_Gain ~ Generation + Method, data = summary_results, 
                FUN = function(x) round(mean(x), 3)))

# values for boundaries
boundary_vals <- all_results %>%
  group_by(Method) %>%
  summarize(
    mean_gain = mean(Genetic_Gain),
    sd_gain = sd(Genetic_Gain),
    upper_bound = mean_gain + sd_gain,
    lower_bound = mean_gain - sd_gain,
    .groups = "drop"
  )
print(boundary_vals)