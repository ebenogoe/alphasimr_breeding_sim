# Load needed libraries
library(AlphaSimR)
library(ggplot2)

# Set Simulation Parameters
i <- 10 # Number of founders
nChr <- 1
segSites <- 100 
species <- "MAIZE"
heritability <- 0.3
qtlPerChr <- 10
selectProp <- 0.1
numOfGens <- 5

# Create founders
founderPop <- runMacs(nInd = i, nChr = nChr, segSites = segSites, species = species)
SP <- SimParam$new(founderPop)

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

#Global selection from now onwards refers to the first method where selection was done
#after combining all progeny into one super cluster and selecting the top 10%
#whereas cluster selection refers to the new method -- selecting the top 10% from the 
#progeny of each cross

# variables for original global selection 
current_pop_pheno_global <- F1
current_pop_geno_global <- F1

# variables for new cluster-based selection
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

# track the neccessary values for all combinations for F1 pop
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

# Continue selection for subsequent generations
for (gen in 2:numOfGens) {
  # global phenotypic selection i.e across all clusters using "pheno"
  current_pop_pheno_global <- self(current_pop_pheno_global)
  nSelect_global <- ceiling(selectProp * current_pop_pheno_global@nInd)
  current_pop_pheno_global <- selectInd(current_pop_pheno_global, nInd = nSelect_global, trait = 1, use = "pheno")
  
  # global genetic selection
  current_pop_geno_global <- self(current_pop_geno_global)
  nSelect_global <- ceiling(selectProp * current_pop_geno_global@nInd)
  current_pop_geno_global <- selectInd(current_pop_geno_global, nInd = nSelect_global, trait = 1, use = "gv")
  
  # cluster phenotypic selection i.e within each cluster using "pheno"
  selected_cluster_pheno <- self(selected_cluster_pheno)
  nSelect_cluster <- ceiling(selectProp * selected_cluster_pheno@nInd)
  selected_cluster_pheno <- selectInd(selected_cluster_pheno, nInd = nSelect_cluster, trait = 1, use = "pheno")
  
  # cluster genetic selection
  selected_cluster_geno <- self(selected_cluster_geno)
  nSelect_cluster <- ceiling(selectProp * selected_cluster_geno@nInd)
  selected_cluster_geno <- selectInd(selected_cluster_geno, nInd = nSelect_cluster, trait = 1, use = "gv")
  
  # track values for all generations defined by loop
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


# Genetic gain comparison
ggplot(results_comparison, aes(x = Generation, y = Genetic_Gain, color = Method)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(title = "Comparison of Genetic Gain - All Combinations",
       y = "Genetic Gain",
       x = "Generation")

#Genetic variance comparison
ggplot(results_comparison, aes(x = Generation, y = Genetic_Variance, color = Method)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(title = "Comparison of Genetic Variance Across Selection Methods",
       y = "Genetic Variance",
       x = "Generation")

# Seems to assign uninformative/non-existing traits as rownames. Remove with the following:
rownames(results_comparison) <- NULL

print(aggregate(Genetic_Gain ~ Generation + Method, data = results_comparison, 
                FUN = function(x) round(mean(x), 3)))
View(results_comparison)
