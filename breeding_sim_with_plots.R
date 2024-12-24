# Load AlphaSimR
library(AlphaSimR)

#Set Simulation Parameters
i <- 30 #number of founders
nChr <- 1
segSites <- 100 
species <- "MAIZE"
heritability <- 0.3
qtlPerChr <- 10
selectProp <- 0.1
numOfGens <- 5

#Make 10 founders
founderPop <- runMacs(nInd = i, nChr = nChr, segSites = segSites, species = species)
SP <- SimParam$new(founderPop)

#create single trait
SP$addTraitA(nQtlPerChr = qtlPerChr, mean = 0, var = 1)

founders <- newPop(founderPop)

#needed to prevent the error about selection trait having missing values
SP$setVarE(h2 = heritability)

#create crossing scheme
#all pairwise crosses
pairs <- t(combn(1:i, 2))

#Initialize F1 population
F1 <- NULL

for (i in 1:nrow(pairs)) {
  #create 100 F1 offspring for each cross
  #cross <- randCross(pop = founders[pairs[i, ]], nCrosses = 1, nProgeny = 100)
  cross <- makeCross(pop = founders, crossPlan = matrix(pairs[i, ], ncol = 2), nProgeny = 100)
  #Combine all crosses into a single population
  if (is.null(F1)) {
    F1 <- cross
  } else {
    F1 <- c(F1, cross)
  }
}

#Start with F1 population
current_pop <- F1
cat("Generation F1", ": Selected all ", current_pop@nInd, " individuals\n", sep = "")

num_individuals <- c(current_pop@nInd)

# Initialize vectors to track genetic/phenotypic gain and variance
genetic_gain <- c(meanG(current_pop))
phenotypic_gain <- c(meanP(current_pop))
genetic_variance <- c(varG(current_pop))
phenotypic_variance <- c(varP(current_pop))

#repeat for F2 though F5 
for (gen in 2:numOfGens) {
  #Self all F1s
  current_pop <- self(current_pop)
  
  #select top 10% of F1
  nSelect <- ceiling(selectProp * current_pop@nInd)
  current_pop <- selectInd(current_pop, nInd = nSelect, trait = 1, use = "pheno")
  
  num_individuals <- c(num_individuals, nSelect)
  
  # Record the mean genetic and phenotypic values and variances
  genetic_gain <- c(genetic_gain, meanG(current_pop))
  phenotypic_gain <- c(phenotypic_gain, meanP(current_pop))
  genetic_variance <- c(genetic_variance, varG(current_pop))
  phenotypic_variance <- c(phenotypic_variance, varP(current_pop))
  accuracy <- cor(gv(current_pop), pheno(current_pop))
  
  # Print progress
  cat("Generation F", gen, ": Selected ", nSelect, " individuals\n", sep = "")
  
}

# Store results in a dataframe
results_df <- data.frame(
  Generation = 1:numOfGens,
  Num_of_Ind = num_individuals,
  Genetic_Gain = genetic_gain,
  Phenotypic_Gain = phenotypic_gain,
  Genetic_Variance = genetic_variance,
  Phenotypic_Variance = phenotypic_variance
)

# Plot genetic and phenotypic gain
plot(1:numOfGens, genetic_gain, type = "o", col = "blue", ylim = range(c(genetic_gain, phenotypic_gain)),
     xlab = "Generation", ylab = "Mean Value", main = "Genetic and Phenotypic Gain Across Generations")
lines(1:numOfGens, phenotypic_gain, type = "o", col = "red")
legend("topleft", cex=0.5, legend = c("Genetic Gain", "Phenotypic Gain"), col = c("blue", "red"), lty = 1, pch = 1)

# Plot for variances (genetic and phenotypic variance)
plot(1:numOfGens, genetic_variance, type = "o", col = "darkgreen", ylim = range(c(genetic_variance, phenotypic_variance)),
     xlab = "Generation", ylab = "Variance", main = "Genetic and Phenotypic Variance")
lines(1:numOfGens, phenotypic_variance, type = "o", col = "orange")
legend("topright", cex=0.5, legend = c("Genetic Variance", "Phenotypic Variance"), col = c("darkgreen", "orange"), lty = 1, pch = 1)