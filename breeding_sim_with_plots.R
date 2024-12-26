# Load AlphaSimR
library(AlphaSimR)

# Set Simulation Parameters
i <- 10 # number of founders
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

#all pairwise crosses
pairs <- t(combn(1:i, 2))

# Initialize F1 population
F1 <- NULL
for (k in 1:nrow(pairs)) {
  cross <- makeCross(pop = founders, crossPlan = matrix(pairs[k, ], ncol = 2), nProgeny = 100)
  F1 <- if (is.null(F1)) cross else c(F1, cross)
}

# Initialize populations for pheno and geno selection methods
current_pop_pheno <- F1
current_pop_geno <- F1

# tracking variables for phenotypic and genetic selection
results_pheno <- data.frame(Generation = 1, Num_of_Ind = F1@nInd,
                            Genetic_Gain = meanG(current_pop_pheno),
                            Phenotypic_Gain = meanP(current_pop_pheno),
                            Genetic_Variance = varG(current_pop_pheno)[1],
                            Phenotypic_Variance = varP(current_pop_pheno)[1])

#NB: same starting point as above since both are starting out with the same F1 pop
results_gv <- results_pheno 

#repeat selection for F2 to F5
for (gen in 2:numOfGens) {
  # self both pops
  current_pop_pheno <- self(current_pop_pheno)
  current_pop_geno <- self(current_pop_geno)
  
  # Pheno selection: select top 10%
  nSelect_pheno <- ceiling(selectProp * current_pop_pheno@nInd)
  current_pop_pheno <- selectInd(current_pop_pheno, nInd = nSelect_pheno, trait = 1, use = "pheno")
  
  # Geno selection: select top 10%
  nSelect_geno <- ceiling(selectProp * current_pop_geno@nInd)
  current_pop_geno <- selectInd(current_pop_geno, nInd = nSelect_geno, trait = 1, use = "gv")
  
  #track phenotypic selection results
  results_pheno <- rbind(results_pheno, data.frame(
    Generation = gen,
    Num_of_Ind = nSelect_pheno,
    Genetic_Gain = meanG(current_pop_pheno),
    Phenotypic_Gain = meanP(current_pop_pheno),
    Genetic_Variance = varG(current_pop_pheno)[1],
    Phenotypic_Variance = varP(current_pop_pheno)[1]
  ))
  
  #track genetic selection results
  results_gv <- rbind(results_gv, data.frame(
    Generation = gen,
    Num_of_Ind = nSelect_geno,
    Genetic_Gain = meanG(current_pop_geno),
    Phenotypic_Gain = meanP(current_pop_geno),
    Genetic_Variance = varG(current_pop_geno)[1],
    Phenotypic_Variance = varP(current_pop_geno)[1]
  ))
}

#Geno and Pheno gain across generations using pheno trait values
plot(results_pheno$Generation, results_pheno$Genetic_Gain, type = "o", col = "red", ylim = range(c(results_pheno$Genetic_Gain, results_pheno$Phenotypic_Gain)),
     xlab = "Generation", ylab = "Mean Value", main = "Genetic and Phenotypic Gain Across Generations")
lines(results_pheno$Generation, results_pheno$Phenotypic_Gain, type = "o", col = "blue")
legend("topleft", cex=0.5, legend = c("Genetic Gain", "Phenotypic Gain"), col = c("red", "blue"), lty = 1, pch = 1)

#Genetic and phenotypic variance across generations also using pheno trait values
plot(results_pheno$Generation, results_pheno$Genetic_Variance, type = "o", col = "red", ylim = range(c(results_pheno$Genetic_Variance, results_pheno$Phenotypic_Variance)),
     xlab = "Generation", ylab = "Variance", main = "Genetic and Phenotypic Variance")
lines(results_pheno$Generation, results_pheno$Phenotypic_Variance, type = "o", col = "blue")
legend("topright", cex=0.5, legend = c("Genetic Variance", "Phenotypic Variance"), col = c("red", "blue"), lty = 1, pch = 1)


#Comparison of genetic gain from pheno and geno selection USING GENETIC VALUES (GV)
##Hypothesis tested: Genetic selection results in faster genetic gain compared
##to phenotypic selection at a stringent selection intensity.
plot(results_gv$Generation, results_gv$Genetic_Gain, type = "o", col = "red",
     ylim = range(c(results_gv$Phenotypic_Gain, results_gv$Genetic_Gain)),
     xlab = "Generation", ylab = "Gain", main = "Geno vs. Pheno Selection")
lines(results_gv$Generation, results_gv$Phenotypic_Gain, type = "o", col = "blue")
legend("topleft", cex=0.5, legend = c("Geno Selection", "Pheno Selection"), col = c("red", "blue"), lty = 1, pch = 1)


# Sombined results for display
results_df <- data.frame(
  Generation = 1:numOfGens,
  Num_of_Ind = results_pheno$Num_of_Ind,
  Genetic_Gain_Pheno = results_pheno$Genetic_Gain,
  Phenotypic_Gain_Pheno = results_pheno$Phenotypic_Gain,
  Genetic_Gain_Geno = results_gv$Genetic_Gain,
  Phenotypic_Gain_Geno = results_gv$Phenotypic_Gain
)

print(results_df)
names(results_df)
#NB: Genetic_gain_geno refers to mean genetic values obtained using the "gv" selection criterion in selectInd()
#Likewise Phenotypic_gain_geno refers to mean trait values obtained using the "pheno" selection criterion
