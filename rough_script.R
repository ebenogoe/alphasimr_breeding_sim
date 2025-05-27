library(AlphaSimR)
library(e1071)

indis = 100
loci = 10
trait_mean = 10
h2 = 0.2

founders = runMacs(indis, nChr = 2, segSites = 100, species = "MAIZE")

SP <- SimParam$new(founders)
SP$addTraitA(loci, trait_mean, var = 5)
SP$setVarE(h2)

Parents = newPop(founders)
F1 = randCross(Parents, 5, 20)
F1 = setPheno(F1)
F2 = self(F1)

F3 <- selectWithinFam(F2, nInd = 10, use = "pheno", selectTop = TRUE)
F3 <- setPheno(F3)

## Rapid cycling
for (i in 1:8) {
  F3Sel <- selectInd(F3, nInd = 6)
  F3SelCross <- randCross(F3Sel, 10)
  F3 <- self(F3SelCross)
}

trainingPopPheno <- pheno(F3)
trainingPopGeno <- pullSegSiteGeno(F3)

trainingData <- as.data.frame(cbind(trainingPopPheno, trainingPopGeno))
colnames(trainingData) <- paste0("ID", 1:(ncol(trainingPopGeno) + ncol(trainingPopPheno)))

model <- svm(ID1 ~ ., data = trainingData, kernel = "radial", scale = FALSE)


genoF3 <- pullSegSiteGeno(F3)
colnames(genoF3) <- paste0("ID", 2:(ncol(genoF3) + 1))
EBVF3 <- as.numeric(predict(model, genoF3))
F3@ebv <- as.matrix(EBVF3)

gv(F3)
ebv(F3)

as.vector(cor(gv(F3), ebv(F3))) * 100

F4 <- selectInd(F3, 5, use = "ebv")
F4 <- self(F4)

genoF4 <- pullSegSiteGeno(F4)
colnames(genoF4) <- paste0("ID", 2:(ncol(genoF4) + 1))
EBVF4 <- as.numeric(predict(model, genoF4))
F4@ebv <- as.matrix(EBVF4)
