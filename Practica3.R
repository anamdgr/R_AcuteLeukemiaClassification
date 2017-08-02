workingDir <- "C:/Users/usuario/Desktop/Curso14-15/Aprendizaje_Computacional/Practica3/ALL-AML_Leukemia/ALL-AML_Leukemia"
setwd(workingDir)
#importamos el fichero de datos para entrenamiento
leukemia.train <- read.csv("AMLALL_train.data", header=FALSE)
#guardamos en un vector la clase (AML o ALL) a la que pertenece cada muestra (estaba en la última columna)
leukemia.classes <- leukemia.train[,7130]
summary(leukemia.classes)
#ALL:27, AML:11
#eliminamos esa última columna de los datos
leukemia.train <- leukemia.train[,-7130]
#transponemos la matriz
leukemia.train <- t(leukemia.train)
#para poner en row.names los nombres de los genes, importamos el fichero y copiamos los datos
gene.names <- read.table("human_genes.txt", header=FALSE)
row.names(leukemia.train) <- gene.names[,1]

#importamos el fichero de datos para el test y guardamos los datos de igual forma que con el de train
leukemia.test <- read.csv("AMLALL_test.data", header=FALSE)
leukemia.test.classes <- leukemia.test[,7130]
summary(leukemia.test.classes)
#ALL:20, AML:14
leukemia.test <- leukemia.test[,-7130]
leukemia.test <- t(leukemia.test)
row.names(leukemia.test) <- gene.names[,1]

#unimos los datos de train+test
leukemia.alldata <- cbind(leukemia.train, leukemia.test)
#unimos también el vector de clases
leukemia.alldata.classes <- as.factor(c(as.character(leukemia.classes), as.character(leukemia.test.classes)))
summary(leukemia.alldata.classes)
#ALL:47, AML:25

#primer filtro sencillo para eliminar los genes que casi no varían en su expresión
library(genefilter)
#(con p-value = 0.1)
filter <- ttest(leukemia.alldata.classes, p=0.1)
leukemia.filtered <- genefilter(leukemia.alldata, filterfun(filter))
#como podemos ver, se han seleccionado 2325 genes con el filtro
sum(leukemia.filtered)
#(con p-value = 0.05)
filter1 <- ttest(leukemia.alldata.classes, p=0.05)
leukemia.filtered1 <- genefilter(leukemia.alldata, filterfun(filter1))
#como podemos ver, se han seleccionado 2325 genes con el filtro
sum(leukemia.filtered1)
#guardamos el conjunto de 2654 genes
leukemia.2654 <- leukemia.alldata[which(leukemia.filtered),]
#guardamos el conjunto de 2046 genes
leukemia.2046 <- leukemia.alldata[which(leukemia.filtered1),]

#selección de genes con algortimos genéticos
library(galgo)
#la siguiente función crea y configura todos los objetos necesarios para un problema de selección de variables
#para clasificación
leukemia.bb <- configBB.VarSel(data=leukemia.2046, classes=leukemia.alldata.classes, 
                                      force.train=c(1:38), force.test=c(39:72), 
                                      classification.method="nearcent", chromosomeSize=5, 
                                      maxSolutions=100, goalFitness=0.90, 
                                      saveVariable="leukemia.bb", saveFrequency=30, 
                                      saveFile="leukemia_bb.Rdata")
blast(leukemia.bb)
loadObject("leukemia_bb.Rdata")
plot(leukemia.bb, type="fitness")
plot(leukemia.bb, type="confusion")

library(caret)
library(doParallel)
#damos formato a los datos para el uso de la librería
leukemia.train.caret <- t(leukemia.766)
ga_ctrl <- gafsControl(functions = rfGA, # Assess fitness with RF
                       method = "cv",    # 10 fold cross validation
                       genParallel = TRUE, # Use parallel programming
                       allowParallel = TRUE,
                       verbose = TRUE)
lev <- levels(leukemia.classes)     # Set the levels
rf_ga3 <- gafs(x = leukemia.train.caret, y = leukemia.classes,
               iters = 10, # 100 generations of algorithm
               popSize = 20, # population size for each generation
               levels = lev,
               gafsControl = ga_ctrl)
rf_ga3
##demasiado tiempo de ejecución!!!!!
#ni utilizando 4 núcleos:
#registerDoParallel(4)
#getDoParWorkers()

#500 genes, 500 iteraciones, 20-100 tamaño de población
#library GA
#diseñar la función de fitness


#escalamos los datos
leukemia.train.scaled <- leukemia.train
#filtro sencillo para eliminar los genes que casi no varían en su expresión
library(genefilter)
#filter0.01 <- ttest(leukemia.classes, p=0.1)
#leukemia.train.filtered <- genefilter(leukemia.train.scaled, filterfun(filter0.01))
#sum(leukemia.train.filtered)
#guardamos el conjunto filtrado, de 774 genes
#leukemia.774 <- leukemia.train.scaled[which(leukemia.train.filtered),]
#leukemia.774 <- t(leukemia.774)
#leukemia.774 <- data.frame(leukemia.774, Class=as.vector(leukemia.classes))


cor.leukemia <- cor(t(leukemia.train.scaled))
hc <- findCorrelation(cor.leukemia, cutoff=0.7)
hc <- sort(hc)
reduced.leukemia <- t(leukemia.train.scaled)[, -c(hc)]

filter0.07 <- ttest(leukemia.classes, p=0.07)
leukemia.t.filtered <- genefilter(t(reduced.leukemia), filterfun(filter0.07))
sum(leukemia.t.filtered)

class.asNumeric <- as.character(leukemia.classes)
class.asNumeric[class.asNumeric == "AML"] <- 1
class.asNumeric[class.asNumeric == "ALL"] <- 0
class.numeric <- as.numeric(class.asNumeric, stringAsFactor = FALSE)

leukemia.654 <- leukemia.train.scaled[which(leukemia.t.filtered),]
#matriz transpuesta
leukemia.654 <- t(leukemia.654)
#escalamos los datos
leukemia.654 <- scale(leukemia.654)
#convertimos en dataframe
leukemia.654 <- data.frame(leukemia.654, Class=class.numeric)
#comprobamos las clases (ALL=0, AML=1)
leukemia.654$Class


library(GA)
library(MASS)
library(mclust)

#definición de la función de fitness basada en regresión logística
fitness.glm <- function(chromosome) {
  forml <- as.formula(paste("Class~", paste(colnames(leukemia.654[,which(chromosome==1)]), collapse="+"),sep=""))
  regLog <- glm(forml, data=leukemia.654, family=binomial)
  summary(regLog)
  pred <- predict(regLog, leukemia.654, type="response")
  pred.th <- pred
  pred.th[pred.th<0.5]<-0
  pred.th[pred.th>=0.5]<-1
  error <- classError(leukemia.654$Class, pred.th)$errorRate
  result <- - error - sum(chromosome==1)/654
  return(result)
}

#para probar la función fitness con un cromosoma aleatorio:
#chromosome <- sample(c(0,1), 654, replace=TRUE, prob=c(0.9, 0.1))
#fitness.glm(chromosome)

#Algoritmo Genético utilizando la función fitness definida
geneticAlg <- ga(type="real-valued", fitness=fitness.glm, nBits=605, maxiter=500)
summary(geneticAlg)
plot(geneticAlg)

#Lista de genes seleccionados:
selectedGenes <- colnames(leukemia.654[,which(geneticAlg@solution[1,]==1)])
selectedGenes
#49 genes


#visualizamos la expresión de los genes seleccionados en un heatmap
library(gplots)
selectedData <- leukemia.654[,which(geneticAlg@solution[1,]==1)]
color.map <- function(class) { if (class=="ALL") "#FF0000" else "#0000FF" }
classcolors <- unlist(lapply(leukemia.classes, color.map))
heatmap(t(as.matrix(selectedData)), Colv=NA, Rowv=NA, col=redgreen(75), ColSideColors=classcolors)


library(affycoretools)
plotPCA(as.matrix(selectedData), groups=leukemia.classes, groupnames=levels(leukemia.classes),legend=FALSE)
