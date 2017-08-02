
#########################
### DATOS DE LEUCEMIA ###
#########################

#directorio de trabajo
workingDir <- "C:/Users/usuario/Desktop/Curso14-15/Aprendizaje_Computacional/Practica3/ALL-AML_Leukemia/ALL-AML_Leukemia"
setwd(workingDir)
#importamos el fichero de datos para entrenamiento
leukemia.train <- read.csv("AMLALL_train.data", header=FALSE)
#guardamos en un vector la clase (AML o ALL) a la que pertenece cada muestra (estaba en la última columna)
leukemia.classes <- leukemia.train[,7130]
#comprobamos el número de pacientes de cada tipo
summary(leukemia.classes)
#ALL:27, AML:11
#eliminamos esa última columna de los datos
leukemia.train <- leukemia.train[,-7130]
#leemos el fichero con el nombre de los genes
gene.names <- read.table("human_genes.txt", header=FALSE)
#ponemos los nombre de los genes como nombres de columnas
colnames(leukemia.train) <- gene.names[,1]


##################################
### NORMALIZACIÓN DE LOS DATOS ###
##################################

leukemia.train.norm <- apply(leukemia.train, 2, Normalize <- function(dataSet){
  datanorm<-( (dataSet-min(dataSet)) / (max(dataSet)-min(dataSet)) )
})
#Guardamos la matriz en formato data frame
leukemia.train.norm.frame <- data.frame(leukemia.train.norm)

################################
### FILTRADO POR CORRELACIÓN ###
################################

#Eliminación genes muy relacionados entre sí (redundancias)
cor.leukemia <- cor(leukemia.train.norm.frame)
hc <- findCorrelation(cor.leukemia, cutoff=0.9)
hc <- sort(hc)
leukemia.train.norm.frame <- leukemia.train.norm.frame[, -c(hc)]

#save(leukemia.train.norm.frame, file="trainNormFrame.csv")


#Conversión de la clase a numérico (ALL=0, AML=1)
leukemia.classes.binary <- (leukemia.classes =="AML")*1

#cálculo de la correlación
correlation <- sapply(leukemia.train.norm.frame, function(i){
  cor(i,leukemia.classes.binary)
})

#Número de genes que queremos seleccionar
gene.number <- 250
#Ordenación por correlación (en valor absoluto) y selección de los 200 primeros
best.correlation <- sort(abs(correlation), decreasing = TRUE)[1:gene.number] 
#Nombres de los genes seleccionados
best.correlation.names <- names(best.correlation)
#Extracción de la matriz de datos de los genes seleccionados
leukemia.best.corr <- leukemia.train.norm.frame[,best.correlation.names]
#Añadimos la clase (en formato numérico) al data frame con los datos filtrados por correlación
leukemia.best.corr$Class <- leukemia.classes.binary


## Visualización la expresión de los genes seleccionados en un heatmap
library(gplots)
color.map <- function(class) { if (class=="ALL") "#FF0000" else "#0000FF" }
classcolors <- unlist(lapply(leukemia.classes, color.map))
heatmap(t(as.matrix(leukemia.best.corr)), Colv=NA, col=redgreen(75), ColSideColors=classcolors)

##########################
### SELECCION DE GENES ###
##########################

## Utilización de algoritmos genéticos para el problema de selección de genes para clasificación

library(GA)
library(MASS)
library(mclust)
library(caret)

## Definición de la FUNCIÓN DE FITNESS basada en regresión logística
fitness.glm <- function(chromosome) {
  forml <- as.formula(paste("Class~", paste(colnames(leukemia.best.corr[,which(chromosome==1)]), collapse="+"),sep=""))
  regLog <- glm(forml, data=leukemia.best.corr, family=binomial("logit"))
  pred <- predict(regLog, leukemia.best.corr, type="response")
  pred.th <- pred
  pred.th[pred.th<0.5]<-0
  pred.th[pred.th>=0.5]<-1
  #calculamos la matriz de confusión y obtenemos la accuracy, la sensibilidad y la especificidad
  confMatrix <- confusionMatrix(pred.th, leukemia.best.corr$Class)
  accuracy <- confMatrix$overall[1]
  #error <- classError(leukemia.best.corr$Class, pred.th)$errorRate
  #result <- - error - sum(chromosome)/200
  result <- accuracy - sum(chromosome)/gene.number
  return(result)
}

## Aplicación del Algoritmo Genético utilizando la función fitness definida
#200 genes, 50 iteraciones, 50 tamaño de población

suggest <- matrix(sample(0:1, 40 * 250, replace = TRUE, prob=c(0.75,0.25)), 40, 250)

suggest[1:10,c(2,5,9,123,175)] <- 1
suggest[11:20,c(6,9,103,108,126)] <- 1
suggest[21:30,c(2,5,112,123,126)] <- 1
suggest[31:40,c(6,103,108,112,175)] <- 1


geneticAlg <- ga(type = "binary", fitness = fitness.glm, nBits = gene.number, 
                 popSize = 40, elitism = base::max(1, round(40*0.05)), maxiter = 40,
                 suggestions = suggest, seed = 12345)
summary(geneticAlg)
sumGA <- summary(geneticAlg)
#save(sumGA, file="summaryGA.txt")
rm(sumGA)
load(file="summaryGA.txt")
dev.off()
plot(geneticAlg)

## Lista de genes seleccionados
selectedGenes <- colnames(leukemia.best.corr[,which(geneticAlg@solution[1,]==1)])
selectedGenes <- colnames(leukemia.best.corr[,which(sumGA$solution[1,]==1)])
selectedGenes
# ~ 45 genes seleccionados

## Datos de los genes seleccionados
selectedData <- leukemia.best.corr[,which(geneticAlg@solution[1,]==1)]



## Visualización la expresión de los genes seleccionados en un heatmap
library(gplots)
color.map <- function(class) { if (class=="ALL") "#FF0000" else "#0000FF" }
classcolors <- unlist(lapply(leukemia.classes, color.map))
heatmap(t(as.matrix(selectedData)), Colv=NA, col=redgreen(75), ColSideColors=classcolors)


## Visualización en detalle de genes diferencialmente expresados 

#LGALS3
par(mfrow=c(1,2))
LGALS3.Galectin3 <- selectedData[,grep("LGALS3", names(selectedData))]
stripchart(LGALS3.Galectin3 ~ leukemia.classes, method="jitter")
boxplot(LGALS3.Galectin3 ~ leukemia.classes)
dev.off()

histALL <- hist(LGALS3.Galectin3[leukemia.classes=="ALL"])
histAML <- hist(LGALS3.Galectin3[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "LGAS3 Expression in AML/ALL", xlab = "LGAS3 Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "LGAS3 Expression in AML/ALL", xlab = "LGAS3 Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)

#MPO
par(mfrow=c(1,2))
MPO.Myeloperoxidase <- selectedData[,grep("^MPO", names(selectedData))]
stripchart(MPO.Myeloperoxidase ~ leukemia.classes, method="jitter")
boxplot(MPO.Myeloperoxidase ~ leukemia.classes)
dev.off()

histALL <- hist(MPO.Myeloperoxidase[leukemia.classes=="ALL"])
histAML <- hist(MPO.Myeloperoxidase[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "MPO Expression in AML/ALL", xlab = "MPO Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "MPO Expression in AML/ALL", xlab = "MPO Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)

#ZNF
par(mfrow=c(1,2))
ZNF <- selectedData[,grep("^ZNF", names(selectedData))]
stripchart(ZNF ~ leukemia.classes, method="jitter")
boxplot(ZNF ~ leukemia.classes)
dev.off()

histALL <- hist(ZNF[leukemia.classes=="ALL"])
histAML <- hist(ZNF[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "Zinc finger protein Expression in AML/ALL", xlab = "ZNF Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "Zinc finger protein Expression in AML/ALL", xlab = "ZNF Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)

#MB-1 (nombrado en el paper)
par(mfrow=c(1,2))
MB.1 <- selectedData[,grep("MB.1", names(selectedData))]
stripchart(MB.1 ~ leukemia.classes, method="jitter")
boxplot(MB.1 ~ leukemia.classes)
dev.off()

histALL <- hist(MB.1[leukemia.classes=="ALL"])
histAML <- hist(MB.1[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "MB-1 Gene Expression in AML/ALL", xlab = "MB-1 Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "MB-1 Gene Expression in AML/ALL", xlab = "MB-1 Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)

#CD19
par(mfrow=c(1,2))
CD19.Antigen <- selectedData[,grep("CD19", names(selectedData))]
stripchart(CD19.Antigen ~ leukemia.classes, method="jitter")
boxplot(CD19.Antigen ~ leukemia.classes)
dev.off()

histALL <- hist(CD19.Antigen[leukemia.classes=="ALL"])
histAML <- hist(CD19.Antigen[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "CD19 Antigen Expression in AML/ALL", xlab = "CD19 Antigen Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "CD19 Antigen Expression in AML/ALL", xlab = "CD19 Antigen Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)

#Zyxin
if(length(grep("Zyxin", names(selectedData))) != 0){
  par(mfrow=c(1,2))
  Zyxin <- selectedData[,grep("Zyxin", names(selectedData))]
  stripchart(Zyxin ~ leukemia.classes, method="jitter")
  boxplot(Zyxin ~ leukemia.classes)
}

dev.off()

histALL <- hist(Zyxin[leukemia.classes=="ALL"])
histAML <- hist(Zyxin[leukemia.classes=="AML"])
plot(histALL, col=rgb(0,0,1,0.5), xlim=c(0,1), main = "Zyxin Expression in AML/ALL", xlab = "Zyxin Expression")
plot(histAML, col=rgb(1,0,0,0.5), xlim=c(0,1), main = "Zyxin Expression in AML/ALL", xlab = "Zyxin Expression", add=T)
legend('topright',c('ALL','AML'), fill = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), bty = 'o', box.col= "white", border = NA)



#####################
### DATOS DE TEST ###
#####################

#importamos el fichero de datos para el test y guardamos los datos de igual forma que con los de train
leukemia.test <- read.csv("AMLALL_test.data", header=FALSE)
leukemia.test.classes <- leukemia.test[,7130]
summary(leukemia.test.classes)
#ALL:20, AML:14
leukemia.test <- leukemia.test[,-7130]
colnames(leukemia.test) <- gene.names[,1]
#Conversión de la clase a numérico (ALL=0, AML=1)
leukemia.test.classes.binary <- (leukemia.test.classes =="AML")*1
#Normalización de los datos de test
leukemia.test.norm <- sapply(leukemia.test, Normalize <- function(dataSet){
  datanorm<-( (dataSet-min(dataSet)) / (max(dataSet)-min(dataSet)) )
})
#Guardamos la matriz en formato data frame con la clase en numérico
leukemia.test.norm.frame <- data.frame(leukemia.test, Class=leukemia.test.classes.binary)

selectedGenes.testExpression <- leukemia.test.norm.frame[,selectedGenes]

classcolorstest <- unlist(lapply(leukemia.test.classes, color.map))
heatmap(t(as.matrix(selectedGenes.testExpression)), Colv=NA, col=redgreen(75), ColSideColors=classcolorstest)


#####################
### CLASIFICACIÓN ###
#####################

## Clasificación mediante regresión logistica de los casos del grupo de TEST 
## utilizando para la predicción los genes seleccionados por el GA
forml <- as.formula(paste("Class~", paste(selectedGenes, collapse="+"),sep=""))
regLog <- glm(forml, data=leukemia.test.norm.frame, family=binomial)
summary(regLog)
pred <- predict(regLog, leukemia.test.norm.frame, type="response")
pred.th <- pred
pred.th[pred.th<0.5]<-0
pred.th[pred.th>=0.5]<-1
confMatrix <- confusionMatrix(pred.th,leukemia.test.norm.frame$Class)
confMatrix$table
accuracy <- confMatrix$overall[1]
accuracy


###########################
### COMPARACIÓN CON SFS ###
###########################


sfsFunction <- function(dataSet){
  
  genes <- colnames(dataSet)[-201]
  usedGenes <- vector()
  best.acc <- 0
  total.acc <- 0
  
  for(i in 1:length(genes)){
    forml <- as.formula(paste("Class~", genes[i], sep=""))
    regLog <- glm(forml, data=dataSet, family=binomial("logit"))
    pred <- predict(regLog, dataSet, type="response")
    pred.th <- pred
    pred.th[pred.th<0.5]<-0
    pred.th[pred.th>=0.5]<-1
    confMatrix <- confusionMatrix(pred.th, dataSet$Class)
    accuracy <- confMatrix$overall[1]
    if(accuracy > best.acc){
      best.acc <- accuracy
      best.gene <- genes[i]
    }
  }
  usedGenes <- c(usedGenes, best.gene)
  genes <- genes[-which(genes==best.gene)]
  
  while(length(genes)>0) {
    best.acc <- 0
    used <- paste(usedGenes, collapse="+")
    for(i in 1:length(genes)){
      forml <- as.formula(paste("Class~", paste(used, genes[i], sep="+"), sep=""))
      regLog <- glm(forml, data=dataSet, family=binomial("logit"))
      pred <- predict(regLog, dataSet, type="response")
      pred.th <- pred
      pred.th[pred.th<0.5]<-0
      pred.th[pred.th>=0.5]<-1
      confMatrix <- confusionMatrix(pred.th, dataSet$Class)
      accuracy <- confMatrix$overall[1]
      if(accuracy > best.acc){
        best.acc <- accuracy
        best.gene <- genes[i]
      }
    }
    usedGenes <- c(usedGenes, best.gene)
    genes <- genes[-which(genes==best.gene)]
    if(best.acc > total.acc){
      total.acc <- best.acc
      best.subset <- usedGenes
    }else{
      break
    }
  }
  res.list <- list(total.acc = total.acc, best.subset=best.subset)
  return(res.list) 
}  

sfsSubset <- sfsFunction(leukemia.best.corr)

sfsSubset


selectedData.class <- selectedData
selectedData.class$Class <- leukemia.best.corr$Class

library(mlbench)
library(caret)
# prepare training scheme
control.train <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Class~., data=leukemia.best.cor, method="knn", trControl=control.train)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
summary(importance)
# plot importance
plot(importance, top=40)

# define the control using a random forest selection function
control.rfe <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results.rfe <- rfe(leukemia.best.corr[1:200], leukemia.classes, sizes=c(1:200), rfeControl=control.rfe)
# summarize the results
print(results.rfe)
# list the chosen features
predictors(results.rfe)
# plot the results
plot(results.rfe, type=c("g", "o"))


library(randomForest)
# prepare training scheme
control.train <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Class~., data=leukemia.best.corr, method="knn", trControl=control.train)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# plot importance
pdf("VarImpPlot3.pdf", width = 13, height = 10)
plot(importance, top=40, mar=c(30,30,30,30))
dev.off()
