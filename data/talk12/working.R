
## The following R packages are required：
## tidyverse, ROCR, randomForestSRC, pROC, progress, tictoc, doSNOW, ggplot2, caret and magrittr

## source rf.function.R
source("./rf.function.PPI.R");

## data
## --------------------------------------------
## Qin
qc.fail.run<- read.csv("./data/qc.fail.run.csv")
## discovery, species
discovery.featdata <- read.csv("./data/PRJEB6337.LC.s.discovery.txt",sep = "\t",header = T,row.names = 1)
discovery.metadata <- read.csv("./data/PRJEB6337.LC.s.discovery.metadata.txt",sep = "\t",row.names = 1)
## external validation, species
validation.featdata<-read.csv("./data/PRJEB6337.LC.s.validation.txt",sep = "\t",header = T,row.names = 1)
validation.metadata<-read.csv("./data/PRJEB6337.LC.s.validation.metadata.txt",sep = "\t",header = T,row.names = 1)
## geuns
genus.data <-read.csv(file = "./data/Qin.genus.data.csv",row.names = 1)

## Iebba
Iebba.genus<-read.csv("./data/V4.rel.all.genus.txt",sep = "\t",header = T,row.names = 1)
Iebba.genus.metadata <- read.csv("./data/samples2run2",sep = "\t",row.names = 1)

## Loomba's features
loomba.features <- read.csv(file = "./data/loomba.features.final.txt")
###
discovery.featdata.loomba <- discovery.featdata %>% dplyr::select(loomba.features$features)
discovery.featdata.loomba <-discovery.featdata.loomba[!rownames(discovery.featdata.loomba) %in% qc.fail.run$run.id,]
validation.featdata.loomba <- validation.featdata %>% dplyr::select(loomba.features$features)
validation.featdata.loomba <-validation.featdata.loomba[!rownames(validation.featdata.loomba) %in% qc.fail.run$run.id,]
qc.dis.metadata <- discovery.metadata[!rownames(discovery.metadata) %in% qc.fail.run$run.id,]
qc.val.metadata <- validation.metadata[!rownames(validation.metadata) %in% qc.fail.run$run.id,]

## ----------------------------------------------------------------------------------------
## modeling
## 19 species
## --------------------------------------------
discovery <- rf.setdata(discovery.featdata.loomba,qc.dis.metadata,grouping = "Group",control = "Healthy")
rf.discovery <- rf.train(discovery, parallel = T, num.cores = 20, class.weights = T)
rf.discovery.eval <- rf.evaluate(rf.discovery)
#Overall performance (AUC):  0.93 (95%CI: 0.89-0.97)
rf.ext.19<-rf.external.validation(rf.discovery,validation.featdata.loomba,qc.val.metadata)
#Overall performance (AUC):  0.92 (95%CI: 0.86-0.98)

## 4 species
## --------------------------------------------
discovery.4.features <- discovery.featdata.loomba %>% dplyr::select(c("Veillonella.atypica","Veillonella.parvula","Streptococcus.parasanguinis","Streptococcus.salivarius"))
discovery.4.data <- rf.setdata(discovery.4.features,qc.dis.metadata,grouping = "Group",control = "Healthy")
rf.discovery.4 <- rf.train(discovery.4.data,parallel = T, num.cores = 20, class.weights = T)
rf.discovery.4.eval <- rf.evaluate(rf.discovery.4)
#Overall performance (AUC):  0.91 (95%CI: 0.87-0.95)
rf.ext.4<-rf.external.validation(rf.discovery.4,validation.featdata.loomba,qc.val.metadata)
#Overall performance (AUC):  0.89 (95%CI: 0.81-0.96)

## 2 genus
## Qin
qc.genus.data <- genus.data[which(!row.names(genus.data) %in% qc.fail.run$run.id),]
qc.2.genus.data <- qc.genus.data %>% dplyr::select(c("Veillonella","Streptococcus"))
all.metadata <- rbind(qc.dis.metadata,qc.val.metadata)
guens.2 <- rf.setdata(qc.2.genus.data,all.metadata,grouping = "Group",control = "Healthy")
#In total, the feature data contains 302 samples, with 2 features.
#Case group: 'LC', w/ 160 samples, 
#Control group: 'Healthy', w/ 142 samples,rf.all.5to19 <- rf.train(all.5to19.data,parallel = T, num.cores = 20, class.weights = T)
rf.guens.2 <- rf.train(guens.2,parallel = T, num.cores = 20, class.weights = T)
rf.guens.2.eval <- rf.evaluate(rf.guens.2)
#Overall performance (AUC): Overall performance (AUC):  0.9 (95%CI: 0.86-0.93) 
## Iebba
Iebba.genus<-Iebba.genus %>% dplyr::select(c("Veillonella","Streptococcus"))
Iebba.genus<-Iebba.genus*100 ## (0-1) to (0-100)
Iebba.genus.metadata$Group<-str_replace_all(Iebba.genus.metadata$Group,"C","Healthy")
Iebba.genus.metadata$Group<-str_replace_all(Iebba.genus.metadata$Group,"D","LC")
Iebba<-rf.setdata(Iebba.genus,Iebba.genus.metadata,grouping = "Group",control = "Healthy")
rf.Iebba<-rf.train(Iebba,parallel = T, num.cores = 20, class.weights = T)
rf.Iebba.eval <- rf.evaluate(rf.Iebba)
#Overall performance (AUC):  0.85 (95%CI: 0.74-0.96)

Iebba.val<-rf.external.validation(rf.guens.2,Iebba.genus,Iebba.genus.metadata)
#Overall performance (AUC):  Overall performance (AUC):  0.94 (95%CI: 0.89-1)

qin.val<-rf.external.validation(rf.Iebba,qc.2.genus.data,all.metadata)
#Overall performance (AUC): Overall performance (AUC):  0.89 (95%CI: 0.85-0.93) 
