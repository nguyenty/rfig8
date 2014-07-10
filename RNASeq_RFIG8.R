pm0 <- proc.time()
# cleaning data
library(edgeR)
library(fields)
library(plyr)
library(reshape)
dir.source <- "U:/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
#dir.source <- "/home/ntyet/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
source(paste(dir.source, "QL.fit.R",sep=""))
source(paste(dir.source, "NBDev.R",sep =""))
source(paste(dir.source, "PoisDev.R",sep =""))
source(paste(dir.source, "QL.results.R",sep =""))


# datdir <- "/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data"
# datdir <- "/home/ntyet/cyfiles/R/RA/Data"
datdir <-  "U:/R/RA/Data"
ndata <- read.csv(
  paste(datdir,
        "/Reanalysis Data/g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv",
        sep =""))


dat <- read.table(paste(datdir,"/RFI_uniq_comb_count_corrected.txt", sep = ""))

meta.data <- read.csv(paste(datdir,"/RIN values RNAseq Obj3a Martine.csv",sep = ""))
nname <- paste("20800",sprintf("%04.0f",meta.data$Sample.Name), sep = "")
meta.data2 <- cbind(ID=as.numeric(nname),meta.data)

nndata <- ndata[ndata[,"idpig"]%in% nname,c("idpig","rfi.ADJUSTED") ]
for (i in 1:length(nname)){
  for(j in 1:length(nname)){
    if (meta.data2$ID[i]==nndata$idpig[j])
      meta.data2$RFI.value[i] <- nndata$rfi.ADJUSTED[j]
  }
}

write.table(meta.data2, paste(datdir,"/Meta.data2.txt",sep = ""), sep = "\t")
write.csv(meta.data2, paste(datdir,"/Meta.data2.csv",sep =""))

meta.data2[,c("ID","Sample.Name","RFI.value")]
meta.data[,c("Sample.Name","RFI.value")]
name <- paste("X",meta.data$Sample.Name, sep = "")
dat2 <- dat[, name]

Lane <- as.factor(meta.data$Lane.RNAseq)
Diet <- as.factor(meta.data$diet)
Line <- as.factor(meta.data$line)
RFI <- meta.data2$RFI.value
RINb <- meta.data$RIN.before.GD
RINa  <- meta.data$RIN.after.GD
Conc <- meta.data$Conc.after.GD.ng.ul.
dateGD <- meta.data$date.GD
dateRNA <- meta.data$date.RNA.extraction
variable_name <- c("Lane", "Diet", "Line", "RFI", 
                   "RINb", "RINa", 
                   "Conc", "dateGD", "dateRNA")


counts <- as.matrix(dat2[rowSums(dat2>0)>3&
                           rowMeans(dat2)>1,])
log.offset <- log(apply(counts, 2, quantile, .75))

# function to compute the constant term in pdf of NegBin with 
# paramters mean \mu, observation y and dispersion parameter 1/disp
log.gamma <- function(counts, disp){
  log.g <- NULL
  n <- length(counts)
  for(i in 1:n){
    log.g[i] <- sum(log(0:(counts[i]-1)+disp)) - sum(log(1:counts[i]) )
  }
  return(log.g)  
}

SAT.LIKE2<-function(count,disp){
  means<-count
  like<-disp*log(disp/(disp+means))
  like[count!=0]<-like[count!=0]+count[count!=0]*log(means[count!=0]/(disp+means[count!=0]))+
    log.gamma(count[count!=0],disp )
  sum(like)
}

## Function to calculate AIC of the QL.fit model 

AIC.QL <- function(counts,QL.fit.object){
  n <- dim(counts)[2]
  m <- dim(counts)[1]
  disp <- 1/QL.fit.object$NB.disp
  den.df <- QL.fit.object$den.df
  phi.hat.dev <- QL.fit.object$phi.hat.dev
  p <- n - den.df
  dev <- phi.hat.dev*den.df
  L0 <- NULL
  for (i in 1:m){
    L0[i] <- SAT.LIKE2(counts[i,],disp[i])
  }
  
  return(dev-2*L0+2*p)
}
# colnames(counts)



## Function do the model list and test.mat 

list_model <- function(full_model){
  n <- dim(full_model)[2]
  variable_name <- colnames(full_model)
  variable_name <- gsub(":", "", variable_name)
  variable_name <- gsub("2", "", variable_name)
  variable_name <- gsub("7", "", variable_name)
  test.mat <- NULL
  design.list <- vector("list", n)
  design.list[[1]] <- full_model
  for (i in 2:n) {
    design.list[[i]] <- as.matrix(full_model[,-i])
    test.mat <- rbind(test.mat, c(1,i))
  }
  
  row.names(test.mat) <-  variable_name[-1]
  
  if (any(variable_name == "dateGD1/13/01")){
    ind_dateGD <- which(variable_name == "dateGD1/13/01")
    design.list <- vector("list", n-3)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(n-4), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 2:(ind_dateGD-1)){
      design.list[[i]] <- as.matrix(full_model[,-i])
      row.names(test.mat)[i-1] <- variable_name[i]
    }
    
    design.list[[ind_dateGD]] <- full_model[,-(ind_dateGD:(ind_dateGD+3))]
    row.names(test.mat)[ind_dateGD -1] <- "dateGD"
    
    if((ind_dateGD+1)<=(n-3)){
      for(i in ((ind_dateGD+1):(n-3))){
        design.list[[i]] <- full_model[,-(i+3)]
        row.names(test.mat)[i-1] <- variable_name[i+3]    
      }
    }
  }
  
  if (any(variable_name == "dateRNA11/14/01")){
    ind_dateRNA <- which(variable_name == "dateRNA11/14/01")
    design.list <- vector("list", n-1)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(n-2), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 2:(ind_dateRNA-1)){
      design.list[[i]] <- as.matrix(full_model[,-i])
      row.names(test.mat)[i-1] <- variable_name[i]
    }
    
    design.list[[ind_dateRNA]] <- full_model[,-(ind_dateRNA:(ind_dateRNA+1))]
    row.names(test.mat)[ind_dateRNA -1] <- "dateRNA"
    
     if ((ind_dateRNA+1)<=(n-1)) {for(i in ((ind_dateRNA+1):(n-1))){
      design.list[[i]] <- full_model[,-(i+1)]
      row.names(test.mat)[i-1] <- variable_name[i+1]    
    }
    }
  }
  if (n ==2) design.list[[2]] <- rep(1, nrow(full_model))
  return(list(design.list = design.list, test.mat = test.mat))
}

  
## Function do all the things with input Full model

fit_model <- function(full_model, model_th){
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, 
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  k <- nrow(test.mat)
  name_model <- NULL
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
  dir.create(model_dir, showWarnings = FALSE)
  save(result, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
}



# Model 1: 
#variable_name
m <- 1
full_model <- model.matrix(~Line*Diet*RFI + Conc + RINa + RINb + Lane)
#colnames(full_model)
# list_model(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 2: 
#variable_name
m <- 2
full_model <- model.matrix(~Line*Diet + RFI + Conc + RINa + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
pm1 <- proc.time()
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1
proc.time() -pm1

# Model 3: 
#variable_name
m <- 3
full_model <- model.matrix(~Diet + Line*RFI + Conc + RINa + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 4: 
#variable_name
m <- 4
full_model <- model.matrix(~Line + Diet*RFI + Conc + RINa + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)


assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 5: 
#variable_name
m <- 5
full_model <- model.matrix(~Line+Diet+RFI + Conc + RINa + RINb)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 6: 
#variable_name
m <- 6
full_model <- model.matrix(~Line+Diet+RFI + Conc + RINa + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 7: 
#variable_name
m <- 7
full_model <- model.matrix(~Line+Diet+RFI + Conc  + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)


assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 8: 
#variable_name
m <- 8
full_model <- model.matrix(~Line+Diet+RFI + RINa+ RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 9: 
#variable_name
m <- 9
full_model <- model.matrix(~Line+Diet+RFI + Conc  + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 10: 
#variable_name
m <- 10
full_model <- model.matrix(~Line+Diet+RFI +  RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 11: 
#variable_name
m <- 11
full_model <- model.matrix(~Line+Diet+RFI + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 12: 
#variable_name
m <- 12
full_model <- model.matrix(~Line+Diet+RFI + RINa+  RINb)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 13: 
#variable_name
m <- 13
full_model <- model.matrix(~Line + RFI + RINa + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 14: 
#variable_name
m <- 14
full_model <- model.matrix(~Line + RFI + RINa + RINb)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 15: 
#variable_name
m <- 15
full_model <- model.matrix(~Line + RFI + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 16: 
#variable_name
m <- 16
full_model <- model.matrix(~Line + RFI)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 17: 
#variable_name
m <- 17
full_model <- model.matrix(~Line)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

file_result <- paste(model_dir,"/Model",model_th, "_result.RData", sep ="")
load(file_result)
# result
# str(result)
# hist(result$Q.values[[3]])
# plot(result$Q.values[[3]], result$P.values[[3]])
# 
# which(result$Q.values[[3]] <.15)
# 
# counts[13985, Line ==1]
# counts[13985, Line ==2]
# 
# counts[916, Line ==1]
# counts[916, Line ==2]
# 
# 
# counts[5621, Line ==1]
# counts[5621, Line ==2]
# 
# counts[5823, Line ==1]
# counts[5823, Line ==2]
# 
# 
# counts[13411, Line ==1]
# counts[13411, Line ==2]

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))


#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 18: 
#variable_name
m <- 18
full_model <- model.matrix(~Line*Diet + RFI + Conc + RINa+RINb+Lane+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)


assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 19: 
#variable_name
m <- 19
full_model <- model.matrix(~Line*Diet + RFI + Conc + RINa+RINb+Lane+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 20: 
#variable_name
m <- 20
full_model <- model.matrix(~Line+ Diet + RFI + Conc + RINa+RINb+Lane+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


#counts <- counts[1:50,]
# Model 21: 
#variable_name
m <- 21
full_model <- model.matrix(~Line+ Diet + RFI + Conc + RINa+RINb+Lane+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 22: 
#variable_name
m <- 22
full_model <- model.matrix(~Line*Diet*RFI)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 23: 
#variable_name
m <- 23
full_model <- model.matrix(~Line*Diet + RFI)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1




# Model 24: 
#variable_name
m <- 24
full_model <- model.matrix(~Line*Diet)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 25: 
#variable_name
m <- 25
full_model <- model.matrix(~Line+Diet + RFI + dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 26: 
#variable_name
m <- 26
full_model <- model.matrix(~Line+Diet + RFI+ dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


#counts <- counts[1:50,]
# Model 27: 
#variable_name
m <- 27
full_model <- model.matrix(~Line+ RFI+ dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 28: 
#variable_name
m <- 28
full_model <- model.matrix(~Line+ RFI+ dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 29: 
#variable_name
m <- 29
full_model <- model.matrix(~Diet + Line+ RFI+ RINb + RINa + dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 30: 
#variable_name
m <- 30
full_model <- model.matrix(~Diet + Line+  dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 31: 
#variable_name
m <- 31
full_model <- model.matrix(~Diet + Line+  RFI)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 32: 
#variable_name
m <- 32
full_model <- model.matrix(~Diet + Line)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 33: 
#variable_name
m <- 33
full_model <- model.matrix(~Line + RINb + RINa)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

proc.time() -pm0

# Model 34: 
#variable_name
m <- 34
full_model <- model.matrix(~Line*Diet + RFI +  RINa + RINb + Lane)
#colnames(full_model)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
pm1 <- proc.time()
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1
proc.time() -pm1

# Model 35: 
#variable_name
m <- 35
full_model <- model.matrix(~Line*Diet + RFI + RINa+RINb+Lane+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 36: 
#variable_name
m <- 36
full_model <- model.matrix(~ Line + Diet + RFI + RINa + RINb + Lane + dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1




# Model 37: 
#variable_name
m <- 37
full_model <- model.matrix(~ Line + RFI + RINa + RINb + Lane + dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 38: 
#variable_name
m <- 38
full_model <- model.matrix(~ Line + RINa + RINb + Lane + dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 39: 
#variable_name
m <- 39
full_model <- model.matrix(~ Line + RINa + RINb + dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 40: 
#variable_name
m <- 40
full_model <- model.matrix(~Line+ Diet + RFI + Conc + RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 41: 
#variable_name
m <- 41
full_model <- model.matrix(~Line+ Diet + RFI +  RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 42: 
#variable_name
m <- 42
full_model <- model.matrix(~Line+ Diet +  Conc + RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 43: 
#variable_name
m <- 43
full_model <- model.matrix(~Line+ Diet +  Conc + RINa+RINb+Lane+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 44: 
#variable_name
m <- 44
full_model <- model.matrix(~Line+ Diet + RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 45: 
#variable_name
m <- 45
full_model <- model.matrix(~Line+ Diet + RINa+RINb)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 46: 
#variable_name
m <- 46
full_model <- model.matrix(~Line + RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 47: 
#variable_name
m <- 47
full_model <- model.matrix(~Diet +Line +RFI+ RINa+RINb+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 48: 
#variable_name
m <- 48
full_model <- model.matrix(~Diet + Line  + RINa+RINb+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


# Model 49: 
#variable_name
m <- 49
full_model <- model.matrix(~Diet + Line  +dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 50: 
#variable_name
m <- 50
full_model <- model.matrix(~Line  +dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



## model 51:
m <- 51
full_model <- model.matrix(~Line + dateGD)
model_th <- m
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))

#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1




## model 52:
m <- 52
full_model <- model.matrix(~Line*Diet + RFI + RINa + RINb)
model_th <- m
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))

#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1

# Model 53: 
#variable_name
m <- 53
full_model <- model.matrix(~Line+ Conc + RINa+RINb+Lane+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 54: 
#variable_name
m <- 54
full_model <- model.matrix(~Line+ Conc + RINa+RINb+dateGD)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1




# Model 55: 
#variable_name
m <- 55
full_model <- model.matrix(~Diet + Line  + RINb+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1



# Model 56: 
#variable_name
m <- 56
full_model <- model.matrix(~Diet + Line  + RINa+dateRNA)
model_th <- m 
list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(datdir, "/Reanalysis Data/result/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
file_fit <- paste(model_dir,"/Model",model_th, "_fit.RData", sep ="")
load(file_fit)

assign(paste("mean", model_th, sep = "_" ),mean(fit$phi.hat.dev))
get(paste("mean", model_th, sep = "_" ))

assign(paste("median", model_th, sep = "_" ),median(fit$phi.hat.dev))
get(paste("median", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),mean(AIC.QL(counts, fit)))
get(paste("AICQL", model_th, sep = "_" ))
#pm1 <- proc.time(); fit_model(full_model, model_th) ;proc.time() -pm1


aicql <- NULL
for (i in 1:56){
  aicql[i] <- get(paste("AICQL", i, sep = "_" ))
}
dput(aicql)
# 
# aicql <-
#   c(254.09743480009, 247.186014414358, 247.076309937566, 246.890982693105, 
#     243.086727310635, 243.520456591509, 243.589321393447, 242.601525323455, 
#     243.589321393447, 241.781142640803, 240.951447828559, 240.924294553224, 
#     240.189614988258, 238.558290989328, 238.635833463567, 237.078056292756, 
#     235.742446580315, 255.639139521706, 251.299567866567, 248.738734362595, 
#     252.953575771078, 248.182652274894, 241.547730330366, 239.846098780661, 
#     244.907680600988, 240.740984907204, 238.535003002783, 242.808201793557, 
#     248.822660679789, 242.706250102576, 239.358860320335, 237.671996687007, 
#     236.707563474895, 236.707563474895, 248.577126366692, 246.125139178975, 
#     243.841882811463, 241.597066985525, 240.059069716067, 251.14902530212, 
#     248.82053699604, 248.973124784279, 250.920989489666, 246.796864880925, 
#     238.905904442673, 244.282239603957, 244.367236657274, 242.353281253212, 
#     238.793592170467, 236.940544460193, 240.857807291276, 240.857807291276, 
#     248.311643543953, 246.391524017405)
model_ind <- c(19, 35:39, 33, 12, 14, 31, 32, 16, 17, 47:50, 18, 21,43,53,54 , 
               46, 25, 30, 28, 51)
aicout <- aicql[model_ind]
names(aicout) <- as.character(model_ind)
aicout[order(aicout)]

meanql <- medianql <- NULL
for (i in 1:51){
  meanql[i] <- get(paste("mean", i, sep = "_" ))
  medianql[i] <- get(paste("median", i, sep = "_" ))
}


model_ind <- c(19, 35:39, 33, 12, 14, 31, 32, 16, 17, 47:50, 18, 21, 40:42, 
               46, 25, 30, 28, 51)


meanout <- meanql[model_ind]
names(meanout) <- as.character(model_ind)
meanout


medianout <- medianql[model_ind]
names(medianout) <- as.character(model_ind)
which.min(medianql)
which.min(meanql)
which.min(aicql)
aicout
dput(aicout)
dput(aicql)
dput(meanql)
# 
# meanql <- 
#   c(0.334721646024284, 0.350116481528789, 0.347672260342115, 0.350809204222422, 
#     0.355401826647718, 0.384166632385809, 0.382662066233458, 0.358729210186931, 
#     0.382662066233458, 0.388845742054197, 0.419958265040179, 0.359476955865531, 
#     0.362729682818776, 0.362981252909053, 0.421142652777239, 0.423970834155595, 
#     0.443845315651025, 0.312302309139007, 0.344739202718954, 0.34978215262864, 
#     0.312548301464419, 0.410276105839358, 0.417996679663526, 0.432821665183204, 
#     0.348773100738718, 0.379381733755932, 0.382940256826451, 0.361051393299416, 
#     0.325998719896125, 0.352603271415917, 0.422296883915852, 0.435883428360209, 
#     0.364865657175356, 0.364865657175356, 0.350388440600861, 0.355869103323954, 
#     0.352657798747701, 0.360096851210259, 0.361209762853011, 0.310444734472288, 
#     0.326065166862968, 0.308868200449161, 0.309287721531872, 0.32208267726487, 
#     0.362043890132158, 0.327403628758292, 0.358413366052591, 0.362068587453515, 
#     0.386244102421621, 0.400043240543351, 0.369759896598516)

dput(medianql)


# medianql <- 
#   c(0.184732474976862, 0.19738492710699, 0.197217148795153, 0.198866129285064, 
#     0.20352322145403, 0.227717320004234, 0.231638327434404, 0.204433092619715, 
#     0.231638327434404, 0.235595167743339, 0.266071892814122, 0.205941615216103, 
#     0.209406611452538, 0.210591250459231, 0.26887664889006, 0.272101072081571, 
#     0.284322986522261, 0.157224262699112, 0.183740908811615, 0.186431038858083, 
#     0.16035592566512, 0.261358129194564, 0.269143360641632, 0.277109676245007, 
#     0.201919673396674, 0.223426043598974, 0.230715334362308, 0.216171008005133, 
#     0.172149155011531, 0.208359188372984, 0.269977011794668, 0.278506558337507, 
#     0.209385547473775, 0.209385547473775, 0.189023737463352, 0.193245488150766, 
#     0.195368976261989, 0.199933911558764, 0.203375854277504, 0.163999927266782, 
#     0.172141869426716, 0.165577396251722, 0.161832444573113, 0.17410724898945, 
#     0.206328487907058, 0.180635965453675, 0.197986710658171, 0.199605804926009, 
#     0.22993927554767, 0.243739713432757, 0.223225395554762)

which.min(medianout)
which.min(meanout)
which.min(aicout)
aicout