require(Matrix)
#library(QuasiSeq)
library(edgeR)
require(reshape)
require(plyr)
library(fields)
library(reshape)
dir.source <- "U:/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
#dir.source <- "/home/ntyet/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
source(paste(dir.source, "QL.fit.R",sep=""))
source(paste(dir.source, "NBDev.R",sep =""))
source(paste(dir.source, "PoisDev.R",sep =""))
source(paste(dir.source, "QL.results.R",sep =""))

resultdir <- "U:/R/RA/Data/Reanalysis Data/resultg8"
scount <- read.table("single end uniquely mapped reads count table for Yet.txt", 
                     header = T)
cbc <- read.table('CBC data for pigs with RNA-seq data avaible.txt',
                  header =T)

metadata <- read.table("meta_data_RNA-seq_G9P2.txt", 
                       header = T)

rfiadj <- read.csv("g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv", 
                   header = T)

##### cleaning data####



cbc <- cbc[order(cbc$ear), ]
metadata <- metadata[order(metadata$idpig), ]
rfiadj <- rfiadj[order(rfiadj$idpig),]

fullidpig <- as.numeric(paste("20900", metadata$idpig, sep = ""))

covset <- cbind(metadata[, -4], rfiadj[rfiadj$idpig %in% fullidpig, c("rfi.ADJUSTED")], 
                cbc[, c("iddam", "idsire", "Neutrophil", "Lymphocyte", "Monocyte",
                        "Eosinophil", "Basophil" )])
colnames(covset) <- c("idpig", "Line", "Diet",  "Block", "Blockorder", "Concb", 
                      "RINb", "Conca", "RINa", "RFI",
                      "iddam", "idsire", "neut",   
                      "lymp", "mono","eosi", "baso")

#####set of covariates considered ####
covset <- cbind(covset, lneut = log(covset$neut), llymp = log(covset$lymp), 
                lmono = log(covset$mono), leosi = log(covset$eosi), 
                lbaso = log(covset$baso)) 
covset$Line <- as.factor(covset$Line)
covset$Diet <- as.factor(covset$Diet)
covset$Block <- as.factor(covset$Block)
covset$Blockorder <- as.factor(covset$Blockorder)
covset$iddam <- as.factor(covset$iddam)
covset$idsire <- as.factor(covset$idsire)
levels(covset$idsire) <- 1:11
levels(covset$iddam) <- 1:20

covset[, c("iddam", "idsire")]
#detach(covset)
attach(covset)

counts <- as.matrix(scount[rowSums(scount[,-1]>0)>3&
                             rowMeans(scount[,-1])>8 & 
                             rowSums(scount[,-1][,Line ==1] > 0) >0 &
                             rowSums(scount[,-1][, Line ==2] >0) >0 ,-1])
dim(counts)
dim(scount)
dim(counts)
log.offset <- log(apply(counts, 2, quantile, .75))
###List of models function ####
### Case 1: no cbc data ####
# dim(covset)
# colnames(covset)
# full_model <- model.matrix(~Line*Diet*RFI + Block+ Blockorder)
# colnames(full_model)
# rankMatrix(full_model)
# list_model(full_model)
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



###############

list_model <- function(full_model){
  n <- dim(full_model)[2]
  variable_name <- colnames(full_model)[-1]
  variable_name <- gsub(":", "", variable_name)
  for (i in 2:8){variable_name <- gsub(i, "", variable_name)}
  test.mat <- NULL
  design.list <- vector("list", n)
  design.list[[1]] <- full_model
  for (i in 2:n) {
    design.list[[i]] <- as.matrix(full_model[,-i])
    test.mat <- rbind(test.mat, c(1,i))
  }
  
  row.names(test.mat) <-  variable_name
  
  if (any(variable_name == "Block") & any(variable_name == "Blockorder")){
    ind_block <- which(variable_name == "Block")[1]
    ind_blockorder <- which(variable_name == "Blockorder")[1] 
    nlist <- n - sum(variable_name == "Block") -
      sum(variable_name == "Blockorder") + 2
    design.list <- vector("list", nlist)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 1:(ind_block-1)){ # i <- ind_block - 1
      design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
      row.names(test.mat)[i] <- variable_name[i]
    }
    
    design.list[[ind_block + 1]] <- full_model[,-((ind_block+1):(ind_block+3))]
    # row.names(test.mat)
    #colnames(design.list[[ind_block+1]])
    design.list[[ind_block +2]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
    row.names(test.mat)[ind_block + 1] <- "Blockorder" 
    # colnames(design.list[[ind_block+2]])
    # colnames(design.list[[1]])
    if(ind_blockorder + 7 < n){ # i <- 15 # colnames(design.list[[i]])
      for(i in ((ind_block+3):(n-8))){
        design.list[[i]] <- full_model[,-(i + 8)] # colnames (full_model)
        row.names(test.mat)[i-1] <- variable_name[i+7]    
      }
    }
  }
  
  if (any(variable_name == "Block") & all(variable_name != "Blockorder")){
    ind_block <- which(variable_name == "Block")[1]
    nlist <- n - sum(variable_name == "Block") + 1
    design.list <- vector("list", nlist)
    design.list[[1]] <- full_model
    test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
    row.names(test.mat) <- variable_name[1:nrow(test.mat)]
    
    for(i in 1:(ind_block-1)){ # i <- ind_block - 1
      design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
      row.names(test.mat)[i] <- variable_name[i]
    }
    
    design.list[[ind_block + 1]] <- full_model[,-((ind_block+1):(ind_block+3))]
    # row.names(test.mat)[ind_block] <- "Block"
    # row.names(test.mat)
    #colnames(design.list[[ind_block+1]])
#     design.list[[ind_block +2]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
#     row.names(test.mat)[ind_block + 1] <- "Blockorder" 
    # colnames(design.list[[ind_block+2]])
    # colnames(design.list[[1]])
    if(ind_block + 3 < n){ # i <- 15 # colnames(design.list[[i]])
      for(i in ((ind_block+2):(n-2))){
        design.list[[i]] <- full_model[,-(i + 2)] # colnames (full_model)
        row.names(test.mat)[i-1] <- variable_name[i+1]    
      }
    }
  }
  
if (all(variable_name != "Block") & any(variable_name == "Blockorder")){
  ind_blockorder <- which(variable_name == "Blockorder")[1] 
  nlist <- n - sum(variable_name == "Blockorder") + 1
  design.list <- vector("list", nlist)
  design.list[[1]] <- full_model
  test.mat <- laply(1:(nlist-1), function(i) c(1,i+1))
  row.names(test.mat) <- variable_name[1:nrow(test.mat)]
  
  for(i in 1:(ind_blockorder-1)){ # i <- ind_blockorder - 1
    design.list[[i+1]] <- as.matrix(full_model[,-(i+1)])
    row.names(test.mat)[i] <- variable_name[i]
  }
  
  design.list[[ind_blockorder + 1]] <- full_model[,-((ind_blockorder+1):(ind_blockorder+7))]
  # row.names(test.mat)[[ind_blockorder]]
  # colnames(design.list[[ind_blockorder+1]])
  row.names(test.mat)[ind_blockorder] <- "Blockorder" 
  # colnames(design.list[[ind_block+2]])
  # colnames(design.list[[1]])
  if(ind_blockorder + 7 < n){ # i <- 15 # colnames(design.list[[i]])
    for(i in ((ind_blockorder+2):(n-6))){
      design.list[[i]] <- full_model[,-(i + 6)] # colnames (full_model)
      row.names(test.mat)[i-1] <- variable_name[i+5]    
    }
  }
}
  if (n ==2) design.list[[2]] <- rep(1, nrow(full_model))
  return(list(design.list = design.list, test.mat = test.mat))
}


## Function do all the things with input Full model

fit_model <- function(full_model, model_th){ # model_th <- 1
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, 
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  pvalue_05 <- apply(result$P.values[[3]][,rownames(test.mat)]<=0.05, 2, sum)
  k <- nrow(test.mat)
  name_model <- NULL 
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir, "/Model",model_th,name_model, sep ="")
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
  
  return(list(mean_model = mean(fit$phi.hat.dev), 
              median_model = median(fit$phi.hat.dev), 
              AIC_model = mean(AIC.QL(counts, fit)),
              pvalue_05 = pvalue_05))
}


# ### check correlation of cbc data####
# #pairs(cbind(lneut, llymp, lmono, leosi, lbaso))
# 
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
# {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y))
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste0(prefix, txt)
#   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor * r)
# }
# pairs(cbind(lneut, llymp, lmono, leosi, lbaso), lower.panel = panel.smooth, upper.panel = panel.cor)
# 


# Model 1
m <- 1
model_th <- m
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 2
m <- 2
model_th <- m
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 3
m <- 3
model_th <- m
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + lbaso + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 305
m <- 305
model_th <- m
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono +  
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 306
m <- 306
model_th <- m
full_model <- model.matrix(~Line + Diet + RFI + Concb + RINb +  RINa + 
                             lneut + llymp + lmono +  lbaso+
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 4
m <- 4
model_th <- m
full_model <- model.matrix(~Line + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + lbaso + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 

# Model 5
m <- 5
model_th <- m
full_model <- model.matrix(~Line + RFI + Concb + RINb + RINa + 
                             lneut + llymp + lmono + lbaso + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
out_model$pvalue_05

assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


# Model 6
m <- 6
model_th <- m

full_model <- model.matrix(~Line + Concb + RINb + RINa + 
                             lneut + llymp + lmono + lbaso + 
                             Block)


#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 7
m <- 7
model_th <- m
full_model <- model.matrix(~Line + Concb + RINa + 
                             lneut + llymp + lmono + lbaso + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



# Model 8
m <- 8
model_th <- m
full_model <- model.matrix(~Line + Concb + RINa + 
                             lneut + llymp + lmono + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 9
m <- 9
model_th <- m
full_model <- model.matrix(~Line + Concb + RINa + 
                             lneut + llymp + 
                             Block)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("pvalue05", model_th, sep = "_" ),out_model$pvalue_05)
get(paste("pvalue05", model_th, sep = "_" ))
assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 


#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 



## Model 100

model_th <- 100
full_model <- model.matrix(~Line*Diet)

design.list <- vector("list",4)
design.list[[1]] <- full_model
Cdiet <- matrix(c(1, 0, 0, 
              0, 1, 0, 
              0, 0, 1,
              0, -2, 0), nrow = 4, byrow = T)
design.list[[2]] <- full_model%*% Cdiet
Cline <- matrix(c(1, 0, 0, 
                           0, 1, 0, 
                           0, 0, 1,
                           0, 0, -2), nrow = 4, byrow = T)

design.list[[3]] <- full_model%*% Cline

design.list[[4]] <- model.matrix(~Line+Diet)

test.mat <- rbind(c(1,2), c(1,3), c(1,4))
rownames(test.mat) <- c("Diet", "Line", "LineDiet")




# Model 101
m <- 101
model_th <- m
full_model <- model.matrix(~Line*Diet*RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)
design.list <- vector("list",2)
design.list[[1]] <- full_model

design.list[[2]] <- model.matrix(~Line*Diet + Line*RFI + Diet*RFI + Concb + RINb + Conca + RINa + 
                                   lneut + llymp + lmono + leosi + lbaso + 
                                   Block + Blockorder)

test.mat <- rbind(c(1,2))
rownames(test.mat) <- c("LineDietRFI")




#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




# Model 102
m <- 102
model_th <- m
full_model <- model.matrix(~Line*Diet + RFI + Concb + RINb + Conca + RINa + 
                             lneut + llymp + lmono + leosi + lbaso + 
                             Block + Blockorder)
design.list <- vector("list",2)
design.list[[1]] <- full_model

design.list[[2]] <- model.matrix(~Line + Diet + RFI + Concb + RINb + Conca + RINa + 
                                   lneut + llymp + lmono + leosi + lbaso + 
                                   Block + Blockorder)

test.mat <- rbind(c(1,2))
rownames(test.mat) <- c("LineDiet")




#colnames(full_model)
#rankMatrix(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)

assign(paste("AICQL", model_th, sep = "_" ),out_model$AIC_model)
get(paste("AICQL", model_th, sep = "_" ))

assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))


assign(paste("mean", model_th, sep = "_" ),out_model$mean_model)
get(paste("mean", model_th, sep = "_" ))

proc.time() -pm1 




list_out <- list_model(full_model)
design.list <- list_out$design.list
test.mat <- list_out$test.mat
k <- nrow(test.mat)
name_model <- NULL
for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
model_dir <- paste(resultdir, "/Model",model_th,name_model, sep ="")
dir.create(model_dir, showWarnings = FALSE)
load(paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
load(paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
str(result)
which(result$Q.values[[3]][,1]<=.05)
sum(abs(log(apply(counts[which(result$Q.values[[3]][,1]<=.05), Line ==1]+1, 1, mean)/apply(counts[which(result$Q.values[[3]][,1]<=.05), Line ==2]+1, 1, mean) )>1))

for (i in 1:11) get(paste("AICQL", i, sep = "_" ))
which.min(c(AICQL_1, AICQL_2, AICQL_3, AICQL_4, AICQL_5, AICQL_6, AICQL_7, AICQL_8, AICQL_9, AICQL_10, AICQL_11))
model_aic <- c(AICQL_1, AICQL_2, AICQL_3, AICQL_4, AICQL_5, AICQL_6, AICQL_7, AICQL_8, AICQL_9, AICQL_10, AICQL_11)

names(model_aic) <- paste("M", 1:11, sep = "")

# dput(model_aic)
# structure(c(318.082732897024, 339.866177153431, 322.237901868796, 
#             320.961821003267, 318.705315044833, 316.637592138859, 314.503342582272, 
#             312.791184638659, 310.914392806018, 309.037439804778, 307.04581821839
# ), .Names = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", 
#               "M9", "M10", "M11"))

write.table(model_aic, file = "model_aic.txt")
save(model_aic, file = "model_aic.RData")

