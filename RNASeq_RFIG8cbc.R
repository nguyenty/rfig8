pm0 <- proc.time()
# cleaning data
library(edgeR)
library(fields)
library(plyr)
library(reshape)
dir.source <- "U:/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
#dir.source <- "/home/ntyet/stevescode/QuasiSeq_1.0-2/QuasiSeq/R/"
source("QL.fit.R")
source("NBDev.R")
source("PoisDev.R")
source("QL.results.R")


# datdir <- "/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data"
# datdir <- "/home/ntyet/cyfiles/R/RA/Data"
datdir <-  "U:/R/RA/Data"
ndata <- read.csv("g8p2_g9p2-rfiadj-FINAL_Jan_2014_rfiadjusted.csv")


# cbc <- read.csv('cbcdata.csv',
#                 header =T)
dat <- read.table("RFI_uniq_comb_count_corrected.txt")
# #dat <- read.table("RFI_uniq_comb_count_2.txt")
# dim(dat)
# which(rownames(dat) %in% "ENSSSCG00000007978")
# which(rownames(dat) %in% "ENSSSCG00000014725")
# dat[which(rownames(dat) %in% "ENSSSCG00000007978"), ]
# dat[which(rownames(dat) %in% "ENSSSCG00000014725"),]

meta.data <- read.csv("RIN values RNAseq Obj3a Martine.csv")
nname <- paste("20800",sprintf("%04.0f",meta.data$Sample.Name), sep = "")
meta.data2 <- cbind(ID=as.numeric(nname),meta.data)

nndata <- ndata[ndata[,"idpig"]%in% nname,c("idpig","rfi.ADJUSTED") ]
for (i in 1:length(nname)){
  for(j in 1:length(nname)){
    if (meta.data2$ID[i]==nndata$idpig[j])
      meta.data2$RFI.value[i] <- nndata$rfi.ADJUSTED[j]
  }
}
# meta.data2
# write.table(meta.data2, "Meta.data2.txt")
# write.csv(meta.data2, paste(datdir,"/Meta.data2.csv",sep =""))
# str(meta.data2)
# str(meta.data)
# meta.data2[,c("ID","Sample.Name","RFI.value")]
# meta.data[,c("Sample.Name","RFI.value")]

# extract cbc data
cbc.data <- read.csv("cbcdata.csv")

which(nname%in%cbc.data$Idpig) # missing cbc for animal 5
nname[5]
idpig_cbc <- laply(1:(length(nname)-1), function(i) which(cbc.data$Idpig==nname[-5][i]))
idpig_cbc 
cbc_cov <- cbc.data[idpig_cbc,]
dim(cbc_cov)
name <- paste("X",meta.data$Sample.Name[-5], sep = "") # meta.data2$Sample.Name
del_row <- which(rownames(dat) %in%c("ENSSSCG00000007978", "ENSSSCG00000014725"))
dat2 <- dat[-del_row, name]
dim(dat2)

Lane <- as.factor(meta.data2[-5,]$Lane.RNAseq)
Diet <- as.factor(meta.data2[-5,]$diet)
Line <- as.factor(meta.data2[-5,]$line)
RFI <- meta.data2[-5,]$RFI.value
RINb <- meta.data2[-5,]$RIN.before.GD
RINa  <- meta.data2[-5,]$RIN.after.GD
Conc <- meta.data2[-5,]$Conc.after.GD.ng.ul.
dateGD <- meta.data2[-5,]$date.GD
dateRNA <- meta.data2[-5,]$date.RNA.extraction
colnames(cbc_cov)

llymp <- log(cbc_cov[,"Lymphocyte"])
lneut <- log(cbc_cov[,"Neutrophil"])
lmono <- log(cbc_cov[,"Monocyte"])
leosi <- log(cbc_cov[,"Eosinophil"])
lbaso <- log(cbc_cov[,"Basophil"])
variable_name <- c("Lane", "Diet", "Line", "RFI", 
                   "RINb", "RINa", 
                   "Conc", "dateGD", "dateRNA")


# 
# colnames(cbc)
# cbc$Idpig
# cbccov <- cbc[cbc$Idpig ]
# nname %in% cbc$Idpig
# nname[5]
counts <- as.matrix(dat2[rowSums(dat2>0)>3&
                           rowMeans(dat2)>8,])
log.offset <- log(apply(counts, 2, quantile, .75))
#dim(counts)

#colnames(counts)


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

library(fdrtool)

g_cdf <- function(z){
  e <- ecdf(z)
  g <- grenander(e)
  g
}
sel_criteria <- function(result){
  dat <- result$P.values[[3]][,colnames(result$P.values[[3]])]
  # Crames Von Miser statistics
  cvm <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2 *
                                        diff(c(0,g_cdf(z)$x.knots))))
  # Kolmogorow Smirnov statistics 
  ks <- apply(dat, 2, function(z)max((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2))
  
  # Anderson-Darling statistics
  ad <- apply(dat, 2, function(z)sum((g_cdf(z)$F.knots - g_cdf(z)$x.knots)^2/
                                       g_cdf(z)$x.knots*(1-g_cdf(z)$x.knots)*
                                       diff(c(0,g_cdf(z)$x.knots))))
  # Proportion of pvalue less than 0.05
  pvalue_05 <- apply(dat<=0.05, 2, sum)
  
  out <- data.frame(pvalue05 = order(pvalue_05),
                    ad = order(ad),
                    cvm = order(cvm),
                    ks = order(ks))
  
  return(out)
}


## Function do the model list and test.mat 

list_model <- function(full_model){
  #colnames(full_model)
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

# colnames(design.list[[13]])
# test.mat
  
## Function do all the things with input Full model

fit_model <- function(full_model, model_th){
  list_out <- list_model(full_model)
  design.list <- list_out$design.list
#   design.list <- vector("list", 2)
#   design.list[[1]] <- full_model
#   
#   design.list[[2]] <- full_model[,-13]
#   test.mat <- matrix(c(1,2), ncol = 2, nrow = 1, byrow = F)
#rownames(test.mat) <- "LineDiet"
   test.mat <- list_out$test.mat
  fit <- QL.fit(counts, design.list, test.mat, 
                log.offset = log.offset, print.progress=FALSE,
                Model = "NegBin")
  result<- QL.results(fit, Plot = FALSE)
  k <- nrow(test.mat)
  name_model <- NULL
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(datdir, "/Reanalysis Data/resultcbc/Model",model_th,name_model, sep ="")
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
  return(sel_criteria(result))
}



# Model 0: ####
m <- 0
model_th <- m
full_model <- model.matrix(~Line*Diet + RFI + RINb + RINa + Conc + Lane + dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 1: ####
m <- 1
model_th <- m
full_model <- model.matrix(~Line+Diet + RFI + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono + leosi + lbaso + dateRNA)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 2: ####
m <- 2
model_th <- m
full_model <- model.matrix(~Line+ RFI + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono + leosi + 
                             lbaso + dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;

# Model 3: ####
m <- 3
model_th <- m
full_model <- model.matrix(~Line+ RFI + RINb + RINa + 
                             Lane + llymp + 
                             lneut + lmono + leosi + 
                             lbaso + dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 4: ####
m <- 4
model_th <- m
full_model <- model.matrix(~Line+ RFI + RINb + RINa + 
                             Lane + llymp + 
                             lneut + lmono + leosi + 
                             dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 5: ####
m <- 5
model_th <- m
full_model <- model.matrix(~Line+ RFI + RINb + RINa + 
                             Lane + llymp + 
                             lneut +  leosi + 
                             dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1



# Model 6: ####
m <- 6
model_th <- m
full_model <- model.matrix(~Line+ RFI + RINb + RINa + 
                             Lane +llymp + 
                             lneut +  
                             dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 7: ####
m <- 7
model_th <- m
full_model <- model.matrix(~Line+  RINb + RINa + 
                             Lane +llymp + 
                             lneut +  
                             dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 8: ####
m <- 8
model_th <- m
full_model <- model.matrix(~Line+  RINb + RINa + 
                             llymp + 
                             lneut +  
                             dateRNA)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
list_model(full_model)$test.mat
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))

proc.time() -pm1
#; fit_model(full_model, model_th) ;

# Model 9: ####
m <- 9
model_th <- m
full_model <- model.matrix(~Line+Diet + RFI + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono + leosi + lbaso + dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;

# Model 10: ####
m <- 10
model_th <- m
full_model <- model.matrix(~Line+Diet + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono + leosi + lbaso + dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 11: ####
m <- 11
model_th <- m
full_model <- model.matrix(~Line+Diet + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono + leosi +  dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;

# Model 12: ####
m <- 12
model_th <- m
full_model <- model.matrix(~Line+Diet + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut + lmono +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;

# Model 13: ####
m <- 13
model_th <- m
full_model <- model.matrix(~Line+Diet + RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut  +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 14: ####
m <- 14
model_th <- m
full_model <- model.matrix(~Line+ RINb + RINa + Conc + 
                             Lane + llymp + 
                             lneut  +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 15: ####
m <- 15
model_th <- m
full_model <- model.matrix(~Line+ RINb + RINa + 
                             Lane + llymp + 
                             lneut  +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;




# Model 15: ####
m <- 15
model_th <- m
full_model <- model.matrix(~Line+ RINb + RINa + 
                             Lane + llymp + 
                             lneut  +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;


# Model 16: ####
m <- 16
model_th <- m
full_model <- model.matrix(~Line+ RINb + RINa + 
                             llymp + 
                             lneut  +   dateGD)
#dim(full_model)
pm1 <- proc.time()
out_model <- fit_model(full_model, model_th)
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
get(paste("ms_criteria", model_th, sep = "_" ))
list_model(full_model)$test.mat
proc.time() -pm1
#; fit_model(full_model, model_th) ;

