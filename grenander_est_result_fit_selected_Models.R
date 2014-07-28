\documentclass{article}
% \usepackage[sc]{mathpazo}
% \usepackage[T1]{fontenc}
\usepackage{geometry}
\usepackage{amsmath, amssymb, mathtools }
\usepackage{enumerate}
\usepackage{array}
\usepackage{fancyhdr}
\usepackage{verbatim}
\usepackage{color}
\usepackage{pstricks}
% \usepackage{longtable}
% \usepackage{fancyvrb}
% \usepackage{fancyhdr}
\usepackage{eqnarray}
%\pagestyle{fancy}
\usepackage{psfrag}
\usepackage{epsfig,epsf}
\usepackage{pstricks}
\geometry{verbose,tmargin=2.5cm,bmargin=2.5cm,lmargin=2.5cm,rmargin=2.5cm}
\setcounter{secnumdepth}{2}
\setcounter{tocdepth}{2}
\usepackage{url}
\usepackage[unicode=true,pdfusetitle,
            bookmarks=true,bookmarksnumbered=true,bookmarksopen=true,bookmarksopenlevel=2,
            breaklinks=false,pdfborder={0 0 1},backref=false,colorlinks=false]
{hyperref}
\hypersetup{
  pdfstartview={XYZ null null 1}}
\usepackage{breakurl}
\begin{document}



<<setup, include=FALSE>>=
  #render_listings()
  #pdf.options(useDingbats = TRUE)
  opts_chunk$set(fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize',message=FALSE,error=FALSE,warning=FALSE)
opts_chunk$set(fig.width=5, fig.height=5, out.width='.6\\linewidth', fig.align='center')

@
  
  <<echo=FALSE,results='hide'>>=
  # some setup
  options(width=60) # make the printing fit on the page
set.seed(1121) # make the results repeatable
require(xtable)
require(PASWR)
@
  
  \title{Analyze RNASeq Data from G8P2 RFI Lines Using QuasiSeq Package (paired end read)}

\author{Yet Nguyen}

\maketitle
<<echo=FALSE>>=
  options(width=90)
@
  
<<>>=
library(reshape)
library(plyr)
library(fdrtool)
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
# extract cbc data
cbc.data <- read.csv("cbcdata.csv")
which(nname%in%cbc.data$Idpig) # missing cbc for animal 5
nname[5]
idpig_cbc <- laply(1:(length(nname)-1), function(i) which(cbc.data$Idpig==nname[-5][i]))
cbc_cov <- cbc.data[idpig_cbc,]
name <- paste("X",meta.data$Sample.Name[-5], sep = "") # meta.data2$Sample.Name
del_row <- which(rownames(dat) %in%c("ENSSSCG00000007978", "ENSSSCG00000014725"))
dat2 <- dat[-del_row, name]
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
counts <- as.matrix(dat2[rowSums(dat2>0)>3&
                           rowMeans(dat2)>8,])

dim(counts)
str(counts)
pval.hist.grenander <- function(p.value){
  grenander.out <- grenander(ecdf(p.value))
  p.brks <- c(0, grenander.out$x.knots)
  b.edf <- c(0, grenander.out$F.knots)
  p.diffs <- diff(p.brks)
  h.cdf <- approx(p.brks, b.edf, xout = p.value)$y  # get the histogram-based CDF estimates from each p-value
  p.hist <- exp(log(diff(b.edf))-log(diff(p.brks))) # get the hight for each histogram bar
  pi0.hat <- min(p.hist)                            # get the pi0 estimate from histogram bar
  h.ebp <- approx(p.brks, pi0.hat/c(p.hist, p.hist[length(p.hist)]), xout = p.value)$y # get the conservative EBP interpolation 
  h.fdr <- exp(log(pi0.hat) + log(p.value) - log(h.cdf))                                     # Get the histogram based FDR estimate
  h.ebp[p.value==0] <- 0
  h.fdr[p.value==0] <- 0
  return(list( p.value = p.value,          # input p-value,
               h.cdf = h.cdf,              # the histogram Grenander based cdf estimate
               h.fdr = h.fdr,              # the histogram Grenander based FDR
               h.ebp = h.ebp,              # the histogram Grenander based EBP
               p.brks = p.brks,            # the p-value break-points of the histogram
               p.hist = p.hist,            # the heights of each histogram bar
               edf.brks = b.edf,           # the breaks points in the EDF of the histogram estimator
               pi0.hat = pi0.hat))         # the histogram Grenander based estimate of the proportion of tests with a true null hypothesis
}

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
  qvalue_dat <- result$Q.values[[3]][,colnames(result$Q.values[[3]])]
  qvalue_20 <- apply(qvalue_dat<=.10, 2, sum)
return( data.frame(pvalue05 = order(pvalue_05),
                   
                   ad = order(ad),
                   cvm = order(cvm),
                   ks = order(ks)))
}

panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}

panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

logmean_gene1 <- log2(apply(counts[,Line ==1] +1, 1, mean))
logmean_gene2 <- log2(apply(counts[,Line ==2] +1, 1, mean))
log2fc <- abs(logmean_gene1-logmean_gene2)
logfc1_gene <- which(log2fc >=1)
gene_list <- rownames(counts)
# function result_sum ####
result_sum <- function(result, threshold){
  pvalue_line <- result$P.values[[3]][,"Line"]
  de_gene <- which(result$Q.values[[3]][,"Line"] <= threshold)
  de_log2fc1_gene <- intersect(logfc1_gene, de_gene)
  id_de_gene <- gene_list[de_gene]
  id_de_log2fc1_gene <- gene_list[de_log2fc1_gene]
  out <- list(de_gene = de_gene, 
              de_log2fc1_gene = de_log2fc1_gene, 
              id_de_gene = id_de_gene, 
              id_de_log2fc1_gene=id_de_log2fc1_gene,
              n_de_gene = length(de_gene), 
              n_de_log2fc1_gene = length(de_log2fc1_gene),
              pvalue_line = pvalue_line)
  return(out)
}

# model 1####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model1.Line.Diet.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model1_result.RData")

sel_criteria(result)

# model 2####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model2.Line.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model2_result.RData")
sel_criteria(result)

# model 3####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model3.Line.RFI.RINb.RINa.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model3_result.RData")

sel_criteria(result)
# model 4####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model4.Line.RFI.RINb.RINa.Lane.llymp.lneut.lmono.leosi.dateRNA/Model4_result.RData")

sel_criteria(result)

# model 5####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model5.Line.RFI.RINb.RINa.Lane.llymp.lneut.leosi.dateRNA/Model5_result.RData")

sel_criteria(result)
# model 6####

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model6.Line.RFI.RINb.RINa.Lane.llymp.lneut.dateRNA/Model6_result.RData")
sel_criteria(result)

pvalue_line6 <- (result$P.values[[3]][,"Line"])
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")


#plot(pvalueline6, qvalueline6)
# model 7####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model7.Line.RINb.RINa.Lane.llymp.lneut.dateRNA/Model7_result.RData")
sel_criteria(result)

pvalue_line7 <-(result$P.values[[3]][,"Line"])
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out


# model 8####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model8.Line.RINb.RINa.llymp.lneut.dateRNA/Model8_result.RData")
sel_criteria(result)
pvalue_line8 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out


pairs(log(cbind(pvalue_line6, pvalue_line7, pvalue_line8)),
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)

# model 9####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model9.Line.Diet.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateGD/Model9_result.RData")

sel_criteria(result)

# model 10####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model10.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateGD/Model10_result.RData")
sel_criteria(result)

# model 11####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model11.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.dateGD/Model11_result.RData")

sel_criteria(result)
# model 12####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model12.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.dateGD/Model12_result.RData")

sel_criteria(result)

pvalue_line12 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out
# model 13####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model13.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.dateGD/Model13_result.RData")

sel_criteria(result)
pvalue_line13 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out

# model 14####

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model14.Line.RINb.RINa.Conc.Lane.llymp.lneut.dateGD/Model14_result.RData")
sel_criteria(result)
pvalue_line14 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out

# model 15####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model15.Line.RINb.RINa.Lane.llymp.lneut.dateGD/Model15_result.RData")
sel_criteria(result)
pvalue_line15 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out


# model 16####
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model16.Line.RINb.RINa.llymp.lneut.dateGD/Model16_result.RData")
sel_criteria(result)
pvalue_line16 <- result$P.values[[3]][,"Line"]
## Total number of Line  DE Genes when FDR is controlled at 0.05, 0.10, 0.15, .20

degene <- c(result_sum(result, 0.05)$n_de_gene,
            result_sum(result, 0.10)$n_de_gene,
            result_sum(result, 0.15)$n_de_gene,
            result_sum(result, 0.2)$n_de_gene)



## Total number of Line DE Genes with log(FC) >=1 when FDR is controlled at 0.05, 0.10, 0.15, .20

lf1 <- c(result_sum(result, 0.05)$n_de_log2fc1_gene,
         result_sum(result, 0.10)$n_de_log2fc1_gene,
         result_sum(result, 0.15)$n_de_log2fc1_gene,
         result_sum(result, 0.2)$n_de_log2fc1_gene)


# Summary table
out <- data.frame(FDR =c(0.05, 0.10, 0.15, .20), degene = degene, lf1 = lf1)
colnames(out) <- c("FDR", "DEGs", "$\\mbox{log(FC)}\\geq 1$")
out


pairs(log(cbind(pvalue_line13, pvalue_line14, pvalue_line15, pvalue_line16)),
      lower.panel=panel.smooth, upper.panel=panel.cor,diag.panel=panel.hist)
