library(fdrtool)

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

# model 1
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model1.Line.Diet.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model1_result.RData")

sel_criteria(result)

# model 2
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model2.Line.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model2_result.RData")
sel_criteria(result)

# model 3
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model3.Line.RFI.RINb.RINa.Lane.llymp.lneut.lmono.leosi.lbaso.dateRNA/Model3_result.RData")

sel_criteria(result)
# model 4
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model4.Line.RFI.RINb.RINa.Lane.llymp.lneut.lmono.leosi.dateRNA/Model4_result.RData")

sel_criteria(result)

# model 5
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model5.Line.RFI.RINb.RINa.Lane.llymp.lneut.leosi.dateRNA/Model5_result.RData")

sel_criteria(result)
# model 6

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model6.Line.RFI.RINb.RINa.Lane.llymp.lneut.dateRNA/Model6_result.RData")
sel_criteria(result)

# model 7
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model7.Line.RINb.RINa.Lane.llymp.lneut.dateRNA/Model7_result.RData")
sel_criteria(result)

# model 8
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model8.Line.RINb.RINa.llymp.lneut.dateRNA/Model8_result.RData")
sel_criteria(result)

# model 9
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model9.Line.Diet.RFI.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateGD/Model9_result.RData")

sel_criteria(result)

# model 10
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model10.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.lbaso.dateGD/Model10_result.RData")
sel_criteria(result)

# model 11
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model11.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.leosi.dateGD/Model11_result.RData")

sel_criteria(result)
# model 12
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model12.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.lmono.dateGD/Model12_result.RData")

sel_criteria(result)

# model 13
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model13.Line.Diet.RINb.RINa.Conc.Lane.llymp.lneut.dateGD/Model13_result.RData")

sel_criteria(result)
# model 14

load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model14.Line.RINb.RINa.Conc.Lane.llymp.lneut.dateGD/Model14_result.RData")
sel_criteria(result)

# model 15
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model15.Line.RINb.RINa.Lane.llymp.lneut.dateGD/Model15_result.RData")
sel_criteria(result)

# model 8
load("/run/user/1000/gvfs/smb-share:server=cyfiles.iastate.edu,share=09/22/ntyet/R/RA/Data/Reanalysis Data/resultcbc/Model16.Line.RINb.RINa.llymp.lneut.dateGD/Model16_result.RData")
sel_criteria(result)
