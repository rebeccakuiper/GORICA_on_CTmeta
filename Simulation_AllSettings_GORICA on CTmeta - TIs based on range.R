
# Simulation based on Example_GORICAonCTmeta.R


# Prelimanary - based on Example_CTmeta.R
if (!require("expm")) install.packages("expm")
library(expm)
if (!require("metafor")) install.packages("metafor")
library(metafor)
if (!require("tsDyn")) install.packages("tsDyn")
library(tsDyn)
if (!require("vars")) install.packages("vars")
library(vars)
#
if (!require("restriktor")) install.packages("restriktor") # install this package first (once)
library(restriktor)


NrSim <- 1000 #1000

q <- 2 # Number of variables
S <- 25 # = length(N)

# Time intervals ('time lags'), using Work engagement paper
TI <- matrix(c(
  12,
  12,
  12,
  4,
  9,
  12,
  12,
  2,
  36,
  12,
  12,
  8,
  24,
  12*(1/365),
  24,
  12,
  2,
  14,
  12,
  10,
  48,
  12,
  12,
  1,
  8), byrow=TRUE)
TI <- TI/12 # now in years
G <- length(unique(TI))

#Hypothesis of interest
H1 <- "Phi12 < Phi21"
# vs complement: Hc <- "Phi12 > Phi21"


# Output: (True) Hypothesis Rates and descr stat of GORICA weights
HypoNames <- c("H1", "complement")
TI_Names <- paste0("TI_", unique(TI)[order(unique(TI))])
ES_Names <- c("large", "medium", "small", "zero")
N_Names <- c("Example", "100", "100and1x500", "100and2x500", "100and1x1000", "100and2x1000", "100and3x1000")
HR <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
            dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
meanGORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                           dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
sdGORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                         dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
minGORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                          dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
maxGORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                          dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
q05GORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                          dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
q95GORICAweights <- array(NA, dim = c(length(HypoNames), G, length(ES_Names), length(N_Names)), 
                          dimnames = list(HypoNames, TI_Names, ES_Names, N_Names))
     
# Simulation over Settings for ES (Phi) and N     
set.seed(123)
for(teller_ES in 1:length(ES_Names)){
  for(teller_N in 1:length(N_Names)){
    # teller_ES = 1; teller_N = 1
    
    print(paste0("teller_ES = ", teller_ES, " out of ", length(ES_Names)))
    print(paste0("teller_N = ", teller_N, " out of ", length(N_Names)))
    
    # ES specification, that is, specification of population values of Phi
    if(teller_ES == 1){    
      Phi_pop <- matrix(c(0.50, 0.15,
                          0.25, 0.40), nrow=q, byrow=TRUE) # population Phi matrix (i.e., lagged effects matrix)
      ESName <- "large"
    }
    if(teller_ES == 2){    
      Phi_pop <- matrix(c(0.50, 0.15,
                          0.20, 0.40), nrow=q, byrow=TRUE) # population Phi matrix (i.e., lagged effects matrix)
      ESName <- "medium"
    }
    if(teller_ES == 3){    
      Phi_pop <- matrix(c(0.50, 0.15,
                          0.17, 0.40), nrow=q, byrow=TRUE) # population Phi matrix (i.e., lagged effects matrix)
      ESName <- "small"
    }
    if(teller_ES == 4){    
      Phi_pop <- matrix(c(0.50, 0.15,
                          0.15, 0.40), nrow=q, byrow=TRUE) # population Phi matrix (i.e., lagged effects matrix)
      ESName <- "zero"
    }
    #
    #
    A_pop <- logm(Phi_pop) # underlying population drift matrix
    #Phi_pop <- expm(A_pop)
    vecPhi_pop <- as.vector(t(Phi_pop)) # vectorized population Phi
    Gamma_pop <- matrix(c(1,    0.3,
                          0.3,    1), nrow=q, byrow=TRUE) # population stationary covariance matrix
    # Since Gamma is now a correlation matrix, Phi is already standardized
    SigmaVAR_pop <- Gamma_pop - Phi_pop %*% Gamma_pop %*% t(Phi_pop) # population residual covariance matrix
    
    #######################################


    # study-specific sample size
    if(teller_N == 1){    
      N <- matrix(c(
        643,
        651,
        473,
        387,
        187,
        209,
        2897,
        160,
        1964,
        848,
        926,
        274,
        433,
        256,
        409,
        926,
        162,
        262,
        247,
        102,
        171,
        201,
        309,
        77,
        67), byrow=TRUE)
      NName <- "Example"
      #sum(N)
    }
    #
    if(teller_N > 1){  # teller_N == 2
      N <- rep(100, 25)
      NName <- "100" 
      #sum(N)
    }
    if(teller_N > 2){ # teller_N == 3
      N[1] <- 500
      NName <- "100and1x500"
      #sum(N)
    }
    if(teller_N > 3){ # teller_N == 4
      N[2] <- 500
      NName <- "100and2x500"
      #sum(N)
    }
    #
    if(teller_N > 4){ # teller_N == 5
      N <- rep(100, 25)
      N[1] <- 1000
      NName <- "100and1x1000"
      #sum(N)
    }
    if(teller_N > 5){ # teller_N == 6
      N[2] <- 1000
      NName <- "100and2x1000"
      #sum(N)
    }
    if(teller_N > 6){ # teller_N == 7
      N[3] <- 1000
      NName <- "100and3x1000"
      #sum(N)
    }

    
    #######################################################################################################################################
    
    ## Simulation based on Example ##
    GORICAweights_g <- array(NA, dim = c(NrSim, 2, G))
    #GORICAweights_g_s <- array(NA, dim = c(NrSim, 2, S, G)) 
    
    
    # Needed in the meta-an
    sub = NULL
    for(i in 1:q){
      sub = c(sub, paste(i, 1:q, sep=""))
    }
    outcome <- rep(sub, S) 
    
    
    # Simulate data and do the meta-analyses
    NotPosDefCovMx <- array(NA, dim=c(NrSim,S))
    #set.seed(123)
    for(sim in 1:NrSim){
      # sim <- 1
    
      # Example
      
      # Sample for each study s in the meta-analysis data based on population values and N_s and DeltaT_s
      S <- length(N)
      vecA <- array(NA, dim = c(S*q*q))
      vecPhi <- array(NA, dim = c(S*q*q))
      CovMx <- array(0, dim = c(S*q*q, S*q*q))
      vecPhi1 <- array(NA, dim = c(S*q*q))
      CovMx1 <- array(0, dim = c(S*q*q, S*q*q))
      #
      G <- length(unique(TI))
      vecPhi_G <- array(NA, dim = c(S*q*q, G))
      CovMx_G <- array(0, dim = c(S*q*q, S*q*q, G))
      Phi_G <- array(NA, dim = c(S*q, q, G))
      SigmaVAR_G <- array(NA, dim = c(S*q, q, G))
      #
      s <- 1
      while(s <= S){
        Y_mu <- VAR.sim(Phi_pop, n = N[s], lag = 1, include = c("none"), starting = NULL, varcov = SigmaVAR_pop)
        Y <- scale(Y_mu, scale = FALSE) # substracts the means (so , it centers, does not standardize now because of scale = F)
        colnames(Y) <- paste0("y", 1:q)
        outcomeVAR <- VAR(Y, p = 1)
        Phi_VARest <- Acoef(outcomeVAR)[[1]]
        if(any(Re(eigen(Phi_VARest)$values) < 0)){
          s <- s # No CTM-equivalent, so do not proceed
        }else{  
          SigmaVAR1_VARest <- cov(resid(outcomeVAR))
          Gamma_VARest <- cov(Y)
          # Standardize parameters! 
          Sxy <- sqrt(diag(diag(Gamma_VARest)))
          Gamma_VARest <- solve(Sxy) %*% Gamma_VARest %*% solve(Sxy) 
          Phi_VARest <- solve(Sxy) %*% Phi_VARest %*% Sxy
          SigmaVAR1_VARest <- solve(Sxy) %*% SigmaVAR1_VARest %*% solve(Sxy) 
          #
          invGamma <- solve(Gamma_VARest)
          B_VARest <- -logm(Phi_VARest)/1
          vecA[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(-B_VARest))
          #
          Phi_VARest <- expm(-B_VARest*TI[s])
          vecPhi[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
          SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
          CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
          if(any( eigen( CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
            s <- s # Cov mx should be pos def
          }else{
            #
            Phi_VARest <- expm(-B_VARest)
            vecPhi1[((s-1)*(q*q)+1):(s*q*q)] <- as.vector(t(Phi_VARest))
            SigmaVAR1_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
            CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] <- kronecker(SigmaVAR1_VARest, invGamma) / (N[s]-q)
            # this equates
            #CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] * kronecker(SigmaVAR1_VARest / SigmaVAR_VARest, matrix(1,q,q))
            if(any( eigen( CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] )$values < 0 )){
              s <- s # Cov mx should be pos def
            }else{
              for(g in 1:G){
                Phi_VARest <- expm(-B_VARest*unique(TI)[g])
                vecPhi_G[((s-1)*(q*q)+1):(s*q*q),g] <- as.vector(t(Phi_VARest))
                Phi_G[((s-1)*(q)+1):(s*q), 1:q, g] <- Phi_VARest
                #
                SigmaVAR_VARest <- Gamma_VARest - Phi_VARest %*% Gamma_VARest %*% t(Phi_VARest)
                CovMx_G[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2),g] = kronecker(SigmaVAR_VARest, invGamma) / (N[s]-q)
                # this equates
                #CovMx1[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] * kronecker(SigmaVAR_VARest / SigmaVAR1_VARest, matrix(1,q,q))
                SigmaVAR_G[((s-1)*(q)+1):(s*q), 1:q, g] <- SigmaVAR_VARest
              }
              if(any( apply(x <- CovMx_G[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2),], 3, function(x){eigen(x)$values}) < 0 )){
                s <- s # Cov mx should be pos def
              }else{
                s <- s+1
              }
            }
          }
        }
      }
      
      ###########################################################################################################
      
      # Multivariate Meta-analyses
      if( any( abs(CovMx[lower.tri(diag(S*q^2))] - CovMx[upper.tri(diag(S*q^2))]) < 0.1^10 ) ){
        CovMx <- (CovMx + t(CovMx))/2
      }
      
      # CTmeta
      G <- length(unique(TI))
      Phi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
      sePhi_Trans <- matrix(NA, ncol=(q^2), nrow=G)
      CovMxPhi_Trans <- array(NA, dim = c(G, q^2, q^2))
      for(g in 1:G){
        # g <- 1
        if( any( abs(CovMx_G[,,g][lower.tri(diag(S*q^2))] - CovMx_G[,,g][upper.tri(diag(S*q^2))]) < 0.1^10 ) ){
          CovMx_G[,,g] <- (CovMx_G[,,g] + t(CovMx_G[,,g]))/2
        }
        metaan <- rma.mv(yi=vecPhi_G[,g], V=CovMx_G[,,g], mods = ~ outcome - 1, method = "FE") 
        Phi_Trans[g,] <- coef(metaan)
        sePhi_Trans[g,] <- metaan$se
        CovMxPhi_Trans[g,,] <- metaan$vb
      }
      
      
      # Apply GORICA
      # On overall Phi
      for(g in 1:G){
        #g <- 1
        est <- Phi_Trans[g,]
        names(est) <- c("Phi11", "Phi12", "Phi21", "Phi22")
        VCOV <- CovMxPhi_Trans[g,,]
        results_g <- goric(est, VCOV = VCOV, H1, comparison = "complement", type = "gorica") 
        #summary(results_g)
        GORICAweights_g[sim, , g] <- results_g$result[,5]
      }
      ## On study-specific Phi; on Phi from each of the primary studies
      #for(g in 1:G){
      #  for(s in 1:S){
      #    #g <- 1; s <- 1
      #    est <- vecPhi_G[((s-1)*(q*q)+1):(s*q*q), g]
      #    names(est) <- c("Phi11", "Phi12", "Phi21", "Phi22")
      #    VCOV <- CovMx_G[((s-1)*(q*q)+1):(s*q*q), ((s-1)*(q*q)+1):(s*q*q), g]
      #    results_g_s <- goric(est, VCOV = VCOV, H1, comparison = "complement", type = "gorica") 
      #    #summary(results_g_s)
      #    GORICAweights_g_s[sim, , s, g] <- results_g_s$result[,5]
      #  }
      #}
    
    } # End simulation CTmeta and GORICA


    gw <- GORICAweights_g[,,order(unique(TI))]
    
    # (True) Hypothesis Rates
    HR[1, , teller_ES, teller_N] <- apply(gw[,1,] > gw[,2,], 2, mean)
    HR[2, , teller_ES, teller_N] <- apply(gw[,2,] > gw[,1,], 2, mean)
    #
    # mean weights
    meanGORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), mean)
    sdGORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), sd)
    minGORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), min)
    maxGORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), max)
    q05GORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), quantile, probs = .05)
    q95GORICAweights[ , , teller_ES, teller_N] <- apply(gw[,,], c(2,3), quantile, probs = .95)
    #
    #
    #HR
    #meanGORICAweights
    #sdGORICAweights
    #minGORICAweights 
    #maxGORICAweights
    #q05GORICAweights 
    #q95GORICAweights
    
    #--- Save and Load (for new data)
    save(list = ls(), file = paste0("Sim_GORICAonCTmeta_ES", ESName, "_N", NName, "_NrSim", NrSim, ".RData"))  # saves all objects including data in a compressed format; use extensions like .rda, .Rda, .RData etc.

  } # end loop N
} # end loop ES


# Plots
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("hrbrthemes")) install.packages("hrbrthemes")
library(hrbrthemes)
if (!require("extrafont")) install.packages("extrafont")
library(extrafont)
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("cowplot")) install.packages("cowplot")
library(cowplot)
if (!require("jtools")) install.packages("jtools")
library(jtools)
if (!require("DataCombine")) install.packages("DataCombine")
library(DataCombine)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("ggpubr")) install.packages("ggpubr")
library(ggpubr)
if (!require("scales")) install.packages("scales")
library(scales)
if (!require("ggsci")) install.packages("ggsci")
library(ggsci)

# TO DO make plots
#HR
#meanGORICAweights
#sdGORICAweights
#minGORICAweights 
#maxGORICAweights
#q05GORICAweights 
#q95GORICAweights

indexHtrue <- 1 # In case ES is zero it is both H1 and H2...
TI_ordered <- unique(TI)[order(unique(TI))]
goricresult <- data.frame(criterion = NA, thr = NA, meanWeight = NA, q05Weight = NA, q95Weight = NA, Htrue = NA, TI = NA, n = NA, ES_pop = NA)
for(teller_ES in 1:length(ES_Names)){
  for(teller_N in 1:length(N_Names)){ 
    for(teller_TI in 1:G){ 
      goricresult[nrow(goricresult) + 1,] = list("goric", HR[indexHtrue, teller_TI, teller_ES, teller_N], meanGORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], q05GORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], q95GORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], indexHtrue, TI_ordered[teller_TI], N_Names[teller_N], ES_Names[teller_ES]) 
      #goricresult[nrow(goricresult) + 1,] = list("goric", HR[indexHtrue, teller_TI, teller_ES, teller_N], meanGORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], q05GORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], q95GORICAweights[indexHtrue, teller_TI, teller_ES, teller_N], indexHtrue, TI_ordered[teller_TI], N_Names[teller_N], (length(ES_Names)-teller_ES))   
    }
  }
}
goricresult <- goricresult[-1,] # Note: delete first row here, because it consists of NAs
#goricresult


# Plot per ES
for(teller_ES in 1:length(ES_Names)){
  #teller_ES <- 1
  #teller_ES <- 2
  #teller_ES <- 3
  #teller_ES <- 4
  #
  # Filter data for that ES
  goricresult_ES <- goricresult[which(goricresult$ES_pop == ES_Names[teller_ES]),]
  #goricresult_ES <- goricresult[which(goricresult$ES_pop == (length(ES_Names)-teller_ES)),] 
  #
  if(teller_ES == 4){ # ES is zero and then H1 and compl both true
    title_HR <- paste0("Hypothesis rates for H1, for a ",  ES_Names[teller_ES], " effect size")
    ylab_HR <- paste0("Hypothesis rate  for H1")
    title_w <- paste0("Mean GORICA weights for H1 with 95% interval, for a ", ES_Names[teller_ES], " effect size")
    ylab_w <- paste0("Mean and 95% range of GORICA weights for H1")
  }else{
    title_HR <- paste0("True hypothesis rates for H1, for a ",  ES_Names[teller_ES], " effect size")
    ylab_HR <- paste0("True hypothesis rate")
    title_w <- paste0("Mean GORICA weights for H1 with 95% interval, for a ", ES_Names[teller_ES], " effect size")
    ylab_w <- paste0("Mean and 95% range of GORICA weights for H1")
  }
  #
  #min <- min(goricresult_ES$thr)
  #max <- max(goricresult_ES$thr)
  min <- 0
  max <- 1
  ggplot(data = goricresult_ES) +
    geom_line(mapping = aes(x = TI, y = thr, colour = n)) +
    theme_apa() +
    geom_point(mapping = aes(x = TI, y = thr, colour = n)) +
    theme_apa() +
    ylim(min, max) +
    theme(legend.position = "bottom") +
    ggtitle(title_HR) +
    scale_color_jco() +
    xlab("Targeted time interval") +
    ylab(ylab_HR) +
    theme(text = element_text(family = "Arial", size=13))
  # save
  ggsave(paste0("Plot_THR_ES", ES_Names[teller_ES], "_NrSim", NrSim, ".png"), plot = last_plot(), device = "png") #voor .png format
  ggsave(paste0("Plot_THR_ES", ES_Names[teller_ES], "_NrSim", NrSim, ".pdf"), plot = last_plot(), device = cairo_pdf) #voor .pdf format
  # Close the pdf and jpg file
  dev.off() 
  #
  #
  ggplot(data = goricresult_ES) +
    geom_line(mapping = aes(x = TI, y = meanWeight, colour = n)) +
    theme_apa() +
    geom_point(mapping = aes(x = TI, y = meanWeight, colour = n)) +
    theme_apa() +
    geom_line(mapping = aes(x = TI, y = q05Weight, colour = n)) +
    theme_apa() +
    geom_line(mapping = aes(x = TI, y = q95Weight, colour = n)) +
    theme_apa() +
    ylim(0, 1) +
    theme(legend.position = "bottom") +
    ggtitle(title_w) +
    scale_color_jco() +
    xlab("Targeted time interval") +
    ylab(ylab_w) +
    theme(text = element_text(family = "Arial", size=13))
  # save
  ggsave(paste0("Plot_weights_ES", ES_Names[teller_ES], "_NrSim", NrSim, ".png"), plot = last_plot(), device = "png") #voor .png format
  ggsave(paste0("Plot_weights_ES", ES_Names[teller_ES], "_NrSim", NrSim, ".pdf"), plot = last_plot(), device = cairo_pdf) #voor .pdf format
  # Close the pdf and jpg file
  dev.off() 
}

