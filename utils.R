install.packages("pacman")
pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra,brms,htmlTable,binom,scales)
remotes::install_github("njtierney/mmcc")
library(mmcc)

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Load and clean data
dat <- read_xlsx("SARS_CoV2_data_3_waves_15dec2021.xlsx",sheet = 2) %>% 
  pivot_longer(CoV2S_1w_m:delta_3w_m) %>% 
  mutate(variant=case_when(str_detect(name,"CoV2")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta"),
         variant=fct_relevel(variant,"WT","Beta","Delta"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name) #%>% filter(value>1)

#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50",mult=1){
  #browser()
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(parameter==!!thresh)
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter%in%c("a","b","c","tm")) %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    mutate(x_pred=exp(-(log((c-y$median*mult)/(y$median*mult-a))/b)+tm)) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975),na.rm=T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(values_from = x,names_from = q)
  
  y %>% bind_cols(x)
}



run_model <- function(data.list){
  
  model <- function(){
    
    #Likelihood
    for(i in 1:n)  
    {
      
      S[i] <- a + (c-a)/(1+exp(-b*(titre[i]-tm)))
      
      Y[i]~dbern(S[i])
      
    }
    
    #Priors
    a~dbeta(1, 1) 
    
    b~dlnorm(1,tau1)
    tau1 <- pow(sigma1,-2)
    sigma1 ~ dunif(0.1,2)
    
    c~dbeta(1, 1)
    
    tm~dlnorm(1,tau2)
    tau2 <- pow(sigma2,-2)
    sigma2 ~ dunif(0.1,2)
    
    exp_tm <- exp(tm)
    
    thresh_80 <- 0.8*(c-a)+a
    thresh_50 <- 0.5*(c-a)+a
  }
  
  inits <-  function(){list(
    a=runif(1,0.5,0.99),
    #b=runif(1,0.01,0.99),
    c=runif(1,0.01,0.5))}
  
  jags(model.file=model,
       data=data.list,
       parameters.to.save=c("S",
                            "b",
                            "a",
                            "c",
                            "tm",
                            "exp_tm",
                            "thresh_50",
                            "thresh_80",
                            "mu",
                            "lambda1",
                            "lambda2"),
       #inits=inits,
       n.iter=1e5,
       n.chains = 2,
       n.burnin = 1e4)
  
}

remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}
