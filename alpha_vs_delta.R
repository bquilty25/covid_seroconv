pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra,brms)
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

dat <- dat %>% 
  filter(variant == "WT" & wave == 2 | variant=="Delta" & wave==3) %>% 
  select(-variant) %>% 
  pivot_wider(values_from = value,names_from = wave) %>% 
  mutate(increase_3_v_2=as.integer(`3`>`2`*1.1)) #increase if >10% increase in titer

dat %>%
  pivot_longer(c(`1`,`2`,`3`)) %>% 
  ggplot(aes(x=name,y=value,group=pid_child))+
  geom_point(alpha=0.2)+
  geom_path(alpha=0.2)+
  #facet_wrap(~name,scale="free_x",labeller = labeller(`fct_rev(age)`=Hmisc::capitalize))+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))+
  scale_colour_brewer(type="qual",guide=F,direction=-1)


# Code from https://github.com/xbouteiller/GompertzFit, 
# modified to include a Bernoulli likelihood rather than Normal

# S = Probability of seroconversion
# D = Maximum probability of seroconversion for a given antibody titre
# tm = threshold titre for seroconversion


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
    
    #b~dnorm(0,lambda1)
    #lambda1~dgamma(0.001,0.001)
    
    b~dlnorm(1,tau1)
    tau1 <- pow(sigma1,-2)
    sigma1 ~ dunif(0.1,2)
    
    c~dbeta(1, 1)
    
    tm~dlnorm(1,tau2)
    tau2 <- pow(sigma2,-2)
    sigma2 ~ dunif(0.1,2)
    
    #tm~dnorm(0,lambda2)
    #lambda2~dgamma(0.001,0.001)
    
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


#Run model for post-wave 3
wave3_dat <- dat %>% 
  #filter(!(`2`<1&`3`<1)) %>% 
  drop_na(`2`,`3`)

(wave3_plot <- wave3_dat %>% 
    pivot_longer(c(`2`,`3`)) %>% 
    filter(name!=1) %>% 
    ggplot(aes(x=name,y=value,group=pid_child,colour=factor(increase_3_v_2)))+
    geom_point(alpha=0.2)+
    geom_path(alpha=0.2)+
    #facet_wrap(~name,scale="free_x",labeller = labeller(`fct_rev(age)`=Hmisc::capitalize))+
    scale_y_log10("Anti-Spike IgG titre")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = NA))
)

pred_t_wave3 <- lseq(from=min(wave3_dat$`2`),to=max(wave3_dat$`2`),length.out = 1000)
wave3_res <- run_model(data.list = list(n=nrow(wave3_dat)+length(pred_t_wave3),
                                        Y=c(as.integer(wave3_dat$increase_3_v_2),rep(NA,length(pred_t_wave3))),
                                        titre=c(log(wave3_dat$`2`),log(pred_t_wave3))))

mcmcplot(wave3_res,random = T)
#saveRDS(wave3_res,"wave3_mcmc.rds")

#wave2_plot+wave3_plot&theme(legend.position = "bottom")&scale_y_log10("Titre",limits=c(NA,10000))&labs(x="After wave x")&scale_colour_brewer("Increase post-wave",type="qual",palette = "Set1",labels=c("No","Yes"))

#ggsave("change_plot.png",width=210,height=150,units="mm",dpi=600,bg="white")

#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50"){
  #browser()
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(str_detect(parameter,thresh))
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter%in%c("a","b","c","tm")) %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    mutate(x_pred=exp(-(log((c-y$median)/(y-a))/b)+tm)) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975),na.rm=T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(values_from = x,names_from = q)
  
  y %>% bind_cols(x)
}

(wave3_plot <- mmcc::tidy(as.mcmc(wave3_res)) %>% 
    filter(str_detect(parameter,"S")) %>% 
    select(median,`2.5%`,`97.5%`) %>% 
    #mutate_all(function(x)1-x) %>% 
    slice(nrow(wave3_dat)+1:n()) %>% 
    bind_cols(pred_t_wave3) %>% 
    ggplot()+
    geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#8856a7", fill="#9ebcda", stat="identity")+
    geom_point(data=wave3_dat,aes(y=increase_3_v_2, x=`2`),pch=124,size=5,alpha=0.5,colour="#8856a7")+
    geom_pointrange(data=extract_ab_thresholds(wave3_res,"thresh_50"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
    geom_pointrange(data=extract_ab_thresholds(wave3_res,"thresh_80"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
    geom_hline(data=mmcc::tidy(as.mcmc(wave3_res)) %>%
                 filter(parameter%in%c("a","c")),
               aes(yintercept=median),
               linetype="dashed")+
    scale_x_log10("Beta S-specific IgG pre-wave (WHO BAU/ml)" )+
    scale_y_continuous("Probability of increased Delta S-specific IgG titres following Delta wave",labels = scales::percent)+
    ggtitle("Delta wave")
)

