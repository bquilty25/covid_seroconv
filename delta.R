pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra)
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
  select(-name)

dat %>%
  ggplot(aes(x=name,y=value,group=pid_child))+
  geom_point(alpha=0.2)+
  geom_path(alpha=0.2)+
  facet_wrap(~variant,scale="free_x",labeller = labeller(`fct_rev(age)`=Hmisc::capitalize))+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))+
  scale_colour_brewer(type="qual",guide=F,direction=-1)

ggsave("change_new.png",width=400,height=150,units="mm",dpi=600,bg="white")

dat %<>% 
  filter(variant == "WT" & wave == 1 | variant == "Beta" & wave == 2 | variant=="Delta" & wave==3) %>% 
  select(-variant) %>% 
  pivot_wider(values_from = value,names_from = wave) %>% 
  mutate(increase_2_v_1=as.integer(`2`>`1`*1.1),
         increase_3_v_2=as.integer(`3`>`2`*1.1)) #increase if >10% increase in titer


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
      
      #S[i] <- a + c * exp(-exp(-b*(t[i])))
      S[i] <- a/(1+exp(-b*(titre[i]-tm)))
      #logit(S[i]) <- a + b*titre[i]
      
      Y[i]~dbern(S[i])
      
    }
    
    
    mu.t ~ dnorm(0, .0000001)
    prec.t ~ dgamma(.01, .01)
    
    
    min_Y <- min(S)
    max_Y <- max(S)
    
    thresh_50 <- -min_Y*0.5+a*0.5+min_Y
    
    thresh_80 <- -min_Y*0.8+a*0.8+min_Y
    
    reduction <- 1-min_Y/a
    
    maximal <- (1-a)*100
    
    minimal <- (1-min_Y)*100
    
    #Priors
    a~dbeta(1, 1) 
    
    
    b~dnorm(0, lambda1)
    lambda1~dgamma(0.0001,0.0001) #uninformative Gamma prior
    
    tm~dnorm(0, lambda2)
    exp_tm <- exp(tm)
    lambda2~dgamma(0.0001,0.0001) #uninformative gamma prior
    
  }
  
  inits <-  function(){list(
    a=runif(1,-5,5),
    b=runif(1,-5,5))}
  
  jags(model.file=model,
       data=data.list,
       parameters.to.save=c("S",
                            "b",
                            "a",
                            "c",
                            "tm",
                            "exp_tm",
                            "reduction",
                            "lambda1",
                            "thresh_50",
                            "thresh_60",
                            "thresh_70",
                            "thresh_80",
                            "thresh_90",
                            "thresh_95",
                            "min_Y",
                            "max_Y",
                            "maximal",
                            "minimal"),
       #inits=inits,
       n.iter=2.5e4,
       n.chains = 2,
       n.burnin = 2.5e3)
  
}

pred_t <- lseq(from=0.55,to=1000,length.out = 1000)


#Run model for post-wave 2
wave2_dat <- dat %>% drop_na(`1`,`2`)
pred_t_wave2 <- lseq(from=min(wave2_dat$`1`),to=max(wave2_dat$`1`),length.out = 1000)
wave2_res <- run_model(data.list = list(n=nrow(wave2_dat)+length(pred_t_wave2),
                                         Y=c(as.integer(wave2_dat$increase_2_v_1),rep(NA,length(pred_t_wave2))),
                                         titre=c(log(wave2_dat$`1`),log(pred_t_wave2))))

saveRDS(wave2_res,"wave2_mcmc.rds")

#Run model for post-wave 3
wave3_dat <- dat %>% drop_na(`2`,`3`)
pred_t_wave3 <- lseq(from=min(wave3_dat$`2`),to=max(wave3_dat$`2`),length.out = 1000)
wave3_res <- run_model(data.list = list(n=nrow(wave3_dat)+length(pred_t_wave3),
                                        Y=c(as.integer(wave3_dat$increase_3_v_2),rep(NA,length(pred_t_wave3))),
                                        titre=c(log(wave3_dat$`2`),log(pred_t_wave3))))

saveRDS(wave3_res,"wave3_mcmc.rds")


#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50"){
  
  browser()
  
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(str_detect(parameter,thresh))
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter%in%c("a","b","tm")) %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    #mutate(x_pred=exp((-a+boot::logit(y$median)/b))) %>% 
    mutate(x_pred=exp(tm-(log((a/y$median)-1)/b))) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975),na.rm=T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(values_from = x,names_from = q)
  
  y %>% bind_cols(x)
}

wave2_thresh_50 <- extract_ab_thresholds(wave2_res,"thresh_50")
wave2_thresh_80 <- extract_ab_thresholds(wave2_res,"thresh_80")

wave2_plot <- mmcc::tidy(as.mcmc(wave2_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(wave2_dat)+1:n()) %>% 
  bind_cols(pred_t_wave2) %>% 
  ggplot()+
  geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_point(data=wave2_dat,aes(y=as.integer(increase_2_v_1), x=`1`),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
  geom_pointrange(data=wave2_thresh_50,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=median))+
  geom_hline(data=mmcc::tidy(as.mcmc(wave2_res)) %>% 
               filter(str_detect(parameter,"maximal")),
             aes(yintercept=median/100),
             linetype="dashed")+
  scale_x_log10("WT S-specific IgG following Wave 2 (WHO BAU/ml)" )+
  scale_y_continuous("Probability of increased WT S-specific IgG titres following Wave 3",labels = scales::percent)+
  ggtitle("Wave 3")


wave3_thresh_50 <- extract_ab_thresholds(wave3_res,"thresh_50")
wave3_thresh_80 <- extract_ab_thresholds(wave3_res,"thresh_80")

(wave3_plot <- mmcc::tidy(as.mcmc(wave3_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(wave3_dat)+1:n()) %>% 
  bind_cols(pred_t_wave3) %>% 
  ggplot()+
  geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_point(data=wave3_dat,aes(y=increase_3_v_2, x=`2`),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
  geom_pointrange(data=wave3_thresh_50,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=median))+
  geom_hline(data=mmcc::tidy(as.mcmc(wave3_res)) %>% 
               filter(str_detect(parameter,"maximal")),
             aes(yintercept=median/100),
             linetype="dashed")+
  scale_x_log10("WT S-specific IgG following Wave 1 (WHO BAU/ml)" )+
  scale_y_continuous("Probability of increased WT S-specific IgG titres following Wave 2",labels = scales::percent)+
  ggtitle("Wave 2")
)


mother_plot+child_plot&
  theme_minimal()&
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))&
  scale_fill_brewer()&
  coord_cartesian(xlim=c(0.5,1000),expand=F)

ggsave("res_logistic.png",width=210,height=150,units="mm",dpi=600)

#Table 1
bind_rows(
  mmcc::tidy(as.mcmc(mother_res)) %>% 
    mutate(age="Mothers"),
  mmcc::tidy(as.mcmc(child_res)) %>% 
    mutate(age="Children")) %>% 
  filter(parameter %in% c("maximal","minimal","b")) %>% 
  bind_rows(extract_ab_thresholds(mother_res,thresh="thresh_50") %>% 
              mutate(age="Mothers") %>% select(parameter,age,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(child_res,thresh="thresh_50") %>% 
              mutate(age="Children") %>% select(parameter,age,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(mother_res,thresh="thresh_80") %>% 
              mutate(age="Mothers") %>% select(parameter,age,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(child_res,thresh="thresh_80") %>% 
              mutate(age="Children") %>% select(parameter,age,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  mutate_at(vars(median,`2.5%`,`97.5%`),~round(.,1)) %>% 
  mutate(estimate=paste0(median, " (",`2.5%`,", ",`97.5%`,")")) %>% 
  select(age,parameter,estimate) %>% 
  pivot_wider(values_from = estimate,names_from = parameter) %>% 
  select(age,minimal,maximal,b,thresh_50,thresh_80) %>% 
  htmlTable::htmlTable()

