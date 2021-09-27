pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra)
remotes::install_github("njtierney/mmcc")
library(mmcc)

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Load and clean data
dat <- read_csv("all_SARS_CoV2_data_03EP2021_updated.csv") %>%
  select(-c(pid_child, mat_pair, chid_pair)) %>%
  rename("child_wave_1_titre"=CoV2S_1wave_child,
         "child_wave_2_titre"=CoV2S_2wave_child,
         "mother_wave_1_titre"=CoV2S_1wave_mat,
         "mother_wave_2_titre"=CoV2S_2wave_mat,
         "child_change"=child_wave1_wave2_incr_decr,
         "mother_change"=mat_w1_w2_incr_decr)

# dat <- read_csv("all_SARS_CoV2_data_17SEP2021_updated.csv") %>% 
#   select(-c(pid_child, mat_pair, chid_pair)) %>% 
#   select(-c(1:6)) %>% 
#   rename("child_wave_1_titre"=SB1351_1wave_child,
#          "child_wave_2_titre"=SB1351_2wave_child,
#          "mother_wave_1_titre"=SB1351_1wave_mat,
#          "mother_wave_2_titre"=SB1351_2wave_mat) %>% 
#   mutate(child_change=ifelse(child_wave_2_titre>child_wave_1_titre,"Increase","Decease or stays the same"),
#          mother_change=ifelse(mother_wave_2_titre>mother_wave_1_titre,"Increase","Decease or stays the same")) %>% 
#   select(sort(names(.)))

dat_child <- dat %>% select(1:3) %>% rename_all(~stringr::str_replace(.,"^child_","")) %>%  mutate(age="child") %>% 
  drop_na(change) %>% mutate(change=as.factor(change),change=fct_rev(change))  #%>% filter(wave_1_titre>1)

dat_mother <- dat %>% select(4:6) %>% rename_all(~stringr::str_replace(.,"^mother_","")) %>% mutate(age="mother")%>% 
  drop_na(change) %>% mutate(change=as.factor(change),change=fct_rev(change)) #%>% filter(wave_1_titre>1)


bind_rows(dat_mother,dat_child) %>% 
  rename("Wave 1"=wave_1_titre,"Wave 2"=wave_2_titre) %>% 
  pivot_longer(cols=c(`Wave 1`,`Wave 2`)) %>% 
  ggplot(aes(x=name,y=value,group=age,colour=age))+
  geom_point(alpha=0.5)+
  geom_path(alpha=0.5)+
  facet_wrap(~fct_rev(age),labeller = labeller(`fct_rev(age)`=Hmisc::capitalize))+
  scale_x_discrete("")+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))+
  scale_colour_brewer(type="qual",guide=F,direction=-1)

ggsave("change.png",width=210,height=150,units="mm",dpi=600)

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
    #S[i] <- a/(1+exp(-b*(t[i]-tm)))
    logit(S[i]) <- a + b*titre[i]
    
    Y[i]~dbern(S[i])
    
  }
  
  
  mu.t ~ dnorm(0, .0000001)
  prec.t ~ dgamma(.01, .01)
  

  min_Y <- min(S)
  max_Y <- max(S)
  reduction <- 1-(max_Y/min_Y)

  thresh_50 <- -min_Y*0.5+max_Y*0.5+min_Y
  
  thresh_60 <- -min_Y*0.6+max_Y*0.6+min_Y
  thresh_70 <- -min_Y*0.7+max_Y*0.7+min_Y
  thresh_80 <- -min_Y*0.8+max_Y*0.8+min_Y
  thresh_90 <- -min_Y*0.9+max_Y*0.9+min_Y
  thresh_95 <- -min_Y*0.95+max_Y*0.95+min_Y
  

  #Priors
  a~dbeta(1, 1) 
  
  #c~dbeta(1,1)
  
  b~dnorm(0, lambda1)
  lambda1~dgamma(0.0001,0.0001) #uninformative Gamma prior
  
   #tm~dnorm(0, lambda2)
   #exp_tm <- exp(tm)
   #lambda2~dgamma(0.0001,0.0001) #uninformative gamma prior
  
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
                                     "max_Y"),
                #inits=inits,
                n.iter=2.5e4,
                n.chains = 2,
                n.burnin = 2.5e3)

}

pred_t <- lseq(from=0.55,to=1000,length.out = 1000)


#Run model for mother

mother_res <- run_model(data.list = list(n=nrow(dat_mother)+length(pred_t),
                                         Y=c(as.integer(dat_mother$change)-1,rep(NA,length(pred_t))),
                                         titre=c(log(dat_mother$wave_1_titre),log(pred_t))))

saveRDS(mother_res,"mother_mcmc.rds")

#Run model for children
child_res <-  run_model(data.list = list(n=nrow(dat_child)+length(pred_t),
                                         Y=c(as.integer(dat_child$change)-1,rep(NA,length(pred_t))),
                                         titre=c(log(dat_child$wave_1_titre),log(pred_t))))
saveRDS(child_res,"child_mcmc.rds")

#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50"){
  
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(str_detect(parameter,thresh))
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter=="a"|parameter=="b") %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    mutate(x_pred=exp((-a+boot::logit(y$median)/b))) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(values_from = x,names_from = q)
  
  y %>% bind_cols(x)
}

mother_thresh_50 <- extract_ab_thresholds(mother_res,"thresh_50")

mother_plot <- mmcc::mcmc_to_dt(as.mcmc(mother_res)) %>% 
  filter(parameter=="a"|parameter=="b") %>% 
  pivot_wider(values_from = "value",names_from = "parameter") %>% 
  crossing(y_vals=seq(0,1,length.out=100)) %>% 
  mutate(x_pred=exp((-a+boot::logit(y_vals)/b))) %>% 
  group_by(y_vals) %>% 
  summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(values_from = x,names_from = q) %>% 
  ggplot()+
  geom_smooth(aes(x=1-y_vals,ymin=`0.025`, ymax=`0.975`,y=`0.5`),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_point(data=dat_mother,aes(x=1/1-as.integer(change)+1, y=wave_1_titre),pch=124,size=4,alpha=0.5,colour="#2ca25f")+
  geom_pointrange(data=mother_thresh_50,aes(ymin=`0.025`,ymax=`0.975`,y=`0.5`,x=1-median))+
  scale_y_log10("WT S-specific IgG following Wave 1 (WHO BAU/ml)" )+
  scale_x_continuous("Probability of increased WT S-specific IgG titres following Wave 2")+
  coord_flip(ylim=c(0.5,1000))+
  ggtitle("Mothers")


#child threshold
child_thresh_50 <- extract_ab_thresholds(child_res,"thresh_50")

child_plot <- mmcc::mcmc_to_dt(as.mcmc(child_res)) %>% 
  filter(parameter=="a"|parameter=="b"|parameter=="tm") %>% 
  pivot_wider(values_from = "value",names_from = "parameter") %>% 
  crossing(y_vals=seq(0,1,length.out=100)) %>% 
  mutate(x_pred=exp((-a+boot::logit(y_vals)/b))) %>% 
  group_by(y_vals) %>% 
  summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
  pivot_wider(values_from = x,names_from = q) %>% 
  ggplot()+
  geom_smooth(aes(x=1-y_vals,ymin=`0.025`, ymax=`0.975`,y=`0.5`),colour="#8856a7", fill="#9ebcda", stat="identity")+
  geom_point(data=dat_child,aes(x=1/1-as.integer(change)+1, y=wave_1_titre),pch=124,size=4,alpha=0.5,colour="#8856a7")+
  geom_pointrange(data=child_thresh_50,aes(ymin=`0.025`,ymax=`0.975`,y=`0.5`,x=1-median))+
  scale_y_log10("WT S-specific IgG following Wave 1 (WHO BAU/ml)")+
  scale_x_continuous("Probability of increased WT S-specific IgG titres following Wave 2")+
  coord_flip(ylim=c(0.5,1000))+
  ggtitle("Children")

mother_plot+child_plot&theme_minimal()&theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(fill = NA))&scale_fill_brewer()

ggsave("res_logistic.png",width=210,height=150,units="mm",dpi=600)
