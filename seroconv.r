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
  
  reduction <- 1-(a/min_Y)

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


#Run model for mother
pred_t_mother <- lseq(from=min(dat_mother$wave_1_titre),to=max(dat_mother$wave_1_titre),length.out = 1000)
mother_res <- run_model(data.list = list(n=nrow(dat_mother)+length(pred_t_mother),
                                         Y=c(as.integer(dat_mother$change)-1,rep(NA,length(pred_t_mother))),
                                         titre=c(log(dat_mother$wave_1_titre),log(pred_t_mother))))

saveRDS(mother_res,"mother_mcmc.rds")

#Run model for children
pred_t_child <- lseq(from=min(dat_child$wave_1_titre),to=max(dat_child$wave_1_titre),length.out = 1000)
child_res <-  run_model(data.list = list(n=nrow(dat_child)+length(pred_t_child),
                                         Y=c(as.integer(dat_child$change)-1,rep(NA,length(pred_t_child))),
                                         titre=c(log(dat_child$wave_1_titre),log(pred_t_child))))
saveRDS(child_res,"child_mcmc.rds")


#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50"){
  
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(str_detect(parameter,thresh))
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter=="a"|parameter=="b"|parameter=="tm") %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    #mutate(x_pred=exp((-a+boot::logit(y$median)/b))) %>% 
    mutate(x_pred=exp(tm-(log((a/y$median)-1)/b))) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975),na.rm=T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(values_from = x,names_from = q)
  
  y %>% bind_cols(x)
}

mother_thresh_50 <- extract_ab_thresholds(mother_res,"thresh_50")
mother_thresh_80 <- extract_ab_thresholds(mother_res,"thresh_80")

mother_plot <- mmcc::tidy(as.mcmc(mother_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(dat_mother)+1:n()) %>% 
  bind_cols(pred_t_mother) %>% 
  ggplot()+
  geom_smooth(aes(ymin=1-`2.5%`, ymax=1-`97.5%`,y=1-median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_point(data=dat_mother,aes(y=1/1-as.integer(change)+1, x=wave_1_titre),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
  geom_pointrange(data=mother_thresh_50,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=1-median))+
  geom_pointrange(data=mother_thresh_80,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=1-median))+
  scale_x_log10("WT S-specific IgG following Wave 1 (WHO BAU/ml)" )+
  scale_y_continuous("Probability of increased WT S-specific IgG titres following Wave 2",labels = scales::percent)+
  ggtitle("Mothers")


#child threshold
child_thresh_50 <- extract_ab_thresholds(child_res,"thresh_50")
child_thresh_80 <- extract_ab_thresholds(child_res,"thresh_80")

child_plot <- mmcc::tidy(as.mcmc(child_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(dat_child)+1:n()) %>% 
  bind_cols(pred_t_child) %>% 
  ggplot()+
  geom_smooth(aes(ymin=1-`2.5%`, ymax=1-`97.5%`,y=1-median,x=...4),colour="#8856a7", fill="#9ebcda", stat="identity")+
  geom_point(data=dat_child,aes(y=1/1-as.integer(change)+1, x=wave_1_titre),pch=124,size=5,alpha=0.5,colour="#8856a7")+
  geom_pointrange(data=child_thresh_50,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=1-median))+
  geom_pointrange(data=child_thresh_80,aes(xmin=`0.025`,xmax=`0.975`,x=`0.5`,y=1-median))+
  scale_x_log10("WT S-specific IgG following Wave 1 (WHO BAU/ml)" )+
  scale_y_continuous("",labels = scales::percent)+
  ggtitle("Children")


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
