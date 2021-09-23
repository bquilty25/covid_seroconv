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

dat_child <- dat %>% select(1:3) %>% rename_all(~stringr::str_replace(.,"^child_","")) %>%  mutate(age="child") %>% 
  drop_na(change) %>% mutate(change=as.factor(change))  %>% filter(wave_1_titre>1)
dat_mother <- dat %>% select(4:6) %>% rename_all(~stringr::str_replace(.,"^mother_","")) %>% mutate(age="mother")%>% 
  drop_na(change) %>% mutate(change=as.factor(change)) %>% filter(wave_1_titre>1)


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

mother_tab <- dat_mother %>% 
  mutate(sero1=wave_1_titre>1.09,sero2=wave_2_titre>1.09) 

mother_tab <- table(mother_tab$sero1,mother_tab$change)
1-epitools::oddsratio.wald(mother_tab)$measure

child_tab <- dat_child %>% 
  mutate(sero1=wave_1_titre>1.09,sero2=wave_2_titre>1.09) 

child_tab <- table(child_tab$sero1,child_tab$change)
1-epitools::oddsratio.wald(child_tab)$measure

#rev cum dens
child_rcd <- Hmisc::Ecdf(dat_child$wave_1_titre,what="1-F")

#titre at attack rate efficacy
data.frame(x=child_rcd$x,y=child_rcd$y) %>% 
  filter(y>0.64,y<0.65)

#rev cum dens
mother_rcd <- Hmisc::Ecdf(dat_mother$wave_1_titre,what="1-F")

#titre at attack rate efficacy
data.frame(x=mother_rcd$x,y=mother_rcd$y) %>% 
  filter(y>0.62,y<0.63)


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
    
    logit(S[i]) <- a + b*t[i]
    
    Y[i]~dbern(S[i])
    
  }
  
  ld50 <- -a/b
  
  #Priors
  a~dnorm(0, lambda1)
  lambda1 ~ dgamma(0.0001,0.0001) #uninformative Gamma prior
  
  b~dnorm(0, lambda2)
  lambda2~dgamma(0.0001,0.0001) #uninformative Gamma prior
  
}

inits <-  function(){list(
  a=runif(1,-5,5),
  b=runif(1,-5,5))}

jags(model.file=model,
                data=data.list,
                parameters.to.save=c("S","b","a","ld50"),
                n.iter=2.5e4,
                n.chains = 2,
                n.burnin = 2.5e3,
                inits = inits)

}

pred_t <- lseq(from=0.55,to=1000,length.out = 1000)

#Run model for mother

mother_res <- run_model(data.list = list(n=nrow(dat_mother)+length(pred_t),
                                         Y=c(as.integer(dat_mother$change)-1,rep(NA,length(pred_t))),
                                         t=c(log(dat_mother$wave_1_titre),log(pred_t))))

saveRDS(mother_res,"mother_mcmc_exc_seroneg.rds")

#Run model for children
child_res <-  run_model(data.list = list(n=nrow(dat_child)+length(pred_t),
                                         Y=c(as.integer(dat_child$change)-1,rep(NA,length(pred_t))),
                                         t=c(log(dat_child$wave_1_titre),log(pred_t))))
saveRDS(child_res,"child_mcmc_exc_seroneg.rds")
