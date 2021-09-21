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
    
    S[i] <- a + c * exp(-exp(-b*(t[i])))
    
    Y[i]~dbern(S[i])
    
  }
  
  max_Y <- max(S)
  reduction <- 1-(a/max_Y)
  
  thresh_50 <- -a*0.5+max_Y*0.5+a
  thresh_95 <- -a*0.05+max_Y*0.05+a
  
  #Priors
  a~dbeta(1, 1) #
  
  c~dbeta(1,1)
  
  b~dnorm(0, lambda1)
  lambda1~dgamma(0.0001,0.0001) #uninformative Gamma prior
  
  # tm~dnorm(0, lambda2)
  # exp_tm <- exp(tm)
  # lambda2~dgamma(0.0001,0.0001) #uninformative gamma prior
  
}

jags(model.file=model,
                data=data.list,
                parameters.to.save=c("S","b","a","c","reduction","lambda1","thresh_50","thresh_95","min_Y","max_Y"),
                n.iter=2.5e5,
                n.chains = 2,
                n.burnin = 2.5e4)

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

mmcc::tidy(as.mcmc(child_res)) %>% 
  mutate(name="child") %>% 
  bind_rows(mmcc::tidy(as.mcmc(mother_res)) %>% 
              mutate(name="mother")) %>% 
  filter(!str_detect(parameter,"S")) %>% 
  arrange(parameter)

#mother threshold
mother_tidy <- mmcc::tidy(as.mcmc(mother_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(dat_mother)+1:n())

mother_thresh_50_y <- mmcc::tidy(as.mcmc(mother_res)) %>% 
  filter(parameter=="thresh_50") %>%  
  pull(median)

mother_thresh_50_med <- exp(approx(y=log(pred_t),x=mother_tidy$median,
       xout = mother_thresh_50_y)$y)

mother_thresh_50_2.5 <- exp(approx(y=log(pred_t),x=mother_tidy$`2.5%`,
       xout = mother_thresh_50_y)$y)

mother_thresh_50_97.5 <- exp(approx(y=log(pred_t),x=mother_tidy$`97.5%`,
       xout = mother_thresh_50_y)$y)

paste0(mother_thresh_50_med,"(",mother_thresh_50_2.5," - ",mother_thresh_50_97.5,")")

mother_plot <-mother_tidy%>% 
  select(median,`2.5%`,`97.5%`) %>% 
  bind_cols(pred_t) %>% 
  ggplot(aes(y=median, x=...4)) + 
  geom_point(data=dat_mother,aes(y=as.integer(change)-1, x=wave_1_titre),pch=124,size=4,alpha=0.5,colour="#2ca25f")+
  geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_pointrange(aes(xmin=mother_thresh_50_2.5,
                     xmax=mother_thresh_50_97.5,
                     x=mother_thresh_50_med,y=mother_thresh_50_y))+
  scale_x_log10("Anti-Spike IgG (Wave 1)")+
  scale_y_continuous("Probability of increased anti-Spike IgG titres (Wave 2)")+
  ggtitle("Mothers")

#child threshold
child_tidy <- mmcc::tidy(as.mcmc(child_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(dat_child)+1:n())

child_thresh_50_y <- mmcc::tidy(as.mcmc(child_res)) %>% 
  filter(parameter=="thresh_50") %>%  
  pull(median)

child_thresh_50_med <- exp(approx(y=log(pred_t),x=child_tidy$median,
                                    xout = child_thresh_50_y)$y)

child_thresh_50_2.5 <- exp(approx(y=log(pred_t),x=child_tidy$`2.5%`,
                                    xout = child_thresh_50_y)$y)

child_thresh_50_97.5 <- exp(approx(y=log(pred_t),x=child_tidy$`97.5%`,
                                     xout = child_thresh_50_y)$y)

child_plot <- child_tidy %>% 
  bind_cols(pred_t) %>% 
  ggplot(aes(y=(median), x=...4)) + 
  geom_point(data=dat_child,aes(y=as.integer(change)-1, x=wave_1_titre),pch=124,size=4,alpha=0.5,colour="#8856a7")+
  geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`),colour="#8856a7", fill="#9ebcda", stat="identity")+
  geom_pointrange(aes(xmin=child_thresh_50_2.5,
                      xmax=child_thresh_50_97.5,
                      x=child_thresh_50_med,y=child_thresh_50_y))+
  scale_x_log10("Anti-Spike IgG (Wave 1)")+
  scale_y_continuous("")+
  ggtitle("Children")

mother_plot+child_plot&theme_minimal()&theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(fill = NA))&coord_fixed(ratio=4/1)&scale_fill_brewer()

ggsave("re.png",width=210,height=150,units="mm",dpi=600)
