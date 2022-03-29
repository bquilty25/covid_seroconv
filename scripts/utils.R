require("pacman")
#remotes::install_github("njtierney/mmcc")
pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra,brms,htmlTable,binom,scales,here,mmcc)


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Load and clean data
dat <- read_xlsx(here("data","SARS_CoV2_data_3_waves_15dec2021.xlsx"),sheet = 2) %>% 
  pivot_longer(CoV2S_1w_m:delta_3w_m) %>% 
  mutate(variant=case_when(str_detect(name,"CoV2")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta"),
         variant=fct_relevel(variant,"WT","Beta","Delta"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name) 

#Load and clean wave 1-4 data
datw4 <- read_xlsx(here("data","data_for_billy_all_4waves_339_29MAR2022_with_vaccine.xlsx"),sheet = 1) %>% 
  filter(is.na(mat_vacc_covid)) %>% #exclude those who ever got vaccinated
  select(pid_child:omicron_4w_m) %>%
  pivot_longer(CoV2S_1w_m:omicron_4w_m) %>% 
  mutate(variant=case_when(str_detect(name,"CoV")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta",
                           str_detect(name,"omicro")~"Omicron"),
         variant=fct_relevel(variant,"WT","Beta","Delta","Omicron"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name) 

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
  
  jags.parallel(model.file=model,
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
       n.iter=2000,
       n.chains = 4)
  
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


###########################################
# define function to produce all the required results
###########################################
calc_wave <- function(dat, wav, preVar, postVar, threshold=.10){
  
  wave_dat <- dat %>% 
    filter( (variant == preVar & wave == wav-1) | (variant == postVar & wave == wav)) %>% 
    select(-variant) %>% 
    mutate(wave=ifelse(wave==wav-1,"pre","post")) %>%
    pivot_wider(values_from = value,names_from = wave) %>% 
    mutate(increase_2_v_1=as.integer(post > (pre*1+threshold))) %>% #serconversion defined as 10% increase
    drop_na(pre,post)
  
  #plot individual level change in titre
  (wave_change_plot <- wave_dat %>% 
      pivot_longer(c(pre,post)) %>% 
      mutate(variant=ifelse(name=="pre",paste0("pre (",preVar,")"),paste0("post (",postVar,")"))) %>% 
      mutate(variant=factor(variant, levels = c(paste0("pre (",preVar,")"),paste0("post (",postVar,")")))) %>%
      ggplot(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)))+
      geom_point(alpha=0.2)+
      geom_path(alpha=0.2)+
      scale_y_log10("S-specific IgG titre (WHO BAU/ml)",limits=c(NA,10000))+
      xlab("")+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA))+
      scale_colour_brewer("Increase\npost-wave",type="qual",palette = "Set1",labels=c("No","Yes"),direction=-1)+
      ggtitle(paste("Wave",wav))
  )
  ggsave(paste0("results/changeplot",wav,preVar,postVar,".png"),wave_change_plot,width=150,height=100,units="mm",dpi=600,bg="white")
  
  #run model
  pred_t_wave <- lseq(from=min(wave_dat$pre),to=max(wave_dat$pre),length.out = 1000)
  wave_res <- run_model(data.list = list(n=nrow(wave_dat)+length(pred_t_wave),
                                         Y=c(as.integer(wave_dat$increase_2_v_1),rep(NA,length(pred_t_wave))),
                                         titre=c(log(wave_dat$pre),log(pred_t_wave))))
  
  #diagnostics
  mcmcplot(wave_res,random = T)
  
  #CoP estimate 
  (wave_plot <- mmcc::tidy(as.mcmc(wave_res)) %>% 
      filter(str_detect(parameter,"S")) %>% 
      select(median,`2.5%`,`97.5%`) %>% 
      slice(nrow(wave_dat)+1:n()) %>% 
      bind_cols(pred_t_wave) %>% 
      ggplot()+
      geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
      geom_point(data=wave_dat,aes(y=as.integer(increase_2_v_1), x=pre),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
      geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_50"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
      geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_80"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
      geom_hline(data=mmcc::tidy(as.mcmc(wave_res)) %>%
                   filter(parameter%in%c("a","c")),
                 aes(yintercept=median),
                 linetype="dashed")+
      scale_x_log10(paste(preVar,"S-specific IgG pre-wave (WHO BAU/ml)"))+
      scale_y_continuous(paste("Probability of increased",postVar,"S-specific\nIgG titres following wave",wav),labels = scales::percent)+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA))+ggtitle(paste("Wave",wav))
  )
  ggsave(paste0("results/waveplot",wav,preVar,postVar,".png"),wave_plot,width=120,height=100,units="mm",dpi=600,bg="white")
  
  #Tabled results
  res <- mmcc::tidy(as.mcmc(wave_res)) %>% 
    filter(parameter %in% c("a","c")) %>% 
    mutate(across(c(median,`2.5%`,`97.5%`),function(x)x*100)) %>% 
    bind_rows(extract_ab_thresholds(wave_res,thresh="a",mult=0.5) %>% 
                mutate(parameter="a_0.5") %>% select(parameter,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
    bind_rows(extract_ab_thresholds(wave_res,thresh="thresh_50") %>% 
                select(parameter,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
    bind_rows(extract_ab_thresholds(wave_res,thresh="thresh_80") %>% 
                select(parameter,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
    mutate(`97.5%`=ifelse(`97.5%`>1000,Inf,`97.5%`),
           across(c(median,`2.5%`,`97.5%`),~formatC(.,digits=1,format = "f"))) %>% 
    mutate(estimate=paste0(median, " (",`2.5%`,", ",`97.5%`,")")) %>% 
    select(parameter,estimate) %>% 
    pivot_wider(values_from = estimate,names_from = parameter) %>% 
    select(a,c,thresh_50,thresh_80,a_0.5) %>% 
    mutate(Wave = wav,
           pre = preVar,
           post = postVar, .before="a") 
  
  #goodness of fit
  (wave_gof <- remove_geom(wave_plot,"GeomPointrange")+
      geom_pointrange(data=wave_dat %>% 
                        mutate(titre_group=cut(pre,breaks=c(-Inf,1.09,5,10,50,100,500,Inf))) %>% 
                        group_by(titre_group) %>% 
                        summarise(N=n(),
                                  avg_titre=median(pre),
                                  x=sum(increase_2_v_1),
                                  p=x/N,
                                  p.lo=binom::binom.confint(x,N,methods = "exact")$lower,
                                  p.hi=binom::binom.confint(x,N,methods = "exact")$upper),
                      aes(x=avg_titre,y=p,ymin=p.lo,ymax=p.hi))
  )
  ggsave(paste0("results/wave_gof",wav,preVar,postVar,".png"),wave_gof,width=120,height=100,units="mm",dpi=600,bg="white")
  
  return(res)
}



