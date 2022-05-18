require("pacman")
#remotes::install_github("njtierney/mmcc")
pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra,brms,htmlTable,binom,scales,here,mmcc,gridExtra,tableHTML,covidregionaldata)


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

#Load and clean data
# dat <- read_xlsx(here("data","SARS_CoV2_data_3_waves_15dec2021.xlsx"),sheet = 2) %>% 
#   pivot_longer(CoV2S_1w_m:delta_3w_m) %>% 
#   mutate(variant=case_when(str_detect(name,"CoV2")~"WT",
#                            str_detect(name,"beta")~"Beta",
#                            str_detect(name,"delta")~"Delta"),
#          variant=fct_relevel(variant,"WT","Beta","Delta"),
#          wave=parse_number(str_sub(name,start=-4))) %>% 
#   select(-name) 

#Load and clean wave 1-4 data
datw4 <- read_excel("data/Covid data for Billy May 2022 (all data).xlsx") %>% 
  rename_all(~stringr::str_replace(.,"^S","")) %>% 
  rename_all(tolower) %>%
  rename("mat_vacc_cov_date_2" =mat_cov_vacc_date_2 ) %>% 
  select(c(pid_child:omicron_4w_m),-contains("barcode"),-contains("cat"),contains("date"),-contains("collection")) %>%
  pivot_longer(cov2s_1w_m:omicron_4w_m,values_to = "igg") %>% 
  mutate(variant=case_when(str_detect(name,"cov")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta",
                           str_detect(name,"omicron")~"Omicron"),
         variant=fct_relevel(variant,"WT","Beta","Delta","Omicron"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name) %>% 
  left_join(
    read_excel("data/Covid data for Billy May 2022 (all data).xlsx") %>% 
      rename_all(~stringr::str_replace(.,"^S","")) %>% 
      rename_all(tolower) %>% 
      select(pid_child,contains("collectiondate")) %>% 
      pivot_longer(collectiondate_1w_m:collectiondate_4w_m, values_to="collection_date") %>% 
      mutate(wave=parse_number(str_sub(name))) %>% 
      select(-name)
  ) %>% 
  mutate(n_doses=case_when(collection_date>mat_vacc_cov_date_1+3&collection_date<mat_vacc_cov_date_2+3|
                             collection_date>mat_vacc_cov_date_1+3&is.na(mat_vacc_cov_date_2+3)~1,
                           collection_date>mat_vacc_cov_date_1+3&collection_date>mat_vacc_cov_date_2+3~2,
                           TRUE~0)) %>% 
  select(-c(mat_vacc_cov_date_1,mat_vacc_cov_date_2)) %>% 
  replace_na(list(vaccinated=F))

#take posterior samples of parameters to estimate values of titre at probability thresholds
extract_ab_thresholds <- function(jags_res,thresh="thresh_50",mult=1,...){
  #browser()
  y <- mmcc::tidy(as.mcmc(jags_res)) %>% 
    filter(parameter==!!thresh)
  
  x <- mmcc::mcmc_to_dt(as.mcmc(jags_res)) %>% 
    filter(parameter%in%c("a","b","c","tm")) %>% 
    pivot_wider(values_from = "value",names_from = "parameter") %>% 
    mutate(x_pred=exp(-(log((c-y$median*mult)/(y$median*mult-a))/b)+tm)) %>% 
    summarise(x = quantile(x_pred, c(0.025, 0.5, 0.975),na.rm=T), q = c(0.025, 0.5, 0.975)) 
  
  y %>% bind_cols(x)
}



run_model <- function(data.list,n_iter){
  
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
       n.iter=n_iter,
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



#### define function to produce all the required results ----

calc_wave <- function(dat, wav, preVar, postVar, threshold=.01, doses_pre=0, doses_post=0, sero_pos_pre=FALSE, n_iter=5000,diag=F, browsing=F){

  if(browsing){browser()}
  
  if(sero_pos_pre){
    min_igg <- 1.09
  } else {
    min_igg <- -Inf
  }
  
  wave_dat <- dat %>% 
    filter( (variant == preVar & wave == wav-1 & n_doses == doses_pre & igg > min_igg) | (variant == postVar & wave == wav & n_doses == doses_post)) %>% 
    select(-c(variant,collection_date, n_doses)) %>% 
    mutate(wave=ifelse(wave==wav-1,"pre","post")) %>%
    pivot_wider(values_from = igg,names_from = wave) %>% 
    mutate(increase_2_v_1=as.integer(post > (pre*1+threshold))) %>% 
    drop_na(pre,post)
  
  #run model
  pred_t_wave <- lseq(from=min(wave_dat$pre),to=max(wave_dat$pre),length.out = 1000)
  wave_res <- run_model(data.list = list(n=nrow(wave_dat)+length(pred_t_wave),
                                         Y=c(as.integer(wave_dat$increase_2_v_1),rep(NA,length(pred_t_wave))),
                                         titre=c(log(wave_dat$pre),log(pred_t_wave))),
                        n_iter=n_iter)
  
  #diagnostics
  if(diag){mcmcplot(wave_res,random = T)}
  
  #CoP estimate 
  (wave_plot <- mmcc::tidy(as.mcmc(wave_res)) %>% 
      filter(str_detect(parameter,"S")) %>% 
      select(median,`2.5%`,`97.5%`) %>% 
      slice(nrow(wave_dat)+1:n()) %>% 
      bind_cols(pred_t_wave) %>% 
      ggplot()+
      geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
      geom_point(data=wave_dat,aes(y=as.integer(increase_2_v_1), x=pre),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
      geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_50") %>% pivot_wider(values_from = x,names_from = q) %>% 
                        select(-c(mean,sd)),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
      geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_80") %>% pivot_wider(values_from = x,names_from = q) %>% 
                        select(-c(mean,sd)),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
      geom_hline(data=mmcc::tidy(as.mcmc(wave_res)) %>%
                   filter(parameter%in%c("a","c")),
                 aes(yintercept=median),
                 linetype="dashed")+
      scale_x_log10(paste(preVar,"S-specific IgG pre-wave (WHO BAU/ml)"))+
      coord_cartesian(xlim=c(min(wave_dat$pre),max(wave_dat$pre)),expand = expansion(mult=0.01))+
      scale_y_continuous(paste("Probability of increased",postVar,"S-specific\nIgG titres following wave",wav),labels = scales::percent)+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA),axis.ticks = element_line())+ggtitle(paste("Wave",wav))
  )
  ggsave(paste0("results/waveplot",wav,preVar,postVar,".png"),wave_plot,width=120,height=100,units="mm",dpi=600,bg="white")
  
  wave_change_plot <- wave_dat %>% 
      pivot_longer(c(pre,post)) %>% 
      mutate(variant=ifelse(name=="pre",paste0("pre (",preVar,")"),paste0("post (",postVar,")"))) %>% 
      mutate(variant=factor(variant, levels = c(paste0("pre (",preVar,")"),paste0("post (",postVar,")")))) %>%
      ggplot()+
      geom_point(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)),alpha=0.2)+
      geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)),alpha=0.2)+
      scale_y_log10("S-specific IgG titre (WHO BAU/ml)",limits=c(NA,10000))+
      xlab("")+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA))+
      scale_colour_manual("Increase\npost-wave",values=c("grey","#e41a1c"),labels=c("No","Yes"))+
      ggtitle(paste0("Wave ",wav,", Doses pre: ",doses_pre,", Doses post: ", doses_post, ", Only seropositives: ", sero_pos_pre))+
    geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_50") %>% pivot_wider(values_from = x,names_from = q) %>% 
                      select(-c(mean,sd)),aes(y=`0.5`,ymin=`0.025`, ymax=`0.975`,x=0.9))+
    geom_pointrange(data=extract_ab_thresholds(wave_res,"thresh_80") %>% pivot_wider(values_from = x,names_from = q) %>% 
                      select(-c(mean,sd)),aes(y=`0.5`,ymin=`0.025`, ymax=`0.975`,x=1.1))+
      geom_text(data=extract_ab_thresholds(wave_res,"thresh_50") %>% pivot_wider(values_from = x,names_from = q) %>% 
                 select(-c(mean,sd)),aes(y=`0.5`,x=0.9, label="50% threshold"),hjust=1)+
    geom_text(data=extract_ab_thresholds(wave_res,"thresh_50") %>% pivot_wider(values_from = x,names_from = q) %>% 
               select(-c(mean,sd)),aes(y=`0.5`,x=1.1, label="80% threshold"),hjust=0)
  
  ggsave(paste0("results/changeplot",wav,preVar,postVar,"_Doses pre",doses_pre,"_Doses post", doses_post, "_Only seropositives",sero_pos_pre,".png"),wave_change_plot,width=150,height=100,units="mm",dpi=600,bg="white")
  
  #Tabled results
  res_est <- mmcc::tidy(as.mcmc(wave_res)) %>%
    filter(parameter %in% c("a", "c", "exp_tm")) %>%
    select(-mean,-sd) %>% 
    pivot_longer(c(median, `2.5%`, `97.5%`)) %>% 
    mutate(value=ifelse(parameter=="a"|parameter=="c",value*100,value))

  res_trans <- res_est%>%
    pivot_wider(names_from = name,values_from = value) %>% 
    bind_rows(
      extract_ab_thresholds(wave_res, thresh = "a", mult = 0.5) %>%
        mutate(parameter = "a_0.5") %>% 
        rowwise() %>% 
        mutate(
          protected=wave_dat %>% 
            filter(pre>1.09) %>% 
            summarise(n=n(),
                      protected_p=sum(pre>x)/n) %>% 
            pull(protected_p)) %>% 
        pivot_wider(values_from = c(x,protected),names_from = q) %>% 
        select(-c(mean,sd))
    ) %>%
    bind_rows(
      extract_ab_thresholds(wave_res, thresh = "thresh_50") %>%
        rowwise() %>% 
        mutate(
          protected=wave_dat %>% 
            filter(pre>1.09) %>% 
            summarise(n=n(),
                      protected_p=sum(pre>x)/n) %>% 
            pull(protected_p)) %>%  
        pivot_wider(values_from = c(x,protected),names_from = q) %>% 
        select(-c(mean,sd))
    ) %>%
    bind_rows(
      extract_ab_thresholds(wave_res, thresh = "thresh_80") %>% 
        rowwise() %>% 
        mutate(
          protected=wave_dat %>% 
            filter(pre>1.09) %>% 
            summarise(n=n(),
                      protected_p=sum(pre>x)/n) %>% 
            pull(protected_p)) %>%  
        pivot_wider(values_from = c(x,protected),names_from = q) %>% 
        select(-c(mean,sd))
    ) %>%
    filter(parameter%!in%c("a", "c", "exp_tm")) %>% 
    select(-c(median, `2.5%`, `97.5%`)) %>% 
    mutate(across(c(x_0.5, x_0.025, x_0.975),  ~ formatC(., digits = 1, format = "f")),
           across(c(protected_0.5, protected_0.025, protected_0.975),~ formatC(.*100, digits = 1, format = "f")),
           x_0.975 = case_when(as.numeric(x_0.975) > max(wave_dat$pre) ~ paste0(">",plyr::round_any(max(wave_dat$pre),500)), 
                               TRUE ~ x_0.975)) %>%
    mutate(estimate = paste0(x_0.5, " (", x_0.975, ", ", x_0.975,")"),
           prop_protected = paste0(protected_0.5, " (", protected_0.025, ", ", protected_0.975,")")
    ) %>% 
    select(parameter,prop_protected, estimate) %>%
    pivot_wider(values_from = c(estimate,prop_protected), names_from = parameter)
  
  res <- res_est %>% 
    pivot_wider(names_from = name,values_from = value) %>% 
    mutate(across(c(median, `2.5%`, `97.5%`),  ~ formatC(., digits = 1, format = "f")),
           `97.5%` = case_when(as.numeric(`97.5%`) > max(wave_dat$pre) ~ paste0(">",plyr::round_any(max(wave_dat$pre),500)), 
                               TRUE ~ `97.5%`)) %>%
    mutate(estimate = paste0(median, " (", `2.5%`, ", ", `97.5%`,")")
    ) %>% 
    select(parameter, estimate) %>%
    pivot_wider(values_from = c(estimate), names_from = parameter) %>% 
    bind_cols(res_trans) %>%
    mutate(
      Wave = wav,
      pre = preVar,
      post = postVar,
      doses_pre=doses_pre,
      doses_post=doses_post,
      only_seropos=only_seropos,
      .before = "a"
    ) 
  
  dat_size<-wave_dat %>% summarise(n=n(),n_increased=sum(increase_2_v_1==T))
  
  res <- res %>% bind_cols(dat_size)
  
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
  
  return(list(res=res,wave_change_plot=wave_change_plot,wave_plot=wave_plot))
}

`%!in%` <- Negate(`%in%`)


