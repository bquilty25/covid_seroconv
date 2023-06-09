require("pacman")
#remotes::install_github("njtierney/mmcc")
pacman::p_load(tidyverse,R2jags,mcmcplots,readxl,bayesplot,patchwork,ggExtra,brms,htmlTable,binom,scales,here,mmcc,gridExtra,tableHTML,covidregionaldata,lmerTest,janitor,sjPlot,fastDummies,ggnewscale,parameters)


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

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
  replace_na(list(vaccinated=F)) %>% 
  mutate(age="adult")

child_dat <- read_excel("data/Covid child 4waves data for Billy August2022.xls") %>% 
  rename_all(~stringr::str_replace(.,"^S","")) %>% 
  rename_all(tolower) %>%
  select(c(pid_child:omicron_4w),-contains("barcode"),-contains("cat"),contains("date"),-contains("collection")) %>%
  pivot_longer(cov2s_1w:omicron_4w,values_to = "igg") %>% 
  mutate(variant=case_when(str_detect(name,"cov")~"WT",
                           str_detect(name,"beta")~"Beta",
                           str_detect(name,"delta")~"Delta",
                           str_detect(name,"omicron")~"Omicron"),
         variant=fct_relevel(variant,"WT","Beta","Delta","Omicron"),
         wave=parse_number(str_sub(name,start=-4))) %>% 
  select(-name) %>% 
  left_join(
    read_excel("data/Covid child 4waves data for Billy August2022.xls") %>% 
      rename_all(~stringr::str_replace(.,"^S","")) %>% 
      rename_all(tolower) %>% 
      select(pid_child,contains("collectiondate")) %>% 
      pivot_longer(collectiondate_1w:collectiondate_4w, values_to="collection_date") %>% 
      mutate(wave=parse_number(str_sub(name))) %>% 
      select(-name)
  ) %>% 
  mutate(age="child",
         n_doses=0)

dat <- bind_rows(datw4,child_dat) %>% 
  drop_na(igg)

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

run_model <- function(data.list,n_iter,vacc_diff=F){
 
  if(vacc_diff){
  model <- function(){

    #Likelihood
    for(i in 1:n)
    {

    
      S[i] <- a + (c-a)/(1+exp(beta1*(titre[i]-tm) )) + xeta*vacc[i]

      Y[i] ~ dbern(S[i])

      
    }
    
      # for(k in 1:n_dat){
      # p_over_thresh <- step(titre[k]-tm)
      # 
      # }
    
    # for(j in 1:n_vacc){
    #   p[j]~dbeta(1,1)
    # }
    

    #Priors
    
    #Lower asymptote
    a ~ dbeta(1, 2)
    
    xeta ~dnorm(0,1e-3)
    beta1 ~dnorm(0,1e-3);T(0,)#dgamma(0.01,0.01)
    #beta2 ~dgamma(0.01,0.01)
    
    #Upper asymptote
    c ~ dbeta(2, 1)
    
    #Inflection point
    tm ~ dnorm(1, 1e-4)# log(exp_tm)
    exp_tm <- exp(tm)
    
    diff <- c-a
    
    
  }
  }else{
    
    model <- function(){
      
      #Likelihood
      for(i in 1:n)
      {
        
        
        S[i] <- a + (c-a)/(1+exp(beta1*(titre[i]-tm))) 
        
        Y[i] ~ dbern(S[i])
      
        
      }
      
      
      #Priors
      
      #Lower asymptote
      a ~ dbeta(1, 2)
      
      xeta ~dnorm(0,1e-3)
      beta1 ~ dnorm(0,1e-3);T(0,)#dgamma(0.01,0.01)
      
      #Upper asymptote
      c ~ dbeta(2, 1)
      
      #Inflection point
      tm ~ dnorm(1, 1e-4)# log(exp_tm)
      exp_tm <- exp(tm)
      
      diff <- c-a
      
      
     p_over_thresh_unvacc   <-  sum(step(unvacc_titres-tm))/length(unvacc_titres)
     #p_over_thresh_one_dose <-  sum(step(one_dose_titres-tm))/length(one_dose_titres)
     #p_over_thresh_two_dose <-  sum(step(two_dose_titres-tm))/length(two_dose_titres)
      
      # for(i in 1:length(unvacc_titres)){
      #   unvacc_titres[i]~dnorm(mu_unvacc,0.001)
      # }
      # 
      # mu_unvacc ~ dunif(0,10)
      # 
      # 
      # for(i in 1:length(vacc_titres)){
      #   vacc_titres[i]~dnorm(mu_vacc,0.001)
      # }
      # 
      # mu_vacc ~ dunif(0,10)
      # 
      # p_over_thresh_unvacc <- equals(mu_unvacc,tm)
      # p_over_thresh_vacc <-  equals(mu_vacc,tm)
      
      #p_ratio <- p_over_thresh_vacc/p_over_thresh_unvacc
      
      #p_value <- step(p_ratio-1)
      
    }
  }

  jags.parallel(model.file=model,
       data=data.list,
       parameters.to.save=c("S",
                            "b",
                            "a",
                            "c",
                            "tm",
                            "xeta",
                            "beta0",
                            "beta1",
                            "beta2",
                            "exp_tm",
                            "diff",
                            "over_thresh_vacc",
                            "p_over_thresh_unvacc",
                            "p_over_thresh_one_dose",
                            "p_over_thresh_two_dose",
                            "p_ratio",
                            "p_value",
                            "mu_unvacc",
                            "mu_vacc"
                            ),
       n.iter=n_iter,
       n.chains = 4)
}


or_model <- function(data.list,n_iter){
    
    model <- function(){
      
      #Likelihood
      for(i in 1:n)  
        {
        
        Y[i]~dbern(S[i])
        logit(S[i]) <- inprod(beta[],X[i,])
        
           }
           #Priors
             for (i in 1:ngroups) {
               beta[i] ~ dnorm(0, 1.0E-6) 
               }
           sigma ~ dunif(0, 100)
           tau <- 1 / (sigma * sigma)
      }
      
    jags.parallel(model.file=model,
                  data=data.list,
                  parameters.to.save=c("S",
                                       "beta",
                                       "sigma"),
                  #inits=inits,
                  n.iter=n_iter,
                  n.chains = 4)
    
  }


#### define function to produce all the required results ----

calc_wave <- function(dat, age, wav, preVar, postVar, threshold=.01, sero_pos_pre=FALSE, vacc_agnostic_thresh=TRUE, n_iter=10000,vacc_diff=F, waning=F, diag=F, browsing=F){

  if(browsing){browser()}
    
  wave_dat <- dat %>% 
      filter( age == age,
        (variant == preVar & wave == wav-1) | (variant == postVar & wave == wav)) 
  
  model_dat <- wave_dat %>% 
    select(-c(variant,collection_date, n_doses)) %>% 
    mutate(wave=ifelse(wave==wav-1,"pre","post")) %>%
    pivot_wider(values_from = igg,names_from = wave) %>% 
    drop_na(pre,post) %>% 
    #add in n_doses
    left_join(wave_dat %>% 
                select(-c(variant,collection_date, igg)) %>% 
                mutate(wave=ifelse(wave==wav-1,"doses_pre","doses_post")) %>% 
                pivot_wider(values_from = n_doses,names_from = c(wave))) %>% 
    #add in time between bloods
    left_join(wave_dat %>% 
                select(-c(variant, igg)) %>%  
                mutate(wave=ifelse(wave==wav-1,"pre","post")) %>%
                pivot_wider(values_from = collection_date,names_from=wave) %>% 
                mutate(t_diff=as.numeric(difftime(units = "weeks",post,pre))) %>% 
                select(-pre,-post) %>% 
                mutate(t_diff=ifelse(is.na(t_diff),mean(t_diff,na.rm=T),t_diff))) %>% 
    #assumed decay = 1%  per week (4% per month) (Israel et al. 2022)
    mutate(pre_waning=pre*(1-0.02)^t_diff) %>% 
    mutate(increase_2_v_1=as.integer(post > (pre*1+threshold)),
           increase_2_v_1_waned=as.integer(post > pre_waning),
           increase_2_v_1_10_thresh=as.integer(post > (pre*1+0.1))) %>% 
    pivot_longer(cols=c(increase_2_v_1,increase_2_v_1_waned,increase_2_v_1_10_thresh),
                 names_to = "assump",
                 values_to = "increase")
  
  model_dat %>% 
    tabyl(assump,increase) %>%
    adorn_totals(where = "col") %>% 
    adorn_percentages(denominator = "row") %>%
    adorn_pct_formatting() %>%
    adorn_ns(position = "front")
  
   if(vacc_agnostic_thresh){
    
  model_dat <- model_dat %>% 
    filter(doses_pre==doses_post) %>% 
    mutate(vacc = as.integer(doses_pre!=0)) %>% 
    arrange(vacc)
  
  } else {
    
    n_pre=0
    n_post=0
    
    model_dat <- model_dat %>% 
      filter(doses_pre==n_pre&doses_post==n_post) %>% 
      mutate(vacc = as.integer(doses_pre!=0)) %>% 
      arrange(vacc)
    
  }

  if(waning){
    model_dat <- model_dat %>% 
      filter(assump=="increase_2_v_1_waned")
  } else {
    model_dat <- model_dat %>% 
      filter(assump=="increase_2_v_1")
  }
  
  #if(sero_pos_pre){
    
    min_igg <- 1.09
    
  #} else {
    
    #min_igg <- -Inf
    
  #}  
 
  #run model
  pred_t_wave <- lseq(from=min(model_dat$pre),to=max(model_dat$pre),length.out = 1000)
  
  #if(vacc_diff){
    data.list <-  list(n=nrow(model_dat)+2*length(pred_t_wave),
                       Y=c(as.integer(model_dat$increase),rep(NA,2*length(pred_t_wave))),
                       titre=c(log(model_dat$pre),log(pred_t_wave),log(pred_t_wave)),
                       vacc=1+c(model_dat$vacc,rep(1,length(pred_t_wave)),rep(0,length(pred_t_wave))),
                       unvacc_titres=log(model_dat %>% filter(pre>min_igg,n_doses==0) %>% pull(pre)),
                       one_dose_titres=log(model_dat %>% filter(pre>min_igg,n_doses==1) %>% pull(pre)),
                       two_dose_titres=log(model_dat %>% filter(pre>min_igg,n_doses==2) %>% pull(pre))
    )
  # } else{
  #   data.list <-  list(n=nrow(model_dat)+length(pred_t_wave),
  #                      Y=c(as.integer(model_dat$increase),rep(NA,length(pred_t_wave))),
  #                      titre=c(log(model_dat$pre),log(pred_t_wave)),
  #                      vacc=NA,
  #                      unvacc_titres=NA,
  #                      one_dose_titres=NA,
  #                      two_dose_titres=NA
  #   )
  # }
  
  wave_res <- run_model(data.list = data.list,
                        n_iter=n_iter,
                        vacc_diff = vacc_diff)
  
  # #alt model
  # or_model_dat <- model_dat %>% 
  #   mutate(doses_pre=factor(doses_pre,levels=c(0,1,2)),
  #                                   log_pre=log(pre)) 
  # 
  # or_res <- or_model(data.list = list(n=nrow(or_model_dat),
  #                                        X=model.matrix(~doses_pre,or_model_dat),
  #                                        ngroups=ncol(model.matrix(~doses_pre,or_model_dat)),
  #                                        Y=c(as.integer(or_model_dat$increase_2_v_1)),
  #                                        titre=c(or_model_dat$log_pre)),
  #                       n_iter=n_iter*10)

  # or_res <- mmcc::tidy(as.mcmc(or_res)) %>% 
  #   filter(str_detect(parameter,"beta")) %>% 
  #   mutate_if(is.numeric,exp) %>% 
  #   mutate(across(c(median, `2.5%`, `97.5%`),  ~ formatC(., digits = 2, format = "f")),
  #          `97.5%` = case_when(as.numeric(`97.5%`) > max(model_dat$pre) ~ paste0(">",plyr::round_any(max(model_dat$pre),500)), 
  #                              TRUE ~ `97.5%`)) %>%
  #   mutate(estimate = paste0(median, " (", `2.5%`, ", ", `97.5%`,")")
  #   ) %>% 
  #   select(parameter, estimate) %>%
  #   pivot_wider(values_from = c(estimate), names_from = parameter) 
  
  #diagnostics
  
  if(diag){mcmcplot(wave_res,random = T)}
  
  # if(vacc_diff==F){
  # bayes_p <- wave_res %>% 
  #   as.mcmc() %>% 
  #   mcmc_to_dt() %>% 
  #   filter(str_detect(parameter,"p_over_thresh")) %>% 
  #   pivot_wider(values_from = value,names_from = parameter) %>% 
  #   mutate(ratio_one_dose=p_over_thresh_one_dose/p_over_thresh_unvacc,
  #          ratio_two_dose=p_over_thresh_two_dose/p_over_thresh_unvacc) 
  # 
  # bayes_p %>% 
  #   pull(ratio_one_dose) %>% 
  #   bayestestR::describe_posterior(.)
  # 
  # bayes_p %>% 
  #   pull(ratio_two_dose) %>% 
  #   bayestestR::describe_posterior(.)
  # }
  
  wave_res %>% 
    as.mcmc() %>% 
    tidy() %>% 
    filter(str_detect(parameter,"p_over_thresh")) 
  
  proportion_protected <- model_dat %>% 
    filter(pre>min_igg) %>% 
    bind_cols(mmcc::tidy(as.mcmc(wave_res)) %>% 
                            filter(str_detect(parameter,"exp_tm"))) %>% 
    mutate(across(c(`2.5%`:`97.5%`),.fns=list(protected=function(x){pre>x}),.names="{fn}_{col}")) %>%
    unite("doses", doses_pre,doses_post) %>% 
    pivot_longer(c(`protected_2.5%`,protected_median,`protected_97.5%`)) %>% 
    tabyl(doses,value,name,increase) %>% 
    adorn_totals("row") %>% 
    adorn_percentages("row")  %>%  
    adorn_pct_formatting(digits = 1)%>% 
    adorn_ns() %>% 
    bind_rows(.id="group") %>%
    select(-`FALSE`) %>% 
    pivot_wider(names_from = group,values_from = `TRUE`) %>% 
    select(doses,protected_median,`protected_97.5%`,everything()) %>% 
    bind_cols(model_dat %>% 
                group_by(doses_pre,doses_post) %>% 
                count() %>% 
                adorn_totals()) %>% 
    select(-c(doses_pre,doses_post)) 
  
  #CoP estimate 
  (wave_plot <- mmcc::tidy(as.mcmc(wave_res)) %>% 
      filter(str_detect(parameter,"S")) %>% 
      select(median,`2.5%`,`97.5%`) %>% 
      slice(nrow(model_dat)+1:n()) %>% 
      bind_cols(data.frame(titre=c(pred_t_wave,pred_t_wave),
                           vacc=as.factor(c(rep(1,length(pred_t_wave)),rep(0,length(pred_t_wave)))))) %>% 
      ggplot()+
      geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=titre, group=vacc,colour= vacc,fill=vacc),stat = "identity")+
      scale_x_log10(paste(preVar,"S-specific IgG pre-wave (WHO BAU/ml)"))+
      coord_cartesian(xlim=c(min(model_dat$pre),max(model_dat$pre)),expand = TRUE)+
      scale_y_continuous(paste("Probability of increased",postVar,"S-specific\nIgG titres following wave",wav),labels = scales::percent,limits=c(0,1))+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA),axis.ticks = element_line())#+
      #ggtitle(paste0("Wave ",wav,", vaccine agnostic threshold: ",vacc_agnostic_thresh, ", only seropositives: ", sero_pos_pre))
  )
  ggsave(paste0("results/waveplot","age",age,"wave",wav,preVar,postVar,"vacc_ag_thresh",vacc_agnostic_thresh,"seropositivesonly",sero_pos_pre,".png"),wave_plot,width=150,height=100,units="mm",dpi=600,bg="white")
  
  (wave_change_plot <- model_dat %>% 
      pivot_longer(c(pre,post)) %>% 
      mutate(variant=ifelse(name=="pre",paste0("pre (",preVar,")"),paste0("post (",postVar,")"))) %>% 
      mutate(variant=factor(variant, levels = c(paste0("pre (",preVar,")"),paste0("post (",postVar,")")))) %>%
      ggplot()+
      geom_point(aes(x=variant,y=value,group=pid_child,colour=factor(increase)),alpha=0.2)+
      geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase)),alpha=0.2)+
      scale_y_log10("S-specific IgG titre (WHO BAU/ml)",limits=c(NA,10000))+
      xlab("")+
      theme_minimal()+
      theme(panel.border = element_rect(fill = NA))+
      #scale_colour_manual("Increase\npost-wave",values=c("grey","#e41a1c"),labels=c("No","Yes"))+
      ggtitle(paste0("Wave ",wav,", vaccine agnostic threshold: ",vacc_agnostic_thresh, ", only seropositives: ", sero_pos_pre))+
    geom_pointrange(data=mmcc::tidy(as.mcmc(wave_res)) %>% 
                      filter(str_detect(parameter,"exp_tm")),aes(y=median,ymin=`2.5%`, ymax=`97.5%`,x=0.9))+
      geom_text(data=mmcc::tidy(as.mcmc(wave_res)) %>% 
                  filter(str_detect(parameter,"exp_tm")),aes(y=median,x=0.9, label="50%\nthreshold"),hjust=1))
  
  ggsave(paste0("results/changeplot","age",age,"wave",wav,preVar,postVar,"vacc_ag_thresh",vacc_agnostic_thresh,"seropositivesonly",sero_pos_pre,".png"),wave_change_plot,width=150,height=100,units="mm",dpi=600,bg="white")
  
  #Tabled results
  res_est <- mmcc::tidy(as.mcmc(wave_res)) %>%
    filter(parameter %in% c("a", "c", "exp_tm","diff")) %>%
    select(-mean,-sd) %>% 
    pivot_longer(c(median, `2.5%`, `97.5%`)) %>% 
    mutate(value=ifelse(parameter=="a"|parameter=="c"|parameter=="diff",value*100,value))
  
  res <- res_est %>% 
    pivot_wider(names_from = name,values_from = value) %>% 
    mutate(across(c(median, `2.5%`, `97.5%`),  ~ formatC(., digits = 1, format = "f"))) %>%
    mutate(estimate = paste0(median, " (", `2.5%`, ", ", `97.5%`,")")
    ) %>% 
    select(parameter, estimate) %>%
    pivot_wider(values_from = c(estimate), names_from = parameter) %>% 
    mutate(
      age=age,
      Wave = wav,
      pre = preVar,
      post = postVar,
      vacc_agnostic_thresh = vacc_agnostic_thresh,
      sero_pos_pre=sero_pos_pre,
      .before = "a"
    ) 
  
  dat_size <- model_dat %>% summarise(n=n(),n_increased=sum(increase)) 
  
  res <- res %>% bind_cols(dat_size) 
  
  res
  
  #goodness of fit
  (wave_gof <- remove_geom(wave_plot,"GeomPointrange")+
      geom_pointrange(data=model_dat %>% 
                        mutate(titre_group=cut(pre,breaks=c(-Inf,1.09,5,10,50,100,500,1000,5000,Inf))) %>% 
                        group_by(titre_group) %>% 
                        summarise(N=n(),
                                  avg_titre=median(pre),
                                  x=sum(increase),
                                  p=x/N,
                                  p.lo=binom::binom.confint(x,N,methods = "exact")$lower,
                                  p.hi=binom::binom.confint(x,N,methods = "exact")$upper),
                      aes(x=avg_titre,y=p,ymin=p.lo,ymax=p.hi))
  )
  ggsave(paste0("results/wave_gof","age",age,"wave",wav,preVar,postVar,"vacc_ag_thresh",vacc_agnostic_thresh,"seropositivesonly",sero_pos_pre,".png"),wave_gof,width=150,height=100,units="mm",dpi=600,bg="white")
  
  return(
    list(
      res = res,
      proportion_protected = proportion_protected,
      wave_change_plot = wave_change_plot,
      wave_plot = wave_plot,
      #or_res = or_res,
      wave_res = wave_res
    )
  )
}

`%!in%` <- Negate(`%in%`)

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
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


