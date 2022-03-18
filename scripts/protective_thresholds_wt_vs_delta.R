source("utils.R")

dat <- dat %>% 
  filter(variant == "WT" & wave == 2 | variant=="Delta" & wave==3) %>% 
  select(-variant) %>% 
  pivot_wider(values_from = value,names_from = wave) %>% 
  mutate(increase_3_v_2=as.integer(`3`>`2`*1.1)) #increase if >10% increase in titer


wave3_wt_dat <- dat %>% 
  drop_na(`2`,`3`)

(wave3_wt_change_plot <- wave3_wt_dat %>% 
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

#Vector of Ag values for prediction (appended to data list)
pred_t_wave3_wt <- lseq(
  from = min(wave3_wt_dat$`2`),
  to = max(wave3_wt_dat$`2`),
  length.out = 1000
)

# Data for model 
model_dat <- list(
  n = nrow(wave3_wt_dat) + length(pred_t_wave3_wt),
  Y = c(as.integer(wave3_wt_dat$increase_3_v_2), rep(NA, length(pred_t_wave3_wt))),
  titre = c(log(wave3_wt_dat$`2`), log(pred_t_wave3_wt))
)

#Run model
wave3_wt_res <- run_model(data.list = model_dat)

#Check convergence
mcmcplot(wave3_wt_res,random = T)

saveRDS(wave3_wt_res,here("model_output","wave3_wt_mcmc.rds"))

#Plot model curve with thresholds
(wave3_wt_plot <- mmcc::tidy(as.mcmc(wave3_wt_res)) %>% 
    filter(str_detect(parameter,"S")) %>% 
    select(median,`2.5%`,`97.5%`) %>% 
    #mutate_all(function(x)1-x) %>% 
    slice(nrow(wave3_wt_dat)+1:n()) %>% 
    bind_cols(pred_t_wave3_wt) %>% 
    ggplot()+
    geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#d95f02", fill="#fc8d62", stat="identity")+
    geom_point(data=wave3_wt_dat,aes(y=increase_3_v_2, x=`2`),pch=124,size=5,alpha=0.5,colour="#fc8d62")+
    geom_pointrange(data=extract_ab_thresholds(wave3_wt_res,"thresh_50"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
    geom_pointrange(data=extract_ab_thresholds(wave3_wt_res,"thresh_80"),aes(x=`0.5`,xmin=`0.025`, xmax=Inf,y=median))+
    geom_hline(data=mmcc::tidy(as.mcmc(wave3_wt_res)) %>%
                 filter(parameter%in%c("a","c")),
               aes(yintercept=median),
               linetype="dashed")+
    scale_x_log10("Ancestral S-specific IgG pre-wave (WHO BAU/ml)" )+
    scale_y_continuous("Probability of increased Delta S-specific\nIgG titres following Delta wave",labels = scales::percent)+
    ggtitle("Ancestral vs. Delta wave")
)

wave3_wt_plot +
theme_minimal()&
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))&
  scale_fill_brewer()&
  coord_cartesian(xlim=c(0.5,NA),expand=F)

ggsave(here("results","ancestral_vs_delta.png"),width=150,height=150,units="mm",dpi=600,bg="white")

#Table 1
bind_rows(
  mmcc::tidy(as.mcmc(wave3_wt_res)) %>% 
    mutate(wave="Wave 3")) %>% 
  filter(parameter %in% c("a","c")) %>% 
  mutate(across(c(median,`2.5%`,`97.5%`),function(x)x*100)) %>% 
  bind_rows(extract_ab_thresholds(wave3_wt_res,thresh="a",mult=0.5) %>% 
              mutate(wave="Wave 3",parameter="a_0.5")%>% 
              select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave3_wt_res,thresh="thresh_50") %>% 
              mutate(wave="Wave 3")%>% 
              select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave3_wt_res,thresh="thresh_80")%>% 
              mutate(wave="Wave 3")%>% 
              select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  mutate(`97.5%`=ifelse(`97.5%`>250,Inf,`97.5%`),
         across(c(median,`2.5%`,`97.5%`),~formatC(.,digits=1,format = "f"))) %>% 
  mutate(estimate=paste0(median, " (",`2.5%`,", ",`97.5%`,")")) %>% 
  select(wave,parameter,estimate) %>% 
  pivot_wider(values_from = estimate,names_from = parameter) %>% 
  select(wave,a,c,thresh_50,thresh_80,a_0.5) %>% 
  htmlTable::htmlTable(rownames = FALSE, header=c("Wave", 
                                                    "Probability of increased titres at minimal pre-wave antibody levels (%, 95% CrI)",
                                                    "Probability of increased titres at maximal pre-wave antibody levels (%, 95% CrI)",
                                                    "50% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                    "80% reduction threshold (WHO BAU/ml, median, 95% CrI)",
                                                    "IgG titres required to reduce probability of seroconversion by 50% (WHO BAU/ml, median, 95% CrI)"))


# Goodness of fit plot
(ancestral_vs_delta_gof <- remove_geom(wave3_wt_plot,"GeomPointrange")+
  geom_pointrange(data=wave3_wt_dat %>% 
                    mutate(titre_group=cut(`2`,breaks=c(-Inf,1.09,10,100,Inf))) %>% 
                    group_by(titre_group) %>% 
                    summarise(N=n(),
                              avg_titre=median(`2`),
                              x=sum(increase_3_v_2),
                              p=x/N,
                              p.lo=binom::binom.confint(x,N,methods = "exact")$lower,
                              p.hi=binom::binom.confint(x,N,methods = "exact")$upper),
                  aes(x=avg_titre,y=p,ymin=p.lo,ymax=p.hi))+
  theme_minimal()&
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))&
  scale_fill_brewer()&
  coord_cartesian(xlim=c(0.5,1000),expand=F))


ggsave(ancestral_vs_delta_gof,here("results","ancestral_vs_delta_gof.png"),width=150,height=150,units="mm",dpi=600,bg="white")

