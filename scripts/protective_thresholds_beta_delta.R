source("scripts/utils.R")

dat <- dat %>% 
  filter(variant == "WT" & wave == 1 | variant == "Beta" & wave == 2 | variant=="Delta" & wave == 3) %>% 
  select(-variant) %>% 
  pivot_wider(values_from = value,names_from = wave) %>% 
  mutate(increase_2_v_1=as.integer(`2`>`1`*1.1),
         increase_3_v_2=as.integer(`3`>`2`*1.1)) #increase if >10% increase in titer

#Run model for post-wave 2
wave2_dat <- dat %>% 
  drop_na(`1`,`2`)

(wave2_change_plot <- wave2_dat %>% 
  pivot_longer(c(`1`,`2`,`3`)) %>% 
  filter(name!=3) %>% 
  mutate(variant=ifelse(name==1,"Ancestral","Beta")) %>% 
  ggplot(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)))+
  geom_point(alpha=0.2)+
  geom_path(alpha=0.2)+
  labs(title="Beta wave")+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))
)

pred_t_wave2 <- lseq(from=min(wave2_dat$`1`),to=max(wave2_dat$`1`),length.out = 1000)
wave2_res <- run_model(data.list = list(n=nrow(wave2_dat)+length(pred_t_wave2),
                                         Y=c(as.integer(wave2_dat$increase_2_v_1),rep(NA,length(pred_t_wave2))),
                                         titre=c(log(wave2_dat$`1`),log(pred_t_wave2))))

mcmcplot(wave2_res,random = T)
saveRDS(wave2_res,here("model_output", "wave2_mcmc.rds"))

#Run model for post-wave 3
wave3_dat <- dat %>% 
  #filter(!(`2`<1&`3`<1)) %>% 
  drop_na(`2`,`3`)

(wave3_change_plot <- wave3_dat %>% 
  pivot_longer(c(`1`,`2`,`3`)) %>% 
  filter(name!=1) %>% 
  mutate(variant=ifelse(name==2,"Beta","Delta")) %>% 
  ggplot(aes(x=variant,y=value,group=pid_child,colour=factor(increase_3_v_2)))+
  geom_point(alpha=0.2)+
  geom_path(alpha=0.2)+
  scale_y_log10("Anti-Spike IgG titre")+
  theme_minimal()+
  labs(title="Delta wave")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))
  )

pred_t_wave3 <- lseq(from=min(wave3_dat$`2`),to=max(wave3_dat$`2`),length.out = 1000)
wave3_res <- run_model(data.list = list(n=nrow(wave3_dat)+length(pred_t_wave3),
                                        Y=c(as.integer(wave3_dat$increase_3_v_2),rep(NA,length(pred_t_wave3))),
                                        titre=c(log(wave3_dat$`2`),log(pred_t_wave3))))

mcmcplot(wave3_res,random = T)
saveRDS(wave3_res,here("model_output","wave3_mcmc.rds"))

plot_a <- wave2_change_plot+wave3_change_plot+plot_layout(guides="collect")&
  theme(legend.position = "bottom")&
  scale_y_log10("S-specific IgG titre (WHO BAU/ml)",limits=c(NA,10000))&
  labs(x="")&
  scale_colour_brewer("Increase post-wave",type="qual",palette = "Set1",labels=c("No","Yes"),direction=-1)

ggsave(here("results", "change_plot.png"),width=210,height=150,units="mm",dpi=600,bg="white")

(wave2_plot <- mmcc::tidy(as.mcmc(wave2_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
  slice(nrow(wave2_dat)+1:n()) %>% 
  bind_cols(pred_t_wave2) %>% 
  ggplot()+
  geom_smooth(aes(ymin=`2.5%`, ymax=`97.5%`,y=median,x=...4),colour="#2ca25f", fill="#99d8c9", stat="identity")+
  geom_point(data=wave2_dat,aes(y=as.integer(increase_2_v_1), x=`1`),pch=124,size=5,alpha=0.5,colour="#2ca25f")+
    geom_pointrange(data=extract_ab_thresholds(wave2_res,"thresh_50"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
    geom_pointrange(data=extract_ab_thresholds(wave2_res,"thresh_80"),aes(x=`0.5`,xmin=`0.025`, xmax=`0.975`,y=median))+
  geom_hline(data=mmcc::tidy(as.mcmc(wave2_res)) %>%
              filter(parameter%in%c("a","c")),
            aes(yintercept=median),
            linetype="dashed")+
  scale_x_log10("WT S-specific IgG pre-wave (WHO BAU/ml)" )+
  scale_y_continuous("Probability of increased Beta S-specific\nIgG titres following Beta wave",labels = scales::percent)+
  ggtitle("Beta wave")
)

(wave3_plot <- mmcc::tidy(as.mcmc(wave3_res)) %>% 
  filter(str_detect(parameter,"S")) %>% 
  select(median,`2.5%`,`97.5%`) %>% 
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
  scale_y_continuous("Probability of increased Delta S-specific\nIgG titres following Delta wave",labels = scales::percent)+
  ggtitle("Delta wave")
)

plot_b <- wave2_plot+wave3_plot&
  #scale_x_log10("IgG pre-wave (WHO BAU/ml)" )&
  #scale_y_continuous("Probability of increased titres following wave",labels = scales::percent)&
  theme_minimal()&
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA),
        legend.position = "bottom")&
  scale_fill_brewer()&
  coord_cartesian(xlim=c(0.5,1000),expand=F)

ggsave(here("results", "res_logistic.png"),width=210,height=150,units="mm",dpi=600,bg="white")

(plot_a/plot_b)+plot_annotation(tag_levels = "A")
ggsave(here("results","combined.png"),width=210,height=210,units="mm",dpi=600,bg="white")


#Table 1
bind_rows(
  mmcc::tidy(as.mcmc(wave2_res)) %>% 
    mutate(wave="Wave 2"),
  mmcc::tidy(as.mcmc(wave3_res)) %>% 
    mutate(wave="Wave 3")) %>% 
  filter(parameter %in% c("a","c")) %>% 
  mutate(across(c(median,`2.5%`,`97.5%`),function(x)x*100)) %>% 
  bind_rows(extract_ab_thresholds(wave2_res,thresh="a",mult=0.5) %>% 
              mutate(wave="Wave 2",parameter="a_0.5") %>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave3_res,thresh="a",mult=0.5) %>% 
              mutate(wave="Wave 3",parameter="a_0.5")%>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave2_res,thresh="thresh_50") %>% 
              mutate(wave="Wave 2") %>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave3_res,thresh="thresh_50") %>% 
              mutate(wave="Wave 3")%>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave2_res,thresh="thresh_80") %>% 
              mutate(wave="Wave 2") %>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  bind_rows(extract_ab_thresholds(wave3_res,thresh="thresh_80")%>% 
              mutate(wave="Wave 3")%>% select(parameter,wave,`2.5%`=`0.025`,median=`0.5`,`97.5%`=`0.975`)) %>% 
  mutate(`97.5%`=ifelse(`97.5%`>1000,Inf,`97.5%`),
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

#goodness of fit
wave2_gof <- remove_geom(wave2_plot,"GeomPointrange")+
  geom_pointrange(data=wave2_dat %>% 
                    mutate(titre_group=cut(`1`,breaks=c(-Inf,1.09,10,100,Inf))) %>% 
                    group_by(titre_group) %>% 
                    summarise(N=n(),
                              avg_titre=median(`1`),
                              x=sum(increase_2_v_1),
                              p=x/N,
                              p.lo=binom::binom.confint(x,N,methods = "exact")$lower,
                              p.hi=binom::binom.confint(x,N,methods = "exact")$upper),
                  aes(x=avg_titre,y=p,ymin=p.lo,ymax=p.hi))

wave3_gof <- remove_geom(wave3_plot,"GeomPointrange")+
  geom_pointrange(data=wave3_dat %>% 
                    mutate(titre_group=cut(`2`,breaks=c(-Inf,1.09,10,100,Inf))) %>% 
                    group_by(titre_group) %>% 
                    summarise(N=n(),
                              avg_titre=median(`2`),
                              x=sum(increase_3_v_2),
                              p=x/N,
                              p.lo=binom::binom.confint(x,N,methods = "exact")$lower,
                              p.hi=binom::binom.confint(x,N,methods = "exact")$upper),
                  aes(x=avg_titre,y=p,ymin=p.lo,ymax=p.hi))

#This requires running the other script also
wave2_gof+wave3_gof+ancestral_vs_delta_gof&
  theme_minimal()&
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = NA))&
  scale_fill_brewer()&
  coord_cartesian(xlim=c(0.5,1000),expand=F)

ggsave(here("results","combined_gof.png"),width=300,height=100,units="mm",dpi=600,bg="white")

