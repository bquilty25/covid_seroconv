
source("scripts/utils.R")

wav=2

decay_est <- dat %>% 
  filter( (variant == "WT" & wave == wav-1) | (variant == "WT" & wave == wav)) %>% 
  select(-c(variant,collection_date, n_doses)) %>% 
  mutate(wave=ifelse(wave==wav-1,"igg_pre","igg_post")) %>%
  pivot_wider(values_from = igg,names_from = wave) %>% 
  drop_na(igg_pre,igg_post) %>% 
  #add in time between bloods
  left_join(dat %>% 
              filter( (variant == "WT" & wave == wav-1) | (variant == "WT" & wave == wav)) %>% 
              select(-c(variant, igg)) %>%  
              mutate(wave=ifelse(wave==wav-1,"pre","post")) %>%
              pivot_wider(values_from = collection_date,names_from=wave) %>% 
              mutate(t_diff=as.numeric(difftime(units = "weeks",post,pre))) %>% 
              rename("date_pre"=pre,"date_post"=post) %>% 
              mutate(t_diff=ifelse(is.na(t_diff),mean(t_diff,na.rm=T),t_diff))) 


mod_dat <- decay_est %>% 
  select(-n_doses) %>% 
  filter(igg_post<igg_pre) %>% 
  mutate(date_post=as.numeric(difftime(date_post,date_pre),units="weeks"),
         date_pre=0) %>% 
  mutate(igg_baseline=igg_pre
  ) %>% 
  pivot_longer(names_to = c(".value","set"),
               names_pattern = "(.+)_(.+)",
               cols=c(igg_pre,igg_post,date_pre,date_post)) 


fit_brm <- mod_dat %>% 
  brms::brm(log(igg)~date+(1|pid_child),data=.)

summary(fit_lmer)
confint(fit_lmer)

parameters::model_parameters(fit_lmer)

fit_lmer %>% plot_model(.,type="pred")#,terms = c("date","igg_baseline"))

ggplot(mod_dat)+
  geom_path(aes(x=date,y=igg,group=pid_child))+
  scale_y_log10()


## look at decay of antibodies
datw4 %>%
  filter(variant == "WT") %>%
  mutate(wave = paste0("wave",wave)) %>%
  select(-c(n_doses,collection_date)) %>% 
  pivot_wider(names_from = wave, values_from = igg) %>%
  filter(!(wave1 <= wave2) | (wave2 <= wave3) | (wave3 <= wave4)) %>%
  #mutate(selected = (wave1 >= wave2) & (wave2 >= wave3) & (wave3 > wave4)) %>%
  pivot_longer(cols = wave1:wave4, names_to = "wave") %>%
  ggplot(aes(x = wave, y = value, group = pid_child))+#, color = selected))+
    geom_point(alpha = .4) +
    geom_line(alpha = .2) +
    scale_y_log10() +
    scale_colour_brewer("Declined",type="qual",palette = "Set1",labels=c("No","Yes"),direction=-1)+
    theme_classic()

## plot correlation of variant specific IgG with wild type
dat %>% filter(age=="child") %>% 
  pivot_wider(names_from = variant, values_from = igg) %>%
  pivot_longer(cols = Beta:Omicron, names_to = "Variant", values_to = "Other") %>%
  ggplot(aes(x=WT, y=Other, color = factor(wave))) +
  geom_point(size = 0.5, alpha =.5, stroke = 1) + 
  geom_abline(slope = 1, intercept = 0, lty = "dashed", lwd = 1, alpha=.5)+
  scale_x_log10("S-IgG titres to ancestral SARS-CoV-2") + scale_y_log10("S-IgG titres to variant SARS-CoV-2") +
  scale_colour_brewer("Titres after wave",type="qual",palette = "Set2", labels=c("WT","Beta", "Delta","Omicron"), direction=1)+
  facet_wrap(.~Variant) +
  theme_classic()+
  theme(strip.background = element_blank())
ggsave("results/titre_cor.png",width=200,height=100,units="mm",dpi=600,bg="white")
ggsave("results/titre_cor.pdf",width=200,height=100,units="mm",dpi=600,bg="white")


# look at sample dates
read_xlsx(here("data","data_for_billy_all_4waves_339_29MAR2022_with_vaccine.xlsx"),sheet = 1) %>% 
  filter(is.na(mat_vacc_covid)) %>%
  select(pid_child, CoV2S_1w_m, CoV2S_2w_m, CoV2S_3w_m, CoV2S_4w_m,
         CollectionDate_1w_m, CollectionDate_2w_m, CollectionDate_3w_m, CollectionDate_4w_m) %>%
  select(pid_child,CollectionDate_1w_m:CollectionDate_4w_m) %>%
  pivot_longer(CollectionDate_1w_m:CollectionDate_4w_m) %>%
  ggplot(aes(x=value)) +
    geom_histogram(bins=60,color = "white") +
    xlab("") + ylab("Number of samples collected") +
    scale_x_datetime(labels = date_format("%Y %b"),
                     date_breaks = "3 month") +
    theme_classic()

# WT IgG by collection date
# IgG = read_xlsx(here("data","data_for_billy_all_4waves_339_29MAR2022_with_vaccine.xlsx"),sheet = 1) %>% 
#   filter(is.na(mat_vacc_covid)) %>%
#   select(pid_child, CoV2S_1w_m, CoV2S_2w_m, CoV2S_3w_m, CoV2S_4w_m) %>%
#   pivot_longer(cols = -1, names_pattern = "(.*)(....)$", names_to = c("strain","wave"), values_to = "IgG")
# wave = read_xlsx(here("data","data_for_billy_all_4waves_339_29MAR2022_with_vaccine.xlsx"),sheet = 1) %>% 
#   filter(is.na(mat_vacc_covid)) %>%
#   select(pid_child, CollectionDate_1w_m, CollectionDate_2w_m, CollectionDate_3w_m, CollectionDate_4w_m) %>%
#   pivot_longer(cols = -1, names_pattern = "(.*)(....)$", names_to = c("cat","wave"), values_to = "Date")
# merge(IgG, wave, by=c("pid_child","wave")) %>%
#   tibble() #%>%
  
sa_dat <- get_national_data(countries = "South Africa") 

sa_plot <- sa_dat%>% 
  mutate(date=as.Date(date)) %>%
  mutate(variant=cut(date,breaks=c(as.Date("2020-01-01"),
                                   as.Date("2020-12-01"), 
                                   as.Date("2021-04-01"), 
                                   as.Date("2021-11-15"), 
                                   lubridate::now()),
                                   labels=c("WT","Beta","Delta","Omicron"))) %>% 
  group_by(variant) %>% 
  mutate(peak=date[which.max(cases_new)]) %>% 
  ggplot(aes(x=date,y=cases_new,fill=variant))+
  geom_col(alpha=0.75)+
  labs(x="",y="Daily reported cases")+
  scale_fill_brewer("Variant wave",palette = "Set2")


sero_time <- datw4 %>% 
  mutate(collection_date=as.Date(collection_date)) %>% 
  group_by(wave) %>% 
  #mutate(collection_date=median(as.Date(collection_date))) %>% 
  filter(variant=="WT") %>% 
  arrange(-n_doses) %>% 
  ggplot() +
  geom_point(aes(x=collection_date, y=igg, group=pid_child,colour=factor(n_doses)),
             alpha=0.25,
             size=1) +
  geom_line(aes(x=collection_date, y=igg, group=pid_child,colour=factor(n_doses)),alpha=0.5,size=0.25) +
  scale_color_manual("Number of vaccine doses received before sampling", values=c("grey",brewer_pal(palette = "Set1",direction=-1)(2)))+
  geom_text(data=datw4 %>% 
              mutate(collection_date=as.Date(collection_date)) %>% 
              group_by(wave) %>% 
              summarise(min_cd=min(collection_date),
                        max_cd=max(collection_date),
                        mean_cd=mean(collection_date)),
            
            aes(x=mean_cd,y=15000,label=wave),
            nudge_y = 0.3
  )+
  ggnewscale::new_scale_colour()+
  geom_segment(data=datw4 %>% 
                 mutate(collection_date=as.Date(collection_date)) %>% 
                 group_by(wave) %>% 
                 summarise(min_cd=min(collection_date),
                           max_cd=max(collection_date),
                           mean_cd=mean(collection_date)),
               
               aes(x=mean_cd-15,xend=mean_cd+15,y=15000,yend=15000))+
  # geom_pointrange(data= map(results3,1) %>% bind_rows(.id="group") %>% left_join(map(results3,2) %>% bind_rows(.id="group"),by="group")  %>% select(-group) %>% separate(doses,into=c("doses_pre","doses_post"),sep="_") %>% filter(doses_pre=="Total",sero_pos_pre==FALSE,vacc_agnostic_thresh==TRUE) %>% 
  #                   separate(exp_tm,into = c("thresh_50","thresh_2.5","thresh_97.5"),sep=" ") %>% 
  #                   mutate(across(contains("thresh_"),~parse_number(.))),
  #                 aes(x=sa_dat%>% 
  #                       mutate(date=as.Date(date)) %>% 
  #                       mutate(variant=cut(date,breaks=c(as.Date("2020-01-01"),
  #                                                        as.Date("2020-12-01"), 
  #                                                        as.Date("2021-04-01"), 
  #                                                        as.Date("2021-11-15"), 
  #                                                        lubridate::now()),
  #                                          labels=c("WT","Beta","Delta","Omicron"))) %>% 
  #                       group_by(variant) %>%
  #                       filter(variant!="WT") %>% 
  #                       summarise(peak=date[which.max(cases_new)]) %>% pull(peak),#datw4 %>% filter(wave<4) %>% group_by(wave) %>% summarise(avg_cd=median(as.Date(collection_date))) %>% pull(avg_cd),
  #                     y=thresh_50,ymin=thresh_2.5,ymax=thresh_97.5,
  #                     colour=c("Beta","Delta","Omicron")),
  #                 pch="-",
  #                 fatten=10,size=0.5)+
 # scale_color_manual("50% protection threshold", values=c(brewer_pal(palette = "Set2",direction=1)(4)[2:4]),guide="none")+
  xlab("") + 
  ylab("WT IgG titre") +
  scale_y_log10()
  #scale_color_manual("Number of doses recieved sampling",values=c("#006d2c","#08519c","#a63603"))

#ggsave(paste0("results/sero_over_time.png"),width=210/1.5,height=297/1.5,units="mm",dpi=600,bg="white")  

b <- map(results3,"wave_change_plot")[[1]] %>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath") +
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[2])),guide="none")+
  ggtitle("Beta wave") + 
  annotate("text",x=c(1,2),y=c(15000),label=c("1","2"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)),size=0.25,alpha=0.5)

d <- map(results3,"wave_change_plot")[[3]]%>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath")+
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[3])),guide="none")+
  labs(y="",title = "Delta wave")+ annotate("text",x=c(1,2),y=c(15000),label=c("2","3"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)),size=0.25,alpha=0.5)

o <- map(results3,"wave_change_plot")[[7]]%>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath")+
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[4])),guide="none")+
  labs(y="",title="Omicron wave")+ 
  annotate("text",x=c(1,2),y=c(15000),label=c("3","4"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase_2_v_1)),size=0.25,alpha=0.5)

((sa_plot/sero_time)&
  scale_x_date(limits = c(as.Date("2020-09-01"),as.Date("2022-03-01")),date_labels = "%b %Y",breaks = "3 months"))/
  ((b+d+o)&scale_x_discrete(labels=c("Pre-wave","Post-wave"))&scale_y_log10("WT IgG titre",limits=c(NA,15000)))&plot_annotation(tag_level = "A")&
  theme_classic()&
  theme(legend.position = "bottom")

ggsave(paste0("results/sero_over_time.png"),width=210,height=297,units="mm",dpi=600,bg="white")  
ggsave(paste0("results/sero_over_time.pdf"),width=210,height=297,units="mm",dpi=600,bg="white")  


#children plot

sero_time <- dat %>% 
  filter(age=="child") %>% 
  mutate(collection_date=as.Date(collection_date)) %>% 
  group_by(wave) %>% 
  #mutate(collection_date=median(as.Date(collection_date))) %>% 
  filter(variant=="WT") %>% 
  ggplot() +
  geom_point(aes(x=collection_date, y=igg, group=pid_child,
                 #colour=factor(wave)
                 ),
             alpha=0.25,
             size=1,
             colour="grey"
             ) +
  geom_line(aes(x=collection_date, y=igg, group=pid_child),alpha=0.25,size=0.25,colour="grey") +
  geom_text(data=dat %>% 
              filter(age=="child") %>% 
              mutate(collection_date=as.Date(collection_date)) %>% 
              group_by(wave) %>% 
              summarise(min_cd=min(collection_date),
                        max_cd=max(collection_date),
                        mean_cd=mean(collection_date)),
            
            aes(x=mean_cd,y=15000,label=wave),
            nudge_y = 0.3
  )+
  #scale_colour_brewer("Variant wave",palette = "Set2")+
  ggnewscale::new_scale_colour()+
  geom_segment(data=dat %>% 
                 filter(age=="child") %>% 
                 mutate(collection_date=as.Date(collection_date)) %>% 
                 group_by(wave) %>% 
                 summarise(min_cd=min(collection_date),
                           max_cd=max(collection_date),
                           mean_cd=mean(collection_date)),
               
               aes(x=mean_cd-15,xend=mean_cd+15,y=15000,yend=15000))+
xlab("") + 
  ylab("WT IgG titre") +
  scale_y_log10()

b <- map(results_children,"wave_change_plot")[[1]] %>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath") +
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[2])),guide="none")+
  ggtitle("Beta wave") + 
  annotate("text",x=c(1,2),y=c(15000),label=c("1","2"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase)),size=0.25,alpha=0.5)

d <- map(results_children,"wave_change_plot")[[2]]%>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath")+
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[3])),guide="none")+
  labs(y="",title = "Delta wave")+ annotate("text",x=c(1,2),y=c(15000),label=c("2","3"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase)),size=0.25,alpha=0.5)

o <- map(results_children,"wave_change_plot")[[3]]%>% 
  remove_geom(.,geom_type = "GeomText") %>% 
  remove_geom(.,geom_type = "GeomPath")+
  scale_color_manual(values=rev(c("grey",brewer_pal(palette="Set2")(4)[4])),guide="none")+
  labs(y="",title="Omicron wave")+ 
  annotate("text",x=c(1,2),y=c(15000),label=c("3","4"))+
  geom_path(aes(x=variant,y=value,group=pid_child,colour=factor(increase)),size=0.25,alpha=0.5)

((sa_plot/sero_time)&
    scale_x_date(limits = c(as.Date("2020-05-01"),as.Date("2022-03-01")),date_labels = "%b %Y",breaks = "3 months"))/((b+d+o)&scale_x_discrete(labels=c("Pre-wave","Post-wave"))&scale_y_log10("WT IgG titre",limits=c(NA,15000)))&
  plot_annotation(tag_level = "A")&
  theme_classic()&
  theme(legend.position = "bottom")

ggsave(paste0("results/sero_over_time_children.png"),width=210,height=297,units="mm",dpi=600,bg="white")  
ggsave(paste0("results/sero_over_time_children.pdf"),width=210,height=297,units="mm",dpi=600,bg="white")  
