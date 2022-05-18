
source("scripts/utils.R")

## look at decay of antibodies
datw4 %>%
  filter(variant == "WT") %>%
  mutate(wave = paste0("wave",wave)) %>%
  pivot_wider(names_from = wave, values_from = value) %>%
  mutate(selected = (wave1 >= wave2) & (wave2 >= wave3) & (wave3 > wave4)) %>%
  pivot_longer(cols = wave1:wave4, names_to = "wave") %>%
  ggplot(aes(x = wave, y = value, group = pid_child, color = selected))+
    geom_point(alpha = .4) +
    geom_line(alpha = .2) +
    scale_y_log10() +
    scale_colour_brewer("Declined",type="qual",palette = "Set1",labels=c("No","Yes"),direction=-1)+
    theme_classic()

## plot correlation of variant specific IgG with wild type
datw4 %>%
  pivot_wider(names_from = variant, values_from = value) %>%
  pivot_longer(cols = Beta:Omicron, names_to = "Variant", values_to = "Other") %>%
  ggplot(aes(x=WT, y=Other, color = factor(wave))) +
  geom_point(size = 2, alpha =.2, stroke = 1) + 
  geom_abline(slope = 1, intercept = 0, lty = "dashed", lwd = 1, alpha=.5)+
  scale_x_log10() + scale_y_log10() +
  scale_colour_brewer("Titres after wave",type="qual",palette = "Set2", labels=c("WT","Beta", "Delta","Omicron"), direction=-1)+
  facet_wrap(.~Variant) +
  theme_classic()
ggsave("results/titre_cor.png",width=200,height=80,units="mm",dpi=600,bg="white")

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
  ggplot(aes(x=date,y=cases_new))+geom_col()+labs(x="Date",y="Daily reported cases")


sero_time <- datw4 %>% 
  mutate(collection_date=as.Date(collection_date)) %>% 
  filter(variant=="WT") %>% 
  arrange(-n_doses) %>% 
  #select(pid_child, wave, IgG, Date) %>%
  ggplot(aes(x=collection_date, y=igg, group=pid_child,colour=factor(n_doses),alpha=factor(n_doses))) +
  geom_point(size=1,
  ) +
  geom_line() +
  xlab("") + ylab("WT IgG titre") +
  scale_y_log10()

(sero_time/sa_plot)&
  scale_x_date(limits = c(as.Date("2020-01-01"),as.Date("2022-05-01")),date_labels = "%b %Y")&
  theme_classic()&
  scale_alpha_manual(values=c(0.1,0.2,0.2),guide="none")&
  scale_color_manual("Number of doses recieved sampling",values=c("grey",brewer_pal(palette = "Set1",direction = -1)(2)))
  