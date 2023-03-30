# Noori Choi
# MS project 2018 analysis
# 2021-08-26

# import packages ---------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(lme4)
library(nlme)
library(egg)
library(hrbrthemes)
library(lmerTest)
library(ggcorrplot)
library(emmeans)
library(multcomp)
library(multcompView)
library(lattice)
library(ggeffects)
library(MuMIn)

# import dataset ----------------------------------------------------------
bout_df <- read.csv("MS_bout+noise.csv", header=T, na.strings=c("NA"))
overlap_df <- read.csv("pianka_idx.csv", header=T, na.strings=c("NA"))
temp_df <- read.csv("temperature.csv", header=T, na.strings=c("NA"))
schiz_df <- read.csv("MS_bout+schiz_interaction.csv", header=T, na.strings=c("NA"))
pitfall_df <- read.csv("pitfall.csv", header=T, na.strings=c("NA"))

bout_df %>% 
  mutate(duration = bout_glob_end - bout_glob_begin) -> bout_df

spider <- c('2', '3', '5')
bird <- c('0', '1', '6', '7', '18', '20', '26','35')
cricket <- c('10')
human <- c('31', '33')
unknown_vi <- c('4', '9', '45')
unknown_ac <- c('13', '26', '39')
date_ls = c('180518', '180520', '180522', '180525', '180527', 
            '180531', '180604', '180607', '180611', '180614', 
            '180618', '180624', '180628', '180703')

bout_df %>% 
  group_by(glob_cat) %>% # filter glob_cat occur less than 100
  filter(n() > 100,
         glob_cat != 'unknown')%>% # mutate type
  mutate(type = if_else(glob_cat %in% spider, '0', 
                        if_else(glob_cat %in% bird, '1',
                                if_else(glob_cat %in% cricket, '2',
                                        if_else(glob_cat %in% human, '3', 
                                                if_else(glob_cat %in% unknown_vi, '4', '5')
                                        ))))) -> bout_df_100 

# Vibratory soundscape (figure 3 & 4) -------------------------------------
## dataframe -------------------------------------------------------------
spider <- c('2', '3', '5')
bird <- c('0', '1', '6', '7', '18', '20', '26','35')
cricket <- c('10')
human <- c('31', '33')
unknown_vi <- c('4', '9', '45')
unknown_ac <- c('13', '26', '39')
date_ls = c('180518', '180520', '180522', '180525', '180527', 
          '180531', '180604', '180607', '180611', '180614', 
          '180618', '180624', '180628', '180703')

bout_df %>% 
  group_by(glob_cat) %>% # filter glob_cat occur less than 100
  filter(n() > 100,
         glob_cat != 'unknown')%>% # mutate type
  mutate(type = if_else(glob_cat %in% spider, '0', 
                        if_else(glob_cat %in% bird, '1',
                                if_else(glob_cat %in% cricket, '2',
                                        if_else(glob_cat %in% human, '3', 
                                                if_else(glob_cat %in% unknown_vi, '4', '5')
                                                ))))) -> bout_df_100 


## temporal partitioning + temperature (figure 3) --------------------------
### generate labels for wrapped plot
date.labs <- c('May 18', 'May 20', 'May 22', 'May 25', 'May 27', 'May 31',
               'June 04', 'June 07', 'June 11', 'June 14', 'June 18', 'June 24',
               'June 29', 'July 03')
names(date.labs) <- date_ls

### temporal partitioning (figure 3a)
bout_df_100 %>% 
  filter(type %in% c('0', '1', '2', '3', '4')) %>%
  mutate(glob_cat = factor(glob_cat, levels = c('0', '1', '6', '7', '18', '20', '26','35', # birds
                                                '10', # cricket
                                                '31', '33', # human noise
                                                '4', '9', '45', # unknown vibratory sounds 
                                                '2', '3', '5')), # spiders
         date = factor(date, levels = date_ls),
         day_time = if_else(bout_glob_begin > 3600*24, bout_glob_begin - 3600*24, bout_glob_begin),
         time_bin = day_time %/% 3600,
         time_bin = factor(time_bin, levels = seq(0, 23))) %>% 
  group_by(glob_cat, date, time_bin, .drop=FALSE) %>% 
  tally() %>% 
  group_by(glob_cat) %>% 
  mutate(prop = n/max(n)) %>% 
  ggplot(aes(x=time_bin, y=glob_cat)) +
  geom_tile(aes(alpha=prop, fill=glob_cat), color='white') +
  scale_alpha(range = c(0 ,1)) +
  scale_x_discrete(breaks = c(0, 5, 10, 15, 20)) +
  scale_y_discrete(labels=c('0' = 'American crow', '1' = 'Blue jay',
                            '6' = 'Red-eyed vireo', '7' = 'Pine warbler',
                            '18' = 'Northern cardinal', '20' = 'Yellow brested chat',
                            '26' = 'Whipper whirl','35' = 'Eastern wood pewee',
                            '10' = 'Cricket',
                            '31' = 'Airplane', '33' = 'unknown noise',
                            "2" = "S.duplex", "3" = "S. stridulans",
                            "5" = "S. uetzi", "4" = "unknown_1",
                            "9" = "unknown_2", "45" = "unknown_3")) +
  facet_grid(. ~date, switch = 'x', labeller = labeller(date = date.labs)) +
  labs(x = "Date & Time", y = "Species") +
  theme_classic() +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "transparent", linetype='blank'),
        strip.text = element_text(size = 18, face = "bold"),
        axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold.italic"),
        axis.text.x = element_text(vjust = 0.5, size = 14, face = "bold"),
        legend.position = "none") -> date_time_fig

### temperature (figure 3b)
temp_df %>% 
  filter(date %in% date_ls) %>% 
  mutate(habitat = if_else(plot %in% c('A', 'C', 'E'), 'Leaf litter', 'Pine litter'),
         day_time = ifelse(rec_time > 3600*24, rec_time - 3600*24, rec_time),
         time_bin = day_time %/% 3600,
         time_bin = factor(time_bin, levels = seq(0, 23))) %>% 
  group_by(date, time_bin) %>% 
  dplyr::summarize(m_temp = mean(temp, na.rm=TRUE)) %>% 
  ggplot(aes(x=time_bin, y=m_temp, group=1)) +
  geom_line() +
  geom_point() +
  scale_x_discrete(breaks = c(0, 5, 10, 15, 20)) +
  facet_grid(. ~date, switch = 'x', labeller = labeller(date = date.labs)) +
  labs(x = "Date & Time", y = "Temperature ('C)") +
  theme_classic() +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "transparent", linetype='blank'),
        strip.text = element_text(size = 18, face = "bold"),
        axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold.italic"),
        axis.text.x = element_text(size = 14, face = "bold"),
        legend.position = "none") -> temp_fig

### combine figure 3a & 3b
ggarrange(date_time_fig, temp_fig, ncol= 1, heights = c(2, 1))

## spatial signal space partitioning (figure 4) ----------------------------
bout_df_100 %>% 
  filter(type %in% c('0', '4')) %>% 
  mutate(habitat = if_else(plot %in% c('A', 'C', 'E'), 'Leaf litter', 'Pine litter')) %>% 
  group_by(glob_cat, habitat) %>% 
  dplyr::summarize(bout_n = n()/n_distinct(plot)) %>% 
  mutate(freq = bout_n/sum(bout_n),
         glob_cat = factor(glob_cat, levels = c('2', '3', '5', # spider
                                                '4', '9', '45'))) %>% # unknown_ac
  ggplot(aes(x=glob_cat, y=freq, fill=habitat)) +
  geom_bar(position = 'fill', stat = 'identity') + 
  labs(x = "Species", y = "Occurrence proportion") +
  scale_fill_discrete(name = "Habitat type",
                      labels = c("Leaf litter", "Pine litter")) +
  scale_x_discrete(labels=c("2" = "S.duplex", "3" = "S. stridulans",
                            "5" = "S. uetzi", "4" = "unknown_1",
                            "9" = "unknown_2", "45" = "unknown_3")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 14, face = "bold"),
        axis.text.x = element_text(vjust = 0.5, size = 14, face = "bold.italic"),
        legend.title=element_text(size=14, face='bold'),
        legend.text=element_text(size=12))

## spectral partitioning ---------------------------------------------------
bout_df_100 %>% 
  filter(type %in% c('0', '1', '2', '3', '4')) %>%
  mutate(glob_cat = factor(glob_cat, levels = c('0', '1', '6', '7', '18', '20', '26','35', # birds
                                                '10', # cricket
                                                '31', '33', # human noise
                                                '4', '9', '45', # unknown vibratory sounds 
                                                '2', '3', '5')) # spiders
         ) %>% 
  ggplot(aes(x=glob_cat, y=med_dom_freq, fill=type)) +
  geom_boxplot() +
  scale_fill_discrete(name = 'Category', 
                      labels=c('0' = 'Birds', '1' = 'Cricket',
                                '2' = 'Antrhopogenic noise', '3' = 'Unknown vibratory sounds',
                                '4' = 'Spiders')) +
  scale_x_discrete(labels=c('0' = 'American crow', '1' = 'Blue jay',
                            '6' = 'Red-eyed vireo', '7' = 'Pine warbler',
                            '18' = 'Northern cardinal', '20' = 'Yellow brested chat',
                            '26' = 'Whipper whirl','35' = 'Eastern wood pewee',
                            '10' = 'Cricket',
                            '31' = 'Airplane', '33' = 'unknown noise',
                            "2" = "S.duplex", "3" = "S. stridulans",
                            "5" = "S. uetzi", "4" = "unknown_1",
                            "9" = "unknown_2", "45" = "unknown_3")) +
  labs(x = "Vibratory sound types", y = "Dominant frequency (Hz)") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 20, face = "bold.italic", angle=90),
        legend.position = "right",
        legend.title = element_text(size=20, face='bold'),
        legend.text = element_text(size=18))

# Correlation with pitfall trapping (figure 5) ----------------------------
## pearson correlation test ------------------------------------------------
## correlation between heatmap and soundscape data
### merge dataframes
#### pitfall data
pitfall_df %>% 
  mutate(species = factor(species, levels = c('duplex', 'stridulans', 'uetzi')),
         date = factor(date),
         time = factor(time, levels = c('0', '8', '16')),
         plot = factor(plot, levels = c('A', 'B', 'C', 'D', 'E', 'F')),
         habitat = ifelse(plot %in% c('A', 'C', 'E'), 'leaf litter', 'pine litter')) %>%
  complete(., species, date, time, habitat,
           fill = list(male = 0)) %>% 
  group_by(species, date, time) %>% 
  dplyr::summarise(n_indiv = sum(male)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(prop = n_indiv/max(n_indiv)) -> pitfall_join

#### soundscape data
bout_df_100 %>% 
  filter(type == 0) %>%
  mutate(glob_cat = factor(glob_cat),
         date = factor(date, levels = date_ls),
         day_time = if_else(bout_glob_begin > 3600*24, bout_glob_begin - 3600*24, bout_glob_begin),
         time_bin = day_time %/% 28800,
         time_bin = factor(time_bin, levels = c(0, 1, 2)), 
         time = ifelse(time_bin == "1", "8", 
                       ifelse(time_bin == "2", "16", "0")),
         species = ifelse(glob_cat=='2', 'duplex',
                          ifelse(glob_cat=='3', 'stridulans', 'uetzi'))) %>% 
  #group_by(species, date, time_bin) %>% 
  complete(., species, date, time,
           fill = list(male = 0)) %>%
  group_by(species, date, time) %>%
  tally() %>%
  group_by(species) %>% 
  mutate(prop = n/max(n)) -> sound_join

### correlation analysis
corr_df <- inner_join(pitfall_join, sound_join, by=c('species', 'date', 'time'))
cor.test(corr_df$prop.x, corr_df$prop.y) # s = 0.560, p < 0.001

## pitfall trap data + correlation plot (figure 5) -------------------------
### temporal pitfall trap data (figure 5a) ----------------------------------
#### generate labels for wrapped plots
date.labs <- c('May 18', 'May 20', 'May 22', 'May 25', 'May 27','May 31',
               'June 04', 'June 07', 'June 11', 'June 14', 'June 18', 'June 24', 'June 29',
               'July 03')
names(date.labs) <- c('180518', '180520', '180522', '180525', '180527', '180531',
                      '180604', '180607', '180611', '180614', '180618', '180624', '180629',
                      '180703')

#### heatmap of pitfall trapping data
pitfall_df %>% 
  mutate(species = factor(species, levels = c('duplex', 'stridulans', 'uetzi')),
         date = factor(date),
         time = factor(time, levels = c('0', '8', '16')),
         plot = factor(plot, levels = c('A', 'B', 'C', 'D', 'E', 'F')),
         habitat = ifelse(plot %in% c('A', 'C', 'E'), 'leaf litter', 'pine litter')) %>%
  complete(., species, date, time, habitat,
           fill = list(male = 0)) %>% 
  group_by(species, date, time) %>% 
  dplyr::summarise(n_indiv = sum(male)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(prop = n_indiv/max(n_indiv),
         n_indiv = ifelse(n_indiv==0, NaN, n_indiv)) %>% 
  ggplot(aes(x=time, y=species)) +
  geom_tile(aes(alpha=prop, fill=species)) +
  geom_text(aes(label=n_indiv), size=8, face='bold') +
  scale_alpha_continuous(range = c(0, 1)) +
  scale_y_discrete(labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  scale_fill_manual(values = c('#f37735', '#00b159', 'blue')) +
  facet_wrap(date ~., ncol = 13, strip.position = 'bottom', 
             labeller = labeller(date = date.labs), scales="free_x") +
  labs(x = "Date", y = "") +
  theme_classic() +
  theme(strip.placement = "outside",                      
        strip.background = element_rect(fill = "transparent", linetype='blank'),
        strip.text = element_text(size = 20, face = "bold"),
        axis.title = element_blank(),
        axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold.italic"),
        legend.position = "None")

# correlation btw soundscape vs. pitfall (figure 5b) ----------------------
corr_df %>% 
  ggplot(aes(x=prop.x, y=prop.y)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, color='red') +
  labs(x = "Proportion of specimen from pitfall trapping", y = "Proportion of detected sounds") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        legend.position = "none")
  
# Acoustic niche partitioning (figure 6) ----------------------------------
## correlation among signal space dimensions -------------------------------
### pearson correlation test
overlap_df %>% 
  filter(cat_1_type == 'vi',
         cat_2_type == 'vi') %>% 
  mutate(schiz = ifelse(cat_1 %in% spider & cat_2 %in% spider, '1', '0'),
         temp_pi = date_pi*time_pi) %>%
  dplyr::select(c(plot_pi, temp_pi, freq_pi)) %>% 
  cor() 

overlap_df %>% 
  filter(cat_1_type == 'vi',
         cat_2_type == 'vi') %>% 
  mutate(schiz = ifelse(cat_1 %in% spider & cat_2 %in% spider, '1', '0'),
         temp_pi = date_pi*time_pi) %>%
  dplyr::select(plot_pi, temp_pi, freq_pi) %>% 
  cor_pmat()

### figure 5
#### spatial vs. temporal (figure 6a)
overlap_df %>% 
  filter(cat_1_type == 'vi'|
         cat_2_type == 'vi') %>% 
  mutate(schiz = ifelse(cat_1 %in% spider & cat_2 %in% spider, '1', '0'),
         temp_pi = date_pi*time_pi) %>% 
  ggplot(aes(x=plot_pi, y=temp_pi)) +
  geom_point(aes(color=schiz), size=5) +
  geom_smooth(method=lm, se=FALSE, color='red') +
  labs(x = "Spatial overlap", y = "Temporal overlap") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        legend.position = "none") -> plot_temp_pi

#### spatial vs. spectral (figure 6b)
overlap_df %>% 
  filter(cat_1_type == 'vi',
         cat_2_type == 'vi') %>% 
  mutate(schiz = ifelse(cat_1 %in% spider & cat_2 %in% spider, '1', '0'),
         st_overlap = plot_pi * time_pi * date_pi) %>% 
  ggplot(aes(x=plot_pi, y=freq_pi)) +
  geom_point(aes(color=schiz), size=5) +
  geom_smooth(method=lm, se=FALSE, color='red') +
  labs(x = "Spatial overlap", y = "Spectral overlap") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        legend.position = "none") -> plot_freq_pi

#### temporal vs. spectral (figure 6c)
overlap_df %>% 
  filter(cat_1_type == 'vi',
         cat_2_type == 'vi') %>% 
  mutate(schiz = ifelse(cat_1 %in% spider & cat_2 %in% spider, '1', '0'),
         temp_pi = time_pi * date_pi) %>% 
  ggplot(aes(x=temp_pi, y=freq_pi)) +
  geom_point(aes(color=schiz), size=5) +
  geom_smooth(method=lm, se=FALSE, color='red') +
  labs(x = "Temporal overlap", y = "Spectral overlap") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        legend.position = "none") -> temp_freq_pi

ggarrange(plot_temp_pi, plot_freq_pi, temp_freq_pi, nrow=1)


# Variation in Schiz signal properties by general noise (figure 7) --------
## linear mixed effects model ----------------------------------------------
### dataframe
bout_df_100 %>% 
  filter(type == "0") %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) -> bout_df_schiz

bout_df_schiz %>% 
  filter(glob_cat == "2") -> dup_df
bout_df_schiz %>% 
  filter(glob_cat == "3") %>% 
  mutate(idle = ifelse(n_idle > 0, 1, 0))-> str_df
bout_df_schiz %>% 
  filter(glob_cat == "5") -> uet_df

### duration, dominant frequency, peak rate
#### duplex
dr_nn_dup <- lmer(log(duration) ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=dup_df, REML=FALSE)
df_nn_dup <- lmer(med_dom_freq ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=dup_df, REML=FALSE)
pr_nn_dup <- lmer(peak_rate ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=dup_df, REML=FALSE)

MuMIn::r.squaredGLMM(dr_nn_dup)
MuMIn::r.squaredGLMM(df_nn_dup)
MuMIn::r.squaredGLMM(pr_nn_dup)

anova(dr_nn_dup)
anova(df_nn_dup)
anova(pr_nn_dup)

#### stridulans
dr_nn_str <- lmer(log(duration) ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=str_df, REML=FALSE)
df_nn_str <- lmer(med_dom_freq ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=str_df, REML=FALSE)
pr_nn_str <- lmer(peak_rate ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=str_df, REML=FALSE)

MuMIn::r.squaredGLMM(dr_nn_str)
MuMIn::r.squaredGLMM(df_nn_str)
MuMIn::r.squaredGLMM(pr_nn_str)

anova(dr_nn_str)
anova(df_nn_str)
anova(pr_nn_str)

#### uetzi
dr_nn_uet <- lmer(log(duration) ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=uet_df, REML=FALSE)
df_nn_uet <- lmer(med_dom_freq ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=uet_df, REML=FALSE)
pr_nn_uet <- lmer(peak_rate ~ n_noise_30m*ent_noise_30m + (1|mean_temp), data=uet_df, REML=FALSE)

anova(dr_nn_uet)
anova(df_nn_uet)
anova(pr_nn_uet)

## figures (figure 7a, 7b, 7c) ------------------------------------------------------
### duration (figure 7a)
bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%
  ggplot(aes(x=n_noise_30m, y=log(duration), color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  labs(x = "Number of noise", y = "Bout duration (log(s))") +
  scale_color_discrete(name = "Species",
                      labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> dr_Nnum

bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%  # spiders
  ggplot(aes(x=ent_noise_30m, y=log(duration), color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  labs(x = "Shannon entropy", y = "Bout duration (log(s))") +
  scale_color_discrete(name = "Species",
                       labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.title=element_text(size = 30, face='bold'),
        legend.text=element_text(size = 28)) -> dr_Nent

ggarrange(dr_Nnum, dr_Nent, ncol= 2)

### dominant frequency (figure 7b)
bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%  # spiders
  ggplot(aes(x=n_noise_30m, y=med_dom_freq, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  labs(x = "Number of noise", y = "Dominant frequency (Hz)") +
  ylim(0, 2500) +
  scale_color_discrete(name = "Species",
                       labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> df_Nnum

bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%  # spiders
  ggplot(aes(x=ent_noise_30m, y=med_dom_freq, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  ylim(0, 2500) +
  labs(x = "Shannon entropy", y = "Dominant frequency (Hz)") +
  scale_color_discrete(name = "Species",
                       labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.title=element_text(size=30, face='bold'),
        legend.text=element_text(size=28)) -> df_Nent

ggarrange(df_Nnum, df_Nent, ncol= 2)

### peak rate (figure 7c)
bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%  # spiders
  ggplot(aes(x=n_noise_30m, y=peak_rate, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  labs(x = "Number of noise", y = "Peak rate (/s)") +
  #ylim(0, 2500) +
  scale_color_discrete(name = "Species",
                       labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> pr_Nnum

bout_df_100 %>% 
  filter(glob_cat %in% spider) %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) %>%  # spiders
  ggplot(aes(x=ent_noise_30m, y=peak_rate, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  #ylim(0, 2500) +
  labs(x = "Shannon entropy", y = "Peak rate (/s)") +
  scale_color_discrete(name = "Species",
                       labels = c("S. duplex", "S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.title=element_text(size = 30, face='bold'),
        legend.text=element_text(size = 28)) -> pr_Nent

ggarrange(pr_Nnum, pr_Nent, ncol= 2)

# Variation in Schiz signal properties by con/heterospecifics (figure 8) ------
## dataframe ---------------------------------------------------------------
schiz_df %>% 
  group_by(glob_cat) %>% 
  mutate(hs_dg = n_hs_30m/max(n_hs_30m),
         cs_dg = n_cs_30m/max(n_cs_30m),
         duration = bout_glob_end - bout_glob_begin) -> schiz_df

schiz_df %>% 
  filter(glob_cat == "3") %>% 
  mutate(idle = ifelse(n_idle > 0, 1, 0)) -> str_df
schiz_df %>% 
  filter(glob_cat == "5") -> uet_df

## linear mixed effects model ----------------------------------------------
### str
dr_hs_str <- lmer(duration ~ hs_dg*cs_dg + (1|mean_temp), data=str_df, REML=FALSE)
df_hs_str <- lmer(med_dom_freq ~ hs_dg*cs_dg + (1|mean_temp), data=str_df, REML=FALSE)
pr_hs_str <- lmer(peak_rate ~ hs_dg*cs_dg + (1|mean_temp), data=str_df, REML=FALSE)

anova(dr_hs_str)
anova(df_hs_str)
anova(pr_hs_str)

### uet
dr_hs_uet <- lmer(duration ~ hs_dg*cs_dg + (1|mean_temp), data=uet_df, REML=FALSE)
df_hs_uet <- lmer(med_dom_freq ~ hs_dg*cs_dg + (1|mean_temp), data=uet_df, REML=FALSE)
pr_hs_uet <- lmer(peak_rate ~ hs_dg*cs_dg + (1|mean_temp), data=uet_df, REML=FALSE)

anova(dr_hs_uet)
anova(df_hs_uet)
anova(pr_hs_uet)

## figures (figure 8) ------------------------------------------------------
### durations
schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=cs_dg, y=log(duration), color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  ylim(1, 6.1) +
  labs(x = "The abudance of conspecific signal", y = "Duration (log(s))") +
  scale_color_discrete(name = "Species",
                       labels = c("S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> dr_cs

schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=hs_dg, y=log(duration), color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  ylim(1, 6.1) +
  labs(x = "The abudance of heterospecific signal", y = "Duration (log(s))") +
  scale_color_discrete(name = "Species",
                       labels = c("S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> dr_hs

ggarrange(dr_cs, dr_hs, ncol=2)

### df
schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=cs_dg, y=med_dom_freq, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  ylim(5, 2000) +
  labs(x = "The abudance of conspecific signal", y = "Dominant frequency (Hz)") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> df_cs

schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=hs_dg, y=med_dom_freq, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  ylim(5, 2000) +
  labs(x = "The abudance of heterospecific signal", y = "Dominant frequency (Hz)") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'none') -> df_hs

ggarrange(df_cs, df_hs, ncol=2)

### peak rate
schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=cs_dg, y=peak_rate, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  scale_y_continuous(breaks=c(5, 10, 15, 20, 25)) +
  labs(x = "The abudance of conspecific signal", y = "Peak rate (/s)") +
  scale_color_discrete(name = "Species",
                       labels = c("S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'None') -> pr_cs

schiz_df %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('3', '5'))) %>%  # spiders
  ggplot(aes(x=hs_dg, y=peak_rate, color=glob_cat)) +
  geom_point(aes(color=glob_cat)) +
  geom_smooth(method=lm, se=TRUE) +
  scale_y_continuous(breaks=c(5, 10, 15, 20, 25)) +
  labs(x = "The abudance of heterospecific signal", y = "Peak rate (/s)") +
  scale_color_discrete(name = "Species",
                       labels = c("S. stridulans", "S. uetzi")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_blank(),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        legend.position = 'None',
        legend.title=element_text(size = 30, face='bold'),
        legend.text=element_text(size = 28)) -> pr_hs

ggarrange(pr_cs, pr_hs, ncol=2)

# stridulans complexity (figure 9) ----------------------------------------
## dataframe ---------------------------------------------------------------
### variation by general noise
bout_df_100 %>% 
  filter(type == "0") %>% 
  mutate(glob_cat = factor(glob_cat, levels = c('2', '3', '5')),
         ent_noise_30m = ent_noise_30m/log(t_noise_30m)) -> bout_df_schiz

bout_df_schiz %>% 
  filter(glob_cat == "3") %>% 
  mutate(idle = ifelse(n_idle > 0, 1, 0))-> noise_str

### variation by heterospecific sounds (S. uetzi)
schiz_df %>% 
  group_by(glob_cat) %>% 
  mutate(hs_dg = n_hs_30m/max(n_hs_30m),
         cs_dg = n_cs_30m/max(n_cs_30m),
         duration = bout_glob_end - bout_glob_begin) -> schiz_df

schiz_df %>% 
  filter(glob_cat == "3") %>% 
  mutate(idle = ifelse(n_idle > 0, 1, 0)) -> hs_str

## generalized linear mixed effects model ----------------------------------
### ~ n_noise*ent_noise
idle_nn_str <- lmer(n_idle ~ n_noise_30m*ent_noise_30m + (1|mean_temp), 
                    data=noise_str)
anova(idle_nn_str)
MuMIn::r.squaredGLMM(idle_nn_str)

### ~ hs_dg
idle_str <- lmer(n_idle ~ cs_dg*hs_dg + (1|mean_temp), data=hs_str)
anova(idle_str)
MuMIn::r.squaredGLMM(idle_str)

## figure 9 ----------------------------------------------------------------
### probability of idle ~ cs_dg
ggpredict(idle_str, "cs_dg [all]") %>% 
  plot() +
  theme_set(theme_ggeffects()) +
  labs(x = "The abundance of conspecific density", y = "Predicted number of idle in a bout") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        plot.title = element_blank(),
        legend.position = 'none')

### probability of idle ~ hs_dg
ggpredict(idle_str, "hs_dg [all]") %>% 
  plot() +
  theme_set(theme_ggeffects()) +
  labs(x = "The abudance of heterospecific signal", y = "Predicted number of idle in a bout") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 15)),
        axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 15)),
        axis.text.y =  element_text(size = 30, face = "bold"),
        axis.text.x = element_text(size = 30, face = "bold"),
        plot.title = element_blank(),
        legend.position = 'none') 




