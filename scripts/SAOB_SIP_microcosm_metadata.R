library(tidyverse)
library(ggpubr)
library(viridis)


## read in metadata from  spreadsheet, binding rows of each worksheet together
  # experimental values file path
  meta_path = "~/Documents/UBC_Research/PostDoc/Daniel Mulat/Results/SIP_Acetate_BONCAT_Sept_2020_data_compilation.xlsx"
  # microcosm metadata file path
  bottle_meta = "~/Documents/UBC_Research/PostDoc/Daniel Mulat/Results/BONCAT-SIP Batch/SIP_Acetate_BONCAT_Sept_2020_metadata.xlsx"
  
  #read files
  meta <- read_xlsx(path = meta_path, sheet = '24h', col_names = TRUE) %>%
    rbind(read_xlsx(path = meta_path, sheet = '72h', col_names = TRUE)) %>% 
    rbind(read_xlsx(path = meta_path, sheet = '144h', col_names = TRUE)) %>% 
    rbind(read_xlsx(path = meta_path, sheet = '240h', col_names = TRUE)) %>% 
    rbind(read_xlsx(path = meta_path, sheet = '312h', col_names = TRUE)) %>% 
    rbind(read_xlsx(path = meta_path, sheet = '408h', col_names = TRUE)) %>%
    left_join(read_xlsx(path = bottle_meta, col_names = TRUE), by = "bottle")
    

# calculate cumulative sum of methane over time in each bottle
  ch4 <- meta %>% group_by(bottle) %>%
    mutate(cum.ch4 =  cumsum(v_ch4)) %>%
    select(bottle, date, time, cum.ch4, treatment, HPG, reactor) %>%
    filter(HPG == "no") %>%
    drop_na(cum.ch4)

# cumulative methane production
ggplot(ch4, 
       aes(x=time, y=cum.ch4, color = treatment, fill = treatment)) +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'ribbon', alpha = 0.1, colour = NA) +
  stat_summary(fun = "mean", geom = "point", size = 2) +  
  stat_summary(fun = "mean", geom = "line") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE,  option = "D") +
  ylab("Cumulative methane production (mL)") + 
  xlab("Time (hrs)") +
  theme_pubr() 

# vfa concentration
ggplot(meta %>% 
         filter(HPG == "no") %>%  
         drop_na(acetate_ppm, treatment), 
       aes(x=time, y=acetate_ppm, color = treatment, fill = treatment)) +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'ribbon', alpha = 0.1, colour = NA) +
  stat_summary(fun = "mean", geom = "point", size = 2) +  
  stat_summary(fun = "mean", geom = "line") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE,  option = "D") +
  ylab("Acetate concentration (mg/L)") + 
  xlab("Time (hrs)") +
  theme_pubr() 





