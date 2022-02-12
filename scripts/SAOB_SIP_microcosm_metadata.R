library(tidyverse)
library(ggpubr)
library(viridis)
library(readxl)
library(ggpubr)
library(cowplot)

## read in metadata from  spreadsheet, binding rows of each worksheet together
  # experimental values file path
  meta_path = "raw_data/SIP_Acetate_BONCAT_Sept_2020_data_compilation.xlsx"
  # microcosm metadata file path
  bottle_meta = "raw_data/SIP_Acetate_BONCAT_Sept_2020_metadata.xlsx"
  
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
ch4_plot <- ggplot(ch4, 
       aes(x=time, y=cum.ch4, color = treatment, fill = treatment)) +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'ribbon', alpha = 0.1, colour = NA) +
  stat_summary(fun = "mean", geom = "point", size = 2) +  
  stat_summary(fun = "mean", geom = "line") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  scale_fill_viridis(discrete = TRUE,  option = "D") +
  ylab("Cumulative methane production (mL)") + 
  xlab("Time (hrs)") + 
  theme_pubr()
ch4_plot

# vfa concentration
vfa_plot <- ggplot(meta %>% 
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

ggsave("figures/cumulative_ch4_plot_sip_timeseries.png", ch4_plot, width=12, height=8, units=c("cm"))

ggsave("figures/vfa_degradation_sip_timeseries.png", vfa_plot, width=12, height=8, units=c("cm"))

chem_grid <- plot_grid(ch4_plot, NULL, vfa_plot, nrow=1, labels = c("A", "", "B"), rel_widths = c(1.5, 0.05, 2))


# now add the title
title <- ggdraw() + 
  draw_label(
    "Cumulative Methane Production and VFA Uptake over the Time-Series",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
saob_grid <- plot_grid(
  title, chem_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1)
)

ggsave("figures/saob_sip_experiment_metadata_grid.png", saob_grid, width=15, height=12, units=c("cm"))

# Arrange with ggpubr
sip_experiment <- ggarrange(ch4_plot, vfa_plot, ncol=2, nrow=1, common.legend=TRUE, legend="bottom", widths=c(1,1.1), labels = c("A", "B"))

ggsave("figures/SIP-experiment-grid-metadata.png", sip_experiment, width=20, height=12, units=c("cm"))
