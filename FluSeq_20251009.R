############################################################################
## This code is to map the geographical heterogeneity in clade frequencies 
## Created on Oct 9, 2025
## Written by Kyu Lee
## Last motified on Nov 12, 2025 
############################################################################

## Data Description:
# Downloaded from GISAID, EpiFlu(https://platform.epicov.org/epi3/frontend#4630fc)
# Search Filter: 1) Type A H3N2 human, 
#                2) Human host
#                3) North America > United States
#                4) Collection date from 2009-08-01 to 2019-07-31
# Seasonal cycle is defined as the period from 08/01 to 07/31 in the following year

library(readxl)
library(tidyverse)
library(usmap)
library(maps)
library(viridis)



setwd("/Users/kyueunlee/Library/CloudStorage/OneDrive-UW/My Drive/Z Drive/Project_all/2025_FluSeq")
# aa substitution weight data
weight <- read.csv("data/Total_weight_summary.csv")
# built-in state data
state_dt <- data.frame(state = state.name, state_abb = state.abb)


for (season in seq(2009,2019)){
  print(season)
  # clade frequency data
  dt <- read_excel(paste0("data/gisaid_epiflu_isolates_",season,".xls"))
  
  dt <- dt %>%
    mutate(state = str_extract(Isolate_Name, "(?<=^[A-Z]/)[^/]+"),
           date = as.Date(Collection_Date),
           month = month(date)) %>%
    left_join(state_dt, by="state") %>%
    filter(!is.na(state_abb))
  
  
  clade_summary <- dt %>%
    group_by(state_abb, Clade) %>%
    summarise(
      n = n()
    ) %>%
    group_by(state_abb) %>%
    mutate(freq = n / sum(n),
           Season=season) %>%
    ungroup()
  
  cal_weight <- clade_summary %>%
    left_join(weight, by=c("Season","Clade"="clade")) %>%
    mutate(perc_tot_weight = freq*tot_weight)%>%
    group_by(state_abb, Season)%>%
    summarize(state_tot_weight = sum(perc_tot_weight))%>%
    rename(state = state_abb)
  
  # Distribution of clade/subclade
  ggplot(clade_summary , aes(x=state_abb, y=freq, fill=Clade))+
    geom_col(position="fill") +
    theme_bw() +
    scale_y_continuous(expand=c(0,0))+
    ggtitle(paste0(season,"-", season+1))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("graph/clade_distribution_",season,".pdf"),width=7,height=5)
  
  # US heatmap with antigenic weights
  plot_usmap(data = cal_weight, values = "state_tot_weight") +
    scale_fill_viridis(
      option = "plasma",
      name = "Total Weight",
      na.value = "grey90"
    ) +
    labs(
      title = "US State Heatmap by Weight (Season 2019)",
      subtitle = "Each state's total weight value"
    ) +
    theme_minimal() +                               # replaces buggy default theme
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA)
    )
  ggsave(paste0("graph/weight_distribution_",season,".pdf"),width=7,height=5)
  
}

