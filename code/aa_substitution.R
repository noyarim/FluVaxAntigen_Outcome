library(tidyverse)

## This code calculates total antigenic advancement weight per a pair of circulating and vaccine virus, 
## given the amino acid substitutions between the two.
## Written by: Kyueun Lee
## Generated Nov 1, 2025
## Last moditifed: Nov 16, 2025


setwd("/Users/kyueunlee/Library/CloudStorage/OneDrive-UW/My Drive/Z Drive/Project_all/2025_FluSeq/data")
# Nextstrain AA data
dt <- readxl::read_xlsx("h3n2_ha_12y_2025-10-01_metadata.xlsx")
# AA substitutions
dt <- dt %>%
  mutate(
    aa_sub1 = str_match(labels, "HA1:\\s*([^;'}]+)")[,2],
    aa_sub2 = str_match(labels, "HA2:\\s*([^;'}]+)")[,2],
    sigpep = str_match(labels, "SigPep:\\s*([^;'}]+)")[,2]
  )

# A function to flat a vector
flat_vector <- function(path_mutations){
  flat <- path_mutations %>%
  str_split(",") %>%     # split each element by comma
  unlist() %>%            # flatten into one vector
  str_trim() %>%
  discard(~ is.na(.x) | .x == "")
}

# A function to get the list of mutations
get_mutations_path <- function(dt, start_node, end_node) {
  path_mutations_ha1 <- c()
  path_mutations_ha2 <- c()
  path_mutations_sigpep <- c()
  current_node <- start_node
  vax_parent_node <- dt$parent_name[dt$name == end_node]
  # accumulate AA substitutions until the vax parent node
  while (current_node != vax_parent_node) {
    row <- dt %>% filter(name == current_node)
    if (nrow(row) == 0 || is.na(row$parent_name)) {
      message("Reached a node with no parent or parent not found: ", current_node)
      break
    }
    
    # extract mutations string
    muts_ha1 <- row$aa_sub1
    muts_ha2 <- row$aa_sub2
    muts_sigpep <- row$sigpep
    
    path_mutations_ha1 <- c(path_mutations_ha1, muts_ha1)
    path_mutations_ha2 <- c(path_mutations_ha2, muts_ha2)
    path_mutations_sigpep <- c(path_mutations_sigpep, muts_sigpep)
    
    #if (!is.na(muts) && muts != "character(0)") {
    #  path_mutations <- c(path_mutations, muts)
    #}
    
    current_node <- row$parent_name
  }
  # Add the AA substitutions from vax parent node to vax node
  last <- dt %>% filter(name == end_node & parent_name == current_node)
  muts_ha1 <- stringi::stri_reverse(last$aa_sub1)
  muts_ha2 <- stringi::stri_reverse(last$aa_sub2)
  muts_sigpep <- stringi::stri_reverse(last$sigpep)
  
  path_mutations_ha1 <- c(path_mutations_ha1, muts_ha1)
  path_mutations_ha2 <- c(path_mutations_ha2, muts_ha2)
  path_mutations_sigpep <- c(path_mutations_sigpep, muts_sigpep)
  
  v_ha1 <- flat_vector(path_mutations_ha1)
  v_ha2 <- flat_vector(path_mutations_ha2)
  v_sigpep <- flat_vector(path_mutations_sigpep)
  
  list(v_ha1, v_ha2, v_sigpep)
}

# List of vax and circulating viruses by season
dt_vax_clade <- readxl::read_xlsx("Clade_vax_AA_antgdistance.xlsx")

# Save aa substitutions between vax and circulating viruses
dt_vax_clade2 <- dt_vax_clade %>%
  mutate(clade = ifelse(clade == "3C.2a1b.2a.1","3C.2a1b.2a",clade))%>%
  rowwise()%>%
  mutate(
    start_node =min(dt$name[dt$clade_membership==clade & substr(dt$name,1,4) == "NODE"]),
    aa_list_ha1 = list(get_mutations_path(dt, start_node, vaccine)[[1]]),
    aa_list_ha2 = list(get_mutations_path(dt, start_node, vaccine)[[2]]),
    aa_list_sigpep = list(get_mutations_path(dt, start_node, vaccine)[[3]])
  )

# Transform the data into a long form
dt_vax_clade2_long <- dt_vax_clade2 %>%
  # Pivot the three list-columns into one named column
  pivot_longer(
    cols = starts_with("aa_list_"),
    names_to = "type",
    values_to = "aa_list"
  ) %>%
  # Clean up the 'type' names (remove the 'aa_list_' prefix)
  mutate(type = str_remove(type, "^aa_list_")) %>%
  # Expand each list into separate rows
  unnest_longer(aa_list)

write.csv(as.data.frame(dt_vax_clade2_long),"Clade_vax_AA_antgdistance_long.csv")


# per substitution dataset
weight <- read.csv("substitution_weights_by_timepoint.csv")
weight <- weight %>%
  mutate(
    type = ifelse(substr(substitution,1,3)=="HA1","ha1", ifelse(substr(substitution,1,3)=="HA2", "ha2", "sigpep")),
    aa_list = str_extract(substitution, "(?<=:)\\S+"),
    Season = ifelse(month(timepoint)==10,year(timepoint),year(timepoint)-1)
  ) %>%
  filter(month(timepoint)==10)

# link weight dataset to the analytic dataset
dt_vax_clade2_long <- read.csv("Clade_vax_AA_antgdistance_long.csv")

dt_vax_clade2_long_wt <- dt_vax_clade2_long %>%
  left_join(weight, by = c("Season","type","aa_list")) %>%
  mutate(weight = ifelse(is.na(weight),0,weight))

write.csv(dt_vax_clade2_long_wt, "Clade_vax_AA_antgdistance_weight.csv")

# calculate the total weight per vaccine-clade pair per season
summary_weight <- dt_vax_clade2_long_wt %>%
  group_by(Season, clade, vaccine) %>%
  summarize(
    tot_weight=sum(weight)
  )
write.csv(summary_weight, "Total_weight_summary(October).csv")

