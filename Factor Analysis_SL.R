rm(list = ls())

## Libraries
library(here) # directory management
library(tidyverse) # data wrangling and plotting
library(readxl) # read Excel files
library(FactoMineR) # Main package for FAMD and MC analyses
library(factoextra) # helper functions for FactoMineR
library(lme4) #this was needed to install missMDA
library(missMDA) # for imputation of missing values
library(forcats) # tools for working with categorical or nominal variables
library(dplyr)


setwd("C://Users/gunas015/Desktop/SL Comp")
data<-read.csv("Data_FAMD_Hpop.csv")

#making factor variable
data$mgtsys<-factor(data$mgtsys)
data$Risk<-factor(data$Risk)
data$endemic<-factor(data$endemic)
data$epiunit<-factor(data$epiunit)

data<- data %>%
  mutate(epiunit = recode(epiunit,
                        "herd" = "Herd",
                        "vsrange" = "VS Range",
                        "farm"= "Farm",
                        "village"= "Village",
                        "communal"= "Communal"
  ))


data1<-data %>%
  rename(
    "Number slaughtered" = SH,
    "Milk collecting centers"= MilkCC,
    "Management type"= mgtsys,
    "Low value animals"= animalLV,
    "FMD risk"= Risk,
    "FMD status"= endemic,
    "Epi unit"= epiunit,
    "New animals"= newanimal,
    "Environment"= enviornment,
    "Communal grazing"= communal,
    "Outside personnel"= people,
    "Disease reporting"= diseasereport,
    "Selling animals"= soldanimals,
    "FMD vaccination"= vaccination,
    "Nothing specific"= nothing,
    "Vehicle_VS office"= vehicle,
    "Illegal transport"= itransport,
    "Eartagged animals"= eartag,
    "Human population"= Hpop
  )



#for factor analysis, continuous variables should be standardized

set.seed(1)
res.famd <- FAMD(data1, sup.var = 1, graph = FALSE)
eig.val <- get_eigenvalue(res.famd)

ggplot(eig.val, aes(x = factor(1:nrow(eig.val)), y = variance.percent)) +
  geom_col(fill = "grey70") +
  geom_line(aes(group = 1), color = "black") +
  geom_point(color = "black") +
  labs(x = "Diminension", 
       y = "Variance explained (%)",
       title = "Scree Plot \n (explained variance)") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), "cm"),
    legend.position = "none",
    axis.title.x = element_text(size = 20, face = "bold", vjust = -2),
    axis.title.y = element_text(size = 20, face = "bold", vjust = 3),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", colour="gray50", hjust = 0.5)
  )


## Variable Contributions  
#Show which variables (survey responses) most contributed to each of the FAMD dimensions.  
#```{r fig.height=10, fig.width=10}
quant_contrib <- as.data.frame(res.famd$quanti.var$contrib)
quant_contrib$var <- rownames(quant_contrib)

qual_contrib <- as.data.frame(res.famd$quali.var$contrib)
qual_contrib$var <- rownames(qual_contrib)

contrib <- bind_rows(quant_contrib, qual_contrib) %>%
  pivot_longer(
    cols = starts_with("Dim"),
    names_to = "dimension",
    values_to = "contribution"
  )


# top variables for each dimension
top_contrib <- contrib %>%
  group_by(dimension) %>%
  slice_max(order_by = contribution, n = 10) %>%
  ungroup()


ggplot(top_contrib, aes(x = reorder(var, contribution), y = contribution)) +
  geom_segment(aes(xend = var, y = 0, yend = contribution), color = "grey70") +
  geom_point(size = 3, color = "steelblue") +
  coord_flip() +
  facet_wrap(~ dimension, scales = "free_y") +
  labs(
    x = NULL,
    y = "Contribution (%)",
    title = "Top Contributing Variables"
  ) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), "cm"),
    legend.position = "none",
    strip.text     = element_text(size = 18, face = "bold", color = "gray40"),
    axis.title.x = element_text(size = 18, face = "bold", vjust = -2),
    axis.title.y = element_text(size = 16, face = "bold", vjust = 3),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 18, face = "bold", colour="gray50", hjust = 0.5)
  )


## HCPC
##Hierarchical clustering on principal components (HCPC).  Let's us determine which rows (aka *individuals*) are most similar, and then, it will assign them to like or similar groups (clusters).

res.hcpc <- HCPC(res.famd, nb.clust = -1, graph = FALSE)

# How many groups were identified?
length(unique(res.hcpc$data.clust$clust))


### Individuals and Groups Factors Map  
#Each row represents a different *individual*.  
#In our case, the individuals also represent geographic locations. 
#Therefore, by identifying groups or clusters, we can determine which geographic locations are most similar.    
#This can be viewed with a scatter plot (below), which shows that three similar groups were identified and these two dimensions explain about 23\% of the variation.  
#Overall, the 2nd group is *tight* and most distinct, but all have some overlap.     


# Positions of individuals
ind <- as.data.frame(get_famd_ind(res.famd)$coord)

# add cluster assignments from HCPC
ind$Cluster <- factor(res.hcpc$data.clust$clust)

# label axes with the percent each dimension explains
eig <- get_eigenvalue(res.famd)
lab_x <- sprintf("Dim 1 (%.1f%%)", eig[1,2])
lab_y <- sprintf("Dim 2 (%.1f%%)", eig[2,2])


ggplot(ind, aes(x = Dim.1, y = Dim.2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(aes(fill = Cluster), geom = "polygon", alpha = 0.2, color = NA) +
  coord_equal() +
  labs(x = lab_x, 
       y = lab_y, 
       color = "Cluster", 
       fill = "Cluster") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), "cm"),
    legend.position = "bottom",
    axis.title.x = element_text(size = 20, face = "bold", vjust = -2),
    axis.title.y = element_text(size = 20, face = "bold", vjust = 3),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 22, face = "bold", colour="gray50", hjust = 0.5)
  )


## Predictor Data frame
#Collect all dimensions and select variables to be used in other models.


# fix axis names
coord_std <- as_tibble(res.famd$ind$coord, rownames = "row_id") %>%
  rename_with(~ gsub("^Dim[\\. ]([0-9]+)$", "Dim_\\1", .x)) %>%
  mutate(row_id = as.integer(row_id))

# extract all axes from FAMD #not working
dim_scores <- get_famd_dim(coord_std,K = 5)

## choose some variables to model directly # 
##From here onward
#____________________________________________________________#
# base_cov <- tibble(
#   row_id = seq_len(nrow(X_complete)),
#   Q3 = X_complete$Q3, # location names
#   Q36_epiunit     = X_complete$Q36 %||% factor(NA),
#   Q4_managemnet = X_complete[["Q4...2"]] %||% factor(NA),
#   Q18_risk  = X_complete$Q18 %||% factor(NA),
#   Q6_endemic   = X_complete$Q32 %||% factor(NA)
# )

base_cov <- tibble(
    row_id =data1$ID,
    District=data1$District,
   )

# Qualtrics lookup key (may not be needed)
# get_q41 <- function(df_num){
#   if ("Q41_3" %in% names(df_num)) return(df_num$Q41_3)
#   q41_cols <- names(df_num)[grepl("^Q41_[2-6]$", names(df_num))]
#   if (length(q41_cols) > 0) rowMeans(df_num[q41_cols], na.rm = TRUE) else rep(NA_real_, nrow(df_num))
# }
# Q41_num <- get_q41(df2 %>% select(where(is.numeric)))
# num_block <- tibble(row_id = seq_len(nrow(df2)), Q41_key = Q41_num) 

# HCPC clusters for reference
clust_df <- res.hcpc$data.clust %>%
  tibble::rownames_to_column("row_id") %>%
  transmute(row_id = as.integer(row_id), cluster = clust)

# data frame with 5 axes and select variables
model_df <- coord_std %>%
  left_join(base_cov,  by = "row_id") %>%
  #left_join(num_block, by = "row_id") %>%
  left_join(clust_df,  by = "row_id") 
#%>%
  #relocate(Q3, .after = row_id) %>%
# mutate(
#   across(c(Q36_epiunit, Q4_managemnet, Q18_risk, Q6_endemic),
#            ~ if (is.factor(.x)) forcats::fct_na_value_to_level(.x, level = "(Missing)")
#            else forcats::fct_na_value_to_level(as.factor(.x), level = "(Missing)"))



#Quick check of result

head(model_df)[1:8]
    
## Write Variables
    
write_csv(model_df, here("risk_model_predictors_2025_10-23.csv"))








