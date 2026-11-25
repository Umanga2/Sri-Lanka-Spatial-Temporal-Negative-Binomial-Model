rm(list = ls())
library(here)
library("devtools")

# install.packages(
#   "INLA",
#   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
#   dep = TRUE)
library(sp)
library(INLA)
library(spdep)
#library(rgeos);
library(ggpubr); library(tidyverse)
library(gridExtra)
library(sf)
#library(ggsn)
#library(GGally)
#install.packages("remotes")
#remotes::install_github("oswaldosantos/INLAOutputs")
library(INLAOutputs)


#______________________________________________________________________________________#
#' Plots neighborhood adjacency
#'
#' This function takes in a polygon object and its corresponding neighborhood object,
#' and plots the polygons and their neighbors on a map.

plot_neighbors <- function(poly_obj, nb_obj) {
  
  if (inherits(poly_obj, "Spatial")) {
    poly_obj <- st_as_sf(poly_obj)
  }
  
  # get coordinates
  coords <- st_coordinates(st_centroid(poly_obj))
  coords_df <- as.data.frame(coords)
  names(coords_df) <- c("long", "lat")
  
  # neighbors to lines as sf
  neighbors_sl <- nb2lines(nb_obj, coords = coords_df, as_sf = TRUE)
  
  # ensure same CRS
  st_crs(neighbors_sl) <- st_crs(poly_obj)
  
  p <- ggplot() +
    geom_sf(data = poly_obj, fill = "gray90", color = "blue") +
    geom_sf(data = neighbors_sl, color = "red") +
    geom_point(data = coords_df, aes(x = long, y = lat), color = "black") +
    theme_minimal() +
    ggtitle("Neighborhood Adjacency") +
    ylab("Northing") +
    xlab("Easting") +
    coord_sf() +
    theme(
      axis.title.x = element_text(size = 22, face = "bold"),
      axis.title.y = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )
  
  plot(p)
}

#__________________________________________________________________________#


rm(list = ls())

#______Not working__________________________________________________________________
# install.packages("INLA",repos=c(getOption("repos"),
#                                 INLA="https://inla.r-inla-download.org/R/stable"), 
#                 dep=TRUE)
#_______________________________________________________________________________________
# INLA installation step working 
#https://stackoverflow.com/questions/66622811/problems-installing-the-inla-package
#
#_______________________________________________________________________________________

# plot adjacency graph function
#source("plot_neighbors.R")

setwd("C://Users/")

data<-read.csv("BYM4_FAMD.csv") 
#%>%
#  select(-c(ID2, ID3)) # dropping these here to make sure they match IDs in "area" object

# spelled different in the polygons, without the "h" in polygons
data$district[data$district == "Mullaithivu"] <- "Mullaitivu"


# number of cattle, buffalo, and goat farmers, 
#source: https://www.statistics.gov.lk/Agriculture/StaticalInformation/LivestockFarmers#gsc.tab=0
# pop_dens <- read.csv("local/pop_farms_density.csv")
# data <- left_join(data, pop_dens, by = "district") %>%
#   mutate(pop_s = as.numeric(scale(pop_2023))
#   )
# This is already included in the data set now 
# which(is.na(data$pop_2023)) # check for good matches

#map
area <- st_read("lka_admbnda_adm2_slsd_20220816.shp")
area$ID<- 1:nrow(area)
area$ID<- str_to_title(area$ID)
str(area$ID)

# # get IDs from area
# match_id <- as.data.frame(area) %>%
#   mutate(district = ADM2_EN) %>%
#   select(ID, district)

# match to data by district
data <- left_join(data, match_id, by="district")
which(is.na(data$ID)) # verify all IDs matched

# check names
ggplot(area) +
  geom_sf(size = .05,color="Grey") +
  geom_sf_label(aes(label =ADM2_EN))+
  scale_fill_distiller(palette="PuRd", direction=1)+
  theme(text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 10),
        legend.key.size = unit(10, "points"))

#______________________________________________________________________________#
# check clusters identified from survey data
#______________________________________________________________________________#

match_var <- data %>%
  select(district, cluster) %>%
  distinct()

area <- area %>%
  left_join(
    match_var,
    by = c("ADM2_EN" = "district")
  )

ggplot(area) +
  geom_sf(aes(fill = as.factor(cluster)), size = 0.05, color = "grey40") +
  scale_fill_brewer(palette = "PuRd", direction = 1, name = "Cluster") +
  theme_minimal() +
  theme(
    text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 10),
    legend.key.size = unit(10, "points"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#_______________________________________________________________________________#

# neighbor matrix
SL_nb <- poly2nb(as(area, "Spatial"), row.names = area$ID) #to create adjacency matrix 
nb2INLA("SL_adj", SL_nb)
SL_adj <- "SL_adj"
summary(SL_nb)

# quick plot check (function sourced above)
plot_neighbors(area, SL_nb)


data<-data%>%
  mutate(id =as.numeric((ID)),
         expected_cases = TotalL*sum(OB)/sum(TotalL),
         incidence_cases = OB / TotalL * 1e5,
         rr_cases = OB / expected_cases,
         intercept = 1) # manually added intercept

#___________________________________________________________________________#
# simple nonspatial model
non_spatial_form <- OB ~ -1 + intercept + 
  Dim_1 + Dim_2 + Dim_3 + Dim_4 + Dim_5 

testGLM <- inla(non_spatial_form,
                family = "poisson",
                data = data,
                E = expected_cases,
                offset = log(data$farms_cattle),
                control.fixed = list(prec = 1, prec.intercept=1),
                control.predictor = list(link = 1, compute = TRUE),
                control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)

summary(testGLM)

## penalized prior
prec <- list(prec = list(prior = "pc.prec", param = c(.3/.31, .01))) #precision
prec_bym2 <- list(phi = list(prior = "pc", param = c(.5, 2/3)),#mixing 
                  prec = list(prior = "pc.prec", param = c(.3/.31, .01)))

# check the mean and variance 
# Poisson is best when mean = variance
# when variance > mean, there may be overdispersion
# should try a zero-inflated or negative-binomial
mean_OB <- mean(data$OB, na.rm=TRUE); mean_OB
var_OB  <- var(data$OB, na.rm=TRUE); var_OB

#______________________________________________________________________________________#

#Spatial Model [model selection done using the spatial model] 
# Compare Poisson, nbinomial, and zeroinflated Poisson
RR01.poi <- inla(OB ~ 1 +
                   f(id,#this is important, cz earlier ID was changed to id str_to_title 
                     model = "bym2",
                     hyper = prec_bym2,
                     graph = SL_adj,
                     scale.model = TRUE,
                     constr = TRUE),
                 family = "poisson",
                 data =data,
                 verbose=FALSE,
                 offset = log(Farm.N),
                 E =expected_cases,
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE))


RR01.nb <- inla(OB ~ 1 +
                  f(id,#this is important, cz earlier ID was changed to id str_to_title 
                    model = "bym2",
                    hyper = prec_bym2,
                    graph = SL_adj,
                    scale.model = TRUE,
                    constr = TRUE),
                family = "nbinomial",
                data =data,
                verbose=FALSE,
                offset = log(Farm.N),
                E =expected_cases,
                control.predictor = list(link = 1, compute = TRUE),
                control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE))


RR01.zip <- inla(OB ~ 1 +
                   f(id,#this is important, cz earlier ID was changed to id str_to_title 
                     model = "bym2",
                     hyper = prec_bym2,
                     graph = SL_adj,
                     scale.model = TRUE,
                     constr = TRUE),
                 family = "zeroinflatedpoisson1",
                 data =data,
                 verbose=FALSE,
                 offset = log(Farm.N),
                 E =expected_cases,
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE))

# compare waic
RR01.poi$waic$waic
RR01.nb$waic$waic
RR01.zip$waic$waic

# compare dic
RR01.poi$dic$dic
RR01.nb$dic$dic
RR01.zip$dic$dic

# looks like the nbinomial may be best

##___________with temporal resolution_____________________________________________________


data<-data %>%
  mutate(id_year = year - min(year) + 1,
         id_year1 = year - min(year) + 1,
         id_year2 = year - min(year) + 1,
         id_year3 = year - min(year) + 1,
         id_id_year = 1:n()) #id_id_year is for interaction effect 



pc_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.05)))
pc_iid <- list(prec = list(prior = "pc.prec", param = c(0.8, 0.05))) 
pc_rw      <- list(prec = list(prior = "pc.prec", param = c(3, 0.1)))


## Type 1 interaction
form_I <- OB ~ -1 + intercept +
  f(ID,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw1",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(id_id_year, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior)

RR_type1 <- inla(form_I,
                 family = "nbinomial",
                 data = data,
                 E = expected_cases,
                 offset = log(farms_cattle),
                 control.fixed = list(prec = 1, prec.intercept=1),
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)

summary(RR_type1)


## Type II interaction
form_II <- OB ~ -1 + intercept +
  f(ID,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw1",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(id,
    model="iid", 
    group=id_year2,
    control.group=list(model="rw1", hyper = pc_rw),
    constr=TRUE,
    hyper=pc_iid) 

RR_type2 <- inla(form_II,
                 family = "nbinomial",
                 data = data,
                 E = expected_cases,
                 offset = log(Farm.N),
                 control.fixed = list(prec = 1, prec.intercept=1),
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)

summary(RR_type2)

## Type III interaction
form_III <- OB ~ -1 + intercept +
  f(ID,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw1",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(id_year2,
    model="iid", 
    group=ID,
    control.group=list(model="besag", graph=SL_adj, hyper=pc_prior),
    constr=TRUE,
    hyper=pc_iid) 

RR_type3 <- inla(form_III,
                 family = "nbinomial",
                 data = data,
                 E = expected_cases,
                 offset = log(Farm.N),
                 control.fixed = list(prec = 1, prec.intercept=1),
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)

summary(RR_type3)

#For model validation different steps 
ExplainedVariance(RR_type3)

data$muRR_type3<- RR_type3$summary.fitted.values[,"mean"]

write_csv(data,("INLA_NBMsummary_2025_10-11.csv"))

## Type IV interaction
form_IV <- OB ~ -1 + intercept +
  f(id,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw2",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(ID.x, 
    model="besag", 
    graph=SL_adj,
    hyper = pc_prior,
    group=id_year3,
    control.group=list(model="rw1", hyper=pc_rw, scale.model=TRUE))

RR_type4 <- inla(form_IV,
                 family = "nbinomial",
                 data = data,
                 E = expected_cases,
                 offset = log(farms_cattle),
                 control.fixed = list(prec = 1, prec.intercept=1),
                 control.predictor = list(link = 1, compute = TRUE),
                 control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)

summary(RR_type4)


# compare waic
RR_type1$waic$waic
RR_type2$waic$waic
RR_type3$waic$waic
RR_type4$waic$waic

# compare dic
RR_type1$dic$dic
RR_type2$dic$dic
RR_type3$dic$dic
RR_type4$dic$dic

### Type 3 looks the best## ________________________________________________

## Modify Type III formula to add covariates
form_IIImod <- OB ~ -1 + intercept +
  f(id,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw1",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(id_year2,
    model="iid", 
    group=ID,
    control.group=list(model="besag", graph=SL_adj, hyper=pc_prior),
    constr=TRUE,
    hyper=pc_iid) + 
  Dim_1 + Dim_3 


RR_type3mod <- inla(form_IIImod,
                    family = "nbinomial",
                     data = data,
                    E = expected_cases,
                    offset = log(Farm.N),
                    control.fixed = list(prec = 1, prec.intercept=1),
                    control.predictor = list(link = 1, compute = TRUE),
                   control.compute = list(dic = TRUE, waic=TRUE), verbose = FALSE)
 
summary(RR_type3mod)


##________________________________________________________________________###

###_______________________________________________________________________####

#model validation, putting best model form III in one INLA command to get CPO and PIT values

#Plotting the calculated mean
spatio_temporal_fitted <- data %>%
  select(District) %>%
  mutate(RR = data$muRR01)


## Type III interaction
form_III <- inla(OB ~ -1 + 
  f(id,
    model = "bym2",
    graph = SL_adj,
    hyper = prec_bym2) +
  f(id_year, 
    model="rw1",
    constr = TRUE,
    scale.model = TRUE,
    hyper=pc_rw) +
  f(id_year1, 
    model="iid",
    constr=TRUE,
    hyper=pc_prior) +
  f(id_year2,
    model="iid", 
    group=ID,
    control.group=list(model="besag",graph=SL_adj, hyper=pc_prior),
    constr=TRUE,
    hyper=pc_iid),
    family = "nbinomial",
    data = data,
    E = expected_cases,
    offset = log(Farm.N),
    control.fixed = list(prec = 1, prec.intercept=1),
    control.predictor = list(link = 1, compute = TRUE),
    control.compute=list(cpo=TRUE),verbose = TRUE)

summary(form_III)

data$mutype3 <- RR_type3$summary.fitted.values[,"mean"]
data$RR<-as.numeric(data$OB/data$expected_cases)
head(cbind(data$mutype3,data$RR))
cor(data$mutype3,data$RR) # 0.98 getting near 1


qplot(data$OB/data$expected_cases,RR_type3$summary.fitted.values$mean) +
  xlab("Observed") +
  ylab("Predicted") +
  geom_smooth(method = "lm", color = "red", linewidth= .5)

a<-form_III$cpo$cpo
hist(a,plot=TRUE) 

#Identify spatial outliers
neg_log_cpo <- -log(form_III$cpo$cpo)
outliers <- which(neg_log_cpo > quantile(neg_log_cpo, 0.95))
outliers

b<-form_III$cpo$pit
hist(b,plot=TRUE)

form_III$cpo$failure
#>0 is bad for failure, didn't get any value above 0,so the model is good,

#compute the adjusted discrete PIT values

pit_raw <- form_III$cpo$pit
cpo_val <- form_III$cpo$cpo
pit_low <- pit_raw - cpo_val

# 1. Mid-PIT (point adjustment)
pit_mid <- pit_raw - 0.5 * cpo_val

# 2. Randomized PIT
set.seed(123)
pit_rand <- pit_raw - runif(length(pit_raw)) * cpo_val

# 3. Czado et al. (2009) non-randomized histogram (10 bins)
bins <- seq(0, 1, length.out = 11)
F_bar <- sapply(bins, function(u) {
  mean(pmin(1, pmax(0, (u - pit_low) / cpo_val)))
})
bin_heights <- diff(F_bar) # Should be ~0.10 per bin under calibration
plot(F_bar)

barplot(bin_heights,
        names.arg = paste0(seq(0, 0.9, 0.1), "-", seq(0.1, 1.0, 0.1)),
        ylim = c(0, max(c(0.15, bin_heights))),
        ylab = "Relative Frequency", xlab = "PIT",
        main = "Non-Randomized Discrete PIT Histogram (Czado et al., 2009)")
abline(h = 0.10, col = "red", lty = 2, lwd = 2) # Ideal uniform height (10% per bin)


#______________________________________________________________________________#

#To get the mean plot of posterior relative risk
#that is already included as column average in the original dataset. 
#data$id<-as.factor(data$id) #maybe not needed.
area$ID<-as.factor(area$ID)
mp<-area%>%
  dplyr::select(ID)%>%
  left_join(data, by="ID")


ggplot(mp, aes(fill=Ave.RR))+
  geom_sf(size=.5)+
  scale_fill_viridis_c(option="C",name = "FMD")+
  theme(strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(10,"points"))


#to get a nice plot for Ave risk only #__________________________________________

R1<-cut(mp$Ave.RR,breaks = c(0.3588887,0.8103706,1.1302319,1.77667262,3.2843127),
        labels = c("0.35-0.81","0.89-1.13","1.18-1.71","1.84-3.18"))

ggplot(mp, aes(fill=R1))+
  geom_sf(size=.5)+
  scale_colour_manual(values=c("#ece2f0","#fff7bc","#fec44f","#225ea8"),aesthetics = "fill")+
  theme_minimal() +
  guides(fill = guide_legend(title ="Relative Risk"))+
  theme(strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.key.size = unit(10,"points"))

##_________________________________________________________________________________##
