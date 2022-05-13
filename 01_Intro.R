if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC, nnet)

# Contains original PSInSAR (deformation) data
PSInSAR_df <- readRDS("data/PSInSAR_df.rds")
PSInSAR_sf <- st_as_sf(PSInSAR_df, coords = c("Longitude", "Latitude"), 
                       agr = "constant", crs=4326)
PSInSAR_sf_n11 <- st_transform(PSInSAR_sf, crs="epsg:32611")

# We will use these variables to split geothermal/non-geothermal
v_mean <- mean(PSInSAR_df$Velocity)
v_sd <- sd(PSInSAR_df$Velocity)
geo_factor <- 1   # We can change this to get more strict on what is or not geothermal

m <- raster::stack("data/brady_ai_stack")
m <- m[[c("Temperature", "Faults", "Slope", 
          "Chalcedony", "Kaolinite", "Gypsum", "Hematite")]]

# Obtain a sf object with the raster data at relevant points
m_sf <- raster::extract(m, PSInSAR_sf_n11) # extract values at each point
PSInSAR_sf_n11 <- cbind(PSInSAR_sf_n11, m_sf)
PSInSAR_sf_n11 <- na.omit(PSInSAR_sf_n11)


separated_coord <- PSInSAR_sf_n11 %>%
  mutate(x = unlist(map(PSInSAR_sf_n11$geometry,1)),
         y = unlist(map(PSInSAR_sf_n11$geometry,2)))
#PSInSAR_df_n11 <- as.data.frame(separated_coord[c("x", "y", "Velocity")])
#PSInSAR_df_n11 <- subset(PSInSAR_df_n11, select=c("x", "y", "Velocity"))
PSInSAR_df_n11 <- subset(as.data.frame(separated_coord), select=-c(Coherence, geometry))

head(PSInSAR_df_n11)

PSInSAR_df_n11$Geothermal <- ifelse(PSInSAR_df_n11$Velocity<(v_mean-geo_factor*v_sd), 
                                    1, 
                                    ifelse(PSInSAR_df_n11$Velocity>(v_mean+geo_factor*v_sd), 
                                           0, NA)
                                    )
Geothermal_df <- subset(PSInSAR_df_n11, select=-Velocity)

rm(df, PSInSAR_sf, PSInSAR_sf_n11, PSInSAR_df, PSInSAR_df_n11, m, m_sf)
gc() # Garbage collector, frees up memory

########       Spatiotemporal sampling from here and on...
new_ai_data_xydf <- na.omit(Geothermal_df)
head(new_ai_data_df)

# Convert to factor
new_ai_data_xydf$Geothermal <- factor(ifelse(new_ai_data_xydf$Geothermal>=0.5, "Yes", "No"),
                                      levels = c("No", "Yes"))

# Creates a new taks to classify using spatio-temporal cross-validation
task_st_classif <- mlr3spatiotempcv::TaskClassifST$new(new_ai_data_xydf, 
                                                       target = "Geothermal",
                                                       positive = "Yes",
                                                       id = "geot_stcv",
                                                       extra_args = list(coords_as_features = FALSE,
                                                                         coordinate_names=c("x", "y"))
                                                       )
# Creates spatiotemporal resampling. 5 folds, repeated 2 times
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 10, repeats = 4)

# The answer to the ultimate question of life, the universe, and everything
set.seed(42)

# Neural network - 7 hidden neuron
nnet_learner <- mlr3::lrn("classif.nnet", size = 7, maxit = 1000) # Neural network, with hidden neurons


print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_st_classif,
                        learner = nnet_learner,
                        store_models = TRUE,
                        resampling = resampling_spcv)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

####################
nnet_14 <- nn_sp$learners[[14]]$model
m <- raster::stack("data/brady_ai_stack")
prediction_raster <- raster::predict(m, nnet_14)
raster::spplot(prediction_raster)
raster::spplot(prediction_raster<0.5)

