if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC, nnet)

### Raster functions to save and do max-min normalization
doe_writeRaster <- function(x, filename, format="raster", overwrite=TRUE, bandorder="BSQ"){
  if(tools::file_ext(filename) != "grd") {
    filename <- tools::file_path_sans_ext(filename)
    filename <- paste(filename, ".grd", sep="")
  }
  f1<-raster::writeRaster(x=x, filename=filename, bandorder=bandorder, 
                          format=format, overwrite=overwrite)
  raster::hdr(f1, "ENVI")
  return(f1)
}

# Max-min normalization
normalize_raster <- function(s){
  mnv <- raster::cellStats(s,'min', na.rm=TRUE)
  mxv <- raster::cellStats(s,'max', na.rm=TRUE)
  x <- (s - mnv) / (mxv - mnv)
  return(x)
}





# Contains original PSInSAR (deformation) data
PSInSAR_df <- readRDS("data/PSInSAR_df.rds")
PSInSAR_sf <- st_as_sf(PSInSAR_df, coords = c("Longitude", "Latitude"), 
                       agr = "constant", crs=4326)
PSInSAR_sf_n11 <- st_transform(PSInSAR_sf, crs="epsg:32611")

# Contains PSInSAR Inverse Distance Weighted (IDW) interpolation
idw_3x3 <- raster::raster("data/deformation_3x3.tif")
idw_30x30 <- raster::raster("data/deformation_30x30.tif")
idw_mean <- raster::cellStats(idw_3x3, mean)
idw_sd <- raster::cellStats(idw_3x3, sd)
hymap <- raster::raster("data/HyMapFull_osp")

# We will use these variables to split geothermal/non-geothermal
v_mean <- mean(PSInSAR_df$Velocity)
v_sd <- sd(PSInSAR_df$Velocity)
geo_factor <- 0.5   # We can change this to get more strict on what is or not geothermal

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

rm(PSInSAR_sf, PSInSAR_sf_n11, PSInSAR_df, PSInSAR_df_n11, m, m_sf)
rm(separated_coord)
gc() # Garbage collector, frees up memory

########       Spatiotemporal sampling from here and on...
new_ai_data_xydf <- na.omit(Geothermal_df)
head(new_ai_data_xydf)

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
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 8, repeats = 5)

# The answer to the ultimate question of life, the universe, and everything
set.seed(42)

# Neural network - 7 hidden neuron
nnet_learner <- mlr3::lrn("classif.nnet", size = 7, maxit = 2000) # Neural network, with hidden neurons


print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_st_classif,
                        learner = nnet_learner,
                        store_models = TRUE,
                        resampling = resampling_spcv)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

####################
nnet_best <- nn_sp$learners[[3]]$model
b <- raster::stack("data/brady_minerals_stack")
d <- raster::stack("data/desert_minerals_stack")
prediction_raster_brady <- raster::predict(b, nnet_best)
prediction_raster_desert <- raster::predict(d, nnet_best)
#raster::spplot(prediction_raster)
raster::spplot(prediction_raster_brady<0.5, main="Brady Model 3")
raster::spplot(prediction_raster_desert<0.5, main="Desert Model 3")

# AUC: Area under the curve. 1 is perfect, 0.5 is equal to fully random
pROC::auc(predictor=values(prediction_raster_brady), response=values(b$Geothermal), na.rm=T)
pROC::auc(predictor=values(prediction_raster_desert), response=values(d$Geothermal), na.rm=T)
# ROC: Receiver Operating Characteristic curve
# a diagonal line is equal to a random response (50% chance of picking the right answer)
# a curve reaching the top left is perfect classification
pROC::roc(predictor=values(prediction_raster_brady), response=values(b$Geothermal), na.rm=T, plot=T)
pROC::roc(predictor=values(prediction_raster_desert), response=values(d$Geothermal), na.rm=T, plot=T)

# Create a prediction and reference data frame
pred <- data.frame(reference = values(b$Geothermal), data = values(prediction_raster_brady))
pred_cutoff <- 0.5
pred$reference <- factor(ifelse(pred$reference>pred_cutoff, "Yes", "No"))
pred$data <- factor(ifelse(pred$data<pred_cutoff, "Yes", "No"))

# Calculate statistics based on Confusion Matrix
caret::confusionMatrix(data=pred$data, reference=pred$reference, positive="Yes")
