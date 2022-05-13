if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC)

PSInSAR_df <- readRDS("data/PSInSAR_df.rds")

PSInSAR_sf <- st_as_sf(PSInSAR_df, coords = c("Longitude", "Latitude"), 
                       agr = "constant", crs=4326)
PSInSAR_sf_n11 <- st_transform(PSInSAR_sf, crs="epsg:32611")

# We will use these variables to split geothermal/non-geothermal
v_mean <- mean(PSInSAR_df$Velocity)
v_sd <- sd(PSInSAR_df$Velocity)
geo_factor <- 1   # We can change this to get more strict on what is or not geothermal

m <- raster::stack("data/brady_ai_stack")
m <- m[[c("Minerals", "Temperature", "Faults", "Slope", 
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
new_ai_data_dfxy <- na.omit(Geothermal_df)
head(new_ai_data_df)


