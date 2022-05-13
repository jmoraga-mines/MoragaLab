if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}

pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC)

# "d:/ThesisLayers/AI_Stacks/FullStacks/desert_full_stack"
# "d:/ThesisLayers/AI_Stacks/FullStacks/brady_full_stack"
# "d:/ThesisLayers/AI_Stacks/MineralStacks/desert_minerals_stack"
# "d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack"
m <- raster::stack("d:/ThesisLayers/AI_Stacks/MineralStacks/brady_minerals_stack")




# Creates a logistic regression learner
logistic_regression <- mlr3::lrn("classif.log_reg", id = "binary") # Logistic regression

# Non-spatiotemporal classification
new_ai_data_df <-  as.data.frame(m, xy=FALSE)
task_classif <- mlr3::as_task_classif(new_ai_data_df,
                                      target = "Geothermal",
                                      positive = "Yes",
                                      id = "geot_cv")
# Creates repeated cv resampling. 5 folds, repeated 2 times
resampling_cv <- mlr3::rsmp('repeated_cv', folds = 5, repeats = 2)
set.seed(42)
rr_lrcv <- mlr3::resample(task = task_classif,
                          learner = logistic_regression,
                          store_models = TRUE,
                          resampling = resampling_cv)  

print(rr_lrcv)
# rr_lrcv$score()
rr_lrcv$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
rr_lrcv$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(rr_lrcv)

########       Spatioemporal sampling from here and on...



new_ai_data_xydf <- as.data.frame(m, xy=TRUE) # We will use coordinates for spatial sampling
new_ai_data_xydf <- na.omit(new_ai_data_xydf)
# new_ai_data_xydf$Geothermal <- as.factor(new_ai_data_xydf$Geothermal)

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
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)

set.seed(42)
rr_lrspcv <- mlr3::resample(task = task_st_classif,
                            learner = logistic_regression,
                            store_models = TRUE,
                            resampling = resampling_spcv)
print(rr_lrspcv)
# rr_lrspcv$score()
rr_lrspcv$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
rr_lrspcv$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(rr_lrspcv)

#### Non-linear learners
# Neural network - 1 hidden neuron
nnet_learner <- mlr3::lrn("classif.nnet", size = 1) # Neural network, with hidden neurons

print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_st_classif,
                        learner = nnet_learner,
                        #                        store_models = TRUE,
                        resampling = resampling_spcv)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))


# Neural network - 7 hidden neuron
nnet_learner <- mlr3::lrn("classif.nnet", size = 7, maxit = 1000) # Neural network, with hidden neurons

print(nnet_learner)
set.seed(42)
nn_sp <- mlr3::resample(task = task_st_classif,
                        learner = nnet_learner,
                        # store_models = TRUE,
                        resampling = resampling_spcv)

nn_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
nn_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(nn_sp)


# Random forest
pacman::p_load(ranger)
rf_learner <- mlr3::lrn("classif.ranger", num.trees=500, num.threads=3) # Random forest learner

print(rf_learner)
set.seed(42)
rf_sp <- mlr3::resample(task = task_st_classif,
                        learner = rf_learner,
                        #                        store_models = TRUE,
                        resampling = resampling_spcv)

rf_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
rf_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(rf_sp)


############# START ### Once again, set up SP CV
if (!require("pacman")){
  install.packages("pacman")
  require("pacman")
}
pacman::p_load(caret, sf, terra, tidyverse, mlr3verse, pROC, mlr3spatiotempcv)
dir1 = "/store03/thesis/ai_stacks/MineralStacks"
dir2 = "d:/ThesisLayers/AI_Stacks/MineralStacks"
if(dir.exists(dir1)) {
  dir_base = dir1
} else if(dir.exists(dir2)){
  dir_base = dir2
} else {
  print("Error: data directories do not exist")
}
input_stack_file = file.path(dir_base, "brady_minerals_stack.gri")
m <- raster::stack(input_stack_file)
########       Spatioemporal sampling from here and on...
new_ai_data_xydf <- as.data.frame(m, xy=TRUE) # We will use coordinates for spatial sampling
new_ai_data_xydf <- na.omit(new_ai_data_xydf)
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
resampling_spcv <- mlr3::rsmp('repeated_spcv_coords', folds = 5, repeats = 2)
############# END  ### Once again, set up SP CV


# Support vector machine (SVM)
pacman::p_load(e1071)
svm_learner <- mlr3::lrn("classif.svm", 
                         cachesize = 1024, 
                         tolerance = 1, 
                         type = 'C-classification',
                         # gamma=0.2, # For kernels: polynomial, radial, sigmoid
                         cost=100,
                         kernel='linear') # SVM

print(svm_learner)
set.seed(42)
svm_sp <- mlr3::resample(task = task_st_classif,
                        learner = svm_learner,
                        #                        store_models = TRUE,
                        resampling = resampling_spcv)

svm_sp$score(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))
svm_sp$aggregate(measures = c(msr("classif.acc"), msr("classif.bacc"), msr("classif.fbeta")))

rm(svm_sp)

