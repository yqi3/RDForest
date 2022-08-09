#### This file contains useful functions for the modeling part of the toy simulation section
set.seed(1234)

#### ------ GBRT is a function that gives an autotuned GBRT ------ ####
GBRT <- function(df){
  # speed up computation
  all_cores <- parallel::detectCores(logical = FALSE)
  doParallel::registerDoParallel(cores = all_cores)
  
  # Model to be tuned
  xgb_model <- boost_tree(
    mode = "regression",
    trees = tune(),
    tree_depth = tune(), min_n = tune(), 
    loss_reduction = tune(),      ## first three: model complexity
    sample_size = tune(), mtry = ncol(df)-1,         ## randomness
    learn_rate = tune(),                         ## step size
  ) %>% 
    set_engine("xgboost", objective = "reg:squarederror")
  
  # Grid of parameters to search
  xgb_grid <- grid_latin_hypercube(
    trees(),
    tree_depth(),
    loss_reduction(),
    min_n(),
    sample_size = sample_prop(),
    learn_rate(),
    size = 30
  )
  
  # Formula for the model
  xgb_wf <- workflow() %>%
    add_formula(Y ~ .) %>%
    add_model(xgb_model)
  
  # Splitting for cross validation
  myfolds <- vfold_cv(df, v=5) # 5-fold
  
  # Tunning
  xgb_res <- tune_grid(
    xgb_wf,
    resamples = myfolds,
    grid = xgb_grid,
    control = control_grid(save_pred = TRUE)
  )
  
  best_auc <- select_best(xgb_res, "rmse") # Best tuned parameters
  final_xgb <- finalize_workflow(xgb_wf, best_auc)
  
  return(final_xgb)
}