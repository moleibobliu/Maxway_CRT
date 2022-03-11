source('function.R')

#### Required packages

library(MASS)
library(glmnet)
library(CompQuadForm)

#### We illustrate the use of our approach in three simulation settings.
#### See Section 5 of the Maxway CRT paper for more details.


#### Setting (I): linear gaussian model for both X and Y

#### Data generation

# Whether to include interaction between X and some variables in Z in this model

interact = T 
interact = F

# Generate semi-supervised data (a few data with labels and a large amount of unlabeled data). 

gen_data <- Gen_data(n = 250, N = 1000, p = 500, r = 0.5, magn_A = 0.15, magn_x = 0.2, 
                     magn_y = 0.25, magn_share = 0.3, s_x = 25, s_y = 25, s = 5,
                     model_x = 'linear', model_y = 'linear', sign = 'random',
                     interact = interact, seed1 = 1000, seed2 = 100)

y_label <- gen_data$y
Z_label <- gen_data$Z_label
A_label <- gen_data$A_label
Z_add <- gen_data$Z_add
A_add <- gen_data$A_add

# Additional generate some labeled data to fit the hold-out training version of the Maxway CRT.
 

gen_data_label <- Gen_data(n = 250, 10, p = 500, r = 0.5, magn_A = 0.15, magn_x = 0.2, 
                           magn_y = 0.25, magn_share = 0.3, s_x = 25, s_y = 25, s = 5,
                           model_x = 'linear', model_y = 'linear', sign = 'random',
                           interact = interact, seed1 = 1000, seed2 = 500)

y_label_add <- gen_data_label$y
Z_label_add <- gen_data_label$Z_label


# Run the model-X and the Maxway CRT:


Maxway_CRT_result <- Maxway_CRT(A_label, Z_label, y_label, A_add = A_add, Z_add = Z_add, 
                                Z_label_add = Z_label_add, y_label_add = y_label_add,
                                model_x = 'Gaussian_lasso', model_y = 'Gaussian_lasso', 
                                M = 1000) # M: the number of randomization samples.

# P-value of the model-X CRT with d0 statistics

Maxway_CRT_result$MX_d0CRT_pvl

# P-value of the model-X CRT with dI statistics

Maxway_CRT_result$MX_dICRT_pvl

# P-value of the Maxway CRT (the in-sample version) with d0 statistics

Maxway_CRT_result$Maxway_d0CRT_pvl


# P-value of the Maxway CRT (the in-sample version) with dI statistics

Maxway_CRT_result$Maxway_dICRT_pvl


# P-value of the Maxway CRT (the out-of-sample/holdout version) with d0 statistics

Maxway_CRT_result$Maxway_out_d0CRT_pvl


# P-value of the Maxway CRT (the out-of-sample/holdout version) with dI statistics

Maxway_CRT_result$Maxway_out_dICRT_pvl





#### Setting (II): linear gaussian model for Y and logistic model for binary X

#### Data generation

# Whether to include interaction between X and some variables in Z in this model

interact = T 
interact = F

# Generate semi-supervised data (a few data with labels and a large amount of unlabeled data). 

gen_data <- Gen_data(n = 250, N = 1000, p = 500, r = 0.5, magn_A = 0.1, magn_x = 0.15, 
                     magn_y = 0.2, magn_share = 0.3, s_x = 25, s_y = 25, s = 5,
                     model_x = 'linear_binary', model_y = 'linear', sign = 'random',
                     interact = interact, seed1 = 1000, seed2 = 100)

y_label <- gen_data$y
Z_label <- gen_data$Z_label
A_label <- gen_data$A_label
Z_add <- gen_data$Z_add
A_add <- gen_data$A_add

# Additional generate some labeled data to fit the hold-out training version of the Maxway CRT.

gen_data_label <- Gen_data(n = 250, 10, p = 500, r = 0.5, magn_A = 0.1, magn_x = 0.15, 
                           magn_y = 0.2, magn_share = 0.3, s_x = 25, s_y = 25, s = 5,
                           model_x = 'linear_binary', model_y = 'linear', sign = 'random',
                           interact = interact, seed1 = 1000, seed2 = 500)

y_label_add <- gen_data_label$y
Z_label_add <- gen_data_label$Z_label


# Run the model-X and the Maxway CRT:


Maxway_CRT_result <- Maxway_CRT(A_label, Z_label, y_label, A_add = A_add, Z_add = Z_add, 
                                Z_label_add = Z_label_add, y_label_add = y_label_add,
                                model_x = 'Binomial_lasso', model_y = 'Gaussian_lasso', 
                                M = 1000) # M: the number of randomization samples.

# P-value of the model-X CRT with d0 statistics

Maxway_CRT_result$MX_d0CRT_pvl

# P-value of the model-X CRT with dI statistics

Maxway_CRT_result$MX_dICRT_pvl

# P-value of the Maxway CRT (the in-sample version) with d0 statistics

Maxway_CRT_result$Maxway_d0CRT_pvl


# P-value of the Maxway CRT (the in-sample version) with dI statistics

Maxway_CRT_result$Maxway_dICRT_pvl


# P-value of the Maxway CRT (the out-of-sample/holdout version) with d0 statistics

Maxway_CRT_result$Maxway_out_d0CRT_pvl


# P-value of the Maxway CRT (the out-of-sample/holdout version) with dI statistics

Maxway_CRT_result$Maxway_out_dICRT_pvl






#### Setting (III): Non-linear (gaussian) model for X and Y

#### Data generation


# Generate semi-supervised data (a few data with labels and a large amount of unlabeled data). 

gen_data <- Gen_data(n = 250, N = 1000, p = 40, r = 0.2, magn_A = 0.1, magn_x = 0.1, 
                     magn_y = 0.1, model_x = 'non_linear', model_y = 'non_linear', seed2 = 100)

y_label <- gen_data$y
Z_label <- gen_data$Z_label
A_label <- gen_data$A_label
Z_add <- gen_data$Z_add
A_add <- gen_data$A_add

# Additional generate some labeled data to fit the hold-out training version of the Maxway CRT.

gen_data_label <- Gen_data(n = 250, 10, p = 40, r = 0.2, magn_A = 0.1, magn_x = 0.1, 
                           magn_y = 0.1, model_x = 'non_linear', model_y ='non_linear', seed2 = 500)

y_label_add <- gen_data_label$y
Z_label_add <- gen_data_label$Z_label


# Run the model-X and the Maxway CRT:

library(randomForest)

Maxway_CRT_result <- Maxway_CRT(A_label, Z_label, y_label, A_add = A_add, Z_add = Z_add, 
                                Z_label_add = Z_label_add, y_label_add = y_label_add,
                                model_x = 'Gaussian_RF', model_y = 'Gaussian_RF', 
                                RF.num.trees = c(200, 200, 200, 20), M = 500)

# P-value of the model-X CRT

Maxway_CRT_result$MX_CRT_pvl


# P-value of the Maxway CRT (the in-sample version) 

Maxway_CRT_result$Maxway_CRT_pvl


# P-value of the Maxway CRT (the out-of-sample/holdout version)

Maxway_CRT_result$Maxway_out_CRT_pvl




