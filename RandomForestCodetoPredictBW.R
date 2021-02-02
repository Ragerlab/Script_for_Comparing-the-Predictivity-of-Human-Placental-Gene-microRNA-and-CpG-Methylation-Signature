#################################################################################################
#################################################################################################
#### Using random forest modeling to identify gene signatures predictive of outcomes
#### Set-up to upload statistical results files from DESeq2 analyses
####
####
#### This is an example using results from a birth weight analysis based on the ELGAN cohort data
#### Birth weight was analyzed as a continuous variable here as the dependent variable, and miRNA expression signatures as predictor variables alongside covariate data
####
####
#### Drafted by Julia Rager, Caroline Ring, Jeliyah Clark, and Vennela Avula
#### Last Update: February 1, 2021
#################################################################################################



#################################################################################################
#################################################################################################
#### Activating appropriate packages
#################################################################################################
#################################################################################################

# Activating appropriate libraries
library(data.table)
library(randomForest)
library(dplyr)
library(tibble)
library(ggplot2)
library(forestFloor)


#################################################################################################
#################################################################################################
#### Set working directory, output folder path, etc
#################################################################################################
#################################################################################################

# Set working directory
setwd("/Users/vavula/IEHS Dropbox/Rager Lab/Vennela_Avula/ELGAN_PredTox/5_Github/Input")


# Create an output folder
Output_Folder <- ("/Users/vavula/IEHS Dropbox/Rager Lab/Vennela_Avula/ELGAN_PredTox/5_Github/Output")


# Setting title of analyses
AnalysisTitle <- "BirthWeight_miRNAs_SigGenes"


#################################################################################################
#################################################################################################
#### Uploading data that are exported from our DESeq2 pipeline
#################################################################################################
#################################################################################################

# Upload subject information / covariate file
subject_info <- read.csv("SubjectInfo_376.csv")


# Upload the normalized expression count data, per subject
counts <- read.csv("VSD_Counts_376subjects.csv", header=TRUE)
colnames(counts)

# Upload the statistical results from DESeq2
stat_results <- read.csv("StatResults_SVA.csv")


#################################################################################################
#################################################################################################
#### Filtering genes for significance
#################################################################################################
#################################################################################################


# Filter the genes within our results based on the statistical thresholds of our choosing (padj<=0.1)
sig_genes <- stat_results[which(stat_results$padj <= 0.1),]


# Pull this list of genes
genes_to_use <- sig_genes[,1] # this is a factor
rm(sig_genes) # removes transient dataframe


# Create a subset data set containing count data for the significant genes we want to evaluate
counts_sig <- counts %>% 
  filter(X %in% genes_to_use) #CLR: need to do library(dplyr) first to get pipe 
rownames(counts_sig) <- NULL
rownames(counts_sig) <- counts_sig[,1]
counts_sig <- counts_sig[,-1]


# Transpose the sig count data, and correct IDs to exclude the "X" that were imported because they were numbers
counts_sig <- as.data.frame(t(counts_sig))
rownames(counts_sig) <- sub("X", "", rownames(counts_sig))


# Make rownames into first column for easier downstream merging
counts_sig <- counts_sig %>% rownames_to_column("ELGANID") #CLR: neeed to do library(tibble) to get this command
class(counts_sig)


###########################################################################################################################
###########################################################################################################################
#### Organizing master dataframes with data needed for RF modeling, by combining count data with appropriate covariate data
###########################################################################################################################
###########################################################################################################################

# Define variables to be evaluated
# This code provides an example of using continuous outcome data (cont_var)
cont_var <- "bw"

# Define covariates that were included in the DESeq2 statistical model, that also need to be included here
cov_list <- c("sex", "mult", "gadays", "smoke_exp")

#Ensure other variables are dropped from subject_info
subject_info <- subject_info %>% 
  select(c("ELGANID",all_of(cov_list),all_of(cont_var)))
head(subject_info)

# Checking classes of covariate data
lapply(subject_info[cov_list], class)

#Ensure both ELGANIDs recognized as numeric
subject_info$ELGANID <- as.numeric(subject_info$ELGANID)
counts_sig$ELGANID <- as.numeric(counts_sig$ELGANID)


# Merging data, focusing on the outcome and model components of interest (e.g., molecular predictors and/or covs, including the outcome covariate)
# Labeled this master dataframe as 'pred_outcome_vars' because it includes all predictors (molecules & covariates) and outcome (phenotype trying to predict)
pred_outcome_vars <- left_join(counts_sig,subject_info, by = "ELGANID")


# Viewing data
head(pred_outcome_vars)
pred_outcome_vars[, cont_var]  # recognizes the generic continuous variable that we designated above
pred_outcome_vars[ , "bw"] # QCing that this is the same as the above step 
all.equal(pred_outcome_vars[, cont_var], pred_outcome_vars[ , "bw"]) # additional QC that this is the same


# Need to make sure cont_var is interpreted as a numeric variable
pred_outcome_vars[, cont_var] <- as.numeric(pred_outcome_vars[,cont_var])
class(pred_outcome_vars[, cont_var])

# Make sure ELGANID is rowname
pred_outcome_vars <- pred_outcome_vars %>%
                column_to_rownames("ELGANID")

# Viewing data in master dataframe
head(pred_outcome_vars)
pred_outcome_vars[1:5,1:5]



#################################################################################################
#################################################################################################
#### Splitting data into training & test sets
#################################################################################################
#################################################################################################
# We will use 2/3 of the data to first build random forest models (training set)
# Then use the remaining 1/3 of the data to test the built models (test set)

# Randomly pulling 66% of the data for the training set, and remaining for test
set.seed(151)
to_train <- sample(1:nrow(pred_outcome_vars), round(nrow(pred_outcome_vars)*0.66, 0), replace = F)
train_vars <- pred_outcome_vars[to_train,]
test_vars <- pred_outcome_vars[-to_train,]


#################################################################################################
#################################################################################################
#### Creating a RF model for the continuous outcome variable
#################################################################################################
#################################################################################################

# Defining the number of trees in each random forest
ntree <- 10001

# Defining the number of permutations used to estimate null distribution for variable importance significance testing. More is more accurate, but takes longer.
n.perms <- 1000

# Create RF model from the training data set, and then testing against the test data set, regressing against the continuous outcome measure
set.seed(151)
RF <- randomForest( x=train_vars[,1:(ncol(train_vars)-1)] , 
                      y=train_vars[ , ncol(train_vars)] , 
                      xtest=test_vars[, 1:(ncol(test_vars)-1)] ,
                      ytest=test_vars[, ncol(test_vars) ], 
                      ntree=ntree,
                      importance=TRUE,
                      proximities=TRUE,
                      keep.inbag = T)

# Viewing results from the built RF model
RF

# Let's pull some of these values to later organize into a master dataframe, with data to report in our overall study findings
# First need to find the names of items in the resulting RF object:
names(RF)
# Note that the definitions of these output items by name are in ?randomForest -- look at "Value" section


# To extract items in RF object and assign them to their own variables, use the code below
# Note that we have to grab the values from the last element, which corresponds to the final tree created (which includes all trees)
MeanOfSqResiduals_RF_TrainSet <- RF$mse[ntree]
PercentVar_RF_TrainSet <- RF$rsq[ntree]
MeanOfSqResiduals_RF_TestSet <- RF$test$mse[ntree]
PercentVar_RF_TestSet <- RF$test$rsq[ntree]


# Calculate RMSE from the built RF model by extracting the last element of RF$mse, which corresponds to the final tree created (which includes all trees), and take its square root
RMSE_RF_TrainSet <- sqrt(RF$mse[ntree])
RMSE_RF_TestSet <- sqrt(RF$test$mse[ntree])


# As a QC, check that these values are the same as mean of sq residuals from the above
all.equal(sqrt(MeanOfSqResiduals_RF_TrainSet), RMSE_RF_TrainSet)
all.equal(sqrt(MeanOfSqResiduals_RF_TestSet), RMSE_RF_TestSet)



# Viewing variable importance plot from the built RF model
varImpPlot(RF)
dev.copy(png,paste0(Output_Folder,"/", "Variable_Imp_Plot.png"))
dev.off()


# Compile variable importance measures into a df
var_imp <- as.data.frame(importance(RF))
#Sort df by %IncMSE
var_imp <- var_imp[order(-var_imp$`%IncMSE`),]

#Export variable importance file
write.csv(var_imp, paste0(Output_Folder,"/", "Variable_Imp_Rankings.csv"))



#################################################################################################
#################################################################################################
#### Comparing the predicted vs observed results within the test data set
#################################################################################################
#################################################################################################

# First pull the values that were predicted across the test dataset
Pred_TestSet <- as.data.frame(RF$test$predicted)
colnames(Pred_TestSet) <- "Predicted_Values"
# Note that you can also pull values that were predicted across the training dataset using: RF$predicted


# Then pull the values that were the originally observations across the test dataset
Obs_TestSet <- as.data.frame(test_vars[,ncol(test_vars)])
rownames(Obs_TestSet) <- rownames(test_vars)
colnames(Obs_TestSet) <- "Observed_Values"


# Next merge the predicted & observed values into one dataframe
pred_table <- merge(Pred_TestSet, Obs_TestSet, by="row.names")
row.names(pred_table) <- pred_table$Row.names
pred_table$Row.names <- NULL


# Create a column that calculates the percentage difference between the observed values and predicted values
pred_table$pct_diff <- (pred_table$Predicted_Values - pred_table$Observed_Values)/ pred_table$Observed_Values * 100 


# Export table
write.csv(pred_table, file = paste0(Output_Folder,"/","PredvsObs_PctDiff_TestSet.csv"))


########
# Pulling model fit parameters from the linear regression fit between the predicted vs observed (e.g., R^2 values)
########


model_lm <- lm(Predicted_Values ~ Observed_Values, data = pred_table) 
(model_summary <- summary(model_lm))# if you surround the assignment with parentheses, the output will still print


# Saving these model fit values to combine in our master results dataframe below
names(model_summary)
# Pulling R^2 of the predicted vs observed linear regression
RSq_PredvsObs_TestSet <- model_summary$r.squared
# Pulling the root mean squared error (RMSE) of the predicted vs observed linear regression
res <- model_summary$residuals #or residuals(model_lm), they are the same   -- both calculate the residual
RMSE_PredvsObs_TestSet <- sqrt(mean(res^2))


########
# View distribution of differences between predicted vs.observed values
########
hist(pred_table$pct_diff)
qplot(
  x = Observed_Values,
  y = pct_diff,
  data = pred_table,
  xlab = "Observed Birth Weight (g)", ylab = "% Difference in Prediction",
  main = ""
)
dev.copy(png,paste0(Output_Folder,"/", "PredvsObs_PctDiff_TestSet.png"))
dev.off()


########
# Plot predicted vs observed values
########
# First derive the axis limits to force axes to be "square" -- same limits on x and y
# Derive these by calculating the overall min of predicted & observed, and overall max of predicted & observed, and concatenating them into a two-element vector (min, max)
# These will later be used in a call to scale_x_continous() and scale_y_continuous()
axis_lims <- c(min(c(pred_table$Observed_Values, pred_table$Predicted_Values)),
               max(c(pred_table$Observed_Values, pred_table$Predicted_Values)))

# Drafting plots               
ggplot(data = pred_table, aes(x=Observed_Values, y=Predicted_Values)) + 
  labs(x="Observed Values", y = "RF-Predicted Values") +
  geom_point() +
  geom_smooth(method="lm") +
  geom_abline(slope =1, intercept = 0) +
  scale_x_continuous(limits = axis_lims) +
  scale_y_continuous(limits = axis_lims)

dev.copy(png,paste0(Output_Folder,"/", "PredvsObs_TestSet_Plot.png"))
dev.off()
# Here, the black line represents the identity line, with slope = 1, intercept = 0
# The blue line represents the linear regression results fitted against the predicted vs. observed values (i.e., shows how far off we are)


#################################################################################################
#################################################################################################
#### Creating the same RF model, but with additional noise variables
#################################################################################################
#################################################################################################

# Creating a new df to add noise variables to
train_vars_noise <- train_vars

set.seed(42) #to make sure noise is reproducible
# Add random noise predictors as an additional method to evaluate model performance, first focusing on molecular predictors
train_vars_noise$noise1 <- sample(train_vars_noise[,3], replace=TRUE)  # Adding a column that contains randomly shuffled values from one of the molecules; sampling with replacement
train_vars_noise$noise2 <- sample(train_vars_noise[,4], replace=TRUE)
train_vars_noise$noise3 <- sample(train_vars_noise[,5], replace=TRUE)
train_vars_noise$noise4 <- sample(train_vars_noise[,6], replace=TRUE)
train_vars_noise$noise5 <- sample(train_vars_noise[,7], replace=TRUE)

# Add random noise predictors for covariates
# First need to name each covariate to pull, to make code more reproducible across datasets
length(cov_list)
# Here, the length is 4 - so we make 4 covariate variables
cov1 <- cov_list[1]
cov2 <- cov_list[2]
cov3 <- cov_list[3]
cov4 <- cov_list[4]

# Then add these to the training dataframe with noise
train_vars_noise$noise_cov1 <- sample(train_vars_noise[,cov1], replace = TRUE)
train_vars_noise$noise_cov2 <- sample(train_vars_noise[,cov2], replace = TRUE)
train_vars_noise$noise_cov3 <- sample(train_vars_noise[,cov3], replace = TRUE)
train_vars_noise$noise_cov4 <- sample(train_vars_noise[,cov4], replace = TRUE)
head(train_vars_noise)


# Re-establishing the training set to include the random noise variables before the outcome variable
temp_num <- 5 + length(cov_list) + 1   # 5 variables from molecular noise, covariable noise variables, plus birth weight variable

# Use this number in the below function, after ncol(train_vars_noise)-10
# also make sure to manually test all the noise variables below
train_vars_noise <- train_vars_noise %>%
  select(colnames(train_vars_noise)[1:(ncol(train_vars_noise)-10)],
         "noise1", "noise2", "noise3", "noise4", "noise5",
         "noise_cov1", "noise_cov2", "noise_cov3", "noise_cov4",
         cont_var)

head(train_vars_noise)



###############
#Re-running RF model with the noise variable added
############### 

set.seed(151)
RF_noise <- randomForest( x=train_vars_noise[,(1:(ncol(train_vars_noise)-1))], 
                          y=train_vars_noise[ , ncol(train_vars_noise)] , 
                          ntree=ntree,
                          importance=TRUE,
                          proximities=TRUE,
                          keep.inbag = T)


# Viewing variable importance plot from the built RF model, including the random noise variable
varImpPlot(RF_noise)
dev.copy(png,paste0(Output_Folder,"/", "Variable_Imp_Rankings_noise_plot.png"))
dev.off()

# Compile variable importance measures into a df
var_imp_noise <- as.data.frame(importance(RF_noise))
#Sort df by %IncMSE
var_imp_noise <- var_imp_noise[order(-var_imp_noise$`%IncMSE`),]

#Export variable importance file, with noise variables
write.csv(var_imp_noise, paste0(Output_Folder,"/", "Variable_Imp_Rankings_noise.csv"))

# Viewing the rank of the noise variables, out of the total number of predictor variables
which(rownames(var_imp_noise)=="noise1") / nrow(var_imp_noise)
# Do this in a loop, looping over names of noise variables
sapply(c(paste0("noise", 1:5),
         "noise_cov1",
         "noise_cov2",
         "noise_cov3",
         "noise_cov4"),
       function(noisevar) which(rownames(var_imp_noise)==noisevar) / nrow(var_imp_noise),
       USE.NAMES = TRUE)

# Not normalized as %
noise_row_num <- sapply(c(paste0("noise", 1:5),
         "noise_cov1",
         "noise_cov2",
         "noise_cov3",
         "noise_cov4"),
       function(noisevar) which(rownames(var_imp_noise)==noisevar),
       USE.NAMES = TRUE)

# Find the highest ranking noise variable
Noise_Highest_Rank <- noise_row_num[which.min(noise_row_num)]

######################################################################################################
######################################################################################################
#### Organizing a 'master results' dataframe, to export all of our pertinent model performance results
######################################################################################################
######################################################################################################

# Pulling number of variable predictors
No_of_Predictors <- ncol(test_vars) - 1 
No_of_Molecular_Predictors <- ncol(counts_sig) - 1


# Creating the summary of all the RF results 
MasterResults_df <- data.frame(No_of_Predictors, 
                               No_of_Molecular_Predictors, 
                               RMSE_RF_TestSet,
                               RSq_PredvsObs_TestSet,
                               RMSE_PredvsObs_TestSet,
                               PercentVar_RF_TestSet,
                               Noise_Highest_Rank,
                               RMSE_RF_TrainSet,
                               MeanOfSqResiduals_RF_TrainSet,
                               PercentVar_RF_TrainSet,
                               MeanOfSqResiduals_RF_TestSet
                               )

row.names(MasterResults_df) <- AnalysisTitle

write.csv(MasterResults_df, paste0(Output_Folder,"/", "RF_ResultsSummary_", AnalysisTitle, ".csv"))
