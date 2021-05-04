install.packages("TH.data")
library("TH.data")
data("bodyfat")
df = bodyfat

# df <- read.csv("bodyfatwomen.csv", header = TRUE)

model <- lm(DEXfat ~ ., data=df)
summary(model)

#Residual analysis #############################################################
#1) print residuals
Rres = rstudent(model)
res = rstandard(model)
#Residuals 27 & 41 have high value -> outliers?
#They're also higher as Rstudent vs Rstandard -> influential point?

#2) plot the residuals
plot(model, which=2)
hist(res, breaks=30, main="RStandard residuals histogram")
plot(res)
text(res, labels=1:71, pos=4)
#Very light tail on the right (positive skew), and a potential outlier at the very left

#3) residual vs fitted value plot
plot(model, which=1)
#Kinda good, but potential problem at right end (higher fitted values)

#4) residuals vs each regressor
par(mfrow=c(2,3))
plot(df$age, Rres, 
       ylab="RStudent Residuals", xlab="Age", 
       main="Full model: residuals vs regressor") 
abline(0, 0)
plot(df$waistcirc, Rres, 
     ylab="RStudent Residuals", xlab="Waist circumference", 
     main="Full model: residuals vs regressor") 
abline(0, 0)
plot(df$waistcirc, Rres, 
     ylab="RStudent Residuals", xlab="Waist circumference", 
     main="Full model: residuals vs regressor") 
abline(0, 0)
plot(df$hipcirc, Rres, 
     ylab="RStudent Residuals", xlab="Hip circumference", 
     main="Full model: residuals vs regressor") 
abline(0, 0)
plot(df$elbowbreadth, Rres, 
     ylab="RStudent Residuals", xlab="Elbow breadth", 
     main="Full model: residuals vs regressor") 
abline(0, 0)
plot(df$kneebreadth, Rres, 
     ylab="RStudent Residuals", xlab="Knee breadth", 
     main="Full model: residuals vs regressor") 
abline(0, 0)
#Some funneling going on? Especially for age -> uneven variance?
par(mfrow=c(1,1))

#5) Partial residual vs regressor plots 
crp(model)
#age, elbowbreadth, anthro3c uninformative (given others)?
#kneebreadth is weird
#rest are good

#6) Cook's distance
plot(model, which=4)
influencePlot(model)
#Residual 25 may be an influential (high leverage) outlier.
#Residuals 27 & 41 & 48 may be outliers.

#7) CovRatio
covr = covratio(model)
plot(covr, xlab="Residual Index")
text(covr, labels=1:71, pos=4)
#A few points seem to degrade precision

#8) DFFITS
dffits = dffits(model)
plot(dffits)
text(dffits, labels=1:71, pos=4)

#Conclusion: remove 27, 41 and 48 
dfr <- df[-c(27,41,48),]
#Transformations ###############################################################

modelr <- lm(DEXfat ~ ., data=dfr)
summary(modelr)
#R^2 increased from 0.9117 to 0.9434

resr = rstandard(modelr)
plot(resr) #seems to trend downwards?
crp(modelr) #use crPlots(.., layout=c(1,1)) to choose one graph at a time (for the pdf)
plot(modelr, which=2) #not very linear :(
hist(resr, breaks=30)

shapiro.test(resr) #it says it's normal

#try to fix age:
dft <- data.frame(dfr)
dft$age <- sapply(dfr$age, function(x) x**3)
modelt <- lm(DEXfat ~ ., data=dft)
rest = rstandard(modelr)
plot(rest) #A-OK
crPlots(modelt)
plot(modelt, which=2)

#Conclusion: transformations upon the regressors don't fix the partial residual plots
#The model behaves well without any transformations
###############################################################################

# Multicolinearity ############################################################
# Full code in a separate file, the conclusion was that there is severe multi-
# colinearity and there's nothing we can do about it - apart from using ridge
# regression.

# Variable selection ##########################################################
# Stepwise regression
m <- step(
        lm(DEXfat ~ ., data=dfr),
        direction="both"
)
# This chooses the model: DEXfat ~ waistcirc + hipcirc + kneebreadth + anthro3b + anthro4
# Statistics: AIC 120.25; Adjusted R-squared:  0.9444064

# Backward selection
m <- step(
        lm(DEXfat ~ ., data=dfr),
        direction="backward"
)
# This produces the same result

# For evaluating whether using a subset of variables is valuable, use ridge
# regression - thank you https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r
library(glmnet)

eval_results <- function(true, predicted, df) {
        SSE <- sum((predicted - true)^2)
        SST <- sum((true - mean(true))^2)
        R_square <- 1 - SSE / SST
        MSE = SSE/nrow(df)
        
        # Model performance metrics
        data.frame(
                MSE = MSE,
                Rsquare = R_square
        )
}

ridge_result <- function(train, test) {
        x = as.matrix(train[, !(names(train) %in% c("DEXfat"))])
        y_train = train$DEXfat
        
        x_test = as.matrix(test[, !(names(test) %in% c("DEXfat"))])
        y_test = test$DEXfat
        
        lambdas <- seq(0, 1, 0.001)
        ridge_reg = glmnet(x, y_train, nlambda = 25, alpha = 0, family = 'gaussian', lambda = lambdas)
        
        cv_ridge <- cv.glmnet(x, y_train, alpha = 0, lambda = lambdas, nfolds = 5)
        optimal_lambda <- cv_ridge$lambda.min
        print(optimal_lambda)
        
        # Prediction and evaluation on train data
        predictions_train <- predict(ridge_reg, s = optimal_lambda, newx = x)
        print(eval_results(y_train, predictions_train, train))
        
        # Prediction and evaluation on test data
        predictions_test <- predict(ridge_reg, s = optimal_lambda, newx = x_test)
        eval_results(y_test, predictions_test, test)
}

# Get a collection of row ids that we will use for testing
test_idxs = sample(nrow(dfr), 20)

# Then perform the regression. First do it on the complete model
test = dfr[test_idxs, ]
train = dfr[-test_idxs, ]
ridge_result(train, test)

# Then do it on variables that were chosen by variable selection
test = dfr[test_idxs, (names(dfr) %in% c("DEXfat", "waistcirc", "hipcirc", "kneebreadth", "anthro3b", "anthro4"))]
train = dfr[-test_idxs, (names(dfr) %in% c("DEXfat", "waistcirc", "hipcirc", "kneebreadth", "anthro3b", "anthro4"))]
ridge_result(train, test)


