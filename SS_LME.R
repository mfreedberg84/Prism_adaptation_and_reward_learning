##Always start by clearing all existing variables.  
rm(list=ls(all=TRUE))

##Load relevant packages
library(nlme)    #Non-linear mixed effects package
library(car)   #Regression package
library(lme4)    #Linear mixed effects package
library(ggplot2)   #Advanced plotting package
library(languageR)  ##Packages containing useful functions for language research 
library(lattice)  ##Plotting package 
library(lmerTest)  ##do NOT load for Variance Components Analysis!!
library(influence.ME)
library(caTools)
library(kableExtra) # To make variable table

######################################## Schintu (2018) Reward Analysis ####################################

# Visuospatial cognition has an inherent lateralized bias. Individual differences 
# in the direction and magnitude of this bias are associated with asymmetrical 
# D2/3 dopamine binding and dopamine system genotypes. Dopamine level affects 
# feedback-based learning and dopamine signaling asymmetry is related to differential 
# learning from reward and punishment. High D2 binding in the left hemisphere is 
# associated with preference for reward. Prism adaptation (PA) is a simple sensorimotor 
# technique, which modulates visuospatial bias according to the direction of the 
# deviation. Left-deviating prism adaptation (LPA) induces rightward bias in healthy 
# subjects. It is therefore possible that the right side of space increases in 
# saliency along with left hemisphere dopaminergic activity. Right-deviating prism 
# adaptation (RPA) has been used mainly as a control condition because it does not
# modulate behavior in healthy individuals. Since LPA induces a rightward visuospatial 
# bias as a result of left hemisphere modulation, and higher dopaminergic activity 
# in the left hemisphere is associated with preference for rewarding events we hypothesized 
# that LPA would increase the preference for learning with reward. Healthy volunteers 
# performed a computer-based probabilistic classification task before and after 
# LPA or RPA. Consistent with our predictions, PA altered the preference for rewarded 
# versus punished learning, with the LPA group exhibiting increased learning from 
# reward. These results suggest that PA modulates dopaminergic activity in a lateralized 
# fashion.

# Task structure: Block 1 = Instructions, Blocks 2-4 = Reward learning
# Note that the three training blocks were divided into 4 sub-blocks to make 12 total.
#
# For each phase (pre and post adaptation) participants divided four stimuli
# in two categories using either the 'z' or 'm' button on a keyboard. Two of the
# stimuli yielded 25 points if performed correctly and 0 points if performed incorrectly.
# The remaining stimuli yielded 0 points if performed correctly and reduced the 
# point total by 25 if performed incorrectly. Participants were instructed to maximize
# their point total by learning to correctly categorize the stimuli. 

# Create table showing the description of variables
data_dict <- data.frame(row1 = c("Variable", 
                                 "Group", 
                                 "Subject",
                                 "Rew",
                                 "Block",
                                 "Phase",
                                 "Opt"
                                 ),
                        
                        row2 = c("Description",
                                 "Left and Right prism adaptation groups coded as -1 and 1, respectively",
                                 "Subject ID",
                                 "Punished and Rewarded stimuli are coded as -1 and 1, respectively",
                                 "Blocks span from 1-12",
                                 "Pre- and post-adaptation periods are coded as -1 and 1, respectively",
                                 "Proportion correct optimal responding (not log-transformed)"
                                 )
                        )

colnames(data_dict) <- NULL      
kable(data_dict, row.names = F) %>% column_spec (1:2, 
                                                 border_left = T, 
                                                 border_right = T
                                                 ) %>% kable_styling()
                 
############################ SETTING UP SHOP, READING IN DATA #################################

##Read in data
pth <- getwd()
df = read.table(paste(pth, "/SS_LME.csv", sep=), header = TRUE, sep = ",")

attach(df)
names(df)
head(df, n = 20)

# Logit transform the data to improve normality
logit_Opt <- logit(Opt, adjust = 0.04)
plot(Opt, logit_Opt)

# Designate categorical factors
df$Group<-as.factor(df$Group)
df$Rew<-as.factor(df$Rew)
df$Phase<-as.factor(df$Phase)

# Define the levels of each categorical variables
df$Rew<- factor(df$Rew,levels = c(-1,1),labels = c("Punished", "Rewarded")) 
df$Group<- factor(df$Group,levels = c(-1,1),labels = c("LPA", "RPA")) 
df$Phase <- factor(df$Phase,levels = c(-1,1),labels = c("Pre", "Post"))
levels(df$Rew)
levels(df$Group)
levels(df$Phase)

# Effect code the three contrasts
options(contrasts=c("contr.sum", "contr.poly"))
df$Rew <- factor(df$Rew,levels = c("Rewarded", "Punished")) 
df$Group <- factor(df$Group,levels = c("RPA", "LPA")) 
df$Phase <- factor(df$Phase,levels = c("Post", "Pre")) 
contrasts(df$Rew)
contrasts(df$Group) 
contrasts(df$Phase)

# Centering variables 
# NOTE: Centering data makes your intercept the halfway point.
df$Block_rev_c <-  df$Block - 5.5

# Print out descriptive summary
summary(df)  

############################ GRAPHING GROWTH CURVES #################################

# NOTE: for all of these plots make sure the x-variable is correct.  I use the term 
# Block_rev_c for the centered time variable and Block_rev for the  beginning of 
# the time course

#plot of individual Subject growth curves (raw data)
g1 <- ggplot(data = df, aes(x = Block, y = logit_Opt, group = Group))
g2 <- g1 + geom_line() + geom_point() + facet_wrap(~Subject)
g3 <- g2 + theme_bw() + scale_x_continuous(name = "Block_rev (minus 1)")
g4 <- g3 + scale_y_continuous(name = "Logit_Percent Correct")
print(g4)

#plot of individual Subject growth curves (linear fit)
g1 <- ggplot(data = df, aes(x = Block_rev_c, y = logit_Opt, group = Subject))
g2 <- g1 + geom_point() + stat_smooth(method = "lm", 
                                      se = FALSE) + facet_wrap(~Subject )
g3 <- g2 + theme_bw() + scale_x_continuous(name = "Block_rev (minus 1)")
g4 <- g3 + scale_y_continuous(name = "Percent Correct")
print(g4)


#plot of individual Subject growth curves(quadratic fit)
g1 <- ggplot(data = df, aes(x = Block_rev_c, y = logit_Opt, group = Subject))
g2 <- g1 + geom_point() + stat_smooth(method = "lm", 
                                      formula = y ~ poly(x,2), 
                                      se = FALSE
                                      ) + facet_wrap(~Subject )

g3 <- g2 + theme_bw() + scale_x_continuous(name = "Block_rev (minus 1)")
g4 <- g3 + scale_y_continuous(name = "Percent Correct")
print(g4)

############################ FUNCTIONAL FORM EVALUATION (Opt) #################################

## full random-effects and fixed-effects structure
## NOTE: The point of this section is to find the model that best fits the raw data.  Model 1 is the FULL model in that it 
## includes all effects (intercept, linear slope, quadratic effect).  From model 1 we peel back the random effects to see if the 
## model fits better without those effects.  To do this we compare model 1 and model 2 (model 2 is doesn't include a random effect for Quadratic slope)
## using the anova function. If there is a significant difference between the model fits,  we 
## defer to the more complex model.  

## HOW TO READ THE MODEL: I'm specifying each part of the code left to right
## 
## model1 <- lmer(Opt ~ Block_rev_c*Group + I(Block_rev_c^2)*Group + (1 + Block_rev_c + I(Block_rev_c^2)|Subject), data = df)
##
## model1			                      = name of the model
## <-					                      = gives
## Logit_Opt	                      = Log-transformed proportion correct (the DV)
## ~ 					                      = Now specify your model parameters starting with fixed effects
## Block_rev_c*Group*Rew*Phase 		  = We want to include a fixed effect of linear slope (block), our three categorical variables and their interactions	
## I(Block_rev_c^2)*Group*Rew*Phase	= We want to include a fixed effect of linear slope (quardratic), our three categorical variables and their interactions	
## +(					                      = Now were specifying our random effects
## 1 +				                      = This is a random intercept for participant
## Block_rev_c				              = Random linear slope for participant
## I(Block_rev_c^2) 			          = Random quadratic slope for participant
## | Subject),				              = This specifies that each of our random effects are for participant.  
## data = df		                    = data set (you can also look at a subset of data, such as when you are examining interactions.

# Run models with different complexitiy of random effects and compare their fits using an anova
Logit_model1 <- lmer(logit_Opt ~ Block_rev_c*Group*Rew*Phase 
                     + I(Block_rev_c^2)*Group*Rew*Phase 
                     + (1 + (Block_rev_c)
                        + I(Block_rev_c^2) | Subject), data = df)

Logit_model2 <- lmer(logit_Opt ~ Block_rev_c*Group*Rew*Phase 
                     + I(Block_rev_c^2)*Group*Rew*Phase 
                     + (1 + (Block_rev_c) | Subject), data = df)

Logit_model3 <- lmer(logit_Opt ~ Block_rev_c*Group*Rew*Phase 
                     + I(Block_rev_c^2)*Group*Rew*Phase 
                     + (1 | Subject), data = df)

anova(Logit_model1, Logit_model2, Logit_model3)
# The most complex model including random effects of intercept and linear and quadratic slope is the most appropriate model
summary(Logit_model1)

# Results of model analysis
#Fixed effects:
#                                         Estimate    Std. Error  df          t-value   Pr(>|t|)    
# (Intercept)                             1.497e+00   1.508e-01   4.407e+01   9.930     8.16e-13 ***
# Block_rev_c                             1.109e-01   1.310e-02   4.401e+01   8.467     8.68e-11 ***
# Group1                                  -2.545e-02  1.508e-01   4.407e+01   -0.169    0.8667    
# Rew1                                    -8.274e-03  4.465e-02   2.096e+03   -0.185    0.8530    
# Phase1                                  2.551e-01   4.465e-02   2.096e+03   5.713     1.27e-08 ***
# I(Block_rev_c^2)                        -6.619e-03  3.425e-03   7.613e+01   -1.933    0.0570 .  
# Block_rev_c:Group1                      -1.755e-02  1.310e-02   4.401e+01   -1.340    0.1871    
# Block_rev_c:Rew1                        2.821e-03   8.572e-03   2.096e+03   0.329     0.7422    
# Group1:Rew1                             -1.795e-01  4.465e-02   2.096e+03   -4.020    6.03e-05 ***
# Block_rev_c:Phase1                      8.203e-03   8.572e-03   2.096e+03   0.957     0.3387    
# Group1:Phase1                           -4.989e-02  4.465e-02   2.096e+03   -1.117    0.2640    
# Rew1:Phase1                             2.807e-02   4.465e-02   2.096e+03   0.629     0.5296    
# Group1:I(Block_rev_c^2)                 2.187e-03   3.425e-03   7.613e+01   0.639     0.5250    
# Rew1:I(Block_rev_c^2)                   2.783e-03   2.806e-03   2.096e+03   0.992     0.3215    
# Phase1:I(Block_rev_c^2)                 -3.393e-03  2.806e-03   2.096e+03   -1.209    0.2267    
# Block_rev_c:Group1:Rew1                 -2.181e-02  8.572e-03   2.096e+03   -2.545    0.0110 *  
# Block_rev_c:Group1:Phase1               -1.875e-02  8.572e-03   2.096e+03   -2.188    0.0288 *  
# Block_rev_c:Rew1:Phase1                 3.808e-03   8.572e-03   2.096e+03   0.444     0.6569    
# Group1:Rew1:Phase1                      -1.145e-01  4.465e-02   2.096e+03   -2.565    0.0104 *  
# Group1:Rew1:I(Block_rev_c^2)            9.410e-04   2.806e-03   2.096e+03   0.335     0.7374    
# Group1:Phase1:I(Block_rev_c^2)          7.747e-04   2.806e-03   2.096e+03   0.276     0.7825    
# Rew1:Phase1:I(Block_rev_c^2)            -6.205e-04  2.806e-03   2.096e+03   -0.221    0.8250    
# Block_rev_c:Group1:Rew1:Phase1          -1.864e-02  8.572e-03   2.096e+03   -2.175    0.0298 *  
# Group1:Rew1:Phase1:I(Block_rev_c^2)     2.288e-03   2.806e-03   2.096e+03   0.815     0.4150   

## NOTE: This part serves two functions.  The 'cor' function will give you a figure 
## that represents how much of the variability in the model is explained byyour fixed effects.  
## The write.table function writes out your fitted values (this is how I produced the figures in excel). 
## Specify where you want the function to write the text file and what to call it.  

FittedFE <- function(x) model.matrix(x) %*% fixef(x)
head(FittedFE(Logit_model1), n = 20)
plotdata <- data.frame(Logit_model1@frame, fitted = FittedFE(Logit_model1))
head(plotdata, n = 20)
cor(y = Logit_model1@frame$logit_Opt, x = FittedFE(Logit_model1)) ^2
write.table(plotdata, paste(pth, "/SS_Data.txt",sep=""), sep="\t")

# Output the model coefficients
coef(Logit_model1)

############################ Unpack four-way interaction #################################

# The four way interaction was significant, suggesting that Reward learning improved
# significantly more for the LPA group over the RPA group. To confirm this, we will
# split the data by Group (LPA vs. RPA) and rerun the analysis.

### Subset the LPA
df_LPA <- df[ which(df$Group == "LPA"), ]

attach(df_LPA)
keep <- c("Subject", "Rew","Block","Opt","Phase")
df_LPA <- df_LPA[keep]
logit_Opt <- logit(Opt, adjust = 0.04)

# Define contrasts
options(contrasts=c("contr.sum", "contr.poly"))
df_LPA$Rew <- factor(df_LPA$Rew,levels = c("Rewarded", "Punished"))  
df_LPA$Phase <- factor(df_LPA$Phase,levels = c("Post", "Pre")) 
contrasts(df_LPA$Rew) 
contrasts(df_LPA$Phase)

# Center time variable
df_LPA$Block_rev_c <-  df_LPA$Block - 5.5
summary(df_LPA)

# Rerun logit model 1 with only LPA group
LPA_model <- lmer(logit_Opt ~ Block_rev_c*Rew*Phase
                  + I(Block_rev_c^2)*Rew*Phase 
                  + (1 + (Block_rev_c)
                     + I(Block_rev_c^2) |Subject), data = df_LPA)
summary(LPA_model)

# Results for LPA 
# Fixed effects:
#                               Estimate   Std. Error df          t value Pr(>|t|)    
# (Intercept)                   1.523e+00  2.019e-01  2.330e+01   7.544   1.06e-07 ***
# Block_rev_c                   1.285e-01  2.002e-02  2.308e+01   6.418   1.48e-06 ***
# Rew1                          1.712e-01  5.843e-02  1.094e+03   2.930   0.00346 ** 
# Phase1                        3.050e-01  5.843e-02  1.094e+03   5.219   2.15e-07 ***
# I(Block_rev_c^2)             -8.806e-03  4.976e-03  3.194e+01  -1.770   0.08635 .  
# Block_rev_c:Rew1              2.464e-02  1.122e-02  1.094e+03   2.196   0.02830 *  
# Block_rev_c:Phase1            2.696e-02  1.122e-02  1.094e+03   2.403   0.01643 *  
# Rew1:Phase1                   1.426e-01  5.843e-02  1.094e+03   2.440   0.01484 *  
# Rew1:I(Block_rev_c^2)         1.842e-03  3.672e-03  1.094e+03   0.502   0.61609    
# Phase1:I(Block_rev_c^2)      -4.168e-03  3.672e-03  1.094e+03  -1.135   0.25664    
# Block_rev_c:Rew1:Phase1       2.245e-02  1.122e-02  1.094e+03   2.001   0.04561 *  
# Rew1:Phase1:I(Block_rev_c^2) -2.908e-03  3.672e-03  1.094e+03  -0.792   0.42859    

### Subset the RPA
df_RPA <- df[ which(df$Group == "RPA"), ]

attach(df_RPA)
keep <- c("Subject", "Rew","Block","Opt","Phase")
df_RPA <- df_RPA[keep]
logit_Opt <- logit(Opt, adjust = 0.04)

# Define contrasts
options(contrasts=c("contr.sum", "contr.poly"))
df_RPA$Rew <- factor(df_RPA$Rew,levels = c("Rewarded", "Punished"))  
df_RPA$Phase <- factor(df_RPA$Phase,levels = c("Post", "Pre")) 
contrasts(df_RPA$Rew) 
contrasts(df_RPA$Phase)

# Center time variable
df_RPA$Block_rev_c <-  df_RPA$Block - 5.5
summary(df_RPA)

# Rerun logit model 1 with only LPA group
RPA_model <- lmer(logit_Opt ~ Block_rev_c*Rew*Phase
                  + I(Block_rev_c^2)*Rew*Phase 
                  + (1 + (Block_rev_c) 
                     + I(Block_rev_c^2) |Subject), data = df_RPA)

summary(RPA_model)

# Results for RPA 
# Fixed effects:
#                               Estimate  Std. Error  df          t value Pr(>|t|)    
# (Intercept)                   1.472e+00  2.256e-01  2.093e+01   6.524   1.86e-06 ***
# Block_rev_c                   9.335e-02  1.681e-02  2.272e+01   5.554   1.25e-05 ***
# Rew1                         -1.878e-01  6.790e-02  1.002e+03  -2.765   0.00579 ** 
# Phase1                        2.052e-01  6.790e-02  1.002e+03   3.022   0.00257 ** 
# I(Block_rev_c^2)             -4.432e-03  5.020e-03  3.666e+01  -0.883   0.38308    
# Block_rev_c:Rew1             -1.899e-02  1.303e-02  1.002e+03  -1.457   0.14538    
# Block_rev_c:Phase1           -1.055e-02  1.303e-02  1.002e+03  -0.809   0.41850    
# Rew1:Phase1                  -8.644e-02  6.790e-02  1.002e+03  -1.273   0.20325    
# Rew1:I(Block_rev_c^2)         3.724e-03  4.267e-03  1.002e+03   0.873   0.38301    
# Phase1:I(Block_rev_c^2)      -2.618e-03  4.267e-03  1.002e+03  -0.614   0.53959    
# Block_rev_c:Rew1:Phase1      -1.483e-02  1.303e-02  1.002e+03  -1.138   0.25538 No interaction between block, reward, and phase
# Rew1:Phase1:I(Block_rev_c^2)  1.667e-03  4.267e-03  1.002e+03   0.391   0.69610 

############################ Unpack three-way interaction for the LPA group #################################

# the LPA, but not RPA, group had a significant three-way interaction between block, reward, and phase,
# suggesting that rewarded stimuli were learned significantly better than punished stimuli
# To determine whether this is true, let's unpack the three-way interaction by reward

### Subset the LPA int rewarded trials
df_LPA_rew <- df_LPA[ which(df_LPA$Rew == "Rewarded"), ]

attach(df_LPA_rew)
keep <- c("Subject", "Rew","Block","Opt","Phase")
df_LPA_rew <- df_LPA_rew[keep]
logit_Opt <- logit(Opt, adjust = 0.04)

# Define contrasts
options(contrasts=c("contr.sum", "contr.poly"))
df_LPA_rew$Phase <- factor(df_LPA_rew$Phase,levels = c("Post", "Pre")) 
contrasts(df_LPA_rew$Phase)

# Center time variable
df_LPA_rew$Block_rev_c <-  df_LPA_rew$Block - 5.5
summary(df_LPA_rew)

# Rerun logit model LPA with only rewrad trials
LPA_model_rew <- lmer(logit_Opt ~ Block_rev_c*Phase
                      + I(Block_rev_c^2)*Phase 
                      + (1 + (Block_rev_c)
                         + I(Block_rev_c^2) | Subject), data = df_LPA_rew)

summary(LPA_model_rew)

# Results
# Fixed effects:
#                           Estimate   Std. Error df          t value   Pr(>|t|)    
# (Intercept)               1.694036   0.262048   22.848720   6.465     1.39e-06 ***
# Block_rev_c               0.153092   0.025950   23.017496   5.900     5.15e-06 ***
# Phase1                    0.447563   0.079103   524.064606  5.658     2.52e-08 ***
# I(Block_rev_c^2)         -0.006964   0.006509   32.955516   -1.070    0.29240    
# Block_rev_c:Phase1        0.049407   0.015186   524.064606  3.253     0.00121 ** Significant evidence of improved learning from reward in LPA group
# Phase1:I(Block_rev_c^2)  -0.007076   0.004971   524.064606  -1.423    0.15521  

### Subset the LPA int punished trials
df_LPA_pun <- df_LPA[ which(df_LPA$Rew == "Punished"), ]

attach(df_LPA_pun)
keep <- c("Subject", "Rew","Block","Opt","Phase")
df_LPA_pun <- df_LPA_pun[keep]
logit_Opt <- logit(Opt, adjust = 0.04)

# Define contrasts
options(contrasts=c("contr.sum", "contr.poly"))
df_LPA_pun$Phase <- factor(df_LPA_pun$Phase,levels = c("Post", "Pre")) 
contrasts(df_LPA_pun$Phase)

# Center time variable
df_LPA_pun$Block_rev_c <-  df_LPA_pun$Block - 5.5
summary(df_LPA_pun)

# Rerun logit model LPA with only rewrad trials
LPA_model_pun <- lmer(logit_Opt ~ Block_rev_c*Phase
                      + I(Block_rev_c^2)*Phase 
                      + (1 + (Block_rev_c)
                         + I(Block_rev_c^2) | Subject), data = df_LPA_pun)

summary(LPA_model_pun)

# Results
#                           Estimate    Std. Error  df          t value   Pr(>|t|)    
# (Intercept)               1.351608    0.187476    22.573206   7.210     2.72e-07 ***
# Block_rev_c               0.103822    0.022314    23.388537   4.653     0.000106 ***
# Phase1                    0.162384    0.077018    523.966839  2.108     0.035471 *  
# I(Block_rev_c^2)         -0.010648    0.005227    56.553370   -2.037    0.046331 *  
# Block_rev_c:Phase1        0.004505    0.014786    523.966839  0.305     0.760726 <-- No strong evidence of LPA adaptation on punished learning
# Phase1:I(Block_rev_c^2)  -0.001260    0.004840    523.966839  -0.260    0.794748   

# Overall: Left, as opposed to right, prism goggles isgnificantly improved 
# reward learning. The improvement in reward learning was significantly greater 
# in the LPA group compared to the RPA group