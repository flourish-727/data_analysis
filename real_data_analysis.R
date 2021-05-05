# This R script is for the real data analysis of the publication
# "Decomposition of the Total Effect for Two Mediators:
#  A Natural Counterfactual Interaction Effect Framework".
# 2013-2014 2015-2016 and 2017-2018 combined data set


setwd("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM")

# use the following packages
library(haven)
library(dplyr)

##############################################################################
##############################################################################
# read the xpt files for 2013-2014 data set
SA_13 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/DEMO_H.XPT")
BMI_13 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BMX_H.XPT")
SBP_13 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BPX_H.XPT")
GGT_13 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BIOPRO_H.XPT")
Alcohol_13 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/ALQ_H.XPT")

# combine the data sets 
my_data_13 <- list(SA_13, BMI_13, SBP_13, GGT_13, Alcohol_13) %>%
  Reduce(function(df1, df2) inner_join(df1, df2, by = "SEQN"), .)


my_data_13 <- as.matrix(my_data_13)
my_data_13 <- as.data.frame(my_data_13)

head(my_data_13)
############################################### data cleaning #################################
############################################### 2013-2014 #####################################
# Outcome Y is the average of BPXSY1, BPXSY2 and BPXSY3 (SBP)
# Mediator M1 is BMXBMI (BMI)
# Mediator M2 is LBXSGTSI (GGT)
# Exposure A is ALQ130 (alcohol drinking)
# Two confounders are C1-Sex and C2-Age, respectively. 

# collect columns
natural_13 <- data.frame(my_data_13$ALQ130, my_data_13$BMXBMI, my_data_13$LBXSGTSI, my_data_13$BPXSY1, my_data_13$BPXSY2,
                         my_data_13$BPXSY3, my_data_13$RIAGENDR, my_data_13$RIDAGEYR)


names(natural_13)<- c("A","M1","M2","Y1","Y2","Y3","C1","C2")

# the outcome Y is the average of three SBP values 
natural_a_13 <- data.frame(natural_13, (natural_13$Y1 + natural_13$Y2 + natural_13$Y3)/3)

head(natural_a_13)
dim(natural_a_13)

# Assign name to the last column which corresponds to outcome
names(natural_a_13)[dim(natural_a_13)[2]] <- c("Y")

# delete missing values 
natural_a_13 <- na.omit(natural_a_13)

###################################################################################
# ALQ130 data cleaning
length(which(natural_a_13$A >= 1 & natural_a_13$A <= 25)) # 3101 obs
length(which(natural_a_13$A == 777)) #none refused
length(which(natural_a_13$A == 999)) # 3 obs do not know

# delete 999 obs
natural_a_13 <- natural_a_13[-which(natural_a_13$A == 999), ]

# For males #
# define never/moderate drinking as 1-2 drinks or less in a day (0 group)
# define heavy drinking as 3 and more drinks in a day (1 group)
natural_a_13$A[which(natural_a_13$C1 == 1 & (natural_a_13$A == 1 | natural_a_13$A == 2))] <- 0
natural_a_13$A[which(natural_a_13$C1 == 1 & natural_a_13$A != 0)] <- 1

# For females #
# define never/moderate drinking as 1 drink or less in a day (0 group)
# define heavy drinking as 2 and more drinks in a day (1 group)
natural_a_13$A[which(natural_a_13$C1 == 2 & natural_a_13$A == 1)] <- 0
natural_a_13$A[which(natural_a_13$C1 == 2 & natural_a_13$A != 0)] <- 1

###################################################################################
# check if there are any missing values labeled with a number other than 1 or 2 
# for Sex variable
length(which(natural_a_13$C1 != 1 & natural_a_13$C1 != 2)) # none 

# remove extreme obs 
length(which(natural_a_13$Y >= 180)) # 16 obs
natural_a_13 <- natural_a_13[-which(natural_a_13$Y >= 180), ]

# Assign females to 0 group and males to 1 group 
natural_a_13$C1[which(natural_a_13$C1 == 2)] <- 0


# log transformation on GGT is suggested since the skewness
natural_a_13$M2 <- log(natural_a_13$M2)

# reset the row names 
row.names(natural_a_13) <- NULL

# double check the data frame
length(which(natural_a_13$A == 0)) + length(which(natural_a_13$A == 1)) == dim(natural_a_13)[1]

length(which(natural_a_13$C1 == 0)) + length(which(natural_a_13$C1 == 1)) == dim(natural_a_13)[1]

is.vector(natural_a_13$M1)
is.numeric(natural_a_13$M1)
length(natural_a_13$M1) == dim(natural_a_13)[1]

is.vector(natural_a_13$M2)
is.numeric(natural_a_13$M2)
length(natural_a_13$M2) == dim(natural_a_13)[1]

is.vector(natural_a_13$C2)
is.numeric(natural_a_13$C2)
length(natural_a_13$C2) == dim(natural_a_13)[1]

is.vector(natural_a_13$Y)
is.numeric(natural_a_13$Y)
length(natural_a_13$Y) == dim(natural_a_13)[1]

dim(natural_a_13) # 3085 x 9
#############################################################################
################# 2013-2014 data cleaning process ends here #################
#############################################################################

##############################################################################
##############################################################################
# read the xpt files for 2015-2016 data set
SA_15 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/DEMO_I.XPT")
BMI_15 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BMX_I.XPT")
SBP_15 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BPX_I.XPT")
GGT_15 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BIOPRO_I.XPT")
Alcohol_15 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/ALQ_I.XPT")

# combine the data sets 
my_data_15 <- list(SA_15, BMI_15, SBP_15, GGT_15, Alcohol_15) %>%
  Reduce(function(df1, df2) inner_join(df1, df2, by = "SEQN"), .)

my_data_15 <- as.matrix(my_data_15)
my_data_15 <- as.data.frame(my_data_15)

head(my_data_15)

############################################### data cleaning #################################
############################################### 2015-2016 #####################################
# Outcome Y is the average of BPXSY1, BPXSY2 and BPXSY3 (SBP)
# Mediator M1 is BMXBMI (BMI)
# Mediator M2 is LBXSGTSI (GGT)
# Exposure A is ALQ130 (alcohol drinking)
# Two confounders are C1-Sex and C2-Age, respectively. 

# collect columns
natural_15 <- data.frame(my_data_15$ALQ130, my_data_15$BMXBMI, my_data_15$LBXSGTSI, my_data_15$BPXSY1, my_data_15$BPXSY2,
                         my_data_15$BPXSY3, my_data_15$RIAGENDR, my_data_15$RIDAGEYR)

names(natural_15)<- c("A","M1","M2","Y1","Y2","Y3","C1","C2")

# the outcome Y is the average of three SBP values 
natural_a_15 <- data.frame(natural_15, (natural_15$Y1 + natural_15$Y2 + natural_15$Y3)/3)

head(natural_a_15)
dim(natural_a_15)

# Assign name to the last column which corresponds to outcome
names(natural_a_15)[dim(natural_a_15)[2]] <- c("Y")

# delete missing values 
natural_a_15 <- na.omit(natural_a_15)

###################################################################################
# ALQ130 data cleaning
length(which(natural_a_15$A >= 1 & natural_a_15$A <= 15)) # 3019 obs
length(which(natural_a_15$A == 777)) #none refused
length(which(natural_a_15$A == 999)) # 3 obs do not know

# delete 999 obs
natural_a_15 <- natural_a_15[-which(natural_a_15$A == 999), ]

# For males #
# define never/moderate drinking as 1-2 drinks or less in a day (0 group)
# define heavy drinking as 3 and more drinks in a day (1 group)
natural_a_15$A[which(natural_a_15$C1 == 1 & (natural_a_15$A == 1 | natural_a_15$A == 2))] <- 0
natural_a_15$A[which(natural_a_15$C1 == 1 & natural_a_15$A != 0)] <- 1

# For females #
# define never/moderate drinking as 1 drink or less in a day (0 group)
# define heavy drinking as 2 and more drinks in a day (1 group)
natural_a_15$A[which(natural_a_15$C1 == 2 & natural_a_15$A == 1)] <- 0
natural_a_15$A[which(natural_a_15$C1 == 2 & natural_a_15$A != 0)] <- 1

###################################################################################
# check if there are any missing values labeled with a number other than 1/2 
# for Sex variable
length(which(natural_a_15$C1 != 1 & natural_a_15$C1 != 2)) # none 

# remove extreme obs 
length(which(natural_a_15$Y >= 180)) # 27 obs
natural_a_15 <- natural_a_15[-which(natural_a_15$Y >= 180), ]

# Assign females to 0 group and males to 1 group 
natural_a_15$C1[which(natural_a_15$C1 == 2)] <- 0

# log transformation on GGT is suggested since the skewness
natural_a_15$M2 <- log(natural_a_15$M2)

# reset the row names 
row.names(natural_a_15) <- NULL

# double check the data frame
length(which(natural_a_15$A == 0)) + length(which(natural_a_15$A == 1)) == dim(natural_a_15)[1]

length(which(natural_a_15$C1 == 0)) + length(which(natural_a_15$C1 == 1)) == dim(natural_a_15)[1]

is.vector(natural_a_15$M1)
is.numeric(natural_a_15$M1)
length(natural_a_15$M1) == dim(natural_a_15)[1]

is.vector(natural_a_15$M2)
is.numeric(natural_a_15$M2)
length(natural_a_15$M2) == dim(natural_a_15)[1]

is.vector(natural_a_15$C2)
is.numeric(natural_a_15$C2)
length(natural_a_15$C2) == dim(natural_a_15)[1]

is.vector(natural_a_15$Y)
is.numeric(natural_a_15$Y)
length(natural_a_15$Y) == dim(natural_a_15)[1]

dim(natural_a_15) # 2992 x 9
#############################################################################
################# 2015-2016 data cleaning process ends here #################
#############################################################################

#############################################################################
#############################################################################
# 2017-2018 data set
# read the xpt files
SA_17 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/DEMO_J.XPT")
BMI_17 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BMX_J.XPT")
SBP_17 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BPX_J.XPT")
GGT_17 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/BIOPRO_J.XPT")
Alcohol_17 <- read_xpt("D:/JCI/Again_JCI/INTREFAM1M2/Year/Combined/FM/ALQ_J.XPT")

# combine the data sets 
my_data_17 <- list(SA_17, BMI_17, SBP_17, GGT_17, Alcohol_17) %>%
  Reduce(function(df1, df2) inner_join(df1, df2, by = "SEQN"), .)


my_data_17 <- as.matrix(my_data_17)
my_data_17 <- as.data.frame(my_data_17)

head(my_data_17)

############################################### data cleaning #################################
############################################### 2017-2018 #####################################
# Outcome Y is the average of BPXSY1, BPXSY2 and BPXSY3 (SBP)
# Mediator M1 is BMXBMI (BMI)
# Mediator M2 is LBXSGTSI (GGT)
# Exposure A is ALQ130 (alcohol drinking)
# Two confounders are C1-Sex and C2-Age, respectively. 

# collect columns
natural_17 <- data.frame(my_data_17$ALQ130, my_data_17$BMXBMI, my_data_17$LBXSGTSI, my_data_17$BPXSY1, my_data_17$BPXSY2,
                         my_data_17$BPXSY3, my_data_17$RIAGENDR, my_data_17$RIDAGEYR)

names(natural_17)<- c("A","M1","M2","Y1","Y2","Y3","C1","C2")

# the outcome Y is the average of three SBP values
natural_a_17 <- data.frame(natural_17, (natural_17$Y1 + natural_17$Y2 + natural_17$Y3)/3)

head(natural_a_17)
dim(natural_a_17)

# Assign name to the last column which corresponds to outcome
names(natural_a_17)[dim(natural_a_17)[2]] <- c("Y")

# delete missing values 
natural_a_17 <- na.omit(natural_a_17)

#############################################################################
# ALQ130 data cleaning
length(which(natural_a_17$A >= 1 & natural_a_17$A <= 15)) # 2875 obs
length(which(natural_a_17$A == 777)) # 1 obs refused
length(which(natural_a_17$A == 999)) # 4 obs do not know

# delete 999 and 777
natural_a_17 <- natural_a_17[-which(natural_a_17$A == 999), ]
natural_a_17 <- natural_a_17[-which(natural_a_17$A == 777), ]

# For males #
# define never/moderate drinking as 1-2 drinks or less in a day (0 group)
# define heavy drinking as 3 and more drinks in a day (1 group)
natural_a_17$A[which(natural_a_17$C1 == 1 & (natural_a_17$A == 1 | natural_a_17$A == 2))] <- 0
natural_a_17$A[which(natural_a_17$C1 == 1 & natural_a_17$A != 0)] <- 1

# For females #
# define never/moderate drinking as 1 drink or less in a day (0 group)
# define heavy drinking as 2 and more drinks in a day (1 group)
natural_a_17$A[which(natural_a_17$C1 == 2 & natural_a_17$A == 1)] <- 0
natural_a_17$A[which(natural_a_17$C1 == 2 & natural_a_17$A != 0)] <- 1

###################################################################################
# check if there are any missing values labeled with a number other than 1/2 
# for Sex variable
length(which(natural_a_17$C1 != 1 & natural_a_17$C1 != 2)) # none


# remove extreme obs
length(which(natural_a_17$Y >= 180)) # 32 obs
natural_a_17 <- natural_a_17[-which(natural_a_17$Y >= 180), ]


# Assign females to 0 group and males to 1 group
natural_a_17$C1[which(natural_a_17$C1 == 2)] <- 0


# log transformation on GGT is suggested since the skewness
natural_a_17$M2 <- log(natural_a_17$M2)

# reset the row names
row.names(natural_a_17) <- NULL

# double check 
length(which(natural_a_17$A == 0)) + length(which(natural_a_17$A == 1)) == dim(natural_a_17)[1]

length(which(natural_a_17$C1 == 0)) + length(which(natural_a_17$C1 == 1)) == dim(natural_a_17)[1]

is.vector(natural_a_17$M1)
is.numeric(natural_a_17$M1)
length(natural_a_17$M1) == dim(natural_a_17)[1]

is.vector(natural_a_17$M2)
is.numeric(natural_a_17$M2)
length(natural_a_17$M2) == dim(natural_a_17)[1]

is.vector(natural_a_17$C2)
is.numeric(natural_a_17$C2)
length(natural_a_17$C2) == dim(natural_a_17)[1]

is.vector(natural_a_17$Y)
is.numeric(natural_a_17$Y)
length(natural_a_17$Y) == dim(natural_a_17)[1]

dim(natural_a_17) # 2843 x 9

#############################################################################
################# 2017-2018 data cleaning process ends here #################
#############################################################################


#############################################################################
# Combine the 2013-2014, 2015-2016 and 2017-2018 data sets 
natural_a <- rbind.data.frame(natural_a_13, natural_a_15)
natural_a <- rbind.data.frame(natural_a, natural_a_17)

dim(natural_a) # 8920 x 9

length(which(natural_a$A == 0)) # 4420 obs 
length(which(natural_a$A == 1)) # 4500 obs
############################################ NOTE !!!!!!!!!!!! #############################################
# In the data frame, the data for M2 is actually log-transformed M2
# We still use the name "M2" for the column name just for simplicity
############################################################################################################

################################################ Data Cleaning Ends ###################################

################################################ Point Estimation ########################################
# Please note that the results for a non-sequential two-mediator scenario will be obtained
# if beta_2 and beta_3 are set to zero in M2 model. If this is the case, then the term 
# pure indirect effect through M2 should be used. 
# It beta_2 and beta_3 are non-zero, then the term seminatural indirect effect through M2
# should be used since it becomes a sequential two-mediator scenario in this setting. 
# In this code, we used PIEM2 all the time just for simplicity. 
#####################################################################################

attach(natural_a)

# outcome model
model_Y <- lm(Y ~ A + M1 + M2 + A*M1 + A*M2 + M1*M2 + A*M1*M2 + C1 + C2, data = natural_a)

# M2 model
model_M2 <- lm(M2 ~ A + M1 + A*M1 + C1 + C2, data = natural_a)

# M1 model
model_M1 <- lm(M1 ~ A + C1 + C2, data = natural_a)

# calculate estimates 

# MSE of M1 model
mse <- sum(model_M1$residuals^2)/(length(model_M1$residuals)-4)

# set the arbitrary reference level of M1 and M2 to be the mean
m1 <- mean(M1)
m2 <- mean(M2)

# C1 is binary so it has two strata
c1m <- 1
c1f <- 0

# use the mean value of C2 for the calculation 
c2 <- mean(C2)

# the reference level of A is 0 and the treatment level of A is 1
a <- 1
a_s <- 0

t0 <- unname(model_Y$coefficients[1])
t1 <- unname(model_Y$coefficients[2])
t2 <- unname(model_Y$coefficients[3])
t3 <- unname(model_Y$coefficients[4])
t4 <- unname(model_Y$coefficients[7])
t5 <- unname(model_Y$coefficients[8])
t6 <- unname(model_Y$coefficients[9])
t7 <- unname(model_Y$coefficients[10])
t8 <- unname(model_Y$coefficients[5])
t9 <- unname(model_Y$coefficients[6])

b0 <- unname(model_M2$coefficients[1])
b1 <- unname(model_M2$coefficients[2])
b2 <- unname(model_M2$coefficients[3])
b3 <- unname(model_M2$coefficients[6])
b4 <- unname(model_M2$coefficients[4])
b5 <- unname(model_M2$coefficients[5])

g0 <- unname(model_M1$coefficients[1])
g1 <- unname(model_M1$coefficients[2])
g2 <- unname(model_M1$coefficients[3])
g3 <- unname(model_M1$coefficients[4])


# CDE(m1,m2)
CDE <- (t1 + t4*m1 + t5*m2 + t7*m1*m2)*(a - a_s)

# Reference Interaction between A and M1
# INT_ref_AM1(m1,m2)
# C1 = Male
INT_ref_AM1_m <- (g0 + g1*a_s + g2*c1m + g3*c2 - m1)*(t4 + t7*m2)*(a - a_s)

# C1 = Female 
INT_ref_AM1_f <- (g0 + g1*a_s + g2*c1f + g3*c2 - m1)*(t4 + t7*m2)*(a - a_s)

# Reference Interaction between A and M2
# INT_ref_AM2(m1,m2)
# C1 = Male
INT_ref_AM2_m <- (t5 + t7*m1)*(b0 + b1*a_s + b2*m1 + b3*a_s*m1 + b4*c1m + b5*c2 - m2)*(a - a_s)

# C1 = Female 
INT_ref_AM2_f <- (t5 + t7*m1)*(b0 + b1*a_s + b2*m1 + b3*a_s*m1 + b4*c1f + b5*c2 - m2)*(a - a_s)


# Reference Interaction between A, M1 and M2
# First need to find the Sum of two reference interaction effects
# INT_ref_AM2+AM1M2(m2)
# C1 = Male 
INT_ref_AM2_AM1M2_m <- (t1 + t5*(b0 + b1*a_s + b4*c1m + b5*c2) + 
                          t7*(b0 + b1*a_s + b4*c1m + b5*c2)*(g0 + g1*a_s + g2*c1m + g3*c2) +
                          t5*(b2 + b3*a_s)*(g0 + g1*a_s + g2*c1m + g3*c2) + 
                          t7*(b2 + b3*a_s)*(mse + (g0 + g1*a_s + g2*c1m + g3*c2)^2) - 
                          (t1 + t5*m2) - t7*m2*(g0 + g1*a_s + g2*c1m + g3*c2))*(a - a_s)

# C1 = Female 
INT_ref_AM2_AM1M2_f <- (t1 + t5*(b0 + b1*a_s + b4*c1f + b5*c2) + 
                          t7*(b0 + b1*a_s + b4*c1f + b5*c2)*(g0 + g1*a_s + g2*c1f + g3*c2) +
                          t5*(b2 + b3*a_s)*(g0 + g1*a_s + g2*c1f + g3*c2) + 
                          t7*(b2 + b3*a_s)*(mse + (g0 + g1*a_s + g2*c1f + g3*c2)^2) - 
                          (t1 + t5*m2) - t7*m2*(g0 + g1*a_s + g2*c1f + g3*c2))*(a - a_s)

# Next we subtract INT_ref_AM2 from INT_ref_AM2_AM1M2

# INT_ref_AM1M2(m1,m2)
# C1 = Male 
INT_ref_AM1M2_m <- INT_ref_AM2_AM1M2_m - INT_ref_AM2_m

# C1 = Female 
INT_ref_AM1M2_f <- INT_ref_AM2_AM1M2_f - INT_ref_AM2_f






# natural counterfactual interaction effects 
# NAT_INT_AM1
# C1 = Male
NAT_INT_AM1_m <- (t4*g1 + t7*g1*(b0 + b1*a_s + b4*c1m + b5*c2) + t5*g1*(b2 + b3*a_s) + 
                    2*t7*g1*(b2 + b3*a_s)*(g0 + g2*c1m +g3*c2) + 
                    t7*g1^2*(b2 + b3*a_s)*(a + a_s))*(a - a_s)^2

# C1 = Female
NAT_INT_AM1_f <- (t4*g1 + t7*g1*(b0 + b1*a_s + b4*c1f + b5*c2) + t5*g1*(b2 + b3*a_s) + 
                    2*t7*g1*(b2 + b3*a_s)*(g0 + g2*c1f +g3*c2) + 
                    t7*g1^2*(b2 + b3*a_s)*(a + a_s))*(a - a_s)^2

# NAT_INT_AM2
# C1 = Male
NAT_INT_AM2_m <- (t5*b1 + t7*b1*(g0 + g1*a_s + g2*c1m + g3*c2) + t5*b3*(g0 + g1*a_s + g2*c1m + g3*c2) + 
                    t7*b3*(mse + (g0 + g1*a_s + g2*c1m + g3*c2)^2))*(a - a_s)^2

# C1 = Female
NAT_INT_AM2_f <- (t5*b1 + t7*b1*(g0 + g1*a_s + g2*c1f + g3*c2) + t5*b3*(g0 + g1*a_s + g2*c1f + g3*c2) + 
                    t7*b3*(mse + (g0 + g1*a_s + g2*c1f + g3*c2)^2))*(a - a_s)^2

# NAT_INT_AM1M2
# C1 = Male
NAT_INT_AM1M2_m <- (t7*b1*g1 + t5*b3*g1 + 2*t7*b3*g1*(g0 + g2*c1m +g3*c2) + 
                      t7*b3*g1^2*(a + a_s))*(a - a_s)^3

# C1 = Female
NAT_INT_AM1M2_f <- (t7*b1*g1 + t5*b3*g1 + 2*t7*b3*g1*(g0 + g2*c1f +g3*c2) + 
                      t7*b3*g1^2*(a + a_s))*(a - a_s)^3


# NAT_INT_M1M2 
# C1 = Male
NAT_INT_M1M2_m <- (b1*g1*(t6 + t7*a_s) + b3*g1*(t3 + t5*a_s) + 2*b3*g1*(t6 + t7*a_s)*(g0 + g2*c1m +g3*c2) + 
                     b3*g1^2*(t6 + t7*a_s)*(a + a_s))*(a - a_s)^2

# C1 = Female
NAT_INT_M1M2_f <- (b1*g1*(t6 + t7*a_s) + b3*g1*(t3 + t5*a_s) + 2*b3*g1*(t6 + t7*a_s)*(g0 + g2*c1f +g3*c2) + 
                     b3*g1^2*(t6 + t7*a_s)*(a + a_s))*(a - a_s)^2


# Pure Indirect Effects 

# PIE_M1
# C1 = Male
PIE_M1_m <- (g1*(t2 + t4*a_s) + g1*(t6 + t7*a_s)*(b0 + b1*a_s + b4*c1m + b5*c2) + 
               g1*(t3 + t5*a_s)*(b2 + b3*a_s) + 2*g1*(t6 + t7*a_s)*(b2 + b3*a_s)*(g0 + g2*c1m +g3*c2) +
               g1^2*(t6 + t7*a_s)*(b2 + b3*a_s)*(a + a_s))*(a - a_s)

# C1 = Female
PIE_M1_f <- (g1*(t2 + t4*a_s) + g1*(t6 + t7*a_s)*(b0 + b1*a_s + b4*c1f + b5*c2) + 
               g1*(t3 + t5*a_s)*(b2 + b3*a_s) + 2*g1*(t6 + t7*a_s)*(b2 + b3*a_s)*(g0 + g2*c1f +g3*c2) +
               g1^2*(t6 + t7*a_s)*(b2 + b3*a_s)*(a + a_s))*(a - a_s)


# PIE_M2
# C1 = Male
PIE_M2_m <- (b1*(t3 + t5*a_s) + b1*(t6 + t7*a_s)*(g0 + g1*a_s + g2*c1m + g3*c2) +
               b3*(t3 + t5*a_s)*(g0 + g1*a_s + g2*c1m + g3*c2) +
               b3*(t6 + t7*a_s)*(mse + (g0 + g1*a_s + g2*c1m + g3*c2)^2))*(a - a_s)

# C1 = Female
PIE_M2_f <- (b1*(t3 + t5*a_s) + b1*(t6 + t7*a_s)*(g0 + g1*a_s + g2*c1f + g3*c2) +
               b3*(t3 + t5*a_s)*(g0 + g1*a_s + g2*c1f + g3*c2) +
               b3*(t6 + t7*a_s)*(mse + (g0 + g1*a_s + g2*c1f + g3*c2)^2))*(a - a_s)

# Total Effect TE
# C1 = Male
TE_m <- (t1 + t5*(b0 + b4*c1m + b5*c2) + b1*t3 + t4*(g0 + g2*c1m + g3*c2) + g1*t2 +
           t7*(b0 + b4*c1m + b5*c2)*(g0 + g2*c1m + g3*c2) + b1*t6*(g0 + g2*c1m + g3*c2) +
           g1*t6*(b0 + b4*c1m + b5*c2) + t5*b2*(g0 + g2*c1m + g3*c2) + t3*b3*(g0 + g2*c1m + g3*c2) +
           t3*b2*g1 + t7*b2*mse + t6*b3*mse + t7*b2*(g0 + g2*c1m + g3*c2)^2 + 
           t6*b3*(g0 + g2*c1m + g3*c2)^2 + 2*g1*t6*b2*(g0 + g2*c1m + g3*c2))*(a - a_s) + 
  (b1*t5 + g1*t4 + b1*t7*(g0 + g2*c1m + g3*c2) + 
     g1*t7*(b0 + b4*c1m + b5*c2) + g1*b1*t6 + t5*b3*(g0 + g2*c1m + g3*c2) + 
     t5*b2*g1 + t3*b3*g1 + t7*b3*mse + t7*b3*(g0 + g2*c1m + g3*c2)^2 +
     2*g1*t7*b2*(g0 + g2*c1m + g3*c2) + 2*g1*t6*b3*(g0 + g2*c1m + g3*c2) + t6*b2*g1^2)*(a^2 - a_s^2) +
  (g1*b1*t7 + t5*b3*g1 + 2*g1*t7*b3*(g0 + g2*c1m + g3*c2) + t7*b2*g1^2 + t6*b3*g1^2)*(a^3 - a_s^3) + 
  t7*b3*g1^2*(a^4 - a_s^4)

# C1 = Female
TE_f <- (t1 + t5*(b0 + b4*c1f + b5*c2) + b1*t3 + t4*(g0 + g2*c1f + g3*c2) + g1*t2 +
           t7*(b0 + b4*c1f + b5*c2)*(g0 + g2*c1f + g3*c2) + b1*t6*(g0 + g2*c1f + g3*c2) +
           g1*t6*(b0 + b4*c1f + b5*c2) + t5*b2*(g0 + g2*c1f + g3*c2) + t3*b3*(g0 + g2*c1f + g3*c2) +
           t3*b2*g1 + t7*b2*mse + t6*b3*mse + t7*b2*(g0 + g2*c1f + g3*c2)^2 + 
           t6*b3*(g0 + g2*c1f + g3*c2)^2 + 2*g1*t6*b2*(g0 + g2*c1f + g3*c2))*(a - a_s) + 
  (b1*t5 + g1*t4 + b1*t7*(g0 + g2*c1f + g3*c2) + 
     g1*t7*(b0 + b4*c1f + b5*c2) + g1*b1*t6 + t5*b3*(g0 + g2*c1f + g3*c2) + 
     t5*b2*g1 + t3*b3*g1 + t7*b3*mse + t7*b3*(g0 + g2*c1f + g3*c2)^2 +
     2*g1*t7*b2*(g0 + g2*c1f + g3*c2) + 2*g1*t6*b3*(g0 + g2*c1f + g3*c2) + t6*b2*g1^2)*(a^2 - a_s^2) +
  (g1*b1*t7 + t5*b3*g1 + 2*g1*t7*b3*(g0 + g2*c1f + g3*c2) + t7*b2*g1^2 + t6*b3*g1^2)*(a^3 - a_s^3) + 
  t7*b3*g1^2*(a^4 - a_s^4)

# The total effect obtained by adding all terms together
# C1 = Male
TTT_m <- CDE + INT_ref_AM1_m + INT_ref_AM2_AM1M2_m + 
  NAT_INT_AM1_m + NAT_INT_AM2_m + NAT_INT_AM1M2_m + NAT_INT_M1M2_m + PIE_M1_m + PIE_M2_m

# C1 = Female
TTT_f <- CDE + INT_ref_AM1_f + INT_ref_AM2_AM1M2_f + 
  NAT_INT_AM1_f + NAT_INT_AM2_f + NAT_INT_AM1M2_f + NAT_INT_M1M2_f + PIE_M1_f + PIE_M2_f

# The pure direct effect obtained by adding the corresponding components together
# C1 = Male
PDE_m <- CDE + INT_ref_AM1_m + INT_ref_AM2_AM1M2_m

# C1 = Female
PDE_f <- CDE + INT_ref_AM1_f + INT_ref_AM2_AM1M2_f


CDE
INT_ref_AM1_m
INT_ref_AM2_m
INT_ref_AM1M2_m

NAT_INT_AM1_m
NAT_INT_AM2_m
NAT_INT_AM1M2_m
NAT_INT_M1M2_m

PIE_M1_m
PIE_M2_m

PDE_m

TE_m
TTT_m


INT_ref_AM1_f
INT_ref_AM2_f
INT_ref_AM1M2_f

NAT_INT_AM1_f
NAT_INT_AM2_f
NAT_INT_AM1M2_f
NAT_INT_M1M2_f

PIE_M1_f
PIE_M2_f

PDE_f

TE_f
TTT_f

####################################################### Point Estimation Ends ################################

####################################################### Interval Estimation ##################################

# Draw confidence intervals by using bootstrap

# define a storage data frame 
# data frame for male 
Frame_m <- data.frame()

# data frame for female
Frame_f <- data.frame()

# number of iterations 
n <- 100000

# The results in the article can be reproduced with the seed below
set.seed(12345)

# This loop performs the bootstrap process
for(i in 1:n) {
  # sample from the original data with replacement
  b_data <- natural_a[sample(nrow(natural_a), nrow(natural_a), replace = TRUE), ]
  
  # model fitting 
  # outcome model
  b_model_Y <- lm(Y ~ A + M1 + M2 + A*M1 + A*M2 + M1*M2 + A*M1*M2 + C1 + C2, data = b_data)
  
  # M2 model
  b_model_M2 <- lm(M2 ~ A + M1 + A*M1 + C1 + C2, data = b_data)
  
  # M1 model
  b_model_M1 <- lm(M1 ~ A + C1 + C2, data = b_data)
  
  
  # assign preliminary values 
  # MSE of M1 model
  b_mse <- sum(b_model_M1$residuals^2)/(length(b_model_M1$residuals) - 4)
  
  # set the arbitrary reference level of M1 and M2 to be the mean
  b_m1 <- 29.2065 
  b_m2 <- 3.093471 
  
  # C1 is binary so it has two strata
  b_c1m <- 1
  b_c1f <- 0
  
  # use the mean value of C2 for the calculation 
  b_c2 <- 45.9565 
  
  # the reference level of A is 0 and the treatment level of A is 1
  b_a <- 1
  b_a_s <- 0
  
  # assign estimates of parameters from linear models 
  b_t0 <- unname(b_model_Y$coefficients[1])
  b_t1 <- unname(b_model_Y$coefficients[2])
  b_t2 <- unname(b_model_Y$coefficients[3])
  b_t3 <- unname(b_model_Y$coefficients[4])
  b_t4 <- unname(b_model_Y$coefficients[7])
  b_t5 <- unname(b_model_Y$coefficients[8])
  b_t6 <- unname(b_model_Y$coefficients[9])
  b_t7 <- unname(b_model_Y$coefficients[10])
  b_t8 <- unname(b_model_Y$coefficients[5])
  b_t9 <- unname(b_model_Y$coefficients[6])
  
  b_b0 <- unname(b_model_M2$coefficients[1])
  b_b1 <- unname(b_model_M2$coefficients[2])
  b_b2 <- unname(b_model_M2$coefficients[3])
  b_b3 <- unname(b_model_M2$coefficients[6])
  b_b4 <- unname(b_model_M2$coefficients[4])
  b_b5 <- unname(b_model_M2$coefficients[5])
  
  b_g0 <- unname(b_model_M1$coefficients[1])
  b_g1 <- unname(b_model_M1$coefficients[2])
  b_g2 <- unname(b_model_M1$coefficients[3])
  b_g3 <- unname(b_model_M1$coefficients[4])
  
  
  # calculate the estimates 
  # CDE(m1,m2)
  b_CDE <- (b_t1 + b_t4*b_m1 + b_t5*b_m2 + b_t7*b_m1*b_m2)*(b_a - b_a_s)
  
  # Referebce Interaction between A and M1
  # INT_ref_AM1(m1,m2)
  # C1 = Male
  b_INT_ref_AM1_m <- (b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2 - b_m1)*(b_t4 + b_t7*b_m2)*(b_a - b_a_s)
  # C1 = Female
  b_INT_ref_AM1_f <- (b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2 - b_m1)*(b_t4 + b_t7*b_m2)*(b_a - b_a_s)
  
  # Reference Interaction between A and M2
  # INT_ref_AM2(m1,m2)
  # C1 = Male
  b_INT_ref_AM2_m <- (b_t5 + b_t7*b_m1)*(b_b0 + b_b1*b_a_s + b_b2*b_m1 + b_b3*b_a_s*b_m1 + b_b4*b_c1m + b_b5*b_c2 - b_m2)*(b_a - b_a_s)
  # C1 = Female
  b_INT_ref_AM2_f <- (b_t5 + b_t7*b_m1)*(b_b0 + b_b1*b_a_s + b_b2*b_m1 + b_b3*b_a_s*b_m1 + b_b4*b_c1f + b_b5*b_c2 - b_m2)*(b_a - b_a_s)
  
  # Reference Interaction between A, M1 and M2
  # First we need to find Sum of two reference interaction effects
  # INT_ref_AM2+AM1M2(m2)
  # C1 = Male 
  b_INT_ref_AM2_AM1M2_m <- (b_t1 + b_t5*(b_b0 + b_b1*b_a_s + b_b4*b_c1m + b_b5*b_c2) + 
                              b_t7*(b_b0 + b_b1*b_a_s + b_b4*b_c1m + b_b5*b_c2)*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) +
                              b_t5*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) + 
                              b_t7*(b_b2 + b_b3*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2)^2) - 
                              (b_t1 + b_t5*b_m2) - b_t7*b_m2*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2))*(b_a - b_a_s)
  
  # C1 = Female
  b_INT_ref_AM2_AM1M2_f <- (b_t1 + b_t5*(b_b0 + b_b1*b_a_s + b_b4*b_c1f + b_b5*b_c2) + 
                              b_t7*(b_b0 + b_b1*b_a_s + b_b4*b_c1f + b_b5*b_c2)*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) +
                              b_t5*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) + 
                              b_t7*(b_b2 + b_b3*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2)^2) - 
                              (b_t1 + b_t5*b_m2) - b_t7*b_m2*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2))*(b_a - b_a_s)
  
  # next we subtract INT_ref_AM2 from INT_ref_AM2+AM1M2
  # INT_ref_AM1M2(m1,m2)
  # C1 = Male
  b_INT_ref_AM1M2_m <- b_INT_ref_AM2_AM1M2_m - b_INT_ref_AM2_m
  # C1 = Female
  b_INT_ref_AM1M2_f <- b_INT_ref_AM2_AM1M2_f - b_INT_ref_AM2_f
  
  # natural counterfactual interaction effects 
  # NAT_INT_AM1
  b_NAT_INT_AM1_m <- (b_t4*b_g1 + b_t7*b_g1*(b_b0 + b_b1*b_a_s + b_b4*b_c1m + b_b5*b_c2) + b_t5*b_g1*(b_b2 + b_b3*b_a_s) + 
                        2*b_t7*b_g1*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_c1m +b_g3*b_c2) + 
                        b_t7*b_g1^2*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
  b_NAT_INT_AM1_f <- (b_t4*b_g1 + b_t7*b_g1*(b_b0 + b_b1*b_a_s + b_b4*b_c1f + b_b5*b_c2) + b_t5*b_g1*(b_b2 + b_b3*b_a_s) + 
                        2*b_t7*b_g1*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_c1f +b_g3*b_c2) + 
                        b_t7*b_g1^2*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
  # NAT_INT_AM2
  b_NAT_INT_AM2_m <- (b_t5*b_b1 + b_t7*b_b1*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) + b_t5*b_b3*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) + 
                        b_t7*b_b3*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2)^2))*(b_a - b_a_s)^2
  
  b_NAT_INT_AM2_f <- (b_t5*b_b1 + b_t7*b_b1*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) + b_t5*b_b3*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) + 
                        b_t7*b_b3*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2)^2))*(b_a - b_a_s)^2
  
  # NAT_INT_AM1M2
  b_NAT_INT_AM1M2_m <- (b_t7*b_b1*b_g1 + b_t5*b_b3*b_g1 + 2*b_t7*b_b3*b_g1*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + 
                          b_t7*b_b3*b_g1^2*(b_a + b_a_s))*(b_a - b_a_s)^3
  
  b_NAT_INT_AM1M2_f <- (b_t7*b_b1*b_g1 + b_t5*b_b3*b_g1 + 2*b_t7*b_b3*b_g1*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + 
                          b_t7*b_b3*b_g1^2*(b_a + b_a_s))*(b_a - b_a_s)^3
  
  
  # NAT_INT_M1M2 
  b_NAT_INT_M1M2_m <- (b_b1*b_g1*(b_t6 + b_t7*b_a_s) + b_b3*b_g1*(b_t3 + b_t5*b_a_s) + 2*b_b3*b_g1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g2*b_c1m +b_g3*b_c2) + 
                         b_b3*b_g1^2*(b_t6 + b_t7*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
  b_NAT_INT_M1M2_f <- (b_b1*b_g1*(b_t6 + b_t7*b_a_s) + b_b3*b_g1*(b_t3 + b_t5*b_a_s) + 2*b_b3*b_g1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g2*b_c1f +b_g3*b_c2) + 
                         b_b3*b_g1^2*(b_t6 + b_t7*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
  # Pure Indirect Effects 
  
  # PIE_M1
  b_PIE_M1_m <- (b_g1*(b_t2 + b_t4*b_a_s) + b_g1*(b_t6 + b_t7*b_a_s)*(b_b0 + b_b1*b_a_s + b_b4*b_c1m + b_b5*b_c2) + 
                   b_g1*(b_t3 + b_t5*b_a_s)*(b_b2 + b_b3*b_a_s) + 2*b_g1*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_c1m + b_g3*b_c2) +
                   b_g1^2*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)
  
  b_PIE_M1_f <- (b_g1*(b_t2 + b_t4*b_a_s) + b_g1*(b_t6 + b_t7*b_a_s)*(b_b0 + b_b1*b_a_s + b_b4*b_c1f + b_b5*b_c2) + 
                   b_g1*(b_t3 + b_t5*b_a_s)*(b_b2 + b_b3*b_a_s) + 2*b_g1*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_c1f + b_g3*b_c2) +
                   b_g1^2*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)
  
  
  
  # PIE_M2
  b_PIE_M2_m <- (b_b1*(b_t3 + b_t5*b_a_s) + b_b1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) +
                   b_b3*(b_t3 + b_t5*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2) +
                   b_b3*(b_t6 + b_t7*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1m + b_g3*b_c2)^2))*(b_a - b_a_s)
  
  b_PIE_M2_f <- (b_b1*(b_t3 + b_t5*b_a_s) + b_b1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) +
                   b_b3*(b_t3 + b_t5*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2) +
                   b_b3*(b_t6 + b_t7*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_c1f + b_g3*b_c2)^2))*(b_a - b_a_s)
  
  
  # Total Effect TE
  
  b_TE_m <- (b_t1 + b_t5*(b_b0 + b_b4*b_c1m + b_b5*b_c2) + b_b1*b_t3 + b_t4*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + b_g1*b_t2 +
               b_t7*(b_b0 + b_b4*b_c1m + b_b5*b_c2)*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + b_b1*b_t6*(b_g0 + b_g2*b_c1m + b_g3*b_c2) +
               b_g1*b_t6*(b_b0 + b_b4*b_c1m + b_b5*b_c2) + b_t5*b_b2*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + b_t3*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2) +
               b_t3*b_b2*b_g1 + b_t7*b_b2*b_mse + b_t6*b_b3*b_mse + b_t7*b_b2*(b_g0 + b_g2*b_c1m + b_g3*b_c2)^2 + 
               b_t6*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2)^2 + 2*b_g1*b_t6*b_b2*(b_g0 + b_g2*b_c1m + b_g3*b_c2))*(b_a - b_a_s) + 
    (b_b1*b_t5 + b_g1*b_t4 + b_b1*b_t7*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + 
       b_g1*b_t7*(b_b0 + b_b4*b_c1m + b_b5*b_c2) + b_g1*b_b1*b_t6 + b_t5*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + 
       b_t5*b_b2*b_g1 + b_t3*b_b3*b_g1 + b_t7*b_b3*b_mse + b_t7*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2)^2 +
       2*b_g1*b_t7*b_b2*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + 2*b_g1*b_t6*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + b_t6*b_b2*b_g1^2)*(b_a^2 - b_a_s^2) +
    (b_g1*b_b1*b_t7 + b_t5*b_b3*b_g1 + 2*b_g1*b_t7*b_b3*(b_g0 + b_g2*b_c1m + b_g3*b_c2) + b_t7*b_b2*b_g1^2 + b_t6*b_b3*b_g1^2)*(b_a^3 - b_a_s^3) + 
    b_t7*b_b3*b_g1^2*(b_a^4 - b_a_s^4)
  
  b_TE_f <- (b_t1 + b_t5*(b_b0 + b_b4*b_c1f + b_b5*b_c2) + b_b1*b_t3 + b_t4*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + b_g1*b_t2 +
               b_t7*(b_b0 + b_b4*b_c1f + b_b5*b_c2)*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + b_b1*b_t6*(b_g0 + b_g2*b_c1f + b_g3*b_c2) +
               b_g1*b_t6*(b_b0 + b_b4*b_c1f + b_b5*b_c2) + b_t5*b_b2*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + b_t3*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2) +
               b_t3*b_b2*b_g1 + b_t7*b_b2*b_mse + b_t6*b_b3*b_mse + b_t7*b_b2*(b_g0 + b_g2*b_c1f + b_g3*b_c2)^2 + 
               b_t6*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2)^2 + 2*b_g1*b_t6*b_b2*(b_g0 + b_g2*b_c1f + b_g3*b_c2))*(b_a - b_a_s) + 
    (b_b1*b_t5 + b_g1*b_t4 + b_b1*b_t7*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + 
       b_g1*b_t7*(b_b0 + b_b4*b_c1f + b_b5*b_c2) + b_g1*b_b1*b_t6 + b_t5*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + 
       b_t5*b_b2*b_g1 + b_t3*b_b3*b_g1 + b_t7*b_b3*b_mse + b_t7*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2)^2 +
       2*b_g1*b_t7*b_b2*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + 2*b_g1*b_t6*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + b_t6*b_b2*b_g1^2)*(b_a^2 - b_a_s^2) +
    (b_g1*b_b1*b_t7 + b_t5*b_b3*b_g1 + 2*b_g1*b_t7*b_b3*(b_g0 + b_g2*b_c1f + b_g3*b_c2) + b_t7*b_b2*b_g1^2 + b_t6*b_b3*b_g1^2)*(b_a^3 - b_a_s^3) + 
    b_t7*b_b3*b_g1^2*(b_a^4 - b_a_s^4)
  
  
  
  b_TTT_m <- b_CDE + b_INT_ref_AM1_m + b_INT_ref_AM2_AM1M2_m + 
    b_NAT_INT_AM1_m + b_NAT_INT_AM2_m + b_NAT_INT_AM1M2_m + b_NAT_INT_M1M2_m + b_PIE_M1_m + b_PIE_M2_m
  
  b_TTT_f <- b_CDE + b_INT_ref_AM1_f + b_INT_ref_AM2_AM1M2_f + 
    b_NAT_INT_AM1_f + b_NAT_INT_AM2_f + b_NAT_INT_AM1M2_f + b_NAT_INT_M1M2_f + b_PIE_M1_f + b_PIE_M2_f
  
  b_PDE_m <- b_CDE + b_INT_ref_AM1_m + b_INT_ref_AM2_AM1M2_m
  
  b_PDE_f <- b_CDE + b_INT_ref_AM1_f + b_INT_ref_AM2_AM1M2_f
  
  
  Vec_m <- c(b_CDE, b_INT_ref_AM1_m, b_INT_ref_AM2_m, b_INT_ref_AM1M2_m,
             b_NAT_INT_AM1_m, b_NAT_INT_AM2_m, b_NAT_INT_AM1M2_m,
             b_NAT_INT_M1M2_m, b_PDE_m, b_PIE_M1_m, b_PIE_M2_m, b_TE_m, b_TTT_m)
  
  
  Vec_f <- c(b_CDE, b_INT_ref_AM1_f, b_INT_ref_AM2_f, b_INT_ref_AM1M2_f,
             b_NAT_INT_AM1_f, b_NAT_INT_AM2_f, b_NAT_INT_AM1M2_f,
             b_NAT_INT_M1M2_f, b_PDE_f, b_PIE_M1_f, b_PIE_M2_f, b_TE_f, b_TTT_f)
  
  Frame_m <- rbind.data.frame(Frame_m, Vec_m)
  
  Frame_f <- rbind.data.frame(Frame_f, Vec_f)
  
  if(i == n){
    names(Frame_m) <- c("CDE", "RefAM1_m", "RefAM2_m", "RefAM1M2_m", "NATAM1_m", "NATAM2_m", "NATAM1M2_m", "NATM1M2_m",
                        "PDE_m", "PIEM1_m", "PIEM2_m", "TEFormula_m", "TESum_m")
    names(Frame_f) <- c("CDE", "RefAM1_f", "RefAM2_f", "RefAM1M2_f", "NATAM1_f", "NATAM2_f", "NATAM1M2_f", "NATM1M2_f",
                        "PDE_f", "PIEM1_f", "PIEM2_f", "TEFormula_f", "TESum_f")
    
  }
  
}

# find the 2.5% and 97.5% quantile
# Males 
quantile(Frame_m$CDE, c(0.025, 0.975))

quantile(Frame_m$RefAM1_m, c(0.025, 0.975))

quantile(Frame_m$RefAM2_m, c(0.025, 0.975))

quantile(Frame_m$RefAM1M2_m, c(0.025, 0.975))

quantile(Frame_m$NATAM1_m, c(0.025, 0.975))

quantile(Frame_m$NATAM2_m, c(0.025, 0.975))

quantile(Frame_m$NATAM1M2_m, c(0.025, 0.975))

quantile(Frame_m$NATM1M2_m, c(0.025, 0.975))

quantile(Frame_m$PDE_m, c(0.025, 0.975))

quantile(Frame_m$PIEM1_m, c(0.025, 0.975))

quantile(Frame_m$PIEM2_m, c(0.025, 0.975))

quantile(Frame_m$TEFormula_m, c(0.025, 0.975))
################################################################
# Females
quantile(Frame_f$CDE, c(0.025, 0.975))

quantile(Frame_f$RefAM1_f, c(0.025, 0.975))

quantile(Frame_f$RefAM2_f, c(0.025, 0.975))

quantile(Frame_f$RefAM1M2_f, c(0.025, 0.975))

quantile(Frame_f$NATAM1_f, c(0.025, 0.975))

quantile(Frame_f$NATAM2_f, c(0.025, 0.975))

quantile(Frame_f$NATAM1M2_f, c(0.025, 0.975))

quantile(Frame_f$NATM1M2_f, c(0.025, 0.975))

quantile(Frame_f$PDE_f, c(0.025, 0.975))

quantile(Frame_f$PIEM1_f, c(0.025, 0.975))

quantile(Frame_f$PIEM2_f, c(0.025, 0.975))

quantile(Frame_f$TEFormula_f, c(0.025, 0.975))
##################################################################
########## Interval Estimation Ends ##############################