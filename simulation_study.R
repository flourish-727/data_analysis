# This R script is for a simulation study 
# to make a comparason between 
# Bellavia's and Valeri's decomposition (BV method) and our proposed decomposition 
# in a non-sequential two-mediator scenario 

setwd("D:/JCI/Again_JCI/New_Sub/Simulation_BV_Our")

#########################################################################
############################### Generate Data for Simulation ############
#########################################################################

# Assign true parameters 

# Outcome model 
true_t0 <- 0.2
true_t1 <- 0.3
true_t2 <- 0.3
true_t3 <- 0.4
true_t4 <- 0.01
true_t5 <- 0.02
true_t6 <- 0.6
true_t7 <- 0.7
true_t8 <- 0.2

# M2 model
true_b0 <- 0.2
true_b1 <- 0.3
true_b4 <- 0.2

# M1 model
true_g0 <- 0.2
true_g1 <- 0.3
true_g2 <- 0.2

# the true sigma 
# all models use the same true sigma value 
true_sig <- 0.5

# intercept of exposure A on C
alpha0 <- 0.3

# mean value of confounding coovariate C
true_C <- 0.2

# number of observations
n <- 1000

# the results can be restored by using the seed below
set.seed(54321)

# generate data set for C
C <- rnorm(n, true_C, true_sig)

# generate data set for exposure A 
# based on the relationship A = alpha0 + 3*C + N(0, 0.5)
A <- rnorm(n, alpha0 + 3*C, true_sig)

# generate data set for M2
M2 <- rnorm(n, true_b0 + true_b1*A + true_b4*C, true_sig)

# generate data set for M1
M1 <- rnorm(n, true_g0 + true_g1*A + true_g2*C, true_sig)

# generate data set for Y 
Y <- rnorm(n, true_t0 + true_t1*A + true_t2*M1 + true_t3*M2 + 
             true_t4*A*M1 + true_t5*A*M2 + true_t6*M1*M2 + true_t7*A*M1*M2 + true_t8*C, true_sig)

# form up the data frame for data analysis 
Data <- data.frame(A, M1, M2, Y, C)

#########################################################################
############################### Our proposed method          ############
#########################################################################

#########################################################################
############################### Calculate the True Effects   ############
#########################################################################

# The formulas are special cases of those for a sequential two-mediator scenario
# by setting beta_2 and Beta_3 to zero

# a and a star
a <- 1
a_s <- 0

# set the arbitrary reference levels of M1 and M3 to be zero 
m1 <- 0
m2 <- 0

# set the beta_2 and beta_3 in the M1 model to be zero 
# since this is a non-sequential two-mediator scenario 
true_b2 <- 0
true_b3 <- 0

# CDE(m1,m2)
true_CDE <- (true_t1 + true_t4*m1 + true_t5*m2 + true_t7*m1*m2)*(a - a_s)

# Reference Interaction between A and M1
# INT_ref_AM1(m1,m2)
true_INT_ref_AM1 <- (true_g0 + true_g1*a_s + true_g2*true_C - m1)*(true_t4 + true_t7*m2)*(a - a_s)

# Reference Interaction between A and M2
# INT_ref_AM2(m1,m2)
true_INT_ref_AM2 <- (true_t5 + true_t7*m1)*
                    (true_b0 + true_b1*a_s 
                     + true_b2*m1 + true_b3*a_s*m1 + true_b4*true_C - m2)*(a - a_s)

# Reference Interaction between A, M1 and M2
# First need to find the Sum of two reference interaction effects
# INT_ref_AM2+AM1M2(m2)
true_INT_ref_AM2_AM1M2 <- (true_t1 + true_t5*(true_b0 + true_b1*a_s + true_b4*true_C) + 
                          true_t7*(true_b0 + true_b1*a_s + true_b4*true_C)*(true_g0 + true_g1*a_s + true_g2*true_C) +
                          true_t5*(true_b2 + true_b3*a_s)*(true_g0 + true_g1*a_s + true_g2*true_C) + 
                          true_t7*(true_b2 + true_b3*a_s)*(true_sig^2 + (true_g0 + true_g1*a_s + true_g2*true_C)^2) - 
                          (true_t1 + true_t5*m2) - true_t7*m2*(true_g0 + true_g1*a_s + true_g2*true_C))*(a - a_s)

# Next we subtract INT_ref_AM2 from INT_ref_AM2_AM1M2

# INT_ref_AM1M2(m1,m2)
true_INT_ref_AM1M2 <- true_INT_ref_AM2_AM1M2 - true_INT_ref_AM2

# natural counterfactual interaction effects 
# NAT_INT_AM1
true_NAT_INT_AM1 <- (true_t4*true_g1 + true_t7*true_g1*(true_b0 + true_b1*a_s + true_b4*true_C) + true_t5*true_g1*(true_b2 + true_b3*a_s) + 
                    2*true_t7*true_g1*(true_b2 + true_b3*a_s)*(true_g0 + true_g2*true_C) + 
                    true_t7*true_g1^2*(true_b2 + true_b3*a_s)*(a + a_s))*(a - a_s)^2

# NAT_INT_AM2
true_NAT_INT_AM2 <- (true_t5*true_b1 + true_t7*true_b1*(true_g0 + true_g1*a_s + true_g2*true_C) + true_t5*true_b3*(true_g0 + true_g1*a_s + true_g2*true_C) + 
                    true_t7*true_b3*(true_sig^2 + (true_g0 + true_g1*a_s + true_g2*true_C)^2))*(a - a_s)^2

# NAT_INT_AM1M2
true_NAT_INT_AM1M2 <- (true_t7*true_b1*true_g1 + true_t5*true_b3*true_g1 + 2*true_t7*true_b3*true_g1*(true_g0 + true_g2*true_C) + 
                      true_t7*true_b3*true_g1^2*(a + a_s))*(a - a_s)^3

# NAT_INT_M1M2 
true_NAT_INT_M1M2 <- (true_b1*true_g1*(true_t6 + true_t7*a_s) + true_b3*true_g1*(true_t3 + true_t5*a_s) + 2*true_b3*true_g1*(true_t6 + true_t7*a_s)*(true_g0 + true_g2*true_C) + 
                      true_b3*true_g1^2*(true_t6 + true_t7*a_s)*(a + a_s))*(a - a_s)^2

# Pure Indirect Effects 

# PIE_M1
true_PIE_M1 <- (true_g1*(true_t2 + true_t4*a_s) + true_g1*(true_t6 + true_t7*a_s)*(true_b0 + true_b1*a_s + true_b4*true_C) + 
                true_g1*(true_t3 + true_t5*a_s)*(true_b2 + true_b3*a_s) + 2*true_g1*(true_t6 + true_t7*a_s)*(true_b2 + true_b3*a_s)*(true_g0 + true_g2*true_C) +
                true_g1^2*(true_t6 + true_t7*a_s)*(true_b2 + true_b3*a_s)*(a + a_s))*(a - a_s)

# PIE_M2
true_PIE_M2 <- (true_b1*(true_t3 + true_t5*a_s) + true_b1*(true_t6 + true_t7*a_s)*(true_g0 + true_g1*a_s + true_g2*true_C) +
                true_b3*(true_t3 + true_t5*a_s)*(true_g0 + true_g1*a_s + true_g2*true_C) +
                true_b3*(true_t6 + true_t7*a_s)*(true_sig^2 + (true_g0 + true_g1*a_s + true_g2*true_C)^2))*(a - a_s)

# Total Effect TE
true_TE <- (true_t1 + true_t5*(true_b0 + true_b4*true_C) + true_b1*true_t3 + true_t4*(true_g0 + true_g2*true_C) + true_g1*true_t2 +
            true_t7*(true_b0 + true_b4*true_C)*(true_g0 + true_g2*true_C) + true_b1*true_t6*(true_g0 + true_g2*true_C) +
            true_g1*true_t6*(true_b0 + true_b4*true_C) + true_t5*true_b2*(true_g0 + true_g2*true_C) + true_t3*true_b3*(true_g0 + true_g2*true_C) +
            true_t3*true_b2*true_g1 + true_t7*true_b2*true_sig^2 + true_t6*true_b3*true_sig^2 + true_t7*true_b2*(true_g0 + true_g2*true_C)^2 + 
            true_t6*true_b3*(true_g0 + true_g2*true_C)^2 + 2*true_g1*true_t6*true_b2*(true_g0 + true_g2*true_C))*(a - a_s) + 
  (true_b1*true_t5 + true_g1*true_t4 + true_b1*true_t7*(true_g0 + true_g2*true_C) + 
     true_g1*true_t7*(true_b0 + true_b4*true_C) + true_g1*true_b1*true_t6 + true_t5*true_b3*(true_g0 + true_g2*true_C) + 
     true_t5*true_b2*true_g1 + true_t3*true_b3*true_g1 + true_t7*true_b3*true_sig^2 + true_t7*true_b3*(true_g0 + true_g2*true_C)^2 +
     2*true_g1*true_t7*true_b2*(true_g0 + true_g2*true_C) + 2*true_g1*true_t6*true_b3*(true_g0 + true_g2*true_C) + true_t6*true_b2*true_g1^2)*(a^2 - a_s^2) +
  (true_g1*true_b1*true_t7 + true_t5*true_b3*true_g1 + 2*true_g1*true_t7*true_b3*(true_g0 + true_g2*true_C) + true_t7*true_b2*true_g1^2 + true_t6*true_b3*true_g1^2)*(a^3 - a_s^3) + 
  true_t7*true_b3*true_g1^2*(a^4 - a_s^4)

true_PDE <- true_CDE + true_INT_ref_AM1 + true_INT_ref_AM2 + true_INT_ref_AM1M2

true_TTT <- true_PDE + true_NAT_INT_AM1 + true_NAT_INT_AM2 + true_NAT_INT_AM1M2 + 
  true_NAT_INT_M1M2 + true_PIE_M1 + true_PIE_M2

true_six_terms <- true_TE - true_PDE

true_CDE
true_INT_ref_AM1
true_INT_ref_AM2
true_INT_ref_AM1M2
true_NAT_INT_AM1
true_NAT_INT_AM2
true_NAT_INT_AM1M2
true_NAT_INT_M1M2
true_PIE_M1
true_PIE_M2

true_PDE
true_TE
true_TTT

true_six_terms

#########################################################################
############################### BV method  ##############################
#########################################################################

#########################################################################
############################### Calculate the True Effects   ############
#########################################################################

# The two methods have the same CDE, INTrefAM1, INTrefAM2 and INTrefAM1M2
# therefore there is no need to calculate these four effects 

# First find term 1-12 

true_term_1 <- (true_t0 + true_t1*a + true_t8*true_C +
             (true_t2 + true_t4*a)*(true_g0 + true_g1*a + true_g2*true_C)) + 
  (true_t3 + true_t5*a + (true_t6 + true_t7*a)*(true_g0 + true_g1*a + true_g2*true_C))*
  (true_b0 + true_b1*a + true_b4*true_C)


true_term_2 <- (true_t0 + true_t1*a + true_t2*m1 + true_t4*a*m1 + true_t8*true_C) +
  (true_t3 + true_t5*a + true_t6*m1 + true_t7*a*m1)*(true_b0 + true_b1*a + true_b4*true_C)


true_term_3 <- (true_t0 + true_t1*a_s + true_t8*true_C +
             (true_t2 + true_t4*a_s)*(true_g0 + true_g1*a + true_g2*true_C)) + 
  (true_t3 + true_t5*a_s + (true_t6 + true_t7*a_s)*(true_g0 + true_g1*a + true_g2*true_C))*
  (true_b0 + true_b1*a + true_b4*true_C)


true_term_4 <- (true_t0 + true_t1*a + true_t3*m2 + true_t5*a*m2 + true_t8*true_C) + 
  (true_t2 + true_t4*a + true_t6*m2 + true_t7*a*m2)*(true_g0 + true_g1*a + true_g2*true_C)


true_term_5 <- (true_t0 + true_t1*a_s + true_t3*m2 + true_t5*a_s*m2 + true_t8*true_C) + 
  (true_t2 + true_t4*a_s + true_t6*m2 + true_t7*a_s*m2)*(true_g0 + true_g1*a + true_g2*true_C)


true_term_6 <- (true_t0 + true_t1*a_s + true_t2*m1 + true_t4*a_s*m1 + true_t8*true_C) +
  (true_t3 + true_t5*a_s + true_t6*m1 + true_t7*a_s*m1)*(true_b0 + true_b1*a + true_b4*true_C)


true_term_7 <- (true_t0 + true_t1*a + true_t8*true_C +
             (true_t2 + true_t4*a)*(true_g0 + true_g1*a_s + true_g2*true_C)) + 
  (true_t3 + true_t5*a + (true_t6 + true_t7*a)*(true_g0 + true_g1*a_s + true_g2*true_C))*
  (true_b0 + true_b1*a_s + true_b4*true_C)



true_term_8 <- (true_t0 + true_t1*a + true_t2*m1 + true_t4*a*m1 + true_t8*true_C) +
  (true_t3 + true_t5*a + true_t6*m1 + true_t7*a*m1)*(true_b0 + true_b1*a_s + true_b4*true_C)


true_term_9 <- (true_t0 + true_t1*a_s + true_t8*true_C +
             (true_t2 + true_t4*a_s)*(true_g0 + true_g1*a_s + true_g2*true_C)) + 
  (true_t3 + true_t5*a_s + (true_t6 + true_t7*a_s)*(true_g0 + true_g1*a_s + true_g2*true_C))*
  (true_b0 + true_b1*a_s + true_b4*true_C)


true_term_10 <- (true_t0 + true_t1*a + true_t3*m2 + true_t5*a*m2 + true_t8*true_C) + 
  (true_t2 + true_t4*a + true_t6*m2 + true_t7*a*m2)*(true_g0 + true_g1*a_s + true_g2*true_C)


true_term_11 <- (true_t0 + true_t1*a_s + true_t3*m2 + true_t5*a_s*m2 + true_t8*true_C) + 
  (true_t2 + true_t4*a_s + true_t6*m2 + true_t7*a_s*m2)*(true_g0 + true_g1*a_s + true_g2*true_C)


true_term_12 <- (true_t0 + true_t1*a_s + true_t2*m1 + true_t4*a_s*m1 + true_t8*true_C) +
  (true_t3 + true_t5*a_s + true_t6*m1 + true_t7*a_s*m1)*(true_b0 + true_b1*a_s + true_b4*true_C)


bv_true_INT_MED_AM1 <- true_term_4 - true_term_5 - true_term_10 + true_term_11 


bv_true_INT_MED_AM2 <- true_term_2 - true_term_6 - true_term_8 + true_term_12


bv_true_INT_MED_AM1M2 <- true_term_1 - true_term_2 - true_term_3 - true_term_4 + true_term_5 + true_term_6 -
  true_term_7 + true_term_8 + true_term_9 + true_term_10 - true_term_11 - true_term_12

bv_true_PNIE_M1 <- true_term_5 - true_term_11

bv_true_PNIE_M2 <- true_term_6 - true_term_12

bv_true_PNIE_M1M2 <- true_term_3 - true_term_5 - true_term_6 - true_term_9 + true_term_11 + true_term_12


bv_true_six_terms <- bv_true_INT_MED_AM1 + bv_true_INT_MED_AM2 + bv_true_INT_MED_AM1M2 + 
  bv_true_PNIE_M1 + bv_true_PNIE_M2 + bv_true_PNIE_M1M2



bv_true_INT_MED_AM1

bv_true_INT_MED_AM2

bv_true_INT_MED_AM1M2

bv_true_PNIE_M1

bv_true_PNIE_M2

bv_true_PNIE_M1M2

bv_true_six_terms


#########################################################################
############################### Our proposed method          ############
#########################################################################

############################ Point estimation ###########################
attach(Data)

# outcome model
model_Y <- lm(Y ~ A + M1 + M2 + A*M1 + A*M2 + M1*M2 + A*M1*M2 + C, data = Data)

# M2 model
model_M2 <- lm(M2 ~ A + C, data = Data)

# M1 model
model_M1 <- lm(M1 ~ A + C, data = Data)

# calculate estimates 

# MSE of M1 model
mse <- sum(model_M1$residuals^2)/(length(model_M1$residuals)-3)

# set the arbitrary reference level of M1 and M2 to be zero
m1 <- 0
m2 <- 0

# use the estimated mean value of C for the calculation 
cc <- mean(C)

# the reference level of A is 0 and the treatment level of A is 1
#a <- 1
#a_s <- 0

# estimates of coefficients from the outcome model
t0 <- unname(model_Y$coefficients[1])
t1 <- unname(model_Y$coefficients[2])
t2 <- unname(model_Y$coefficients[3])
t3 <- unname(model_Y$coefficients[4])
t4 <- unname(model_Y$coefficients[6])
t5 <- unname(model_Y$coefficients[7])
t6 <- unname(model_Y$coefficients[8])
t7 <- unname(model_Y$coefficients[9])
t8 <- unname(model_Y$coefficients[5])

# estimates of coefficients from the M2 model 
# b2 and b3 are zero
b0 <- unname(model_M2$coefficients[1])
b1 <- unname(model_M2$coefficients[2])
b2 <- 0
b3 <- 0
b4 <- unname(model_M2$coefficients[3])

# estimates of coefficients from the M1 model 
g0 <- unname(model_M1$coefficients[1])
g1 <- unname(model_M1$coefficients[2])
g2 <- unname(model_M1$coefficients[3])

# CDE(m1,m2)
CDE <- (t1 + t4*m1 + t5*m2 + t7*m1*m2)*(a - a_s)

# Reference Interaction between A and M1
# INT_ref_AM1(m1,m2)
INT_ref_AM1 <- (g0 + g1*a_s + g2*cc - m1)*(t4 + t7*m2)*(a - a_s)

# Reference Interaction between A and M2
# INT_ref_AM2(m1,m2)
INT_ref_AM2 <- (t5 + t7*m1)*(b0 + b1*a_s + b2*m1 + b3*a_s*m1 + b4*cc - m2)*(a - a_s)

# Reference Interaction between A, M1 and M2
# First need to find the Sum of two reference interaction effects
# INT_ref_AM2+AM1M2(m2)
INT_ref_AM2_AM1M2 <- (t1 + t5*(b0 + b1*a_s + b4*cc) + 
                          t7*(b0 + b1*a_s + b4*cc)*(g0 + g1*a_s + g2*cc) +
                          t5*(b2 + b3*a_s)*(g0 + g1*a_s + g2*cc) + 
                          t7*(b2 + b3*a_s)*(mse + (g0 + g1*a_s + g2*cc)^2) - 
                          (t1 + t5*m2) - t7*m2*(g0 + g1*a_s + g2*cc))*(a - a_s)


# Next we subtract INT_ref_AM2 from INT_ref_AM2_AM1M2

# INT_ref_AM1M2(m1,m2)
INT_ref_AM1M2 <- INT_ref_AM2_AM1M2 - INT_ref_AM2


# natural counterfactual interaction effects 
# NAT_INT_AM1
NAT_INT_AM1 <- (t4*g1 + t7*g1*(b0 + b1*a_s + b4*cc) + t5*g1*(b2 + b3*a_s) + 
                    2*t7*g1*(b2 + b3*a_s)*(g0 + g2*cc) + 
                    t7*g1^2*(b2 + b3*a_s)*(a + a_s))*(a - a_s)^2

# NAT_INT_AM2
NAT_INT_AM2 <- (t5*b1 + t7*b1*(g0 + g1*a_s + g2*cc) + t5*b3*(g0 + g1*a_s + g2*cc) + 
                    t7*b3*(mse + (g0 + g1*a_s + g2*cc)^2))*(a - a_s)^2

# NAT_INT_AM1M2
NAT_INT_AM1M2 <- (t7*b1*g1 + t5*b3*g1 + 2*t7*b3*g1*(g0 + g2*cc) + 
                      t7*b3*g1^2*(a + a_s))*(a - a_s)^3

# NAT_INT_M1M2 
NAT_INT_M1M2 <- (b1*g1*(t6 + t7*a_s) + b3*g1*(t3 + t5*a_s) + 2*b3*g1*(t6 + t7*a_s)*(g0 + g2*cc) + 
                     b3*g1^2*(t6 + t7*a_s)*(a + a_s))*(a - a_s)^2


# Pure Indirect Effects 

# PIE_M1
PIE_M1 <- (g1*(t2 + t4*a_s) + g1*(t6 + t7*a_s)*(b0 + b1*a_s + b4*cc) + 
               g1*(t3 + t5*a_s)*(b2 + b3*a_s) + 2*g1*(t6 + t7*a_s)*(b2 + b3*a_s)*(g0 + g2*cc) +
               g1^2*(t6 + t7*a_s)*(b2 + b3*a_s)*(a + a_s))*(a - a_s)



# PIE_M2
PIE_M2 <- (b1*(t3 + t5*a_s) + b1*(t6 + t7*a_s)*(g0 + g1*a_s + g2*cc) +
               b3*(t3 + t5*a_s)*(g0 + g1*a_s + g2*cc) +
               b3*(t6 + t7*a_s)*(mse + (g0 + g1*a_s + g2*cc)^2))*(a - a_s)


# Total Effect TE
TE <- (t1 + t5*(b0 + b4*cc) + b1*t3 + t4*(g0 + g2*cc) + g1*t2 +
           t7*(b0 + b4*cc)*(g0 + g2*cc) + b1*t6*(g0 + g2*cc) +
           g1*t6*(b0 + b4*cc) + t5*b2*(g0 + g2*cc) + t3*b3*(g0 + g2*cc) +
           t3*b2*g1 + t7*b2*mse + t6*b3*mse + t7*b2*(g0 + g2*cc)^2 + 
           t6*b3*(g0 + g2*cc)^2 + 2*g1*t6*b2*(g0 + g2*cc))*(a - a_s) + 
  (b1*t5 + g1*t4 + b1*t7*(g0 + g2*cc) + 
     g1*t7*(b0 + b4*cc) + g1*b1*t6 + t5*b3*(g0 + g2*cc) + 
     t5*b2*g1 + t3*b3*g1 + t7*b3*mse + t7*b3*(g0 + g2*cc)^2 +
     2*g1*t7*b2*(g0 + g2*cc) + 2*g1*t6*b3*(g0 + g2*cc) + t6*b2*g1^2)*(a^2 - a_s^2) +
  (g1*b1*t7 + t5*b3*g1 + 2*g1*t7*b3*(g0 + g2*cc) + t7*b2*g1^2 + t6*b3*g1^2)*(a^3 - a_s^3) + 
  t7*b3*g1^2*(a^4 - a_s^4)

# The pure direct effect obtained by adding the corresponding components together
PDE <- CDE + INT_ref_AM1 + INT_ref_AM2 + INT_ref_AM1M2


# The total effect obtained by adding ten components together
TTT <- CDE + INT_ref_AM1 + INT_ref_AM2 + INT_ref_AM1M2 +
  NAT_INT_AM1 + NAT_INT_AM2 + NAT_INT_AM1M2 + NAT_INT_M1M2 + PIE_M1 + PIE_M2

six_terms <- TE - PDE


CDE
INT_ref_AM1
INT_ref_AM2
INT_ref_AM1M2
NAT_INT_AM1
NAT_INT_AM2
NAT_INT_AM1M2
NAT_INT_M1M2
PIE_M1
PIE_M2

PDE
TE
TTT

six_terms

#########################################################################
############################### BV method  ##############################
#########################################################################

############################ Point estimation ###########################

# The two methods have the same CDE, INTrefAM1, INTrefAM2 and INTrefAM1M2
# therefore there is no need to calculate these four effects 

# First find term 1-12 

term_1 <- (t0 + t1*a + t8*cc +
             (t2 + t4*a)*(g0 + g1*a + g2*cc)) + 
  (t3 + t5*a + (t6 + t7*a)*(g0 + g1*a + g2*cc))*
  (b0 + b1*a + b4*cc)


term_2 <- (t0 + t1*a + t2*m1 + t4*a*m1 + t8*cc) +
  (t3 + t5*a + t6*m1 + t7*a*m1)*(b0 + b1*a + b4*cc)


term_3 <- (t0 + t1*a_s + t8*cc +
             (t2 + t4*a_s)*(g0 + g1*a + g2*cc)) + 
  (t3 + t5*a_s + (t6 + t7*a_s)*(g0 + g1*a + g2*cc))*
  (b0 + b1*a + b4*cc)


term_4 <- (t0 + t1*a + t3*m2 + t5*a*m2 + t8*cc) + 
  (t2 + t4*a + t6*m2 + t7*a*m2)*(g0 + g1*a + g2*cc)


term_5 <- (t0 + t1*a_s + t3*m2 + t5*a_s*m2 + t8*cc) + 
  (t2 + t4*a_s + t6*m2 + t7*a_s*m2)*(g0 + g1*a + g2*cc)


term_6 <- (t0 + t1*a_s + t2*m1 + t4*a_s*m1 + t8*cc) +
  (t3 + t5*a_s + t6*m1 + t7*a_s*m1)*(b0 + b1*a + b4*cc)


term_7 <- (t0 + t1*a + t8*cc +
             (t2 + t4*a)*(g0 + g1*a_s + g2*cc)) + 
  (t3 + t5*a + (t6 + t7*a)*(g0 + g1*a_s + g2*cc))*
  (b0 + b1*a_s + b4*cc)



term_8 <- (t0 + t1*a + t2*m1 + t4*a*m1 + t8*cc) +
  (t3 + t5*a + t6*m1 + t7*a*m1)*(b0 + b1*a_s + b4*cc)


term_9 <- (t0 + t1*a_s + t8*cc +
             (t2 + t4*a_s)*(g0 + g1*a_s + g2*cc)) + 
  (t3 + t5*a_s + (t6 + t7*a_s)*(g0 + g1*a_s + g2*cc))*
  (b0 + b1*a_s + b4*cc)


term_10 <- (t0 + t1*a + t3*m2 + t5*a*m2 + t8*cc) + 
  (t2 + t4*a + t6*m2 + t7*a*m2)*(g0 + g1*a_s + g2*cc)


term_11 <- (t0 + t1*a_s + t3*m2 + t5*a_s*m2 + t8*cc) + 
  (t2 + t4*a_s + t6*m2 + t7*a_s*m2)*(g0 + g1*a_s + g2*cc)


term_12 <- (t0 + t1*a_s + t2*m1 + t4*a_s*m1 + t8*cc) +
  (t3 + t5*a_s + t6*m1 + t7*a_s*m1)*(b0 + b1*a_s + b4*cc)


bv_INT_MED_AM1 <- term_4 - term_5 - term_10 + term_11 


bv_INT_MED_AM2 <- term_2 - term_6 - term_8 + term_12


bv_INT_MED_AM1M2 <- term_1 - term_2 - term_3 - term_4 + term_5 + term_6 -
  term_7 + term_8 + term_9 + term_10 - term_11 - term_12

bv_PNIE_M1 <- term_5 - term_11

bv_PNIE_M2 <- term_6 - term_12

bv_PNIE_M1M2 <- term_3 - term_5 - term_6 - term_9 + term_11 + term_12

bv_six_terms <- bv_INT_MED_AM1 + bv_INT_MED_AM2 + bv_INT_MED_AM1M2 + 
  bv_PNIE_M1 + bv_PNIE_M2 + bv_PNIE_M1M2


bv_INT_MED_AM1

bv_INT_MED_AM2

bv_INT_MED_AM1M2

bv_PNIE_M1

bv_PNIE_M2

bv_PNIE_M1M2

bv_six_terms



#########################################################################
############################### Our proposed method          ############
#########################################################################
#################################### and ################################
#########################################################################
############################### BV method  ##############################
#########################################################################

############################ Interval estimation ########################

# Draw confidence intervals by using bootstrap

# define two storage data frames 
Frame_our <- data.frame()
Frame_bv <- data.frame()

# number of iterations 
n <- 100000

# The results in the article can be reproduced with the seed below
set.seed(12345)


# This loop performs the bootstrap process
for(i in 1:n) {
  # sample from the original data with replacement
  b_data <- Data[sample(nrow(Data), nrow(Data), replace = TRUE), ]
  
  # model fitting 
  # outcome model
  b_model_Y <- lm(Y ~ A + M1 + M2 + A*M1 + A*M2 + M1*M2 + A*M1*M2 + C, data = b_data)
  
  # M2 model
  b_model_M2 <- lm(M2 ~ A + C, data = b_data)
  
  # M1 model
  b_model_M1 <- lm(M1 ~ A + C, data = b_data)
  
  
  # assign preliminary values 
  # MSE of M1 model
  b_mse <- sum(b_model_M1$residuals^2)/(length(b_model_M1$residuals) - 3)
  
  # set the arbitrary reference level of M1 and M2 to be zero
  b_m1 <- 0
  b_m2 <- 0
  
  
  # use the estimated mean value of C for the calculation 
  b_cc <- 0.1700977 
  
  # the reference level of A is 0 and the treatment level of A is 1
  b_a <- 1
  b_a_s <- 0
  
  # assign estimates of parameters from linear models 
  b_t0 <- unname(b_model_Y$coefficients[1])
  b_t1 <- unname(b_model_Y$coefficients[2])
  b_t2 <- unname(b_model_Y$coefficients[3])
  b_t3 <- unname(b_model_Y$coefficients[4])
  b_t4 <- unname(b_model_Y$coefficients[6])
  b_t5 <- unname(b_model_Y$coefficients[7])
  b_t6 <- unname(b_model_Y$coefficients[8])
  b_t7 <- unname(b_model_Y$coefficients[9])
  b_t8 <- unname(b_model_Y$coefficients[5])

  
  b_b0 <- unname(b_model_M2$coefficients[1])
  b_b1 <- unname(b_model_M2$coefficients[2])
  b_b2 <- 0
  b_b3 <- 0
  b_b4 <- unname(b_model_M2$coefficients[3])

  
  b_g0 <- unname(b_model_M1$coefficients[1])
  b_g1 <- unname(b_model_M1$coefficients[2])
  b_g2 <- unname(b_model_M1$coefficients[3])

  
  
  # calculate the estimates for our decomposition
  # CDE(m1,m2)
  b_CDE <- (b_t1 + b_t4*b_m1 + b_t5*b_m2 + b_t7*b_m1*b_m2)*(b_a - b_a_s)
  
  # Referebce Interaction between A and M1
  # INT_ref_AM1(m1,m2)
  b_INT_ref_AM1 <- (b_g0 + b_g1*b_a_s + b_g2*b_cc - b_m1)*(b_t4 + b_t7*b_m2)*(b_a - b_a_s)
  
  
  # Reference Interaction between A and M2
  # INT_ref_AM2(m1,m2)
  b_INT_ref_AM2 <- (b_t5 + b_t7*b_m1)*(b_b0 + b_b1*b_a_s + b_b2*b_m1 + b_b3*b_a_s*b_m1 + b_b4*b_cc - b_m2)*(b_a - b_a_s)
  
  
  # Reference Interaction between A, M1 and M2
  # First we need to find Sum of two reference interaction effects
  # INT_ref_AM2+AM1M2(m2)
  b_INT_ref_AM2_AM1M2 <- (b_t1 + b_t5*(b_b0 + b_b1*b_a_s + b_b4*b_cc) + 
                              b_t7*(b_b0 + b_b1*b_a_s + b_b4*b_cc)*(b_g0 + b_g1*b_a_s + b_g2*b_cc) +
                              b_t5*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_cc) + 
                              b_t7*(b_b2 + b_b3*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_cc)^2) - 
                              (b_t1 + b_t5*b_m2) - b_t7*b_m2*(b_g0 + b_g1*b_a_s + b_g2*b_cc))*(b_a - b_a_s)
  
  
  # next we subtract INT_ref_AM2 from INT_ref_AM2+AM1M2
  # INT_ref_AM1M2(m1,m2)
  b_INT_ref_AM1M2 <- b_INT_ref_AM2_AM1M2 - b_INT_ref_AM2
  
  
  # natural counterfactual interaction effects 
  # NAT_INT_AM1
  b_NAT_INT_AM1 <- (b_t4*b_g1 + b_t7*b_g1*(b_b0 + b_b1*b_a_s + b_b4*b_cc) + b_t5*b_g1*(b_b2 + b_b3*b_a_s) + 
                        2*b_t7*b_g1*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_cc) + 
                        b_t7*b_g1^2*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
 
  # NAT_INT_AM2
  b_NAT_INT_AM2 <- (b_t5*b_b1 + b_t7*b_b1*(b_g0 + b_g1*b_a_s + b_g2*b_cc) + b_t5*b_b3*(b_g0 + b_g1*b_a_s + b_g2*b_cc) + 
                        b_t7*b_b3*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_cc)^2))*(b_a - b_a_s)^2
  
  
  # NAT_INT_AM1M2
  b_NAT_INT_AM1M2 <- (b_t7*b_b1*b_g1 + b_t5*b_b3*b_g1 + 2*b_t7*b_b3*b_g1*(b_g0 + b_g2*b_cc) + 
                          b_t7*b_b3*b_g1^2*(b_a + b_a_s))*(b_a - b_a_s)^3
  
  
  # NAT_INT_M1M2 
  b_NAT_INT_M1M2 <- (b_b1*b_g1*(b_t6 + b_t7*b_a_s) + b_b3*b_g1*(b_t3 + b_t5*b_a_s) + 2*b_b3*b_g1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g2*b_cc) + 
                         b_b3*b_g1^2*(b_t6 + b_t7*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)^2
  
  
  # Pure Indirect Effects 
  
  # PIE_M1
  b_PIE_M1 <- (b_g1*(b_t2 + b_t4*b_a_s) + b_g1*(b_t6 + b_t7*b_a_s)*(b_b0 + b_b1*b_a_s + b_b4*b_cc) + 
                   b_g1*(b_t3 + b_t5*b_a_s)*(b_b2 + b_b3*b_a_s) + 2*b_g1*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_g0 + b_g2*b_cc) +
                   b_g1^2*(b_t6 + b_t7*b_a_s)*(b_b2 + b_b3*b_a_s)*(b_a + b_a_s))*(b_a - b_a_s)
  
  
  # PIE_M2
  b_PIE_M2 <- (b_b1*(b_t3 + b_t5*b_a_s) + b_b1*(b_t6 + b_t7*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_cc) +
                   b_b3*(b_t3 + b_t5*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_cc) +
                   b_b3*(b_t6 + b_t7*b_a_s)*(b_mse + (b_g0 + b_g1*b_a_s + b_g2*b_cc)^2))*(b_a - b_a_s)
  
  
  # Total Effect TE
  b_TE <- (b_t1 + b_t5*(b_b0 + b_b4*b_cc) + b_b1*b_t3 + b_t4*(b_g0 + b_g2*b_cc) + b_g1*b_t2 +
               b_t7*(b_b0 + b_b4*b_cc)*(b_g0 + b_g2*b_cc) + b_b1*b_t6*(b_g0 + b_g2*b_cc) +
               b_g1*b_t6*(b_b0 + b_b4*b_cc) + b_t5*b_b2*(b_g0 + b_g2*b_cc) + b_t3*b_b3*(b_g0 + b_g2*b_cc) +
               b_t3*b_b2*b_g1 + b_t7*b_b2*b_mse + b_t6*b_b3*b_mse + b_t7*b_b2*(b_g0 + b_g2*b_cc)^2 + 
               b_t6*b_b3*(b_g0 + b_g2*b_cc)^2 + 2*b_g1*b_t6*b_b2*(b_g0 + b_g2*b_cc))*(b_a - b_a_s) + 
    (b_b1*b_t5 + b_g1*b_t4 + b_b1*b_t7*(b_g0 + b_g2*b_cc) + 
       b_g1*b_t7*(b_b0 + b_b4*b_cc) + b_g1*b_b1*b_t6 + b_t5*b_b3*(b_g0 + b_g2*b_cc) + 
       b_t5*b_b2*b_g1 + b_t3*b_b3*b_g1 + b_t7*b_b3*b_mse + b_t7*b_b3*(b_g0 + b_g2*b_cc)^2 +
       2*b_g1*b_t7*b_b2*(b_g0 + b_g2*b_cc) + 2*b_g1*b_t6*b_b3*(b_g0 + b_g2*b_cc) + b_t6*b_b2*b_g1^2)*(b_a^2 - b_a_s^2) +
    (b_g1*b_b1*b_t7 + b_t5*b_b3*b_g1 + 2*b_g1*b_t7*b_b3*(b_g0 + b_g2*b_cc) + b_t7*b_b2*b_g1^2 + b_t6*b_b3*b_g1^2)*(b_a^3 - b_a_s^3) + 
    b_t7*b_b3*b_g1^2*(b_a^4 - b_a_s^4)
  
 
  b_PDE <- b_CDE + b_INT_ref_AM1 + b_INT_ref_AM2 + b_INT_ref_AM1M2
  
  
  b_TTT <- b_CDE + b_INT_ref_AM1 + b_INT_ref_AM2 + b_INT_ref_AM1M2 +
    b_NAT_INT_AM1 + b_NAT_INT_AM2 + b_NAT_INT_AM1M2 + b_NAT_INT_M1M2 + b_PIE_M1 + b_PIE_M2
  
  
  b_six_terms <- b_TE - b_PDE
  
  #######################################################################################
  # calculate the effects for BV method
  # The two methods have the same CDE, INTrefAM1, INTrefAM2 and INTrefAM1M2
  # therefore there is no need to calculate these four effects 
  
  # calculate the 12 terms
  b_term_1 <- (b_t0 + b_t1*b_a + b_t8*b_cc +
               (b_t2 + b_t4*b_a)*(b_g0 + b_g1*b_a + b_g2*b_cc)) + 
    (b_t3 + b_t5*b_a + (b_t6 + b_t7*b_a)*(b_g0 + b_g1*b_a + b_g2*b_cc))*
    (b_b0 + b_b1*b_a + b_b4*b_cc)
  
  
  b_term_2 <- (b_t0 + b_t1*b_a + b_t2*b_m1 + b_t4*b_a*b_m1 + b_t8*b_cc) +
    (b_t3 + b_t5*b_a + b_t6*b_m1 + b_t7*b_a*b_m1)*(b_b0 + b_b1*b_a + b_b4*b_cc)
  
  
  b_term_3 <- (b_t0 + b_t1*b_a_s + b_t8*b_cc +
               (b_t2 + b_t4*b_a_s)*(b_g0 + b_g1*b_a + b_g2*b_cc)) + 
    (b_t3 + b_t5*b_a_s + (b_t6 + b_t7*b_a_s)*(b_g0 + b_g1*b_a + b_g2*b_cc))*
    (b_b0 + b_b1*b_a + b_b4*b_cc)
  
  
  b_term_4 <- (b_t0 + b_t1*b_a + b_t3*b_m2 + b_t5*b_a*b_m2 + b_t8*b_cc) + 
    (b_t2 + b_t4*b_a + b_t6*b_m2 + b_t7*b_a*b_m2)*(b_g0 + b_g1*b_a + b_g2*b_cc)
  
  
  b_term_5 <- (b_t0 + b_t1*b_a_s + b_t3*b_m2 + b_t5*b_a_s*b_m2 + b_t8*b_cc) + 
    (b_t2 + b_t4*b_a_s + b_t6*b_m2 + b_t7*b_a_s*b_m2)*(b_g0 + b_g1*b_a + b_g2*b_cc)
  
  
  b_term_6 <- (b_t0 + b_t1*b_a_s + b_t2*b_m1 + b_t4*b_a_s*b_m1 + b_t8*b_cc) +
    (b_t3 + b_t5*b_a_s + b_t6*b_m1 + b_t7*b_a_s*b_m1)*(b_b0 + b_b1*b_a + b_b4*b_cc)
  
  
  b_term_7 <- (b_t0 + b_t1*b_a + b_t8*b_cc +
               (b_t2 + b_t4*b_a)*(b_g0 + b_g1*b_a_s + b_g2*b_cc)) + 
    (b_t3 + b_t5*b_a + (b_t6 + b_t7*b_a)*(b_g0 + b_g1*b_a_s + b_g2*b_cc))*
    (b_b0 + b_b1*b_a_s + b_b4*b_cc)
  
  
  
  b_term_8 <- (b_t0 + b_t1*b_a + b_t2*b_m1 + b_t4*b_a*b_m1 + b_t8*b_cc) +
    (b_t3 + b_t5*b_a + b_t6*b_m1 + b_t7*b_a*b_m1)*(b_b0 + b_b1*b_a_s + b_b4*b_cc)
  
  
  b_term_9 <- (b_t0 + b_t1*b_a_s + b_t8*b_cc +
               (b_t2 + b_t4*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_cc)) + 
    (b_t3 + b_t5*b_a_s + (b_t6 + b_t7*b_a_s)*(b_g0 + b_g1*b_a_s + b_g2*b_cc))*
    (b_b0 + b_b1*b_a_s + b_b4*b_cc)
  
  
  b_term_10 <- (b_t0 + b_t1*b_a + b_t3*b_m2 + b_t5*b_a*b_m2 + b_t8*b_cc) + 
    (b_t2 + b_t4*b_a + b_t6*b_m2 + b_t7*b_a*b_m2)*(b_g0 + b_g1*b_a_s + b_g2*b_cc)
  
  
  b_term_11 <- (b_t0 + b_t1*b_a_s + b_t3*b_m2 + b_t5*b_a_s*b_m2 + b_t8*b_cc) + 
    (b_t2 + b_t4*b_a_s + b_t6*b_m2 + b_t7*b_a_s*b_m2)*(b_g0 + b_g1*b_a_s + b_g2*b_cc)
  
  
  b_term_12 <- (b_t0 + b_t1*b_a_s + b_t2*b_m1 + b_t4*b_a_s*b_m1 + b_t8*b_cc) +
    (b_t3 + b_t5*b_a_s + b_t6*b_m1 + b_t7*b_a_s*b_m1)*(b_b0 + b_b1*b_a_s + b_b4*b_cc)
  
  
  b_bv_INT_MED_AM1 <- b_term_4 - b_term_5 - b_term_10 + b_term_11 

  
  b_bv_INT_MED_AM2 <- b_term_2 - b_term_6 - b_term_8 + b_term_12
  
  
  b_bv_INT_MED_AM1M2 <- b_term_1 - b_term_2 - b_term_3 - b_term_4 + b_term_5 + b_term_6 -
    b_term_7 + b_term_8 + b_term_9 + b_term_10 - b_term_11 - b_term_12
  
  b_bv_PNIE_M1 <- b_term_5 - b_term_11
  
  b_bv_PNIE_M2 <- b_term_6 - b_term_12
  
  b_bv_PNIE_M1M2 <- b_term_3 - b_term_5 - b_term_6 - b_term_9 + b_term_11 + b_term_12
  
  b_bv_six_terms <- b_bv_INT_MED_AM1 + b_bv_INT_MED_AM2 + b_bv_INT_MED_AM1M2 + 
    b_bv_PNIE_M1 + b_bv_PNIE_M2 + b_bv_PNIE_M1M2
  
  
  
  ####################################################################################
  ############### store the results ##################################################
  ####################################################################################
  
  
  Vec_our <- c(b_CDE, b_INT_ref_AM1, b_INT_ref_AM2, b_INT_ref_AM1M2,
             b_NAT_INT_AM1, b_NAT_INT_AM2, b_NAT_INT_AM1M2,
             b_NAT_INT_M1M2, b_PIE_M1, b_PIE_M2, b_PDE, b_TE, b_TTT, b_six_terms)
  
  Vec_bv <- c(b_bv_INT_MED_AM1, b_bv_INT_MED_AM2,
              b_bv_INT_MED_AM1M2, b_bv_PNIE_M1, b_bv_PNIE_M2, b_bv_PNIE_M1M2, b_bv_six_terms)
  
  
  Frame_our <- rbind.data.frame(Frame_our, Vec_our)
  
  Frame_bv <- rbind.data.frame(Frame_bv, Vec_bv) 
  
  if(i == n){
    names(Frame_our) <- c("CDE", "RefAM1", "RefAM2", "RefAM1M2", "NATAM1", "NATAM2", "NATAM1M2", "NATM1M2",
                        "PIEM1", "PIEM2", "PDE", "TEFormula", "TESum", "Sum_6_our")
    
    names(Frame_bv) <- c("INTmed_AM1", "INTmed_AM2",
                         "INTmed_AM1M2", "PNIE_M1", "PNIE_M2", "PNIE_M1M2", "Sum_6_bv")
    
  }
  
}

# find the 2.5% and 97.5% quantile for our decomposition
quantile(Frame_our$CDE, c(0.025, 0.975))

quantile(Frame_our$RefAM1, c(0.025, 0.975))

quantile(Frame_our$RefAM2, c(0.025, 0.975))

quantile(Frame_our$RefAM1M2, c(0.025, 0.975))

quantile(Frame_our$NATAM1, c(0.025, 0.975))

quantile(Frame_our$NATAM2, c(0.025, 0.975))

quantile(Frame_our$NATAM1M2, c(0.025, 0.975))

quantile(Frame_our$NATM1M2, c(0.025, 0.975))

quantile(Frame_our$PIEM1, c(0.025, 0.975))

quantile(Frame_our$PIEM2, c(0.025, 0.975))

quantile(Frame_our$PDE, c(0.025, 0.975))

quantile(Frame_our$TEFormula, c(0.025, 0.975))

quantile(Frame_our$TESum, c(0.025, 0.975))

quantile(Frame_our$Sum_6_our, c(0.025, 0.975))


# find the 2.5% and 97.5% quantile for bv method 
quantile(Frame_bv$INTmed_AM1, c(0.025, 0.975))

quantile(Frame_bv$INTmed_AM2, c(0.025, 0.975))

quantile(Frame_bv$INTmed_AM1M2, c(0.025, 0.975))

quantile(Frame_bv$PNIE_M1, c(0.025, 0.975))

quantile(Frame_bv$PNIE_M2, c(0.025, 0.975))

quantile(Frame_bv$PNIE_M1M2, c(0.025, 0.975))

quantile(Frame_bv$Sum_6_bv, c(0.025, 0.975))




