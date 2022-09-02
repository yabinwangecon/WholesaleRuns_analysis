### SOLVING FOR RUN THRESHOLD ###


# preliminaries -----------------------------------------------------------
library(here)

# set parameters ----------------------------------------------------------

rho=0.049
delta=1/0.101
phi=1/5.8
mu=0.049
theta=0.449
sigma=0.043
rbar=0.00598
alpha=0.92
xstar_predicted=1/92
l=alpha*phi/(rho+phi-mu)
xstar_guess=
