### SOLVING FOR RUN THRESHOLD ###


# preliminaries -----------------------------------------------------------
library(here)
library(deSolve)

# set model parameters ----------------------------------------------------------

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

# set parameters of the solution algorithm -----------------------

xstar_guess=0.1
xbig=2


# compute value function below threshold using analytical solution --------
# note: using only expressions needed to evaluate function at threshold
xstar=xstar_guess
a1=mu+delta-delta*rbar
a2=sigma^2/2
a3=-(phi+rho+theta*delta+delta)
a4=theta*delta*l+phi*(xstar<=1)
a5=delta+phi*(xstar>=1)
eta=1/(2*a2)*(a2-a1+sqrt((a2-a1)^2-4*a3*a2))
gamma=-1/(2*a2)*(a2-a1-sqrt((a2-a1)^2-4*a3*a2))
if(xstar<=1){
  A1=(1/rbar+a5/a3)*xstar^(-eta)+a4/(a3+a1)*xstar^(1-eta)
  FirstDerWxstar=A1*eta*xstar^(eta-1)-a4/(a3+a1)
  SecondDerWxstar=A1*eta*(eta-1)*xstar^(eta-2)
}
if(xstar>1){
  B2=phi/(gamma+eta)*((1-eta)/(a3+a1)+eta/a3)
  A2=(1/rbar+a5/a3)*xstar^(-eta)+a4/(a3+a1)*xstar^(1-eta)-B2*xstar^(-gamma-eta)
  B1=A2+phi/(gamma+eta)*(gamma/a3-(gamma+1)/(a3+a1))
  FirstDerWxstar=B1*eta*xstar^(eta-1)+B2*(-gamma)*xstar^(-gamma-1)-a4/(a3+a1)
  SecondDerWxstar=B1*eta*(eta-1)*xstar^(eta-2)+B2*(-gamma)*(-gamma-1)*xstar^(-gamma-2)
}


# solve for value function at large x numerically -------------------------

parameters <- c(mu=mu,
                delta=delta,
                sigma=sigma,
                phi=phi,
                rho=rho)
state<- c(W=1/rbar,
          Z=FirstDerWxstar)
derivatives <- function(x,state,parameters){
  with(as.list(c(state,parameters)), {
    dW<-Z
    dZ<- -2*(mu+delta)/sigma^2*Z/x+2*delta/sigma^2*Z/(x*W)-2*phi/sigma^2*min(1,x)/x^2+2*(rho+phi+delta)/sigma^2*W/x^2-2*delta/sigma^2*1/x^2
    list(c(dW,dZ))
  })
}


ode(y=state, times=c(xstar,xbig), func=derivatives, parms=parameters)

#compare compute W at xbig with analytical limit
