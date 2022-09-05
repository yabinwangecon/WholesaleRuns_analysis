### SOLVING FOR RUN THRESHOLD ###

# a problem in the parameters values?
# W at xstar seems too close to Wlimit: issue with Rbar ?

# preliminaries -----------------------------------------------------------
library(deSolve)
library(ggplot2)


# set model parameters ----------------------------------------------------------

rho=0.049
delta=1/0.101
phi=1/5.8
mu=0.049
theta=0.449
sigma=0.043
rbar=0.00598
Rbar=(rbar+rho)/(phi+delta)+1
alpha=0.92
xstar_predicted=1/92
l=alpha*phi/(rho+phi-mu)
Wlimit=(phi+delta)/(phi+delta+rho)
  

# set parameters of the solution algorithm -----------------------

xstar_guess=0.1
xbig=20


# compute value function below threshold using analytical solution --------
# note: using only expressions needed to evaluate function at threshold
xstar=xstar_guess
a1=mu+delta-delta*Rbar
a2=sigma^2/2
a3=-(phi+rho+theta*delta+delta)
a4=theta*delta*l+phi*(xstar<=1)
a5=delta+phi*(xstar>=1)
eta=1/(2*a2)*(a2-a1+sqrt((a2-a1)^2-4*a3*a2))
gamma=-1/(2*a2)*(a2-a1-sqrt((a2-a1)^2-4*a3*a2))
if(xstar<=1){
  A1=(1/Rbar+a5/a3)*xstar^(-eta)+a4/(a3+a1)*xstar^(1-eta)
  FirstDerWxstar=A1*eta*xstar^(eta-1)-a4/(a3+a1)
  SecondDerWxstar=A1*eta*(eta-1)*xstar^(eta-2)
}
if(xstar>1){
  B2=phi/(gamma+eta)*((1-eta)/(a3+a1)+eta/a3)
  A2=(1/Rbar+a5/a3)*xstar^(-eta)+a4/(a3+a1)*xstar^(1-eta)-B2*xstar^(-gamma-eta)
  B1=A2+phi/(gamma+eta)*(gamma/a3-(gamma+1)/(a3+a1))
  FirstDerWxstar=B1*eta*xstar^(eta-1)+B2*(-gamma)*xstar^(-gamma-1)-a4/(a3+a1)
  SecondDerWxstar=B1*eta*(eta-1)*xstar^(eta-2)+B2*(-gamma)*(-gamma-1)*xstar^(-gamma-2)
}


# compute value function at large x numerically -------------------------

parameters <- c(mu=mu,
                delta=delta,
                sigma=sigma,
                phi=phi,
                rho=rho)
state<- c(W=1/Rbar,
          Z=FirstDerWxstar)
derivatives <- function(x,state,parameters){
  with(as.list(c(state,parameters)), {
    dW<-Z
    dZ<- -(2*(mu+delta)*Z)/(sigma^2*x)+(2*delta*Z)/(sigma^2*x*W)-(2*phi*min(1,x))/(sigma^2*x^2)+(2*(rho+phi+delta)*W)/(sigma^2*x^2)-(2*delta)/(sigma^2*x^2)
    list(c(dW,dZ))
  })
}

points=seq(xstar, xbig, by=0.001 )
#points=c(xstar,xbig)
results<-ode(y=state, times=points, func=derivatives, parms=parameters)
#diagnostics(results)
#Wbig<-results[2,2]
results<-data.frame(results)
ggplot(results, aes(x=time, y=W)) +
  geom_line()
#compare compute W at xbig with analytical limit
