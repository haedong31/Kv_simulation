library(tidyverse)

V <- seq(-70, 50, by = 1)

ass_fn <- function(V) {return(1/(1+exp(-(V+22.5)/7.7)))}
ass <- ass_fn(V)
qplot(V, ass)

iss_fn <- function(V) {return(1/(1+exp((V+45.2)/5.7)))}
iss <- iss_fn(V)
qplot(V, iss)

tau_tis_fn <- function(V) {return(270+1050/(1+exp((V+45.2)/5.7)))}
tau_tis <- tau_tis_fn(V)
qplot(V, tau_tis)

tau_iur_fn <- function(V) {return(1200-170/(1+exp((V+45.2)/5.7)))}
tau_iur <- tau_iur_fn(V)
qplot(V, tau_iur)

tau_tas_fn <- function(V) {return(0.493*exp(-0.0628*V)+2.058)}
tau_tas <- tau_tas_fn(V)
qplot(V, tau_tas)

alpha_a_fn <- function(V) {return(0.18064*exp(0.03577*(V+70)))}
alpha_a <- alpha_a_fn(V)
qplot(V, alpha_a)

beta_a_fn <- function(V) {return(0.3956*exp(-0.06237*(V+30)))}
beta_a <- beta_a_fn(V)
qplot(V, beta_a)

alpha_i_fn <- function(V) {return((0.000152*exp(-(V+13.5)/7.0))/(0.067083*exp(-(V+33.5)/7.0)+1))}
alpha_i <- alpha_i_fn(V)
qplot(V, alpha_i)

beta_i_fn <- function(V) {return((0.00095*exp((V+70.5)/7.0))/(0.051335*exp((V+70.5)/7.0)+1))}
beta_i <- beta_i_fn(V)
qplot(V, beta_i)
