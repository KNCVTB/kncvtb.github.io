require(phaseR)

dy = function(t,y,parameters){
  k = parameters[1]
  k * y
}
expdecay  <- flowField(dy,xlim = c(0,2), ylim = c(0, 2),
                       parameters = -1,
                       points     = 21,
                       system     = "one.dim",
                       add        = FALSE, col = "blue")



### deSolve examples
require(deSolve)

parameters <- c(a = -8/3,
                b = -10,
                c =  28)


state <- c(X = 1,
           Y = 1,
           Z = 1)


Lorenz<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X + Y*Z
    dY <- b * (Y-Z)
    dZ <- -X*Y + c*Y - Z
    
    # return the rate of change
    list(c(dX, dY, dZ))
  })   # end with(as.list ...
}

times = seq(0,30,0.01)
out   = as.data.frame(ode(y = state, times = times, func = Lorenz, parms = parameters))
plot(x=times, y=out$X, type="l",col="blue",ylim=c(-50,50),ylab="x, y or z",xlab="t")
lines(x=times, y=out$Y, type="l",col="red",ylim=c(-50,50))
lines(x=times, y=out$Z, type="l",col="darkgreen",ylim=c(-50,50))

plot(x=out$X, y=out$Y, type="l",col="red",ylim=c(-50,50))
plot(x=out$X, y=out$Z, type="l",col="darkgreen",ylim=c(-50,50))

rm(list=ls())

### deSolve examples
require(deSolve)

parameters <- c(k11 = 0.01,
                k12 = 0.01)

state <- c(C1 = 1,
           C2 = 0)

two.compartments<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dC1 = -k11*C1 + k12*(C2-C1)
    dC2 = -k12*(C2-C1)
    # return the rate of change
    list(c(dC1, dC2))
  })   # end with(as.list ...
}

times = seq(0,500,1)
out   = as.data.frame(ode(y = state, times = times, func = two.compartments, parms = parameters))
plot(x=times, y=out$C1, type="l",col="blue",ylim=c(0,1),ylab="C1 or C2",xlab="t")
lines(x=times, y=out$C2, type="l",col="red",ylim=c(0,1))

# phase plot
plot(x=out$C1,y=out$C2, type="l",col="blue",ylim=c(0,1),ylab="C2",xlab="C1")

rm(list=ls())
require("deSolve")

# time unit : days
parameters <- c(b =  0.1, # 10% probability of transmission per contact
                c =  5,   #  5  (close) contacts per day
                d =  5)   #  5 days infectious


state <- c(S = 1000,
           I = 1,
           R = 0)


SIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    N = S+I+R
    dS = -b*c*I/N*S
    dI =  b*c*I/N*S - I/d
    dR =  I/d
    # return the rate of change
    list(c(dS, dI, dR))
  })   # end with(as.list ...
}

times = seq(0,60,0.02)
out   = as.data.frame(ode(y = state, times = times, func = SIR, parms = parameters))
plot(x=times, y=out$S, type="l",col="blue",ylim=c(0,1000),ylab="S, I or R",xlab="time in days")
lines(x=times, y=out$I, type="l",col="red",ylim=c(0,1000))
lines(x=times, y=out$R, type="l",col="darkgreen",ylim=c(0,1000))

plot(x=out$S, y=out$I, type="l",col="red",lwd=2, ylim=c(0,1000),ylab="I",xlab="S")

# after vaccination start with 50% in R and start with I=50
state <- c(S = 500,
           I = 50,
           R = 500)

times = seq(0,60,0.02)
out   = as.data.frame(ode(y = state, times = times, func = SIR, parms = parameters))
plot(x=out$S, y=out$I, type="l",col="red",lwd=2, ylim=c(0,1000),ylab="I",xlab="S")

# above: no epidemic !!

# no what if the virus was more contageous ...
parameters <- c(b =  0.5, # 50% probability of transmission per contact
                c =  5,   #  5  (close) contacts per day
                d =  5)   #  5 days infectious
state <- c(S = 500,
           I = 50,
           R = 500)

times = seq(0,60,0.02)
out   = as.data.frame(ode(y = state, times = times, func = SIR, parms = parameters))
plot(x=out$S, y=out$I, type="l",col="red",lwd=2, ylim=c(0,1000),ylab="I",xlab="S")

# obviously, we need to increase the fraction vaccinated...