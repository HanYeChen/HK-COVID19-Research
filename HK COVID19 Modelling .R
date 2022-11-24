# packages
install.packages("deSolve")
install.packages("ggplot2")
install.packages("reshape2")
library(deSolve)
library(ggplot2)
library(reshape2)

# working directory & data file
setwd("C:/hkcovid")
mydata<-read.csv("covid19.csv", header=T)
if (FALSE){
   "The first row in mydata1 is extra data for the 0-th day (2022/05/31), 
    which is only used in numerical approximations and should not appear in model calibration.
    So the no. of days is the no. of rows -1"
}
nd<-nrow(mydata)-1
dcases<-mydata$daily_cases[2:nrow(mydata)]
daily_cases<-data.frame(date=1:nd,
                        dcases=dcases)

# function to compute fitted daily cases
pred<-function(sol,par,ran){
  rou<-par
  p<-rep(NA, length(ran))
  if (min(ran)==1){
    for (i in ran){
    p[i]<- (rou/(4*N))*(sol$S[sol$time==i]*sol$I[sol$time==i]
                        +2*sol$S[sol$time==i+1/2]*sol$I[sol$time==i+1/2]
                        +sol$S[sol$time==i+1]*sol$I[sol$time==i+1])
    }
  }else{
    for (i in ran){
      p[i-(min(ran)-1)]<- (rou/(4*N))*(sol$S[sol$time==i]*sol$I[sol$time==i]
                          +2*sol$S[sol$time==i+1/2]*sol$I[sol$time==i+1/2]
                          +sol$S[sol$time==i+1]*sol$I[sol$time==i+1])  
    }}
  
  return(p)
}

# Preliminaries
# vaccination rate (numerical approximations)
v<-data.frame(v=mydata$v)
# an artificial time-series
signal<-data.frame(times=seq(from=1, to=nrow(mydata), by=1.0), v)
# interpolating function
input<-approxfun(signal, rule=2)
# modified SIR model (with time-dependent vaccination impact)
mSIR_model<-function(time, variables, parameters){
  with(as.list(c(variables, parameters)),{
    v<-input(time)
    dSdt<- -rou*I*S/N-theta*v
    dIdt<- rou*I*S/N-gamma*I
    dRdt<- gamma*I+theta*v
    return(list(c(dSdt,dIdt,dRdt)))
  })
}
# conventional SIR model
cSIR_model<-function(time, variables, parameters){
  with(as.list(c(variables, parameters)),{
    dSdt<- -rou*I*S/N
    dIdt<- rou*I*S/N-gamma*I
    dRdt<- gamma*I
    return(list(c(dSdt,dIdt,dRdt)))
  })
}


## One-compartment model
# Calibration 
# Least squares estimation
N<-7413070
times0<-seq(from=1, to=max(mydata$time), by=0.5)
SIR_SSE0<-function(para){
  rou<-para[1]
  gamma<-para[2]
  theta<-para[3]
  R0<-(1212699-10357)+theta*6337282
  I0<-10357
  S0<-N-I0-R0
  solutions<-as.data.frame(ode(y=c(S=S0,I=I0,R=R0), 
                               times=times0,
                               func=mSIR_model,
                               parms=c(rou=rou, gamma=gamma, theta=theta)))
  pcases<-pred(solutions, rou, 1:nd)
  SSE<-sum((pcases-dcases)^2)
  return(SSE)
}
op0<-optim(par=c(0.2,0.1,0.5),
           fn=SIR_SSE0,
           method="Nelder-Mead",
           control=list(maxit=1E5, trace=2))
op0

# Model Fitness
SSTa<-sum((dcases-mean(dcases))^2)
SSEa<-SIR_SSE0(op0$par)
r2a=1-SSEa/SSTa

# Fitted Model
rou0=op0$par[1]
gamma0=op0$par[2]
theta0=op0$par[3]
R0<-(1212699-10357)+theta0*6337282
I0<-10357
S0<-N-I0-R0
initials0<-c(S=S0, I=I0, R=R0)
solutions0<-as.data.frame(ode(y=initials0, 
                              times=times0,
                              func=mSIR_model,
                              parms=c(rou=rou0, gamma=gamma0, theta=theta0)))
pdaily_cases0<-data.frame(date=1:nd,
                          pcases=pred(solutions0, rou0, 1:nd))

# Plot fitted daily cases against reported daily cases
ggplot()+
  geom_point(data=daily_cases,
             mapping=aes(x=date, y=dcases),colour="black")+
  geom_line(data=pdaily_cases0, 
            mapping=aes(x=date, y=pcases), colour="red")+
  labs(title="Fitted daily cases v.s. Reported daily cases",
       x="date",
       y="cases")

# Plot infection dynamics
ggplot()+
  geom_line(data=solutions0, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions0, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions0, mapping=aes(x=time, y=R), colour="green")+
  labs(title="Infection Dynamics",
       x="time",
       y="no. of people")


## Two-compartment model
#COMPARTMENT I
times1a<-seq(from=1, to=125, by=0.5)
dcases1<-dcases[1:124]
theta<-theta0
I0<-10357
R0<-(1212699-10357)+theta*6337282
S0<-N-I0-R0
initials1a<-c(S=S0, I=I0, R=R0)

SIR_SSE1<-function(para){
    rou<-para[1]
    gamma<-para[2]
    solutions<-as.data.frame(ode(y=initials1a, 
                                 times=times1a,
                                 func=mSIR_model,
                                 parms=c(rou=rou, gamma=gamma, theta=theta)))
    pcases1<-pred(solutions, rou, 1:124)
    SSE<-sum((pcases1-dcases1)^2)
    return(SSE)
}
op1<-optim(par=c(0.2,0.1),
           fn=SIR_SSE1,
           method="Nelder-Mead",
           control=list(maxit=1E5, trace=2))
          
op1

SSTb<-sum((dcases1-mean(dcases1))^2)
SSEb<-SIR_SSE1(op1$par)
r2b=1-SSEb/SSTb

# Fitted Model
rou1=op1$par[1]
gamma1=op1$par[2]
solutions1a<-as.data.frame(ode(y=initials1a, 
                              times=times1a,
                              func=mSIR_model,
                              parms=c(rou=rou1, gamma=gamma1, theta=theta)))
pdaily_cases1a<-data.frame(date=1:124,
                          pcases=pred(solutions1a, rou1, 1:124))

#COMPARTMENT II
dcases2<-dcases[125:nd]
times1b<-seq(from=125, to=max(mydata$time), by=0.5)
initials1b<-c(S=solutions1a$S[solutions1a$time==125],
              I=solutions1a$I[solutions1a$time==125],
              R=solutions1a$R[solutions1a$time==125])

SIR_SSE2<-function(para){
  rou<-para[1]
  gamma<-para[2]
  solutions<-as.data.frame(ode(y=initials1b, 
                               times=times1b,
                               func=mSIR_model,
                               parms=c(rou=rou, gamma=gamma, theta=theta)))
  pcases2<-pred(solutions, rou, 125:nd)
  SSE<-sum((pcases2-dcases2)^2)
  return(SSE)
}
op2<-optim(par=c(0.2,0.1),
           fn=SIR_SSE2,
           method="Nelder-Mead",
           control=list(maxit=1E5, trace=2))
op2

SSTc<-sum((dcases2-mean(dcases2))^2)
SSEc<-SIR_SSE2(op2$par)
r2c=1-SSEc/SSTc

rou2=op2$par[1]
gamma2=op2$par[2]
solutions1b<-as.data.frame(ode(y=initials1b, 
                               times=times1b,
                               func=mSIR_model,
                               parms=c(rou=rou2, gamma=gamma2, theta=theta)))
pdaily_cases1b<-data.frame(date=125:nd,
                           pcases=pred(solutions1b, rou2, 125:nd))

# Plot fitted daily cases against reported daily cases
ggplot()+
  geom_point(data=daily_cases,
             mapping=aes(x=date, y=dcases),colour="black")+
  geom_line(data=pdaily_cases1a, 
            mapping=aes(x=date, y=pcases), colour="red")+
  geom_line(data=pdaily_cases1b, 
            mapping=aes(x=date, y=pcases), colour="blue")+
  labs(title="Fitted daily cases v.s. Reported daily cases",
       x="date",
       y="cases")

# Plot Infection Dynamics (S(t), I(t), R(t))
ggplot()+
  geom_line(data=solutions1a, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=R), colour="green")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=R), colour="green")+
  labs(title="Infection Dynamics",
       x="time",
       y="no. of people")

## Projections
# Projection
times2<-seq(from=max(mydata$time), to=270, by=0.5)
initials2<-c(S=solutions1b$S[solutions1b$time==max(mydata$time)],
             I=solutions1b$I[solutions1b$time==max(mydata$time)],
             R=solutions1b$R[solutions1b$time==max(mydata$time)])
solutions2<-as.data.frame(ode(y=initials2, 
                             times=times2,
                             func=mSIR_model,
                             parms=c(rou=rou2, gamma=gamma2, theta=theta)))
pdaily_cases2<-data.frame(date=(nd+1):269,
                          pcases=pred(solutions2, rou2, (nd+1):269))

# Plot projected infection dynamics
ggplot()+
  geom_line(data=solutions1a, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=R), colour="green")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=R), colour="green")+
  geom_line(data=solutions2, mapping=aes(x=time, y=S), colour="purple")+
  geom_line(data=solutions2, mapping=aes(x=time, y=I), colour="purple")+
  geom_line(data=solutions2, mapping=aes(x=time, y=R), colour="purple")+
  labs(title="Infection Dynamics with Projections", x="time", y="no. of people")

# Plot projected daily cases
solutions1a_extd<-as.data.frame(ode(y=initials1b, 
                               times=seq(from=125, to=270, by=0.5),
                               func=mSIR_model,
                               parms=c(rou=rou1, gamma=gamma1, theta=theta)))
pdaily_cases1a_extd<-data.frame(date=125:259,
                           pcases=pred(solutions1a_extd, rou1, 125:259))
ggplot()+
  geom_point(data=daily_cases,
             mapping=aes(x=date, y=dcases),colour="black")+
  geom_line(data=pdaily_cases1a, 
            mapping=aes(x=date, y=pcases), colour="red")+
  geom_line(data=pdaily_cases1a_extd, 
            mapping=aes(x=date, y=pcases), colour="purple")+
  geom_line(data=pdaily_cases1b, 
            mapping=aes(x=date, y=pcases), colour="blue")+
  geom_line(data=pdaily_cases2,
            mapping=aes(x=date, y=pcases), colour="green")+
  labs(title="Daily Cases with Projections",
       x="date",
       y="cases")

# Plot effective reproduction number (starting from 2022/10/01)
re1<-data.frame(time=solutions1b$time,
                Re=rou2*solutions1b$S/(gamma2*N))
re2<-data.frame(time=solutions2$time,
                Re=rou2*solutions2$S/(gamma2*N))
ggplot()+
  geom_line(data=re1, mapping=aes(x=time, y=Re), colour="blue")+
  geom_line(data=re2, mapping=aes(x=time, y=Re), colour="green")+
  scale_y_continuous(limits=c(0,3))+
  labs(title="Real-time effective reproduction number Re(t)")

# Plot projected prevalence
preva1a<-data.frame(time=times1a, 
                    prevalence=(solutions1a$I/N)*100)
preva1b<-data.frame(time=times1b, 
                    prevalence=(solutions1b$I/N)*100)
preva2<-data.frame(time=times2,
                   prevalence=(solutions2$I/N)*100)
ggplot()+
  geom_line(data=preva1a, mapping=aes(x=time, y=prevalence), colour="red")+
  geom_line(data=preva1b, mapping=aes(x=time, y=prevalence), colour="blue")+
  geom_line(data=preva2, mapping=aes(x=time, y=prevalence), colour="green")+
  labs(title="Real-time Prevalence with Projections (%)",
       x="time",
       y="prevalence(%)")

## Vaccination Impacts
# initial values
# population
N<-7413070
# currently infected (on Jun 1, 2022)
I0<-10357
# removed (deaths & recovered)
R0<-(1212699-10357)
# susceptible
S0<-N-I0-R0
initials3<-c(S=S0, I=I0, R=R0)
times3<-times0

solutions3<-as.data.frame(ode(y=initials3, 
                              times=times3,
                              func=cSIR_model,
                              parms=c(rou=rou0, gamma=gamma0)))
pdaily_cases3<-data.frame(date=1:nd,
                          pcases=pred(solutions3, rou0, 1:nd))

# Plot simulated daily cases
ggplot()+
  geom_point(data=daily_cases,
             mapping=aes(x=date, y=dcases),colour="black")+
  geom_line(data=pdaily_cases1a, 
            mapping=aes(x=date, y=pcases), colour="red")+
  geom_line(data=pdaily_cases1b, 
            mapping=aes(x=date, y=pcases), colour="blue")+
  geom_line(data=pdaily_cases3, mapping=aes(x=date, y=pcases), colour="purple")+
  labs(title="Fitted Daily Cases with v.s. without Vaccination",
       x="date",
       y="cases")

# Plot simulated infection dynamics
ggplot()+
  geom_line(data=solutions3, mapping=aes(x=time, y=S), colour="purple")+
  geom_line(data=solutions3, mapping=aes(x=time, y=I), colour="purple")+
  geom_line(data=solutions3, mapping=aes(x=time, y=R), colour="purple")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=S), colour="red")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=I), colour="red")+
  geom_line(data=solutions1a, mapping=aes(x=time, y=R), colour="red")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=S), colour="blue")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=I), colour="blue")+
  geom_line(data=solutions1b, mapping=aes(x=time, y=R), colour="blue")+
  labs(title="Infection Dynamics with v.s. without Vaccination",
       x="time",
       y="no. of people")

# Plot simulated real-time prevalence
preva3<-data.frame(time=times0,
                   prevalence=(solutions3$I/N)*100)
ggplot()+
    geom_line(data=preva1a, mapping=aes(x=time, y=prevalence), colour="red")+
    geom_line(data=preva1b, mapping=aes(x=time, y=prevalence), colour="blue")+
    geom_line(data=preva3, mapping=aes(x=time, y=prevalence), colour="purple")+
    labs(title="Real-time Prevalence with v.s. without Vaccination (%)",
         x="time",
         y="prevalence(%)")






























