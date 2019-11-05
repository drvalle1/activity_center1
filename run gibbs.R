rm(list=ls(all=TRUE))
set.seed(11)

#read important functions
setwd('U:\\GIT_models\\activity_center1')
library('Rcpp')
sourceCpp('aux1.cpp')
source('gibbs sampler.R')
source('gibbs functions.R')

#get data
dat=read.csv('fake data.csv',as.is=T)
grid.coord=read.csv('fake data grid.csv',as.is=T)

#basic setup
ngibbs=1000
nburn=ngibbs/2
n.ac=10
gamma1=0.1

#initial coordinates for activity centers (define this based on data instead of uninformative as below)
ac.coord.init=data.frame(x=runif(n.ac,min=0,max=100),
                         y=runif(n.ac,min=0,max=100))

#run gibbs sampler
options(warn=2)
res=gibbs.activity.center(dat=dat,grid.coord=grid.coord,n.ac=n.ac,
                          ac.coord.init=ac.coord.init,gamma1=gamma1)
