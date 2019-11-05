#calculates part of the marginal loglikelihood after integrating our the z's
get.calc.mloglik=function(dist.mat.sel,phi,log.theta,
                          n.tsegm,n.ac,n.grid,
                          dat){
  prob=exp(-phi*dist.mat.sel)
  lprob=log(prob/rowSums(prob))
  
  res=matrix(NA,n.tsegm,n.ac)
  for (i in 1:n.ac){
    lprob.mat=matrix(lprob[i,],n.tsegm,n.grid,byrow=T)
    res[,i]=rowSums(dat*lprob.mat)+log.theta[i]
  }
  res
}
#----------------------------------------------------
sample.coord=function(ac.coord,log.theta,phi,dist.mat,
                      n.ac,n.grid,n.tsegm,
                      jump,dat,grid.coord){
  dist.mat.old=dist.mat.new=dist.mat
  ac.coord.orig=ac.coord.old=ac.coord
  tmp=rnorm(n.ac*2,mean=ac.coord.old,sd=jump)
  ac.coord.prop=matrix(tmp,n.ac,2) #proposed coordinates
  
  for (i in 1:n.ac){
    for (j in 1:2){ #to sample each coordinate separately
      dist.mat.new=dist.mat.old
      ac.coord.new=ac.coord.old
      ac.coord.new[i,j]=ac.coord.prop[i,j]
        
      #get distances
      x2=(grid.coord[,'x']-ac.coord.old[i,'x'])^2
      y2=(grid.coord[,'y']-ac.coord.old[i,'y'])^2
      dist.mat.old[i,]=sqrt(x2+y2)
        
      x2=(grid.coord[,'x']-ac.coord.new[i,'x'])^2
      y2=(grid.coord[,'y']-ac.coord.new[i,'y'])^2
      dist.mat.new[i,]=sqrt(x2+y2)
          
      #get multinomial probabilities
      tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat.old,phi=phi,log.theta=log.theta,
                               n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,
                               dat=dat)
      tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat.new,phi=phi,log.theta=log.theta,
                               n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,
                               dat=dat)
      
      max1=-apply(cbind(tmp.old,tmp.new),1,max)
      max1=matrix(max1,n.tsegm,n.ac)
      tmp.old1=tmp.old+max1
      tmp.new1=tmp.new+max1
      pold=sum(log(rowSums(exp(tmp.old1))))
      pnew=sum(log(rowSums(exp(tmp.new1))))
      
      #accept or reject MH
      k=acceptMH(p0=pold,p1=pnew,x0=1,x1=2,BLOCK=F)
      if (k$x==2){
        ac.coord.old[i,j]=ac.coord.new[i,j]
        dist.mat.old=dist.mat.new
      }
    }
  }
  list(ac.coord=ac.coord.old,accept=ac.coord.orig!=ac.coord.old,dist.mat=dist.mat)
}
#-----------------------------------
sample.phi=function(dist.mat,phi,log.theta,
                    n.grid,n.ac,n.tsegm,
                    jump,dat){
  old=phi
  new=abs(rnorm(1,mean=old,sd=jump)) #reflection proposal around zero
  
  #get marginal loglikel
  tmp.old=get.calc.mloglik(dist.mat.sel=dist.mat,phi=old,log.theta=log.theta,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,
                           dat=dat)
  tmp.new=get.calc.mloglik(dist.mat.sel=dist.mat,phi=new,log.theta=log.theta,
                           n.tsegm=n.tsegm,n.ac=n.ac,n.grid=n.grid,
                           dat=dat)
  max1=apply(cbind(tmp.old,tmp.new),1,max)
  max1=matrix(max1,n.tsegm,n.ac)
  tmp.old1=tmp.old-max1
  tmp.new1=tmp.new-max1
  pold=sum(log(rowSums(exp(tmp.old1))))
  pnew=sum(log(rowSums(exp(tmp.new1))))
  
  #accept or reject MH
  k=acceptMH(p0=pold,p1=pnew,x0=old,x1=new,BLOCK=F)  
  logl=ifelse(k$accept==1,pnew,pold)
  list(phi=k$x,accept=k$accept)
}
#-----------------------------------
sample.z=function(phi,dist.mat,log.theta,
                  n.grid,n.ac,n.tsegm,
                  dat){
  #get distance
  prob=exp(-phi*dist.mat)
  prob=prob/rowSums(prob)
  lprob=log(prob)
  
  #get loglikel
  logl=matrix(NA,n.tsegm,n.ac)
  for (i in 1:n.ac){
    lprob1=matrix(lprob[i,],n.tsegm,n.grid,byrow=T)
    logl[,i]=rowSums(dat*lprob1)+log.theta[i]    
  }
  maximo=matrix(apply(logl,1,max),n.tsegm,n.ac)
  logl=logl-maximo
  tmp=exp(logl)
  soma=matrix(rowSums(tmp),n.tsegm,n.ac)
  prob=tmp/soma
  
  #sample from multinomial
  z=rmultinom1(prob=prob,randu=runif(n.tsegm))
  z+1
}
#-----------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#----------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.01
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
sample.v=function(z,n.ac,gamma1){
  tmp=table(z)
  n=rep(0,n.ac)
  n[as.numeric(names(tmp))]=tmp
  
  seq1=n.ac:1
  tmp=cumsum(n[seq1])
  n.ge=tmp[seq1]
  n.ge1=n.ge[-1]
  v=rbeta(n.ac-1,n[-n.ac]+1,n.ge1+gamma1)
  c(v,1)
}
#-----------------------------
get.loglikel=function(z,dist.mat,phi,
                      n.grid,
                      dat){
  #get log-probabilities for multinomial
  prob=exp(-phi*dist.mat)
  prob=prob/rowSums(prob)
  lprob=log(prob)
  
  #calculate conditional loglikel
  uni.z=unique(z)
  res=0
  for(i in uni.z){
    cond=z==i
    dat1=dat[cond,]
    n=sum(cond)
    lprob1=matrix(lprob[i,],n,n.grid,byrow=T)
    res=res+sum(dat1*lprob1)
  }
  res
}