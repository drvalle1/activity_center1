gibbs.activity.center=function(dat,grid.coord,n.ac,ac.coord.init,gamma1){
  #basic setup
  n.tsegm=nrow(dat)
  n.grid=nrow(grid.coord)
  grid.coord=data.matrix(grid.coord)
  
  #initial values
  ac.coord=data.matrix(ac.coord.init)
  z=sample(1:n.ac,size=n.tsegm,replace=T) #cluster membership
  phi=0.0001 #distance decay parameter
  theta=rep(1/n.ac,n.ac)
  dist.mat=GetDistance(AcCoord=ac.coord,GridCoord=grid.coord,Ngrid=n.grid,Nac=n.ac) 
    
  #matrices to store results
  store.coord=matrix(NA,ngibbs,n.ac*2)
  store.z=matrix(NA,ngibbs,n.tsegm)
  store.param=matrix(NA,ngibbs,1) #to store phi
  store.logl=matrix(NA,ngibbs,1)
  store.theta=matrix(NA,ngibbs,n.ac)
  
  #MH stuff
  adaptMH=50
  jump1=list(coord=matrix(10,n.ac,2),phi=0.2)
  accept1=list(coord=matrix(0,n.ac,2),phi=0)
  
  #gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    
    #sample coordinates
    tmp=sample.coord(ac.coord=ac.coord,log.theta=log(theta),phi=phi,dist.mat=dist.mat,
                     n.ac=n.ac,n.grid=n.grid,n.tsegm=n.tsegm,
                     jump=jump1$coord,dat=dat,grid.coord=grid.coord)
    ac.coord=tmp$ac.coord
    accept1$coord=accept1$coord+tmp$accept
    dist.mat=tmp$dist.mat
    # ac.coord=ac.coord.true
    
    #sample phi
    tmp=sample.phi(n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,
                   phi=phi,dist.mat=dist.mat,log.theta=log(theta),
                   jump=jump1$phi,dat=dat)
    phi=tmp$phi
    accept1$phi=accept1$phi+tmp$accept

    #sample z
    z=sample.z(dist.mat=dist.mat,phi=phi,log.theta=log(theta),
               n.grid=n.grid,n.ac=n.ac,n.tsegm=n.tsegm,
               dat=dat)
    # z=z.true
    
    #sample theta
    v=sample.v(z=z,n.ac=n.ac,gamma1=gamma1)
    theta=rep(NA,n.ac)
    theta[1]=v[1]
    tmp=(1-v[1])
    for (j in 2:n.ac){
      theta[j]=v[j]*tmp
      tmp=tmp*(1-v[j])
    } 
    
    #calculate loglikelihood
    logl=get.loglikel(z=z,dist.mat=dist.mat,phi=phi,
                      n.grid=n.grid,
                      dat=dat)

    if (i<nburn & i%%adaptMH==0){
      #adapt MH
      tmp=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=adaptMH)
      jump1=tmp$jump1
      accept1=tmp$accept1
      
      #re-order data from time to time according to theta (largest to smallest)
      ordem=order(theta,decreasing=T)
      theta=theta[ordem]
      ac.coord=ac.coord[ordem,]
      dist.mat=GetDistance(AcCoord=ac.coord,GridCoord=grid.coord,Ngrid=n.grid,Nac=n.ac) 
      
      znew=rep(NA,n.tsegm)
      for (j in 1:n.ac){
        cond=z==ordem[j]
        if (sum(cond)>0) znew[cond]=j
      }
      z=znew 
    }
    
    #store results
    store.coord[i,]=unlist(ac.coord)
    store.z[i,]=z
    store.param[i,]=phi
    store.logl[i,]=logl
    store.theta[i,]=theta
  }
  list(coord=store.coord,z=store.z,phi=store.param,logl=store.logl,theta=store.theta)  
}
