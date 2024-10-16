###########################################################################################################
############################################## import functions ###########################################
###########################################################################################################
rm(list=ls())
generate_data<-function(A,G,T_max, traj, dt,dt_EM,kill=FALSE, X0dist="const",X0mean=0,X0sd=0){
  dim= sqrt(length(A))
  N_max=round(T_max/dt_EM)
  X_vec=array(0, dim=c(N_max, dim, traj))
  
  if(kill){
    X_vec_cf=array(0, dim=c(N_max, dim, traj*N_max))
    W_vec_cf=X_vec_cf
    W_vec_cf[1,,]<- 0
    if(X0dist== "multi"){
      whichX0col<-sample(1:ncol(X0mean),traj*N_max,replace = T)
      for(d in 1:dim){
        X_vec_cf[1,d,]<-X0mean[d,whichX0col]
      }
    }
    if(X0dist=="const") X_vec_cf[1,,]<-X0mean
    if(X0dist=="normal") X_vec_cf[1,,]<- rnorm(dim*traj*N_max,mean = X0mean,sd=X0sd)
    for(i in 2:N_max){
      W_vec_cf[i,,]= W_vec_cf[i-1,,]+rnorm(dim*traj*N_max,mean = 0,sd=sqrt(dt_EM))
      X_vec_cf[i,,]= X_vec_cf[i-1,,] + A%*%X_vec_cf[i-1,,]*dt_EM + G%*%(W_vec_cf[i,,]- W_vec_cf[i-1,,])
    }
    
    counter=1
    for(i in 1:N_max){
      for(j in 1:traj){
        X_vec[i,,j]= X_vec_cf[i,,counter]
        counter=counter+1
      }
    }
    return(list(givenData=X_vec[seq(from=1,by=round(dt/dt_EM),length.out=round(T_max/dt)),,,drop=F],Counterfactuals=X_vec_cf))
  } else{
    W_vec=X_vec
    W_vec[1,,]<- 0
    if(X0dist== "multi"){
      whichX0col<-sample(1:ncol(X0mean),traj,replace = T)
      for(d in 1:dim){
        X_vec[1,d,]<-X0mean[d,whichX0col]
      }
    }
    if(X0dist=="const") X_vec[1,,]<-X0mean
    if(X0dist=="normal") X_vec[1,,]<- rnorm(dim*traj,mean = X0mean,sd=X0sd)
    for(i in 2:N_max){
      W_vec[i,,]= W_vec[i-1,,]+rnorm(dim*traj,mean = 0,sd=sqrt(dt_EM))
      X_vec[i,,]= X_vec[i-1,,] + A%*%X_vec[i-1,,]*dt_EM + G%*%(W_vec[i,,]- W_vec[i-1,,])
    }
    return(list(givenData=X_vec[seq(from=1,by=round(dt/dt_EM),length.out=round(T_max/dt)),,,drop=F],Counterfactuals=X_vec))
  }
  
}

library(mvnfast)
library(transport)
library(Rfast)
library(T4transport)
library(pracma)
library(Matrix)
library(foreach)
library(doParallel)
library(mvtnorm)

mysinkho2<-function(K, maxiter,u_thresh){
  a=matrix(rep(1/nrow(K),nrow(K)),ncol = 1)
  b=matrix(rep(1/ncol(K),ncol(K)),ncol=1)
  u=matrix(rep(1,nrow(K)),ncol = 1)
  v=matrix(rep(1,ncol(K)),ncol=1)
  i_in_sink=1
  diff<-1
  u_thresh<-u_thresh/nrow(K)
  cur_out<-diag(as.numeric(u))%*%K%*%diag(as.numeric(v))
  while(i_in_sink<maxiter & diff>u_thresh){
    cur_out_old<-cur_out
    u= a/(K%*%v)
    v = b/(t(K)%*%u) # t(G)
    cur_out= diag(as.numeric(u))%*%K%*%diag(as.numeric(v))
    diff<-max(abs(cur_out_old - cur_out))
    i_in_sink=i_in_sink+1
  }
  print(i_in_sink)
  return(cur_out)
}

mysinkho<-function(cost, lam, maxiter){
  G = exp(-cost/lam)
  a=matrix(rep(1/nrow(cost),nrow(cost)),ncol = 1)
  b=matrix(rep(1/ncol(cost),ncol(cost)),ncol=1)
  uold=matrix(rep(1,nrow(cost)),ncol = 1)
  vold=matrix(rep(1,ncol(cost)),ncol=1)
  for (i in 1:maxiter){
    unew = a/(G%*%vold)
    vnew = b/(t(G)%*%unew) # t(G)
    uold     = unew
    vold     = vnew
  }
  return(diag(as.numeric(unew))%*%G%*%diag(as.numeric(vnew)))
}




###########################################################################################################
############################################## experiment 1 ###############################################
###########################################################################################################
Alist=list(matrix(c(-1),nrow=1),matrix(c(-10),nrow=1)) #SDE3
Glist=list(matrix(c(1),nrow=1),matrix(c(sqrt(10)),nrow=1))
X0=matrix(c(5),nrow=1)
traj=100
N_max=100
dt=0.01
Tmax=dt*N_max

iter=30
rep=10
multi_sinkho_A<-array(0,dim = c(rep,iter,length(Alist)))
multi_sinkho_G<-array(0,dim= c(rep,iter,length(Alist)))

for(s in 1:length(Alist)){
  
  A=Alist[[s]]
  G=Glist[[s]]
  dim=nrow(A)
  dat<-generate_data(A,G,T_max=dt*N_max,traj = traj,dt=dt,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)
  
  for(r in 1:rep){
    A_est= -diag(dim)*0
    t_meanG<-mean(diag(G%*%t(G)))
    GGT_est=diag(dim)*runif(1,min=t_meanG*0.1,max=t_meanG*10)
    
    for(it in 1:iter){
      #GGT_est= G%*%t(G)*1.06
      #A_est= A
      X_vec_OT=dat$givenData
      
      tran_list<-list()
      linest<-(diag(dim)+A_est*dt)
      inv_cov<- -0.5*pinv(GGT_est*dt)
      K=matrix(0,nrow=traj,ncol=traj)
      time<-Sys.time()
      registerDoParallel(cl <- makeCluster(1))
      tran_list=foreach(t=2:N_max) %do% {
        for(i in 1:nrow(K)){
          dX=matrix(X_vec_OT[t,,],nrow = dim)- as.numeric(linest %*%X_vec_OT[t-1,,i])
          Ki1<-rowSums((t(dX)%*%inv_cov)*t(dX))
          Ki<-exp(Ki1)
          if(sum(Ki==0)>0) print(paste("the number of 0s in K is: ",sum(Ki==0)))
          Ki[Ki==0]<-.Machine$double.xmin*rank(Ki1[Ki==0])
          K[i,]<-Ki
          
        }
        mysinkho2(K,maxiter = 600)*traj
      }
      stopCluster(cl)
      Sys.time()- time
      
      # tran_list<-list()
      #    linest<-(diag(dim)+A_est*dt)
      #    inv_cov<- -0.5*pinv(GGT_est*dt)
      #    time<-Sys.time()
      #    for(t in 2:N_max){
      #        K=matrix(0,nrow=traj,ncol=traj)
      #        for(i in 1:nrow(K)){
      #            for(j in 1:ncol(K)){
      #               dX=matrix(X_vec_OT[t,,j]- linest %*%X_vec_OT[t-1,,i],ncol=1)
      #               K[i,j]=max(exp(t(dX)%*%inv_cov%*%dX),.Machine$double.xmin)
      #              }
      #          }
      #        tran_list[[t-1]]<-mysinkho2(K,maxiter = 2000)*traj
      #    }
      #    Sys.time()-time
      
      n_sim_traj=5*traj #create ROT trajectories
      pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      for(i in 1:dim(pos_vec_ROT)[3]){
        for(j in 1:dim(pos_vec_ROT)[1]){
          if(j==1){
            pos_vec_ROT[j,,i]<-sample(1:traj,1)
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          } else {
            pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          }
        }
      }
      
      #X_vec_ROT<-X_vec_OT
      
      dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
      
      sume_dxx=0
      sume_xx=0
      for(i in 2:N_max){
        sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
        if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      }
      A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
      
      X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
      for(i in 1:(N_max-1)){
        X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
      }
      
      dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
      cur_G=0
      for(t in 1:n_sim_traj){
        cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
      }
      GGT_est<- cur_G/n_sim_traj
      
      multi_sinkho_G[r,it,s]<-mean(abs(GGT_est- G%*%t(G)))/mean(abs( G%*%t(G)))
      multi_sinkho_A[r,it,s]<-mean(abs(A_est- A))/mean(abs(A))
      print(it)
      if(it %% 5 ==0){
        print(GGT_est)
        print(A_est)
      }
    }
    
  }
}

multi_sinkho_A<-multi_sinkho_A*100
multi_sinkho_G<-multi_sinkho_G*100

plot(colmeans(multi_sinkho_A[,,1]),type="b",ylim=c(0,max(multi_sinkho_A,multi_sinkho_G)),ylab="MAPE",xlab="iteration")
lines(colmeans(multi_sinkho_G[,,1]),type="b",col="red",pch=19)
lines(colmeans(multi_sinkho_G[,,2]),type="b",col="red",pch = 18)
lines(colmeans(multi_sinkho_A[,,2]),type="b",col="black",pch = 23)
arrows(1:iter, colMins(multi_sinkho_A[,,1],value = T), 1:iter, colMaxs(multi_sinkho_A[,,1],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_A[,,2],value = T), 1:iter, colMaxs(multi_sinkho_A[,,2],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_G[,,1],value = T), 1:iter, colMaxs(multi_sinkho_G[,,1],value = T), length=0.05, angle=90, code=3,col="red")
arrows(1:iter, colMins(multi_sinkho_G[,,2],value = T), 1:iter, colMaxs(multi_sinkho_G[,,2],value = T), length=0.05, angle=90, code=3,col="red")
legend("topright", legend=c("A1", "G1", "A2", "G2"), col=c("black","red","black","red"), cex=0.8,pch=c(1,19,23,18))




###########################################################################################################
############################################## experiment 2 ###############################################
###########################################################################################################
Alist=list(matrix(c(0,0,0,0),nrow=2),matrix(c(0,-1,1,0),nrow=2)) #SDE3
Glist=list(matrix(c(1,0,0,1),nrow=2),matrix(c(1,0,0,1),nrow=2))
X0=matrix(c(3,2,-3,1),nrow=2)
traj=100
N_max=50
dt=0.02
Tmax=dt*N_max

iter=50
rep=10
multi_sinkho_A<-array(0,dim = c(rep,iter,length(Alist)))
multi_sinkho_G<-array(0,dim= c(rep,iter,length(Alist)))



for(s in 1:length(Alist)){
  
  A=Alist[[s]]
  G=Glist[[s]]
  dim=nrow(A)
  dat<-generate_data(A,G,T_max==dt*N_max,traj = traj,dt=dt,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)
  
  for(r in 1:rep){
    A_est= -diag(dim)*0
    t_meanG<-mean(diag(G%*%t(G)))
    GGT_est=diag(dim)*runif(1,min=t_meanG*0.1,max=t_meanG*10)
    
    for(it in 1:iter){
      #GGT_est= G%*%t(G)*1.06
      #A_est= A
      X_vec_OT=dat$givenData
      
      tran_list<-list()
      linest<-(diag(dim)+A_est*dt)
      inv_cov<- -0.5*pinv(GGT_est*dt)
      K=matrix(0,nrow=traj,ncol=traj)
      time<-Sys.time()
      registerDoParallel(cl <- makeCluster(1))
      tran_list=foreach(t=2:N_max) %do% {
        for(i in 1:nrow(K)){
          dX=matrix(X_vec_OT[t,,],nrow = dim)- as.numeric(linest %*%X_vec_OT[t-1,,i])
          Ki1<-rowSums((t(dX)%*%inv_cov)*t(dX))
          Ki<-exp(Ki1)
          if(sum(Ki==0)==length(Ki) & i==N_max) print(paste("the number of 0s in K is: ",sum(Ki==0)))
          Ki[Ki==0]<-.Machine$double.xmin*rank(Ki1[Ki==0])
          K[i,]<-Ki
          
        }
        mysinkho2(K,maxiter = 600)*traj
      }
      stopCluster(cl)
      Sys.time()- time
      
      # tran_list<-list()
      #    linest<-(diag(dim)+A_est*dt)
      #    inv_cov<- -0.5*pinv(GGT_est*dt)
      #    time<-Sys.time()
      #    for(t in 2:N_max){
      #        K=matrix(0,nrow=traj,ncol=traj)
      #        for(i in 1:nrow(K)){
      #            for(j in 1:ncol(K)){
      #               dX=matrix(X_vec_OT[t,,j]- linest %*%X_vec_OT[t-1,,i],ncol=1)
      #               K[i,j]=max(exp(t(dX)%*%inv_cov%*%dX),.Machine$double.xmin)
      #              }
      #          }
      #        tran_list[[t-1]]<-mysinkho2(K,maxiter = 2000)*traj
      #    }
      #    Sys.time()-time
      
      n_sim_traj=5*traj #create ROT trajectories
      pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      for(i in 1:dim(pos_vec_ROT)[3]){
        for(j in 1:dim(pos_vec_ROT)[1]){
          if(j==1){
            pos_vec_ROT[j,,i]<-sample(1:traj,1)
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          } else {
            pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          }
        }
      }
      
      #X_vec_ROT<-X_vec_OT
      
      dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
      
      sume_dxx=0
      sume_xx=0
      for(i in 2:N_max){
        sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
        if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      }
      A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
      
      X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
      for(i in 1:(N_max-1)){
        X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
      }
      
      dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
      cur_G=0
      for(t in 1:n_sim_traj){
        cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
      }
      GGT_est<- cur_G/n_sim_traj
      
      multi_sinkho_G[r,it,s]<-mean(abs(GGT_est- G%*%t(G)))/mean(abs( G%*%t(G)))
      multi_sinkho_A[r,it,s]<-mean(abs(A_est- A))/mean(abs(A))
      print(it)
      if(it %% 5 ==0){
        print(GGT_est)
        print(A_est)
      }
    }
    
  }
}

multi_sinkho_A<-multi_sinkho_A*100
multi_sinkho_G<-multi_sinkho_G*100


plot(colmeans(multi_sinkho_A[,,1]),type="b",ylim=c(0,max(multi_sinkho_A,multi_sinkho_G)),ylab="MAPE",xlab="iteration")
lines(colmeans(multi_sinkho_G[,,1]),type="b",col="red",pch=19)
lines(colmeans(multi_sinkho_G[,,2]),type="b",col="red",pch = 18)
lines(colmeans(multi_sinkho_A[,,2]),type="b",col="black",pch = 23)
arrows(1:iter, colMins(multi_sinkho_A[,,1],value = T), 1:iter, colMaxs(multi_sinkho_A[,,1],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_A[,,2],value = T), 1:iter, colMaxs(multi_sinkho_A[,,2],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_G[,,1],value = T), 1:iter, colMaxs(multi_sinkho_G[,,1],value = T), length=0.05, angle=90, code=3,col="red")
arrows(1:iter, colMins(multi_sinkho_G[,,2],value = T), 1:iter, colMaxs(multi_sinkho_G[,,2],value = T), length=0.05, angle=90, code=3,col="red")
legend("topright", legend=c("A1", "G1", "A2", "G2"), col=c("black","red","black","red"), cex=0.8,pch=c(1,19,23,18))




###########################################################################################################
############################################## experiment 3 ###############################################
###########################################################################################################
Alist=list(matrix(c(1,1,2,0),nrow=2),matrix(c(1/3,2/3,4/3,-1/3),nrow=2)) #SDE3
Glist=list(matrix(c(1,-1,2,-2),nrow=2),matrix(c(1,-1,2,-2),nrow=2))
X0=matrix(c(5.56,-1.7266,-0.93,-1.464),nrow=2)
traj=100
N_max=50
dt=0.02
dt_EM=0.01
Tmax=dt*N_max

iter=50
rep=10
multi_sinkho_A<-array(0,dim = c(rep,iter,length(Alist)))
multi_sinkho_G<-array(0,dim= c(rep,iter,length(Alist)))



for(s in 1:length(Alist)){
  
  A=Alist[[s]]
  G=Glist[[s]]
  dim=nrow(A)
  dat<-generate_data(A,G,T_max=dt*N_max,traj = traj,dt_EM=dt_EM,dt=dt,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)
  
  for(r in 1:rep){
    A_est= -diag(dim)*0
    t_meanG<-mean(diag(G%*%t(G)))
    GGT_est=diag(dim)*runif(1,min=t_meanG*0.1,max=t_meanG*10)
    
    for(it in 1:iter){
      #GGT_est= G%*%t(G)*1.06
      #A_est= A
      X_vec_OT=dat$givenData
      
      tran_list<-list()
      linest<-(diag(dim)+A_est*dt)
      inv_cov<- -0.5*pinv(GGT_est*dt)
      K=matrix(0,nrow=traj,ncol=traj)
      time<-Sys.time()
      registerDoParallel(cl <- makeCluster(1))
      tran_list=foreach(t=2:N_max) %do% {
        for(i in 1:nrow(K)){
          dX=matrix(X_vec_OT[t,,],nrow = dim)- as.numeric(linest %*%matrix(X_vec_OT[t-1,,i],ncol = 1))
          Ki1<-rowSums((t(dX)%*%inv_cov)*t(dX))
          Ki<-exp(Ki1)
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          if(sum(Ki)==0) GGT_est<-GGT_est+diag(dim)*1e-8
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          K[i,]<-Ki
        }
        mysinkho2(K,maxiter = 1000)*traj
      }
      stopCluster(cl)
      Sys.time()- time
      
      # tran_list<-list()
      #    linest<-(diag(dim)+A_est*dt)
      #    inv_cov<- -0.5*pinv(GGT_est*dt)
      #    time<-Sys.time()
      #    for(t in 2:N_max){
      #        K=matrix(0,nrow=traj,ncol=traj)
      #        for(i in 1:nrow(K)){
      #            for(j in 1:ncol(K)){
      #               dX=matrix(X_vec_OT[t,,j]- linest %*%X_vec_OT[t-1,,i],ncol=1)
      #               K[i,j]=max(exp(t(dX)%*%inv_cov%*%dX),.Machine$double.xmin)
      #              }
      #          }
      #        tran_list[[t-1]]<-mysinkho2(K,maxiter = 2000)*traj
      #    }
      #    Sys.time()-time
      
      n_sim_traj=5*traj #create ROT trajectories
      pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      for(i in 1:dim(pos_vec_ROT)[3]){
        for(j in 1:dim(pos_vec_ROT)[1]){
          if(j==1){
            pos_vec_ROT[j,,i]<-sample(1:traj,1)
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          } else {
            pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          }
        }
      }
      
      #X_vec_ROT<-X_vec_OT
      
      dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
      
      sume_dxx=0
      sume_xx=0
      for(i in 2:N_max){
        sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
        if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      }
      A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
      
      X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
      for(i in 1:(N_max-1)){
        X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
      }
      
      dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
      cur_G=0
      for(t in 1:n_sim_traj){
        cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
        #print(t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt))
      }
      GGT_est<- cur_G/n_sim_traj
      
      multi_sinkho_G[r,it,s]<-mean(abs(GGT_est- G%*%t(G)))/mean(abs( G%*%t(G)))
      multi_sinkho_A[r,it,s]<-mean(abs(A_est- A))/mean(abs(A))
      print(it)
      if(it %% 5 ==0){
        print(GGT_est)
        print(A_est)
      }
    }
    
  }
}

multi_sinkho_A<-multi_sinkho_A*100
multi_sinkho_G<-multi_sinkho_G*100

plot(colmeans(multi_sinkho_A[,,1]),type="b",ylim=c(0,max(multi_sinkho_A,multi_sinkho_G)),ylab="MAPE",xlab="iteration")
lines(colmeans(multi_sinkho_G[,,1]),type="b",col="red",pch=19)
lines(colmeans(multi_sinkho_G[,,2]),type="b",col="red",pch = 18)
lines(colmeans(multi_sinkho_A[,,2]),type="b",col="black",pch = 23)
arrows(1:iter, colMins(multi_sinkho_A[,,1],value = T), 1:iter, colMaxs(multi_sinkho_A[,,1],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_A[,,2],value = T), 1:iter, colMaxs(multi_sinkho_A[,,2],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_G[,,1],value = T), 1:iter, colMaxs(multi_sinkho_G[,,1],value = T), length=0.05, angle=90, code=3,col="red")
arrows(1:iter, colMins(multi_sinkho_G[,,2],value = T), 1:iter, colMaxs(multi_sinkho_G[,,2],value = T), length=0.05, angle=90, code=3,col="red")
legend("topright", legend=c("A1", "G1", "A2", "G2"), col=c("black","red","black","red"), cex=0.8,pch=c(1,19,23,18))





###########################################################################################################
######################################### experiment random ###############################################
###########################################################################################################
versions<-3
Alist<-list()
Glist<-list()
for(i in 1:versions){
  cur_d<-sample(x=3:10,size=1)
  curmatrix<-matrix(runif(cur_d^2,min=-10,max=10),nrow=cur_d)
  while(max(Real(eig(curmatrix)))>1){
    curmatrix<-matrix(runif(cur_d^2,min=-10,max=10),nrow=cur_d)
    print(paste("fail",i))
  }
  Alist[[i]]<-curmatrix
  Glist[[i]]<-matrix(runif(cur_d^2,min=-5,max=5),nrow=cur_d)
}

X0=matrix(c(5.56,-1.7266,-0.93,-1.464),nrow=4)
traj=100
N_max=50
dt=0.02
Tmax=dt*N_max

iter=50
rep=10
multi_sinkho_A<-array(0,dim = c(rep,iter,length(Alist)))
multi_sinkho_G<-array(0,dim= c(rep,iter,length(Alist)))



for(s in 1:length(Alist)){
  
  A=Alist[[s]]
  G=Glist[[s]]
  dim=nrow(A)
  dat<-generate_data(A,G,T_max=dt*N_max,traj = traj,dt=dt,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)
  
  for(r in 1:rep){
    A_est= -diag(dim)*0
    t_meanG<-mean(diag(G%*%t(G)))
    GGT_est=diag(dim)*runif(1,min=t_meanG*0.1,max=t_meanG*10)
    
    for(it in 1:iter){
      #GGT_est= G%*%t(G)*1.06
      #A_est= A
      X_vec_OT=dat$givenData
      
      tran_list<-list()
      linest<-(diag(dim)+A_est*dt)
      inv_cov<- -0.5*pinv(GGT_est*dt)
      K=matrix(0,nrow=traj,ncol=traj)
      time<-Sys.time()
      registerDoParallel(cl <- makeCluster(1))
      tran_list=foreach(t=2:N_max) %do% {
        for(i in 1:nrow(K)){
          dX=matrix(X_vec_OT[t,,],nrow = dim)- as.numeric(linest %*%matrix(X_vec_OT[t-1,,i],ncol = 1))
          Ki1<-rowSums((t(dX)%*%inv_cov)*t(dX))
          Ki<-exp(Ki1)
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          if(sum(Ki)==0) GGT_est<-GGT_est+diag(dim)*1e-8
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          K[i,]<-Ki
        }
        mysinkho2(K,maxiter = 1000)*traj
      }
      stopCluster(cl)
      Sys.time()- time
      
      # tran_list<-list()
      #    linest<-(diag(dim)+A_est*dt)
      #    inv_cov<- -0.5*pinv(GGT_est*dt)
      #    time<-Sys.time()
      #    for(t in 2:N_max){
      #        K=matrix(0,nrow=traj,ncol=traj)
      #        for(i in 1:nrow(K)){
      #            for(j in 1:ncol(K)){
      #               dX=matrix(X_vec_OT[t,,j]- linest %*%X_vec_OT[t-1,,i],ncol=1)
      #               K[i,j]=max(exp(t(dX)%*%inv_cov%*%dX),.Machine$double.xmin)
      #              }
      #          }
      #        tran_list[[t-1]]<-mysinkho2(K,maxiter = 2000)*traj
      #    }
      #    Sys.time()-time
      
      n_sim_traj=5*traj #create ROT trajectories
      pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      for(i in 1:dim(pos_vec_ROT)[3]){
        for(j in 1:dim(pos_vec_ROT)[1]){
          if(j==1){
            pos_vec_ROT[j,,i]<-sample(1:traj,1)
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          } else {
            pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          }
        }
      }
      
      #X_vec_ROT<-X_vec_OT
      
      dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
      
      sume_dxx=0
      sume_xx=0
      for(i in 2:N_max){
        sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
        if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      }
      A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
      
      X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
      for(i in 1:(N_max-1)){
        X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
      }
      
      dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
      cur_G=0
      for(t in 1:n_sim_traj){
        cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
        #print(t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt))
      }
      GGT_est<- cur_G/n_sim_traj
      
      multi_sinkho_G[r,it,s]<-mean(abs(GGT_est- G%*%t(G)))/mean(abs( G%*%t(G)))
      multi_sinkho_A[r,it,s]<-mean(abs(A_est- A))/mean(abs(A))
      print(it)
      if(it %% 5 ==0){
        print(GGT_est)
        print(A_est)
      }
    }
    
  }
}

multi_sinkho_A<-multi_sinkho_A*100
multi_sinkho_G<-multi_sinkho_G*100

plot(colmeans(multi_sinkho_A[,,1]),type="b",ylim=c(0,max(multi_sinkho_A,multi_sinkho_G)),ylab="MAPE",xlab="iteration")
lines(colmeans(multi_sinkho_G[,,1]),type="b",col="red",pch=19)
lines(colmeans(multi_sinkho_G[,,2]),type="b",col="red",pch = 18)
lines(colmeans(multi_sinkho_A[,,2]),type="b",col="black",pch = 23)
arrows(1:iter, colMins(multi_sinkho_A[,,1],value = T), 1:iter, colMaxs(multi_sinkho_A[,,1],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_A[,,2],value = T), 1:iter, colMaxs(multi_sinkho_A[,,2],value = T), length=0.05, angle=90, code=3)
arrows(1:iter, colMins(multi_sinkho_G[,,1],value = T), 1:iter, colMaxs(multi_sinkho_G[,,1],value = T), length=0.05, angle=90, code=3,col="red")
arrows(1:iter, colMins(multi_sinkho_G[,,2],value = T), 1:iter, colMaxs(multi_sinkho_G[,,2],value = T), length=0.05, angle=90, code=3,col="red")
legend("topright", legend=c("A1", "G1", "A2", "G2"), col=c("black","red","black","red"), cex=0.8,pch=c(1,19,23,18))



###########################################################################################################
######################################### consistency/convergence experiment ##############################
###########################################################################################################

versions<-70
Alist<-list()
Glist<-list()
for(i in 1:versions){
  cur_d<-12 #sample(x=10:12,size=1)
  curmatrix<-matrix(runif(cur_d^2,min=-2,max=2),nrow=cur_d)
  while(max(Real(eig(curmatrix)))>1){
    curmatrix<-matrix(runif(cur_d^2,min=-2,max=2),nrow=cur_d)
    #print(paste("fail",i))
  }
  #print(paste("success",i))
  Alist[[i]]<-curmatrix
  Glist[[i]]<-matrix(runif(cur_d^2,min=-1,max=1),nrow=cur_d)
}


traj_vec<-c(250,400)
N_max=20
dt=0.05
dt_EM=0.01
Tmax=dt*N_max

iter=30
multi_sinkho_A<-array(0,dim = c(length(traj_vec),iter,length(Alist)))
multi_sinkho_G<-array(0,dim= c(length(traj_vec),iter,length(Alist)))



for(s in 1:length(Alist)){
  
  A=Alist[[s]]
  G=Glist[[s]]
  dim=nrow(A)
  X0=diag(dim)*3
  dat<-generate_data(A,G,T_max=dt*N_max,traj = max(traj_vec),dt=dt,dt_EM=dt_EM,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)
  
  for(r in 1:length(traj_vec)){
    A_est= A*0
    traj<-traj_vec[r]
    t_meanG<-max(diag(G%*%t(G)))
    GGT_est= diag(dim)*t_meanG
    
    
    for(it in 1:iter){
      #GGT_est= G%*%t(G)*1.06
      #A_est= A
      X_vec_OT=dat$givenData[,,1:traj]
      
      tran_list<-list()
      linest<-(diag(dim)+A_est*dt)
      inv_cov<- -0.5*pinv(GGT_est*dt)
      K=matrix(0,nrow=traj,ncol=traj)
      time<-Sys.time()
      registerDoParallel(cl <- makeCluster(1))
      tran_list=foreach(t=2:N_max) %do% {
        for(i in 1:nrow(K)){
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          if(sum(Ki)==0) GGT_est<-GGT_est+diag(dim)*1e-8
          Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
          K[i,]<-Ki
        }
        K[which(rowSums(K)==0),]<-rnorm(length(as.numeric(K[which(rowSums(K)==0),])),sd=1e-8)
        K[,which(colSums(K)==0)]<-.Machine$double.xmin
        mysinkho2(K,maxiter = 10000,u_thresh=1e-4)*traj
      }
      stopCluster(cl)
      Sys.time()- time
      print("step 1 complete")


      n_sim_traj=5*traj #create ROT trajectories
      pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
      for(i in 1:dim(pos_vec_ROT)[3]){
        for(j in 1:dim(pos_vec_ROT)[1]){
          if(j==1){
            pos_vec_ROT[j,,i]<-sample(1:traj,1)
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          } else {
            pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
            X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
          }
        }
      }
      print("step 2 complete")
      
      #X_vec_ROT<-X_vec_OT
      n_sim_traj<-dim(X_vec_ROT)[3]
      
      dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
      
      sume_dxx=0
      sume_xx=0
      for(i in 2:N_max){
        sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
        if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      }
      A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
      
      X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
      for(i in 1:(N_max-1)){
        X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
      }
      
      dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
      cur_G=0
      for(t in 1:n_sim_traj){
        cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
        #print(t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt))
      }
      GGT_est<- cur_G/n_sim_traj
      
      multi_sinkho_G[r,it,s]<-mean(abs(GGT_est- G%*%t(G))^2)
      multi_sinkho_A[r,it,s]<-mean(abs(A_est- A)^2)
      print(it)
      if(it %% 5 ==0){
        print(GGT_est)
        print(A_est)
      }
    }
    write.csv(multi_sinkho_G[,iter,],"C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/consistencyExp_G_GIMPI_12d_big.csv",col.names = F)
    write.csv(multi_sinkho_A[,iter,],"C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/consistencyExp_A_GIMPI_12d_big.csv",col.names = F)
    write.csv(multi_sinkho_G[,1,],"C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/consistencyExp_G_WOT_12d_big.csv",col.names = F)
    write.csv(multi_sinkho_A[,1,],"C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/consistencyExp_A_WOT_12d_big.csv",col.names = F)
  }
  print(paste("Done with SDE:",s))
}

GIMPI_A<-rowMeans(multi_sinkho_A[,30,1:69])
WOT_A<-rowMeans(multi_sinkho_A[,1,1:69])
GIMPI_G<-rowMeans(multi_sinkho_G[,30,1:69])
WOT_G<-rowMeans(multi_sinkho_G[,1,1:69])

jpeg(paste0("C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/Figures/GIMPI_convergenceRate_12d.jpeg"),width = 6, height = 5, units = 'in',res = 400)
par(mar = c(5, 5, 2, 4.5))
plot(traj_vec,GIMPI_A,log='xy',type="b",ylab="A MSE",xlab="Observations per marginal")
lines(traj_vec,traj_vec^(-1)*5,lty=3)
par(new = TRUE)
plot(traj_vec,GIMPI_G,type="b",col="steelblue2",pch=19,log='xy',ylab = "",xlab = "",axes = FALSE, bty = "n")
axis(side=4, at = pretty(range(GIMPI_G)),col="steelblue2",col.axis="steelblue2")
mtext("G MSE", side=4, line=3,col = "steelblue2")
dev.off()

plot(traj_vec,GIMPI_G,type="l")
plot(traj_vec,GIMPI_A,type="l")

jpeg(paste0("C:/Users/joeja/Desktop/researchMasters/ParamEstMarginals/Figures/WOT_convergenceRate_12d.jpeg"),width = 6, height = 5, units = 'in',res = 400)
par(mar = c(5, 5, 2, 4.5))
plot(traj_vec,WOT_A,log='xy',type="b",ylab="A MSE",xlab="Observations per marginal")
par(new = TRUE)
plot(traj_vec,WOT_G,type="b",col="steelblue2",pch=19,log='xy',ylab = "",xlab = "",axes = FALSE, bty = "n")
axis(side=4, at = pretty(range(WOT_G)),col="steelblue2",col.axis="steelblue2")
mtext("G MSE", side=4, line=3,col = "steelblue2")
dev.off()





###########################################################################################################
######################################### multi-start experiment ##############################
###########################################################################################################
N_max=20
dt=0.05
dt_EM=0.01
Tmax=dt*N_max
iter=30
traj=100

cur_d<-5 #sample(x=10:12,size=1)
A<-matrix(runif(cur_d^2,min=-5,max=5),nrow=cur_d)
while(max(Real(eigen(A)$values))>2){
  A<-matrix(runif(cur_d^2,min=-5,max=5),nrow=cur_d)
}
G<-matrix(runif(cur_d^2,min=-2,max=2),nrow=cur_d)

dim=nrow(A)
X0=diag(dim)*3
dat<-generate_data(A,G,T_max=dt*N_max,traj = traj,dt=dt,dt_EM=dt_EM,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)



reps<-100
Alist<-list()
Glist<-list()


for(r in 1:reps){
  A_est= matrix(runif(cur_d^2,min=-5,max=5),nrow=cur_d)
  Gguess<-matrix(runif(cur_d^2,min=-3,max=3),nrow=cur_d)+diag(cur_d)*20
  GGT_est= Gguess%*%t(Gguess)
  
  
  for(it in 1:iter){
    #GGT_est= G%*%t(G)*1.06
    #A_est= A
    X_vec_OT=dat$givenData[,,1:traj]
    
    tran_list<-list()
    linest<-(diag(dim)+A_est*dt)
    inv_cov<- -0.5*pinv(GGT_est*dt)
    K=matrix(0,nrow=traj,ncol=traj)
    time<-Sys.time()
    registerDoParallel(cl <- makeCluster(1))
    tran_list=foreach(t=2:N_max) %do% {
      for(i in 1:nrow(K)){
        Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
        if(sum(Ki)==0) GGT_est<-GGT_est+diag(dim)*1e-8
        Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
        K[i,]<-Ki
      }
      K[which(rowSums(K)==0),]<-rnorm(length(as.numeric(K[which(rowSums(K)==0),])),sd=1e-8)
      K[,which(colSums(K)==0)]<-.Machine$double.xmin
      mysinkho2(K,maxiter = 10000,u_thresh=1e-6)*traj
    }
    stopCluster(cl)
    Sys.time()- time
    print("step 1 complete")
    
    
    n_sim_traj=5*traj #create ROT trajectories
    pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
    X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
    for(i in 1:dim(pos_vec_ROT)[3]){
      for(j in 1:dim(pos_vec_ROT)[1]){
        if(j==1){
          pos_vec_ROT[j,,i]<-sample(1:traj,1)
          X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
        } else {
          pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
          X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
        }
      }
    }
    print("step 2 complete")
    
    #X_vec_ROT<-X_vec_OT
    n_sim_traj<-dim(X_vec_ROT)[3]
    
    dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]
    
    sume_dxx=0
    sume_xx=0
    for(i in 2:N_max){
      sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
      if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
    }
    A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)
    
    X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
    for(i in 1:(N_max-1)){
      X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
    }
    
    dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
    cur_G=0
    for(t in 1:n_sim_traj){
      cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
      #print(t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt))
    }
    GGT_est<- cur_G/n_sim_traj
    
    print(it)
    if(it %% 10 ==0){
      print(GGT_est)
      print(A_est)
    }
  }
  Alist[[r]]<-A_est
  Glist[[r]]<-GGT_est
}





###########################################################################################################
######################################### single step optimization         ##############################
###########################################################################################################
N_max=20
dt=0.05
dt_EM=0.01
Tmax=dt*N_max
iter=30
traj=200

cur_d<-3 #sample(x=10:12,size=1)
A<-matrix(runif(cur_d^2,min=-2,max=2),nrow=cur_d)
while(max(Real(eigen(A)$values))>1){
  A<-matrix(runif(cur_d^2,min=-2,max=2),nrow=cur_d)
}
G<-matrix(runif(cur_d^2,min=-1,max=1),nrow=cur_d)

dim=nrow(A)
X0=diag(dim)*3
dat<-generate_data(A,G,T_max=dt*N_max,traj = traj,dt=dt,dt_EM=dt_EM,kill = F,X0dist = "multi",X0mean = X0,X0sd = 0)


A_est= matrix(runif(cur_d^2,min=-5,max=5),nrow=cur_d)
Gguess<-matrix(runif(cur_d^2,min=-3,max=3),nrow=cur_d)+diag(cur_d)*20
GGT_est= Gguess%*%t(Gguess)
tran_list<-list()
for(it in 1:iter){
  #GGT_est= G%*%t(G)*1.06
  #A_est= A
  X_vec_OT=dat$givenData[,,1:traj]
  
  linest<-(diag(dim)+A_est*dt)
  K=matrix(0,nrow=traj,ncol=traj)
  for(t in 2:N_max){
    if(it==1){
      for(i in 1:nrow(K)){
        Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
        if(sum(Ki)==0) GGT_est<-GGT_est+diag(dim)*1e-8
        Ki=dmvnorm(t(matrix(X_vec_OT[t,,],nrow = dim)),mean = as.numeric(linest %*%X_vec_OT[t-1,,i]),sigma = GGT_est*dt)
        K[i,]<-Ki
      }
    } else{
      K<-tran_list[[t-1]]
    }
    
    K[which(rowSums(K)==0),]<-rnorm(length(as.numeric(K[which(rowSums(K)==0),])),sd=1e-8)
    K[,which(colSums(K)==0)]<-.Machine$double.xmin
    tran_list[[t-1]]<-mysinkho2(K,maxiter = 10000,u_thresh=1e-6)*traj
  }
  print("step 1 complete")
  
}


n_sim_traj=5*traj #create ROT trajectories
pos_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
X_vec_ROT<-array(0,dim=c(N_max, dim, n_sim_traj))
for(i in 1:dim(pos_vec_ROT)[3]){
  for(j in 1:dim(pos_vec_ROT)[1]){
    if(j==1){
      pos_vec_ROT[j,,i]<-sample(1:traj,1)
      X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
    } else {
      pos_vec_ROT[j,,i]<-sample(1:traj,1,prob = tran_list[[j-1]][pos_vec_ROT[j-1,1,i],] )
      X_vec_ROT[j,,i]<-X_vec_OT[j,,pos_vec_ROT[j,1,i]]
    }
  }
}
print("step 2 complete")

#X_vec_ROT<-X_vec_OT
n_sim_traj<-dim(X_vec_ROT)[3]

dX_vec_ROT= X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]-X_vec_ROT[1:(nrow(X_vec_ROT)-1),,,drop=F]

sume_dxx=0
sume_xx=0
for(i in 2:N_max){
  sume_xx<-sume_xx+ matrix(X_vec_ROT[i-1,,],nrow=dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
  if(i>1) sume_dxx<-sume_dxx+ matrix(dX_vec_ROT[i-1,,],nrow = dim) %*% t(matrix(X_vec_ROT[i-1,,],nrow=dim))
}
A_est<- (1/dt)* sume_dxx %*% solve(sume_xx)

X_vec_ROT_mod<-array(NA, dim=dim(X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]))
for(i in 1:(N_max-1)){
  X_vec_ROT_mod[i,,]<-(A_est*dt+diag(dim)) %*% X_vec_ROT[i,,]
}

dX_vec_ROT=X_vec_ROT[2:nrow(X_vec_ROT),,,drop=F]- X_vec_ROT_mod
cur_G=0
for(t in 1:n_sim_traj){
  cur_G= cur_G+ t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt)
  #print(t(matrix(dX_vec_ROT[,,t],ncol=dim))%*% matrix(dX_vec_ROT[,,t],ncol = dim)/(Tmax-dt))
}
GGT_est<- cur_G/n_sim_traj

