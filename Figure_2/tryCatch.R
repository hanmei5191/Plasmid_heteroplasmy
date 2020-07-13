w_max=m_max=40

nt_init=matrix(0,nrow=w_max+1,ncol=m_max+1)
w_init=20
m_init=0
nt_init[w_init+1,m_init+1]=3.12e9
rownames(nt_init)=paste0('w',0:w_max)
colnames(nt_init)=paste0('m',0:m_max)

cell_division <- function(nt,k_lethal){
  nt_div=nt*0
  for(a in 0:(k_lethal-1)){ # wild-type 
    for(b in 0:(k_lethal-1-a)){ # mutant
      p_nt=dbinom(0:(a+b),a+b,0.5)*nt[a+1,b+1] 
      for(jpk in 0:(a+b)){ 
        j_all=max(0,jpk-b):min(a,jpk) 
        new_all=p_nt[jpk+1]*dhyper(j_all,a,b,jpk) 
        nt_div[(jpk-j_all)*(w_max+1)+j_all+1]=nt_div[(jpk-j_all)*(w_max+1)+j_all+1]+new_all 
        nt_div[(b-jpk+j_all)*(w_max+1)+a-j_all+1]=nt_div[(b-jpk+j_all)*(w_max+1)+a-j_all+1]+new_all 
      }
    }
  }
  return (nt_div)
}

plasmid_replication_mutation <- function(nt_div,u,r,R,k_lethal){
  nt_growth=nt_div*0
  for(a in 0:(k_lethal-1)){ # wild-type
    a_plus1=floor(a*r)
    new_wt=a_plus1-a
    w_to_m=0:new_wt
    a_plus2=a_plus1-w_to_m
    a_plus2=ifelse(a_plus2<=w_max,a_plus2,w_max) 
    
    p_w_to_m=dbinom(w_to_m,new_wt,u)
    
    for(b in 0:(k_lethal-1-a)){ # mutant
      b_plus1=floor(b*R)
      b_plus2=b_plus1+w_to_m
      new=nt_div[a+1,b+1]*p_w_to_m
      for (i in w_to_m[b_plus2<=m_max]+1){
          nt_growth[a_plus2[i]+1,b_plus2[i]+1]=nt_growth[a_plus2[i]+1,b_plus2[i]+1]+new[i]
      }
    }       
  }
  return (nt_growth)   
}

cell_g0_selection <-function(nt_growth, S, k_optimal, k_lethal){
  nt_selection=nt_growth*0
  
  for(a in 0:(k_lethal-1)){ # wild-type
    for(b in 0:(k_lethal-1-a)){ # mutant
      if ((a+b)==0){
        nt_selection[a+1,b+1]=0
      } else if ((a+b)<=k_optimal){
        nt_selection[a+1,b+1]=nt_growth[a+1,b+1]*(S+(1-S)*(a+b-1)/(k_optimal-1))
      } else {
        nt_selection[a+1,b+1]=nt_growth[a+1,b+1]*((k_lethal-a-b)/(k_lethal-k_optimal))
      }
    }
  }

  if (sum(nt_selection)<1){
    return (nt_selection*0)
  } else if (sum(nt_selection)>3.12e9){
    return (nt_selection*3.12e9/sum(nt_selection))
  }  else {
    return (nt_selection)
  }
}

cell_g1_selection<-function(nt_growth, k_lethal){
  rotate=function(mat){ t(mat[nrow(mat):1,,drop=FALSE]) }
  nt_selection=nt_growth*0
  nt_selection[seq_len(k_lethal),seq_len(k_lethal)][rotate(lower.tri(matrix(NA,nrow=k_lethal,ncol=k_lethal),diag=TRUE))]=
    nt_growth[seq_len(k_lethal),seq_len(k_lethal)][rotate(lower.tri(matrix(NA,nrow=k_lethal,ncol=k_lethal),diag=TRUE))]
  nt_selection[1,1]=0 # a+b=0
  
  if (sum(nt_selection)<1){
    return (nt_selection*0)
  } else if (sum(nt_selection)>3.12e9){
    return (nt_selection*3.12e9/sum(nt_selection))
  }  else {
    return (nt_selection)
  }
}

fcell_2_fplasmid <- function(nt_g0_selection,nt_g1_selection){
  nt_selection=nt_g0_selection+nt_g1_selection
  wt_plasmid=sum(nt_selection*(0:w_max))
  mutant_plasmid=sum(t(nt_selection)*(0:m_max))  
  return (mutant_plasmid/(wt_plasmid+mutant_plasmid))    
}

nt_g0_init=nt_init
nt_g1_init=nt_init*0

generation=318
mutant_plasmid_freq=rep(-1,generation)
nt_g0=nt_g0_init
nt_g1=nt_g1_init

emperical_data=read.delim("F3_data1.tab", header=TRUE)
emperical_generations=emperical_data$generations

cost<-function(theta){
  gamma=theta[1]
  u=theta[2]
  r=theta[3]
  R=theta[4]
  S=theta[5]  
  k_optimal=round(theta[6])  
  k_lethal=round(theta[7])
  if(gamma<=0 | gamma>=1 | u<=0 | u>=1 | r<1 | R<1 | r>=R | S<=0.5 | S>1 | k_optimal<=0 | k_lethal<=k_optimal | k_lethal>m_max){
        return(+Inf)
  }
  for (i in 1:generation){
    nt_g0_div=cell_division(nt_g0,k_lethal)
    nt_g1_div=cell_division(nt_g1,k_lethal)
    
    nt_g1_div=nt_g1_div+nt_g0_div*gamma
    nt_g0_div=nt_g0_div-nt_g0_div*gamma
    
    nt_g0_growth=plasmid_replication_mutation(nt_g0_div,u,r,R,k_lethal)
    nt_g1_growth=plasmid_replication_mutation(nt_g1_div,u,r,R,k_lethal)
    
    nt_g0_selection=cell_g0_selection(nt_g0_growth, S, k_optimal, k_lethal)
    nt_g1_selection=cell_g1_selection(nt_g1_growth, k_lethal) 
    
    if (sum(nt_g0_selection)<1 | sum(nt_g1_selection)<1)
      return (+Inf)
    
    mutant_plasmid_freq[i]=fcell_2_fplasmid(nt_g0_selection,nt_g1_selection)
      
    nt_g0=nt_g0_selection
    nt_g1=nt_g1_selection
  }
  diff=mutant_plasmid_freq[emperical_generations]-emperical_data$frequency
  return (sum(diff^2))
}

draw_parameters=function(){
  gamma_log=runif(1, -8, -2)
  gamma=10^gamma_log
  u_log=runif(1, -8, -2)
  u=10^u_log
  r=runif(1, 1, 10)
  R=runif(1, r, 10)
  S=runif(1, 0.5, 1)
  k_optimal=sample(m_max-1,1)
  k_lethal=sample(max(w_init+m_init+1,k_optimal+1):m_max,1)
  theta0=c(gamma,u, r, R, S, k_optimal, k_lethal)
  return(theta0)
}

run_optim=function(trial){
  err=simpleError("Fake error to start")
  counter=1
  max_tries=100
  while(inherits(err, "error") && counter < max_tries){
    counter = counter + 1
    tryCatch(
       {theta0=draw_parameters()
        res=optim(theta0,cost,control=c(maxit=1e6))
        err="Run successfully"
        return(res)}, 
        error=function(e) {
        err=simpleError("error during loop")
        }
    )
  }
}

library(parallel)
numCores=detectCores()
trials=seq(1, 100)
ptm=proc.time()
results <- mclapply(trials, run_optim, mc.cores = numCores)

sink ("Plasmid_v3_optim.txt")
print(proc.time()-ptm)

cost_dist=NA
min_loss=1
min_loss_par=rep(-1,7)
num_unconv=0
num_optim=0
for (k in results){
    if (!(is.null(k))){
        num_optim=num_optim+1
        if (k$convergence == 1){
            num_unconv=num_unconv+1
        }
        
        else{
            cost_dist=append(cost_dist,k$value)
            if (k$value < min_loss){
                min_loss=k$value
                min_loss_par=k$par
            }
        }
    }
}

cat ("min_loss:", min_loss, "\n")
cat ("min_loss_par:", min_loss_par, "\n") 
cat ("Number of optimization performed:", num_optim, "\n")
cat ("Number of unconverged optimization:", num_unconv, "\n")
sink()

pdf("cost_dist.pdf")
hist(cost_dist)
dev.off()

save(results,file="results.RData")
