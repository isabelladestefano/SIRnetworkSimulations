library(tidyverse)
m=c(50,100, 150) #number of edges added at each timestep
m0=5 #number of vertecies at timestep 0
N = 5000 #number of verticies in the graph

probs =  seq(0.0005,.005,0.0005) #probability of transmission along an edge
t_contageous = 20 # time contageous

aplha = 2

params = expand.grid(m = m, probs = probs) %>%
  mutate( xm = m*sqrt((N-m0)/N), 
          mu = round(xm^alpha*alpha/ xm^(alpha-1)),
          beta = probs*mu,
          v = 1/t_contageous,
          R0 = t_contageous*probs*mu )

dS = function(beta,s,i,v){
  -beta*s*i
}

dI = function(beta,s,i,v){
  beta*s*i - v*i
}

dR = function(beta,s,i,v){
  v*i
}



SIR_dynamics = function(beta, v, mu , probs){
  s =.999
  i = .001
  r = 0
  
  props = data.frame(
    s=s, i=i,r=r)
  
  while(!isTRUE(all.equal(i, 0)) & !isTRUE(all.equal(i,1)) ){
    ds =  dS(beta,s,i,v)
    di = dI(beta,s,i,v)
    dr = dR(beta,s,i,v)
    s = s + ds
    i = i + di
    r = r + dr
    
    total = s+i+r
    props = bind_rows(props, data.frame(
      s=s/total, i=i/total,r=r/total))
  }
  print("done")
  props$beta = beta
  props$v = v
  props$mu = mu
  props$probs = probs
  return(props)
}

dynamical_SIR = mapply(SIR_dynamics, params$beta,params$v, params$mu, params$probs)




for(n in 1:30){
  print(n)
  if (n == 1){
    dynamical_SIR_nums = data.frame(dynamical_SIR[,n])
  } else{
    dynamical_SIR_nums = bind_rows(dynamical_SIR_nums,
      data.frame(dynamical_SIR[,n]))
  }
}

write.csv(dynamical_SIR_nums, "dynamical_systems/SIR/dynamical_system_N1000_t20_m05.csv")
