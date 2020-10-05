library(igraph)
library(tidyverse)

source("igraphSIR_functions.R")
#parameters
m=c(10,25,50) #number of edges added at each timestep
m0=5 #number of vertecies at timestep 0
N = 1000 #number of verticies in the graph

# parameters for desease
probs =  seq(0.005,0.025,0.005)
t_contageous = 15



path = "simulations/10/"
params = expand.grid(m = m, probs = probs) 
B=1000

apply(params, MARGIN = 1, FUN= function(params){
  sapply(x=1:B,
    function(x){
      ba_network = generateBANetwork(N=N,m0 = m0,m = params[1])
      
      runPASim(ba_network, p_transmit = params[2], t_contageous = t_contageous, x, path)
    })
})



