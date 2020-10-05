library(igraph)
library(tidyverse)
library(foreach)
library(doParallel)



# parameters for desease
probs =  seq(0.005,0.025,0.005)



#parameters

params = expand.grid(m = m, probs = probs) 




apply(params, MARGIN = 1, FUN= function(params){
  
  registerDoParallel()
  #fixed_params
  m=c(10,25,50) #number of edges added at each timestep
  m0=5 #number of vertecies at timestep 0
  N = 1000 #number of verticies in the graph
  t_contageous = 15
  path = "simulations/10/"
  B=1000
  #define function
  runBASim = function(N,m0, m, p_transmit, t_contageous, simnum, path){
    #make network using BA algorithm
    start_graph = make_empty_graph(n = m0,directed = F)
    network = sample_pa(n = N, m=m, start.graph = start_graph , directed = F)
    deg = degree(network) # degrees of the verticies
    
    init_infected = sample(which(deg == round(mean(deg))),1) # pick initial infected to have mean number of connected edges
    
    network = set.vertex.attribute(network, name = "status",value = "susceptable")
    network = set.vertex.attribute(network, name = "status", index =  init_infected, value = "infected")
    network = set.vertex.attribute(network, name = "time_infected", value = 0)
    
    network = set.edge.attribute(network, name = "transmission", value = 0)
    
    t=0
    
    which_infected = which(vertex.attributes(network)$status == "infected")
    start =Sys.time()
    print(paste0("start: ", start))
    while(
      length(which_infected) > 0 & length(which_infected)<N){
      
      which_infected = which(vertex.attributes(network)$status == "infected")
      t=t+1
      
      network = set.vertex.attribute(network, name = "time_infected", index = which_infected, value = vertex.attributes(network, which_infected)$time_infected+1)
      network = set.edge.attribute(network, name = "transmission", value = 0)
      
      sapply(1:length(which_infected), function(i){
        which_infect = which_infected[i]
        possible_infect= adjacent_vertices(network,which_infect)[[1]]
        new_infect = possible_infect[1==rbinom(length(possible_infect),1,p_transmit)]
        if(length(new_infect)>0){
          edges <- c()
          sapply(new_infect, function(x){edges <<- c(edges, which_infect, x)})
          network <<- set.edge.attribute(network, name = 'transmission', index =edges, value = 1)
          network <<- set.vertex.attribute(network, name = 'status', index = new_infect, value = 'infected')
        }
      })
      
      
      new_recover = V(network)[vertex.attributes(network)$time_infected > t_contageous]
      network = set.vertex.attribute(network, name = 'status', index = new_recover, value = 'recovered')
      
      
      save(network, file = paste0(path,"duration",t,"_N",N,"_m0",m0,"_m",m,"_p", p_transmit, "_t", t_contageous, "_",simnum,".Rdata"))
    }
    print(paste0("end: ", Sys.time() - start))
    
  }

  foreach(x=1:B, .packages = c('tidyverse','igraph'), .verbose = T)%dopar%
      runBASim(N=N,m0=m0, m =params[1], p_transmit = params[2], t_contageous = t_contageous, x, path)
})


runBASim(N=N,m0=m0, m =m[1], p_transmit = probs[2], t_contageous = t_contageous, 1, path)

