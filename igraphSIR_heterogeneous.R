library(igraph)
library(tidyverse)


# parameters for simulations----
m=c(50,100, 150) #number of edges added at each timestep
m0=5 #number of vertecies at timestep 0
N = 5000 #number of verticies in the graph

probs =  seq(0.001,.005,0.001) #probability of transmission along an edge
t_contageous = 20 # time contageous
params = expand.grid(m = m, probs = probs) 




# simulation function ----
runSIRsimBAnetwrok = function(N,m0, m, p_transmit, t_contageous, simnum, path){
  
  #make network using BA algorithm
  start_graph = make_empty_graph(n = m0,directed = F) #algorithm begins with an empty graph with m0 nodes
  network = sample_pa(n = N, m=m, start.graph = start_graph , directed = F) #igraph builed prefferential attatchment using barbassi-albert model
  
  deg = degree(network) # degrees of the verticies in network
  
  #calculate theoretical degree mean from pareto distribution
  xm = m*sqrt((N-m0)/N) 
  alpha = 2
  theo_mu = round(xm^alpha*alpha/ xm^(alpha-1)) #theoretical mean of the degree distribution 
  
  
  init_infected = sample(which(abs(deg-theo_mu)==min(abs(deg-theo_mu))),1) #pick initial infected to have degree equal to theoretical mean degree
  init_degree = deg[init_infected]
  
  
  #set up atrributes for network 
  network = set.vertex.attribute(network, name = "status",value = "susceptable") #SIR status
  network = set.vertex.attribute(network, name = "time_infected", value = 0) #time infected with desease once status = "infected"
  network = set.edge.attribute(network, name = "transmission", value = 0) #this tracks along which edge a new infection is transmitted at any given time point
  network = set.vertex.attribute(network, name = "status", index =  init_infected, value = "infected") #SIR status for initial infected
  
  t=0 #start desease clock
  
  
  #start timing simulation for logging purposes
  start =Sys.time() 
  print(paste0("start: ", start))
  
  
  which_infected = which(vertex.attributes(network)$status == "infected") #get vertex idecies of initially infected individuals
  all_SIR_counts <- data.frame(infects = c(), 
                               susceptable = c(),
                               recovered = c(),
                               simnum = c(),
                               ID = c() )
  #continue time stepping until either no-one or everyone is infected
  
  num_infected = length(which_infected)
  while(num_infected > 0 & num_infected<N  ){
    which_infected = which(vertex.attributes(network)$status == "infected") #get vertex idecies of infected individuals
    t=t+1
    
    
    network = set.edge.attribute(network, name = "transmission", value = 0) #reset transmission over edges for this time step
    
    #for every infected individual
    sapply(1:length(which_infected), function(i){
      which_infect = which_infected[i]
      possible_infect= adjacent_vertices(network,which_infect)[[1]] #find immediate neightbors
      new_infect = possible_infect[1==rbinom(length(possible_infect),1,p_transmit)] #independently for each neighbor a bernouli draw decides whether the infection spreads across a particular edge
      if(length(new_infect)>0){ #forthe new infections caused by the individual which_infect
        edges <- c()
        sapply(new_infect, function(x){edges <<- c(edges, which_infect, x)}) #get the edges along which desease was transmitted
        network <<- set.edge.attribute(network, name = 'transmission', index =edges, value = 1) #update edge transmission for this time-step
        network <<- set.vertex.attribute(network, name = 'status', index = new_infect, value = 'infected') #change SIR status of newly infected individuals
      }
    })
    
    network = set.vertex.attribute(network, name = "time_infected", index = which_infected, value = vertex.attributes(network, which_infected)$time_infected+1) #increment time infected
    
    new_recover = V(network)[vertex.attributes(network)$time_infected > t_contageous] #which vertecies are infected longer than time_contageous (recovered)
    network = set.vertex.attribute(network, name = 'status', index = new_recover, value = 'recovered') #update status of recovered individuals 
    
    #summarise and formate network counts to save for results
    SIR_counts  = data.frame(status = vertex.attributes(network)$status)   %>%
      group_by(status)%>%
      summarise(N=n())   %>%
      spread(status,N) %>%
      mutate(t= t)
    
    if(!"infected" %in% names(SIR_counts)){
      SIR_counts$infected <- 0
    } 
    if(!"susceptable" %in% names(SIR_counts)){
      SIR_counts$susceptable <- 0
    } 
    if(!"recovered" %in% names(SIR_counts)){
      SIR_counts$recovered <- 0
    } 
    
    
    #add timestep counts to simulations counts
    
    SIR_counts$simnum = simnum
    SIR_counts$ID = paste("p",p_transmit, "_m",m)
    all_SIR_counts = bind_rows(all_SIR_counts, SIR_counts)
    num_infected = SIR_counts$infected[1]  
    
    #save(network, file = paste0(path,"duration",t,"_N",N,"_m0",m0,"_m",m,"_p", p_transmit, "_t", t_contageous, "_",simnum,".Rdata"))
    print(num_infected)
  }
  
  
  
  
  print(paste0("end: ", Sys.time() - start))
  
  #save parameters in different file for theoretical comparison using SIR dynamical system
  parameters = data.frame(
    init_degree = init_degree,
    N = N, 
    m0 = m0,
    m = m,
    p_transmit = p_transmit,
    t_contageous = t_contageous,
    simnum = simnum,
    ID = paste("p",p_transmit, "_m",m)
  )
  
  
  write.csv(all_SIR_counts, file = paste0(path,"population_counts/",simnum,"_p",p_transmit*10^3,"E-3_m",m, "_N",N,".csv"))
  write.csv(all_SIR_counts, file = paste0(path,"parameters/",simnum,"_p",p_transmit*10^3,"E-3_m",m,"_N",N, ".csv"))
}



# run 1000 simulations for each parameter combination ----
path = "simulations/heterogeneous/"
#B=1000

apply(params, MARGIN = 1, FUN= function(params){
  sapply(1:100,
         FUN = function(x){
           runSIRsimBAnetwrok(N=N,m0=m0, m =params[1], p_transmit = params[2], t_contageous = t_contageous, x, path)
         })
})

runSIRsimBAnetwrok(N=N,m0=m0, m =50, p_transmit = 0.005, t_contageous = t_contageous, 1, path)

