library(tidyverse)

sim_files = list.files("./simulations/heterogeneous/population_counts")
print(sim_files)
simulations = data.frame(infected = c(),
                         suceptable = c(),
                         t = c(),
                         recivered = c(),
                         simnum = c(),
                         ID = c())

for( f in sim_files){
  sim = read_csv(paste0("./simulations/heterogeneous/population_counts/",f))
  simulations = bind_rows(simulations, sim)
}


m0=5
N = 1000
alpha = 2
t_contageous = 20


simulations$p = unlist(lapply(str_split(simulations$ID, " _"), function(x){
   p = str_split(x[1], " ")
   as.numeric(p[[1]][2])
}))

simulations$m = unlist(lapply(str_split(simulations$ID, " _"), function(x){
  m = str_split(x[2], " ")
  as.numeric(m[[1]][2])
}))


write.csv(simulations, file = "simulations/heterogeneous/simulations_N1000_t20_m05.csv")
