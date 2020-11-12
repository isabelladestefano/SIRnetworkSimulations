library(tidyverse)
library(igraph)

path = "C:/Users/isabe/OneDrive/Documents/2020/SIR_network/simulations/10"
files =list.files("C:/Users/isabe/OneDrive/Documents/2020/SIR_network/simulations/10")


all_SIR_counts = data.frame(infected=c(), recovered = c(), susceptable = c())

for(f in 1:length(files)){
  filename = files[f]
  load(paste0(path,"/",filename))
  SIR_counts = data.frame(status = vertex.attributes(network)$status)%>%
  mutate(status = if_else(status == 'suseptable', 'susceptable', status))%>%
  group_by(status)%>%
  summarise(N=n())   %>%
  spread(status,N)%>%
  mutate(id = files[f]) 
  
  if(!"infected" %in% names(SIR_counts)){
    SIR_counts$infected <- 0
  } 
  if(!"susceptable" %in% names(SIR_counts)){
    SIR_counts$susceptable <- 0
  } 
  if(!"recovered" %in% names(SIR_counts)){
    SIR_counts$recovered <- 0
  } 

  
  all_SIR_counts = bind_rows(all_SIR_counts, SIR_counts)
}


all_SIR_counts = all_SIR_counts %>%
  mutate(prop_susceptable = susceptable/sum(susceptable,infected,recovered), 
         new_infected = infected - lag(infected), 
         R = new_infected/infected)

string = strsplit('uration101_N1000_m05_m10_p0.001_t15_1.Rdata', split = "_")[[1]][1] 
grep("d+",string)

