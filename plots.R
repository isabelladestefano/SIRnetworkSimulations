library(tidyverse)

simulations = read_csv(file = "simulations/heterogeneous/simulations_N1000_t20_m05.csv")
dynamical_system = read_csv(file = "dynamical_systems/SIR/dynamical_system_N1000_t20_m05.csv")


m0=5 #number of vertecies at timestep 0
N = 1000 #number of verticies in the graph
t_contageous = 20 # time contageous

alpha = 2


simulations = simulations %>% 
  mutate(xm = m*sqrt((N-m0)/N),
         mu = round(xm^alpha*alpha/ xm^(alpha-1),-1),
         R0 = mu*p*t_contageous,
         s = susceptable/N)
  
dynamical_system = dynamical_system %>% 
  group_by(beta,v,mu,probs)%>%
  mutate(t= 1:length(s),
         p=probs)
sim_plot_data = simulations %>%
  filter(p < 0.002)
  

dynam_plot_data =dynamical_system %>%
  filter(p < 0.002)

sim_plot_data %>%
  ggplot()+
  geom_point(aes(x=t,y=log(1-s)), col= "lightblue", alpha = 0.3, size =0.5)+
  geom_point(data = dynam_plot_data, aes(x=t, y = log(1-s)), col ="red", size = 0.5)+
  facet_grid(mu~p)


# plot log(1-s) vs t
# plot di as a function of s or (1-s) or log(1-s)

sim_plot_data = sim_plot_data %>%
  group_by(p,mu)%>%
  mutate(R = infected/N - lag(infected)/N) 

dynam_plot_data = dynam_plot_data %>%
  group_by(p,mu)%>%
  mutate(R = i- lag(i))

sim_plot_data %>%
  ggplot()+
  geom_point(aes(x=log(1-s),y=R), col= "lightblue", alpha = 0.3, size =0.5)+
  geom_point(data = dynam_plot_data, aes(x=log(1-s), y = R), col ="red", size = 0.5)+
  facet_grid(mu~p)+
  labs(y="new_infected")


simulations %>% group_by(simnum, ID, R0) %>% summarise(S = min(susceptable)) %>% 
  #mutate(round(R0,3))%>%
  ggplot(aes(x = R0, y= 1-S, col = ID, group = simnum))+
  geom_point()+
  xlim(0,5)




simulations %>%
  ggplot(aes(x = t, y = susceptable, col = p, group = simnum))+
  geom_point(size = 0.5, alpha = 0.5)+
  facet_wrap(~mu)+
  labs(y = "Proportion Susceptable", x = "Time since start of epidemic")

simulations %>%
  group_by(p,mu,simnum) %>%
  summarise(t=max(t))%>%
  ggplot(aes(x = t, col = p))+
  geom_histogram()+
  facet_grid(p~mu)+
  labs(x = "Epidemic duration")

simulations %>%
  group_by(p,mu,simnum) %>%
  summarise(S=min(susceptable))%>%
  ggplot(aes(x = S, col = p))+
  geom_histogram()+
  facet_grid(p~mu)+
  labs(x = "Final size of epidemic")


simulations %>%  
  ggplot(aes(x=susceptable, y = ratio, col = ID, group = simnum))+
  geom_point(size = 0.5)