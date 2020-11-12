library(tidyverse)
files =list.files("C:/Users/isabe/OneDrive/Documents/2020/SIR_network/8_27/network")


path = getwd()
degree = list()
power_fit = list()
networks = list()

for(f in 1:length(files)){
  filename = files[f]
  load(paste0(path,"/8_27/network/",filename))
  networks[[f]] = adj_matrix
  degree[[f]] = apply(adj_matrix, 1, sum)
  power_fit[[f]] = power.law.fit(degree[[f]])
  print(f)
}


files =list.files("C:/Users/isabe/OneDrive/Documents/2020/SIR_network/8_27/SIR")
path = getwd()

for(f in 1:length(files)){
  filename = files[f]
  load(paste0(path,"/8_27/SIR/",filename))
  if(f == 1){
    sim_restults = pop_counts
  } else{
    sim_restults = bind_rows(sim_restults,pop_counts)
  }
  print(f)
}


sim_restults %>%
  mutate(R0 = max_degree*0.5*14)%>%
  ggplot(aes(x=num_suseptable, y = R))+
  geom_point()

x=seq(-3,3,.01)

x[200:300]
df = data.frame(x=x, y= dnorm(x)) 
df2 =df[400:500,]

df %>%
  ggplot(aes(x=x,y=y))+
  geom_density(stat = "identity",color= "lightblue", fill = "lightblue")+
  geom_density(data = df2, stat = "identity",color= "darkblue", fill = "darkblue")+
  theme_minimal()+
  scale_x_continuous(breaks = seq(-3,3,1), labels = seq(-3,3,1))+
  labs(x='',y='')+ 
  theme(axis.text = element_text(size = 20))
  

x=0:5

df = data.frame(x=x, y= dbinom(x,5,.5)) 
df2 =df[2:4,]

df %>%
  ggplot(aes(x=x,y=y))+
  geom_histogram(stat = "identity",color= "lightblue", fill = "lightblue")+
  geom_histogram(data = df2, stat = "identity",color= "darkblue", fill = "darkblue")+
  theme_minimal()+
  scale_x_continuous(breaks = seq(0,5,1), labels = seq(0,5,1))+
  labs(x='',y='')+ 
  theme(axis.text = element_text(size = 20))
