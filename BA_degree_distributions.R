library(tidyverse)
library(igraph)
library(Pareto)
##decide on parameters for Barbassi-Albert model heterogeneous graph. Has exponent 3 for degree distribution


#xm = m*sqrt((N-m0)/N) #parameter of theoretical pareto distribution
#alpha = 2
#theoretical mean and meadian for degree distribution
#theo_mu = round(xm^alpha*alpha/ xm^(alpha-1))
#theo_median = round(xm* 2^(1/alpha))

#potential interpreations of R0
#expand.grid(theo_mu = theo_mu, probs = probs) %>%
#mutate(R0_mean = theo_mu*probs*t)

#expand.grid(theo_median = theo_median, probs = probs) %>%
#mutate(R0_med = theo_median*probs*t)

m=10
m0=5
N = 1000
xm = m*sqrt((N-m0)/N)
alpha = 2

#sampling from theoretical distribution
theo_mu = xm^alpha*alpha/ xm^(alpha-1)
theo_median = xm* 2^(1/alpha)

theo_degree =rPareto(10000,xm,alpha)


#simulating
pa_connect = function(N,m,m0){
  a = make_empty_graph(m0,directed = F)
  pa_network = sample_pa(n = N, m=m, start.graph = a , directed = F)
  return(pa_network)
}
pa_network = pa_connect(N,m,m0)
power.law.fit(degree(pa_network))

#compare
compare_degree_distributions = data.frame(
  #sim_degree = degree(pa_network),
  theo_degree = theo_degree
) %>%
  gather(key = 'source', value = "degree")

ggplot(compare_degree_distributions)+
  geom_histogram(aes(x = degree, fill = source), binwidth = 1, alpha = .5)+
  scale_x_continuous(limits = c(0,500), breaks = seq(0,500,20), labels = seq(0,500,20))

compare_degree_distributions %>%
  group_by(source)%>%
  summarise(min = min(degree), max = max(degree), mu = mean(degree), med = median(degree))

## generate network 


pa_connect(N,m,m0)

