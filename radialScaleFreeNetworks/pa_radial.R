# m in pref att is the number of edges added for each node
# does the probablility of being infected proportional to the number of exposures (i.e. those connected to it that are infected)
# m in the prefernential attachment step will be the proportion of people within the radius that become infected if we want a physical analouge
# m could correspond to the effect of masks, and proabability of infection along any given connection would correspond to the effect of masks

alpha = 2.5
d=2
beta = d*(alpha - 1)+1


L =100
V= L^2

draw_R = function(beta){
  R=seq(1.001, L+0.01, .001)
  R_dist = R^(-beta)*2*pi*R
  R_dist = R_dist/sum(R_dist)
  sample(R, 1, prob = R_dist)
}

#R_dist = data.frame(R = replicate(1000, draw_R(beta)) )

# R_dist %>%
#   ggplot(aes(x=R))+
#   geom_histogram()

pa_connect = function(neighbors, adj_matrix, m){
   degree = sapply(neighbors, function(x){sum(adj_matrix[x,])})
   degree = degree/sum(degree)
   sample(neighbors, min(m, length(neighbors)), prob = degree)
}



check_R = function(R, vertex){
  apply(location_matrix, c(1,2), function(v){
    vertex_loc = which(location_matrix==vertex, arr.ind=TRUE)
    v_loc = which = which(location_matrix==v, arr.ind=TRUE)
    sum((vertex_loc - v_loc)^2) < R 
  })
}


adj_matrix = diag(V)
location_matrix = matrix(1:V, nrow= L, ncol = L)
m = 3


vertex_order = sample(1:V, V)

startTime = Sys.time()
sapply(vertex_order, function(vertex){
  possible_neighbors = location_matrix[which(check_R(draw_R(beta), vertex), arr.ind = T)]
  #connected_neighbors = location_matrix[which(check_R(draw_R(beta), vertex), arr.ind = T)]
  connected_neighbors = pa_connect(possible_neighbors, adj_matrix, m)
  sapply(connected_neighbors, function(x){
      adj_matrix[vertex, x] <<-  1
      adj_matrix[x, vertex] <<-  1
  })
})
Sys.time() - startTime
#save(adj_matrix, file = "adj_matrix_alpha2.5_L25_m3.RData")
#load("adj_matrix_alpha2.5_L100_m3.RData")

deg_dist = data.frame(deg = apply(adj_matrix, 1, sum)) %>%
  mutate(deg = (deg-1))%>%
  group_by(deg)%>%
  summarise(dist = n())

deg_dist %>%
  ggplot()+
  geom_histogram(aes(x=deg, y=dist), stat = "identity", binwidth = 1)+
  geom_smooth(formula = y~x^-alpha, aes(x=deg, y=dist))


#library(igraph)
deg_dist %>%
  #filter(deg>5)%>%
  pull(dist)%>%
  power.law.fit()


# plot adjacency matrix
data.frame(adj_matrix) %>%
  mutate(node2 = 1:V)%>%
  gather(key = "node1", val = "connect", - node2)%>%
  mutate(node1 = sub("X", "", node1)) %>%
  mutate_all(as.numeric)%>%
  ggplot(aes(x=node1, y = node2, fill = connect))+
  geom_tile()
