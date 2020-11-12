

N=1000

RMSE = function(beta, data){
  v = 1/20
  SE = 0
  s_hat = .999
  i_hat = 0.001
  r_hat = 0
  tmax = nrow(data)
  for(t in 1:tmax){
    ds_hat =   -beta*s_hat*i_hat
    di_hat = beta*s_hat*i_hat -v*i_hat
    dr_hat = v*i_hat
    s_hat = s_hat + ds_hat
    i_hat = i_hat + di_hat
    r_hat = r_hat + dr_hat
    
    SE = SE + (s_hat - data$s[t])^2+
      (i_hat - data$i[t])^2 +
      (r_hat - data$r[t])^2
  }
  return(sqrt(SE/tmax))
}

fit_data = simulations %>%
  select(s = susceptable, i = infected, r = recovered, t = t, p = p, mu = mu)%>%
  mutate(s = s/N, i = i/N, r=r/N)

sim_pars = fit_data %>%
  select(mu,p) %>%
  unique()

SIR_fits = mapply(function(p,mu){
   data = fit_data %>% filter(mu == mu, p == p)
   
   
   RMSE = function(beta){
     v = 1/20
     SE = 0
     s_hat = .999
     i_hat = 0.001
     r_hat = 0
     tmax = nrow(data)
     for(t in 1:tmax){
       ds_hat =   -beta*s_hat*i_hat
       di_hat = beta*s_hat*i_hat -v*i_hat
       dr_hat = v*i_hat
       s_hat = s_hat + ds_hat
       i_hat = i_hat + di_hat
       r_hat = r_hat + dr_hat
       
       SE = SE + (s_hat - data$s[t])^2+
         (i_hat - data$i[t])^2 +
         (r_hat - data$r[t])^2
     }
     return(sqrt(SE/tmax))
   }
   fit = optim(par = c(beta = mu*p), fn = RMSE)
   return(c(fit = fit,mu = mu,p = p))
  },
  sim_pars$p, sim_pars$mu)


