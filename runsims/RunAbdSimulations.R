
S <- 25

mu_R <- 1.5
sigma_R <- 0

r <- rnorm(S, mean = mu_R, sd = sigma_R)

A <- BuildA(S = S, mu_A = 0, sigma_A = 1, rho_A = 0)
diag(A) <- 0
A[A != 0] <- A[A != 0] / sd(A[A != 0]) / sqrt(S)
A[A != 0] <- A[A != 0] - mean(A[A != 0])

B <- array(rnorm(S*S*S, mean = 0, sd = 1), c(S, S, S))
B <- NoSelfHOIs(B)

B[B != 0] <- B[B != 0] / sd(B[B != 0]) / S
B[B != 0] <- B[B != 0] - mean(B[B != 0])

mu_len <- 2
sigma_len <- 20

mu_As <- seq(-1, -4, length.out = mu_len)
sigma_As <- seq(0.01, 0.6, length.out = sigma_len) # Code breaks if sigma = 0

mu_Bs <- seq(-1, -3, length.out = mu_len)
sigma_Bs <- seq(0.01, 0.6, length.out = sigma_len)  # Code breaks if sigma = 0

end_time <- 500
time_step <- end_time

out_data <- data.frame(Mu = c(), Sigma = c(), SpID = c(), Abd = c(), Int = c())

for(mu in 1:mu_len) {
  print(mu)
  for(sigma in 1:sigma_len) {
    print(sigma)
    
    curB <- B
    curB[curB != 0] <- curB[curB != 0] * sigma_Bs[sigma]
    curB[curB != 0] <- curB[curB != 0] + mu_Bs[mu] / S^2
    #print(mean(curB[curB != 0]));print(sd(curB[curB != 0]))
    
    curA <- A
    curA[curA != 0] <- curA[curA != 0] * sigma_As[sigma]
    curA[curA != 0] <- curA[curA != 0] + mu_As[mu] / S
    #print(mean(curA[curA != 0]));print(sd(curA[curA != 0]))
    #print(sum(curA[curA != 0])); print(sum(curB[curB != 0]))
    
    hoi_A <- A
    hoi_A[hoi_A != 0] <- hoi_A[hoi_A != 0] * 0.1
    hoi_A[hoi_A != 0] <- hoi_A[hoi_A != 0] + -1 / S
    pars <- list(S = S,
                 r = r,
                 d = rep(1, times = S),
                 A = matrix(0, nrow = S, ncol = S),
                 B = curB,
                 h = 0)
    ini_state <- runif(S)
    
    times <- seq(0, end_time, by = time_step)
    out_dyn <- as.data.frame(ode(ini_state, times, Dynamics, pars))  
    eq_abds <- as.numeric(out_dyn[nrow(out_dyn),-1])
    hoi_data <- cbind(data.frame(Mu = mu_Bs[mu], Sigma = sigma_Bs[sigma]), data.frame(SpID = 1:S, Abd = eq_abds, Frac = sum(curB > 0) / (S^3 - S), Int = "HOI"))
    
    pars <- list(S = S,
                 r = r,
                 d = rep(1, times = S),
                 A = curA,
                 B = array(0, c(S, S, S)),
                 h = 0)
    ini_state <- runif(S)
    out_dyn <- IntegrateDynamics(inistate = ini_state,
                                 pars = pars,
                                 endtime = end_time,
                                 timestep = time_step,
                                 fn = Dynamics)
    eq_abds <- as.numeric(out_dyn[nrow(out_dyn),-1])
    pw_data <- cbind(data.frame(Mu = mu_As[mu], Sigma = sigma_As[sigma]), data.frame(SpID = 1:S, Abd = eq_abds, Frac = sum(curA > 0) / (S^2 - S), Int = "PW"))
    
    out_data <- rbind(out_data, pw_data, hoi_data)
    
  }
}
