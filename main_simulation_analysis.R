library(lubridate) # date formatting
library(quantmod) # download VIX data
library(data.table) # must have for large data
library(plyr) # some data manipulation using dlply
library(parallel) # parallel processing using mclapply
library(tictoc) # for time measurement
library(stargazer) # for TeX output
library(plm) # for OLS regression with FE
library(qs) # for efficient storage

rm(list = ls())

############################################################################################

##############################
##### GET THE DATA ###########

# load data? 
BUILD_DATA <- FALSE
# run sim for mu and Sigma?
RUN_R_sim <- FALSE
RUN_MARKET_sim <- FALSE

RUSSEL_DATA_FILE <- "russel_1000.csv"

if(BUILD_DATA) {
  
  sym_data <- "days_russ_20240101.csv"
  sym <- read.csv(sym_data,stringsAsFactors = FALSE)
  sym <- sym[order(sym$symbol),]
  
  # daily CRSP file
  crsp_file <- "CRSP_1960_2023_d.csv"
  select.var <- c("date","PERMCO","TICKER","RET")
  DT  <- fread(crsp_file,select = select.var)
  
  DT <- unique(DT[DT$TICKER %in% sym$symbol,])
  DT$date <- ymd(DT$date)
  DT$RET <- as.numeric(DT$RET)
  DT <- na.omit(DT)
  gc()
  
  # check the ticker and permco
  perm_tick <- unique(DT[,list(PERMCO,TICKER),])
  # we have 1002 tickers but almost 1800 PERMCO.. same goes for PERMNO
  uniqueN(perm_tick$TICKER)
  
  # let's check whether the data matches with the one in the sym in terms of dates
  MIN_date <- DT[,lapply(.SD,min), by = "TICKER",.SDcols = "date"]
  MAX_date <- DT[,lapply(.SD,max), by = "TICKER",.SDcols = "date"]
  min_max_date <- merge(MIN_date,MAX_date,by = "TICKER")
  names(min_max_date) <- c("symbol","start2","end2")
  
  sym2 <- merge(sym,min_max_date)
  sym2$start <- ymd(sym2$start)
  sym2$end <- ymd(sym2$end)
  
  # keep symbols with 500 days
  sym2$days2 <- sym2$end2 - sym2$start2
  drop_sym1 <- sym2[sym2$days2 <= 500,"symbol"]
  # WELL is not evident here it's still trading today but in in CRSP
  drop_sym2 <- sym2[sym2$end - sym2$end2 > 0 ,"symbol"]
  mean(sym2$start2 <= sym2$start)
  
  drop_sym <- c(drop_sym1,drop_sym2)
  # let's aggregate the data at date-ticker to make sure it's unique
  DT <- DT[!DT$TICKER %in% drop_sym,]
  DT <- DT[,lapply(.SD,mean,na.rm = TRUE), by = c("TICKER","date"),.SDcols = "RET"]
  DT[,N:= .N, by = "TICKER"]
  DT <- DT[DT$N > 500,]
  fwrite(DT,RUSSEL_DATA_FILE)
  
}

if(!BUILD_DATA) {
  DT <- fread(RUSSEL_DATA_FILE)
}


valid_syms <- unique(DT$TICKER)


##########################################
######## RANDOM MONTHLY RETURNS ##########
##########################################

N <- 100
T_size <- 240 # number of months
month_index <- sort(rep(1:T_size,21))
D <- diag(N)
e <- rep(1,N)

random_sym_ret <- function(N, seed) {
  set.seed(seed)
  symbols_sample <- sample(valid_syms,N,replace = FALSE)
  DT_sub <- data.frame(DT[DT$TICKER %in% symbols_sample,])
  
  DT_sub_list <- dlply(DT_sub,"TICKER",data.frame)
  
  sample_month_fun <- function(x) {
    # set.seed(seed)
    R <- sample(x$RET,21*T_size,replace = TRUE)
    y <- data.frame(Symbol = unique(x$TICKER),RET = R)
    y$index <- month_index
    y <- data.table(y)
    y <- y[,list(RET = prod(RET+1) - 1),by = "index"]
    y$TICKER <- unique(x$TICKER)
    return(y)
  }
  
  DT_sub_list <- lapply(DT_sub_list,sample_month_fun) 
  DT_sub <- rbindlist(DT_sub_list)
  DT_sub <- dcast.data.table(DT_sub,index~TICKER, value.var = "RET")
  R_sim <- DT_sub[,-1] # <---- for a single simulation, we have 2.5K investors
  Sigma <- var(R_sim)
  Mu <- apply(R_sim,2,mean)
  list(Mu,Sigma)
}

RDS_sim_file <- "R_sim_list.RDS"
# we only need the mu and sigma from the data. I am not sure why need to bother more
if(RUN_R_sim) {
  R_sim_list <- mclapply(1:1000,function(seed) 
    random_sym_ret(N,seed),mc.cores = detectCores())
  saveRDS(R_sim_list,RDS_sim_file)
}

if(!RUN_R_sim) {
  R_sim_list <- readRDS(RDS_sim_file)
}

gc()
##############################################################################################


################################
######## USEFUL FUNCTION #######

fun_compliant <- function(x) {
  X_E <- x$Agents_Port
  x_F <- X_E[1,]
  X_E <- X_E[-1,]
  X_L <- apply(X_E,2,function(w) sum(w[w>0]) )
  X_S <- apply(X_E,2,function(w) -sum(w[w<0]) )
  RULE <- (X_L + 1.5*X_S + x_F) <= 2
  return(mean(RULE))
} 

fun_compliant_market <- function(x) {
  w <- x$Market_Eq_Port
  x_F <- x$Market_RF
  X_L <- sum(w[w>0])
  X_S <- -sum(w[w<0])
  RULE <- (X_L + 1.5*X_S + x_F) <= 2
  return(RULE*1)
}

my_sum_fun <- function(x) {
  N <- length(x)
  Min <- min(x)
  Q_L <- quantile(x,0.05)
  Q_M <- quantile(x,0.5)
  Mean <- mean(x)
  Q_H <- quantile(x,0.95)
  Max <- max(x)
  result <- c(N,Min,Q_L,Q_M,Mean,Q_H,Max)
  return(result)
}

#############################
######## INPUTS #############

# keep r_F fixed for now
r_F_exo <- 0.05/12 # make sure it's monthly since the data is
AG_list1 <- list(A = c(1,10),G = c(0,10))
AG_list2 <- list(A = c(1,10),G = c(0,1)) # with shorter support for G


MARKET_SIM <- function(sim_i,RF_endog = TRUE,AG_list,COMPLIANT_ONLY = FALSE) {
  
  seed <- sim_i
  Mu <- R_sim_list[[sim_i]][[1]]
  Sigma <- R_sim_list[[sim_i]][[2]]
  
  A <- AG_list$A
  G <- AG_list$G
  set.seed(seed)
  {
    # the following generates the main inputs for investors
    As <- sort(runif(50,A[1],A[2]))
    Gs <- sort(runif(50,G[1],G[2]))
    
    AGs <- expand.grid(As,Gs)
    names(AGs) <- c("A","G")
    AGs$w <- runif(nrow(AGs),10^4,10^7)
    w_vec <- (AGs$w) # we have 2.5K investors in total
  }
  M <- nrow(AGs)
  # simulate gammas
  gamma_list <- lapply(1:M, function(seed) {set.seed(seed); return(rep(sample(c(0,1/N),1),N))}   )
  
  # need Omega^i inverse, which takes a couple of seconds to solve for all investors
  Omega_inv_comp <- lapply(1:M,function(i) solve(AGs$A[i]*Sigma +AGs$G[i]*D)  )
  # then we need
  # w^{i}e^{\top}\left[\Omega^{i}\right]^{-1}, which we compute as
  Omega_inv_w_e <- lapply(1:M, 
                          function(i) w_vec[i]*t(e)%*%Omega_inv_comp[[i]] )
  
  num1 <-  sum(sapply(1:M, function(i) Omega_inv_w_e[[i]]%*%(Mu +  AGs$G[i]*gamma_list[[i]]) ))
  num2 <- sum(w_vec)
  denom <- sum(sapply(1:M, function(i) Omega_inv_w_e[[i]]%*%(e) ))
  r_F_end <- (num1 - num2)/denom
  
  if(RF_endog)
    r_F <- r_F_end
  
  if(!RF_endog)
    r_F <- r_F_exo
  
  X_E_num <- lapply(1:M, 
                    function(i) Omega_inv_comp[[i]]%*%(Mu +  AGs$G[i]*gamma_list[[i]]  - r_F*e ) )
  X_E_denom <- sapply(X_E_num, sum)
  x_F <- 1 - X_E_denom
  X_E_list <- lapply(1:M,function(i) X_E_num[[i]]/X_E_denom[i]  )
  
  # then the market portfolio of the is scalled with respect to wealth
  X_E_list_adj <-  sapply(1:M, 
                          function(i) X_E_list[[i]]*(1-x_F[i])  )
  
  # aggregate how much demand there is for each asset
  X_E_list_adj <- rbind(x_F,X_E_list_adj)
  
  # we need to rule out non-compliant portfolios
  if(COMPLIANT_ONLY) {
    X_E_c <- X_E_list_adj[-1,]
    x_F_c <- X_E_list_adj[1,]
    X_L <- apply(X_E_c,2,function(w) sum(w[w>0]) )
    X_S <- apply(X_E_c,2,function(w) -sum(w[w<0]) )
    RULE <- (X_L + 1.5*X_S + x_F_c) <= 2
    X_E_list_adj <- X_E_list_adj[,RULE]
    gamma_list <- gamma_list[RULE]
    w_vec <- w_vec[RULE]
  }
    
  X_M <- (X_E_list_adj%*%w_vec)/(sum(w_vec))
  X_E <- as.matrix(X_M[-1,])
  x_F_eq <- X_M[1,]
  list(Market_Eq_Port = X_E, Market_RF =  x_F_eq, 
       r_F_end =  r_F, PARM = AGs, 
       Agents_Port = X_E_list_adj,
       gamma_list = gamma_list)
}

RDS_market_sim_file <- "market_sim_list.qs"

# BASELINE
if(RUN_MARKET_sim) {
  tic("my function")
  MARKET_LIST_SIM <- mclapply(1:10^3,function(seed) 
    MARKET_SIM(seed,RF_endog = FALSE,AG_list = AG_list1,
               COMPLIANT_ONLY = TRUE),
    mc.cores = detectCores() - 1)
  toc()
  # qsave(MARKET_LIST_SIM,RDS_market_sim_file)
}


if(!RUN_MARKET_sim) {
  MARKET_LIST_SIM <- qread(RDS_market_sim_file)
}

######################################################################################################################

#########################################################################
######## TABLE 1 and 2: SUMMARY STATS AND COMPLIANT PORTFOLIOS ##########

SUM_STATS_TABLE <- function(MARKET_LIST_SIM) {
  # prop of compliant portfolios across investors
  prop_compliant <- sapply(MARKET_LIST_SIM,fun_compliant)
  # market portfolio
  market_compliant <- sapply(MARKET_LIST_SIM,fun_compliant_market)
  
  summary(prop_compliant)
  r_F_eq <- sapply(MARKET_LIST_SIM,function(x) x$r_F_end )*12
  
  # look at the x_F position
  x_F_min <- sapply(MARKET_LIST_SIM,function(x) min(x$Agents_Port[1,]))
  x_F_mean <- sapply(MARKET_LIST_SIM,function(x) mean(x$Agents_Port[1,]))
  x_F_max <- sapply(MARKET_LIST_SIM,function(x) max(x$Agents_Port[1,]))
  x_F_sd <- sapply(MARKET_LIST_SIM,function(x) sd(x$Agents_Port[1,]))
  
  ###### look at the r_F as a function of A
  A_mean <- sapply(MARKET_LIST_SIM,function(x) max(x$PARM$A))
  lm_rf1 <- lm(r_F_eq ~ x_F_sd)
  
  # look at long positions
  Total_Long <- lapply(MARKET_LIST_SIM,
                       function(x) (apply(x$Agents_Port[-1,],2,function(x) sum(x[x>0]) )) ) 
  Total_Long_min <- sapply(Total_Long,min)
  Total_Long_mean <-  sapply(Total_Long,mean)
  Total_Long_max <- sapply(Total_Long,max)
  
  # look at short positions
  Total_Short <- lapply(MARKET_LIST_SIM,
                        function(x) (apply(x$Agents_Port[-1,],2,function(x) sum(x[x<0]) )) ) 
  Total_Short_min <- sapply(Total_Short, min)
  Total_Short_mean <- sapply(Total_Short, mean)
  Total_Short_max <- sapply(Total_Short, max)
  
  # let's stack altogether in a data.frame and report summary stats 
  ds_sum <- data.frame(prop_compliant,market_compliant,
                       r_F_eq,
                       x_F_min,x_F_mean,x_F_max,
                       Total_Long_min,Total_Long_mean,Total_Long_max,
                       Total_Short_min,Total_Short_mean,Total_Short_max)
  
  X_labels <- c("\\% Compliant (investor)", "\\% Compliant (market)",
                "\\$r_F\\$",
                "\\$x_F\\$ (min)", "\\$x_F\\$ (mean)", "\\$x_F\\$ (max)",
                "Total Long (min)", "Total Long (mean)", "Total Long (max)",
                "Total Short (min)", "Total Short (mean)", "Total Short (max)")
  X <- t(apply(ds_sum,2,my_sum_fun))
  rownames(X) <- X_labels
  
  return(X)
}
Y_labels <- c("N.","Min","P5","Median","Mean","P95","Max")

X <- SUM_STATS_TABLE(MARKET_LIST_SIM)
stargazer(X, covariate.labels = c("",Y_labels),
          summary = FALSE)

################################################
#### TABLE 3: Correlation Among Factors ########

# let's compute each factor based on the simulation

compute_factors_fun <- function(sim_i,MARKET_LIST_SIM) {
  Mu <- R_sim_list[[sim_i]][[1]]
  Sigma <- R_sim_list[[sim_i]][[2]]
  SIM_I <- MARKET_LIST_SIM[[sim_i]] 
  r_F <- SIM_I$r_F_end
  Mu_r_F_vec <- Mu - r_F
  X_M <- SIM_I$Market_Eq_Port
  Sigma_X <- Sigma%*%X_M
  sigma_sq_M <- as.numeric(t(X_M)%*%Sigma%*%X_M)
  beta_vec <- Sigma_X/sigma_sq_M
  
  AG_i <-  SIM_I$PARM
  w_vec <- AG_i$w
  A_vec <- AG_i$A
  G_vec <- AG_i$G
  theta_vec <- G_vec/A_vec
  
  X_E_i <- SIM_I$Agents_Port[-1,]
  
  gamma_list_i <- SIM_I$gamma_list
  gamma_mat <- do.call(cbind, gamma_list_i)
  delta_vec1 <- X_E_i - gamma_mat
  # I did confirm sweep manually
  delta_vec2 <- sweep(delta_vec1,2,w_vec*theta_vec,"*")
  delta_vec <- rowSums(delta_vec2)/sum(w_vec)
  
  # gamma_m
  gamma_m_vec <- sweep(gamma_mat,2,w_vec,"*")
  gamma_m_vec <- rowSums(gamma_m_vec)/sum(w_vec)
  delta_m_vec <- X_M - gamma_m_vec
  
  # \Delta vec
  Delta_vec <- sweep(delta_vec1,2,w_vec*(theta_vec - 1),"*")
  Delta_vec <- rowSums(Delta_vec)/sum(w_vec)
  ds_i <- data.frame(sim_i,Mu_r_F_vec,beta_vec,delta_vec,delta_m_vec,Delta_vec)
  return(ds_i)
}

ds_factor_list <- mclapply(1:length(MARKET_LIST_SIM),
                           function(i) compute_factors_fun(i,MARKET_LIST_SIM),mc.cores = detectCores()-1)

DS_factor <- rbindlist(ds_factor_list)
cor(DS_factor[,c("beta_vec","delta_m_vec","Delta_vec"),with = FALSE])

#########################################
#### TABLE 5: REGRESSION RESULTS ########

lm1 <- lm(Mu_r_F_vec ~ beta_vec,data = DS_factor)
lm2 <- lm(Mu_r_F_vec ~ delta_m_vec,data = DS_factor)
lm3 <- lm(Mu_r_F_vec ~ beta_vec + delta_m_vec,data = DS_factor)

lm_list1 <- list(lm1,lm2,lm3)
stargazer(lm_list1)


lm1 <- plm(Mu_r_F_vec ~ beta_vec , index = c("sim_i"),model = "within",data = DS_factor)
lm2 <- plm(Mu_r_F_vec ~ delta_m_vec,index = c("sim_i"),model = "within",data = DS_factor)
lm3 <- plm(Mu_r_F_vec ~ beta_vec + delta_m_vec,index = c("sim_i"),model = "within",data = DS_factor)
lm_list2 <- list(lm1,lm2,lm3)
stargazer(lm_list2,
          covariate.labels = c("$\\beta$","$\\delta_m$"))

# report the average RSQ and p value from all regressions
lm_fun_i <- function(ds_i) {
  lm1 <- lm(Mu_r_F_vec ~ beta_vec,data = ds_i)
  lm2 <- lm(Mu_r_F_vec ~ delta_m_vec,data = ds_i)
  lm3 <- lm(Mu_r_F_vec ~ beta_vec + delta_m_vec,data = ds_i)
  list(lm1,lm2,lm3)
}


lm_list_all <- lapply(ds_factor_list, lm_fun_i)
lm_list_all1 <- lapply(lm_list_all,function(x) x[[1]])
lm_list_all2 <- lapply(lm_list_all,function(x) x[[2]])
lm_list_all3 <- lapply(lm_list_all,function(x) x[[3]])

extract_beta_Rsq <- function(lm_i) {
  lm_i_sum <- summary(lm_i)
  rsq <- lm_i_sum$adj.r.squared
  B <- lm_i_sum$coefficients[-1,c(1,4)]
  list(rsq = rsq,beta = B)
}

LMS1 <- lapply(lm_list_all1,extract_beta_Rsq)
RSQ1 <- sapply(LMS1, function(x) x$rsq)
B1 <- t(sapply(LMS1, function(x) x$beta))
R1 <- apply(data.frame(RSQ1,B1),2,mean)
SE1 <- apply(data.frame(RSQ1,B1),2,sd)
round(R1,4)
round(SE1/sqrt(1000),4)


LMS1 <- lapply(lm_list_all2,extract_beta_Rsq)
RSQ1 <- sapply(LMS1, function(x) x$rsq)
B1 <- t(sapply(LMS1, function(x) x$beta))
R2 <- apply(data.frame(RSQ1,B1),2,mean)
SE2 <- apply(data.frame(RSQ1,B1),2,sd)
round(R2,4)
round(SE2/sqrt(1000),4)


LMS1 <- lapply(lm_list_all3,extract_beta_Rsq)
RSQ1 <- sapply(LMS1, function(x) x$rsq)
B1 <- t(sapply(LMS1, function(x) x$beta))
R3 <- apply(data.frame(RSQ1,B1),2,mean)
SE3 <- apply(data.frame(RSQ1,B1),2,sd)

round(R3,4)
round(SE3/sqrt(1000),4)


