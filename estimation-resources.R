
library(gtools)
library(isotone)
library(tidyverse)
library(ggplot2)
library(rankrate)
library(arrangements)
library(matrixStats)
library(reshape2) # data reformatting
library(pander)   # creating nice tables

componentwise_dwmb <- function(M, alpha, p_q, theta_r, pi_r, pi_0_r, x_q, log) {
  J <- length(p_q)
  R <- J
  
  if(log) {
    C_psi_comp_log <- alpha*mallows_norm(theta_r, J, R, log=T) - mallows_norm(theta_r*alpha, J, R, log=T)
    
    C_prod_comp_log_part1 <- M*(alpha-1)*sum(log(sapply(p_q, function(p_j){1-p_j})))
    C_prod_comp_log_part2 <- sum(log(sapply(p_q, function(p_j){sum((choose(M, (0:M))*(p_j/(1-p_j))^(0:M))^(1-alpha))})))
    C_prod_comp_log <- C_prod_comp_log_part1 - C_prod_comp_log_part2
    
    original_mallows_log <- alpha*dmall(pi_r, pi_0_r, theta_r,log=T)
    original_binomial_log <- (1-alpha)*(sum(mapply(function(x_j, p_j){dbinom(x_j, M, p_j, log=T)}, x_q, p_q)))
    C_log <- C_psi_comp_log + C_prod_comp_log
    C_log <- ifelse(alpha==0, log(1/factorial(J)), C_log)
    
    result <- C_log + original_mallows_log + original_binomial_log
    return(result)
  } else {
    C_psi_comp <- ((mallows_norm(theta_r, J, R, log=F))^(alpha))/(mallows_norm(theta_r*alpha, J, R, log=F))
    C_prod_comp <- prod(sapply(p_q, function(p_j){(((1-p_j)^(M*(alpha-1)))/(sum((choose(M, (0:M))*(p_j/(1-p_j))^(0:M))^(1-alpha))))} ))

    original_mallows <- dmall(pi_r, pi_0_r, theta_r,log=F)
    original_binomial <- (prod(mapply(function(x_j, p_j){dbinom(x_j, M, p_j, log=F)}, x_q, p_q)))
    original_mallows_alpha <- original_mallows^(alpha)
    original_binomial_alpha <- original_binomial^(1-alpha)
    
    C <- C_psi_comp*C_prod_comp
    C <- ifelse(alpha==0, 1/factorial(J), C)
    
    result <- C*original_mallows_alpha*original_binomial_alpha
    
    ranking_component <- C_psi_comp*original_mallows_alpha
    rating_component <- C_prod_comp*original_binomial_alpha
    result_updated <- ranking_component*rating_component
    
    return(list(ranking_component, rating_component, result_updated))
  }
}
componentwise_dwmb <- Vectorize(componentwise_dwmb)

mallows_norm <- function(gamma, J, R, log) {
  if(log) {
    result <- sum(log(sapply(1:R, function(j){(1-exp(-gamma*(J-j+1)))/(1-exp(-gamma))})))
    result <- ifelse(gamma==0, 0, result) # Set result to zero for input of zero
    return(result)
  } else {
    result <- prod(sapply(1:R, function(j){(1-exp(-gamma*(J-j+1)))/(1-exp(-gamma))}))
    result <- ifelse(gamma==0, 1, result) # Set result to one for input of zero
    return(result)
  }
}
mallows_norm <- Vectorize(mallows_norm)

ASTAR_wmb <- function(rankings,ratings,M,alpha,keep_nodes=FALSE,ties_enabled=TRUE,use_log=FALSE){
  
  # calculate constants
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings) != c(I,J))){stop("rankings and ratings must have same dimension")}
  if(is.null(attr(rankings,"assignments"))){
    attr(rankings,"assignments") <- matrix(TRUE,nrow=I,ncol=J)
  }
  
  # calculations for updating theta
  Ji <- apply(attr(rankings,"assignments"),1,sum) #each judge's total number of assignments
  Ri <- apply(rankings,1,function(pi){sum(!is.na(pi))}) #each judge's size of ranking
  i_pi <- which(Ri>0 & Ji>0) #which judges provided rankings
  
  t1 <- function(theta,D){theta*D/I} #first log probability term
  t2 <- function(theta){sum(psi(theta,Ji[i_pi],Ri[i_pi],log=T))/I} #second log probability term
  
  get_theta <- function(D){ # function used to optimize for theta conditional on some order
    # D is the minimum summed kendall distance between some central ordering and the observations.
    # As D increases, theta decreases monotonically and the cost increases monotonically.
    opt_theta <- function(theta){t1(theta,D=D)+t2(theta)}
    result <- optimize(f=opt_theta,interval=c(1e-8,10^8))
    return(list(thetahat = result$minimum,objective = result$objective))
  }
  
  # calculations for updating p
  c1 <- apply(ratings,2,function(x){sum(x,na.rm=T)})/I #set of constants in log binomial density
  c2 <- apply(M-ratings,2,function(x){sum(x,na.rm=T)})/I #set of constants in log binomial density
  
  t3_new <- function(p){1}
  t3_gradient_new <- function(p){1}
  
  if(use_log) {
      # Use logSumExp for better accuracy
      t3_new <- function(p){
      p <- ifelse(p<=0, 1e-10, p)
      p <- ifelse(p>=1, 1-1e-10, p)
      result <- sum((alpha-1)*c1*( log(p) - log(1-p) )) + 
        sum(unlist( lapply(p, function(p){logSumExp( (1-alpha)*( log(choose(M,0:M)) + (0:M)*( log(p) - log(1-p) ) ) )}) )) 
      return(result)
    }
    
    t3_gradient_new <- function(p) {
      p <- ifelse(p<=0, 1e-10, p)
      p <- ifelse(p>=1, 1-1e-10, p)
      
      part1 <- log((alpha-1)*c1)-(log(p)+log(1-p))
      
      # Use logSumExp for better accuracy
      part21 <- Vectorize(function(pj){ -1*logSumExp( (1-alpha)*( log(choose(M,0:M)) + (0:M)*( log(p) - log(1-p) ) ) ) })
      part22 <- Vectorize(function(pj){ ( log(1-alpha) + logSumExp( log(0:M)+ (1-alpha)*log(choose(M,0:M)) + ((0:M)*(1-alpha) - 1)*log(pj) - ((0:M)*(1-alpha) + 1)*log(1-pj) ) ) })
      result <- exp(part1) + exp(part21(p))*exp(part22(p))
      
      return(result)
    }
    
  } else {
    t3_new <- function(p){
      p <- ifelse(p<=0, 1e-10, p)
      p <- ifelse(p>=1, 1-1e-10, p)
      result <- sum((alpha-1)*c1*log(p/(1-p))) + 
        sum(unlist(lapply(p, function(p){ log(sum( (choose(M,0:M)*(p/(1-p))^(0:M))^(1-alpha))) }))) 
      return(result)
    }
    
    t3_gradient_new <- function(p) {
      p <- ifelse(p<=0, 1e-10, p)
      p <- ifelse(p>=1, 1-1e-10, p)
      part1 <- (alpha-1)*c1/(p*(1-p))
      part21 <- Vectorize(function(pj){ ( (sum( (choose(M,0:M)*(pj/(1-pj))^(0:M))^(1-alpha))) )^(-1) })
      part22 <- Vectorize(function(pj){ (1-alpha)*sum( (0:M)*(choose(M,0:M)^(1-alpha))*( ((pj)^((0:M)*(1-alpha) - 1) )/((1-pj)^((0:M)*(1-alpha) + 1) ) ) ) })
      result <- part1 + part21(p)*part22(p)
      return(result)
    }
  }
  
  get_p <- function(order){ # constrained optimization function
    # order is the top-R portion of a central ordering (i.e., if order = c(5,4), then object 5 must be in first place
    # and object 4 must be in second place, the remaining objects must be in places 3:J (any order among those) ).
    
    # initializing the order
    R <- length(order)
    if(R<J){unordered <- setdiff(1:J,order)}else{unordered <- c()}
    order_and_unordered <- c(order,unordered)
    
    # calculation of linear constraints
    ui <- rbind(diag(1,nrow=J),diag(-1,nrow=J))
    ci <- rep(c(0,-1),each=J)
    if(R>1){for(place in 2:R){
      ui_addition <- rep(0,J)
      ui_addition[c(order[place],order[place-1])] <- c(1,-1)
      ui <- rbind(ui,ui_addition)
      ci <- c(ci,0)
    }}
    if(R<J){for(place in (R+1):J){
      ui_addition <- rep(0,J)
      ui_addition[c(order_and_unordered[place],order[R])] <- c(1,-1)
      ui <- rbind(ui,ui_addition)
      ci <- c(ci,0)
    }}
    
    # get intelligent starting values via isotonic regression
    meanratings <- c1/(c1+c2)
    order_meanratings <- c(order,setdiff(order(meanratings),order))
    start <- gpava(1:J,meanratings[order_meanratings])
    start <- start$x[order(order_meanratings)]
    start[start==0] <- 1e-7
    start[start==1] <- 1-1e-7
    
    # optimize for p conditional on an ordering of parameters
    result <- suppressWarnings(constrOptim(theta = start,f = t3_new, grad = t3_gradient_new, ui=ui,ci=ci-1e-8, mu=1e-8,control = list(reltol = 1e-10)))
    
    return(list(phat = result$par,objective = result$value))
  }
  
  # get D (necessary for updating theta)
  Q <- getQ(rankings,I,J)*I
  get_D <- function(order){
    S <- setdiff(1:J,order)
    D <- 0
    if(length(S)>=2){
      D <- D +
        sum(apply(gtools::combinations(length(S),2,S),1,function(uv){
          min(Q[uv[1],uv[2]],Q[uv[2],uv[1]])
        }))
    }
    for(i in 1:length(order)){D <- D + sum(Q[setdiff(1:J,order[1:i]),order[i]])}
    return(D)
  }
  
  # create matrix of open nodes to continue exploring
  open <- matrix(data=c(1:J,rep(NA,J*(J-1))),nrow=J,ncol=J)
  open[,J] <- unlist(lapply(1:J,function(order){
    round(get_theta(D=get_D(order))$objective + get_p(order)$objective,4)
  }))
  num_nodes <- J
  
  # create matrix to store complete paths
  result <- matrix(NA,nrow=0,ncol=J)
  
  # search through the current tree
  while(nrow(open)>0){
    which_current <- which.min(open[,J])
    current <- c(na.exclude(open[which_current,1:(J-1)]))
    if(length(current)==(J-1)){
      # print("found a solution!")
      curr_min <- open[which_current,J]
      which_mins <- which(open[,J] == curr_min)
      result <- rbind(result,c(current,curr_min))
      open <- open[-which_current,]
      if(length(which_mins)==1){break() #only stop if not more possible solutions
      }else{#print("still going!")
      }
    }else{
      neighbors <- setdiff(1:J,current)
      for(neighbor in neighbors){
        order <- c(current,neighbor)
        totalcostheuristic <- round(get_theta(D=get_D(order))$objective + get_p(order)$objective,4)
        open <- rbind(open,c(order,rep(NA,J-1-length(order)),totalcostheuristic))
      }
      num_nodes <- num_nodes + length(neighbors)
      open <- open[-which_current,]
    }
  }
  
  # be sure all solutions are equal in their totalcostheuristic (in case of multiple solutions)
  result <- result[which(result[,J]==min(result[,J])),,drop=FALSE]
  result[,J] <- apply(result[,1:(J-1),drop=FALSE],1,function(pi){setdiff(1:J,pi)}) #complete the rankings
  
  # Check whether ties are enabled
  if((nrow(result)>1) && (ties_enabled)){
    message("There's a tie! Results are shown as a matrix to give multiple solutions.")
    results <- t(apply(result,1,function(curr_node){
      c(get_p(curr_node)$phat,get_theta(D=get_D(curr_node))$thetahat)
    }))
    if(keep_nodes){
      return(list(pi0=result,
                  p=results[,1:J],
                  theta=results[,J+1,drop=FALSE]/alpha,
                  num_nodes=num_nodes,
                  nodes = open))
    }else{
      return(list(pi0=result,
                  p=results[,1:J],
                  theta=results[,J+1,drop=FALSE]/alpha,
                  num_nodes=num_nodes))
    }
  }else{
    curr_node <- result[1,]
    results <- c(get_p(curr_node)$phat,get_theta(D=get_D(curr_node))$thetahat)
    if(keep_nodes){
      return(list(pi0=curr_node,
                  p=results[1:J],
                  theta=results[J+1]/alpha,
                  num_nodes=num_nodes,
                  nodes = open))
    }else{
      return(list(pi0=curr_node,
                  p=results[1:J],
                  theta=results[J+1]/alpha,
                  num_nodes=num_nodes))
    }
  }
  
}
ASTAR_wmb_alpha_vectorized <- Vectorize(FUN=ASTAR_wmb, vectorize.args="alpha")

ASTAR_plot <- function(rankings,ratings,M,n,keep_nodes=FALSE,parameter="p",low=0,high=1,use_log=FALSE,smooth=F,return_plot=F) {
  # Determine number of observations and number of objects
  I <- nrow(rankings)
  J <- ncol(rankings)
  
  # Create sequence 0, 1/n, 2/n, ..., n
  alpha_grid <- seq(from=low, to=high, by=(high-low)/n)
  
  # Label these alpha values with indices i=0,1,...,n (NOTE: 0-indexed)
  i <- seq(from=0, to=n, by=1)
  
  # Vectorize over grid of alpha values, ignoring ties to prevent errors later
  fit_list <- ASTAR_wmb_alpha_vectorized(rankings, ratings, M, 
                                         alpha=alpha_grid, 
                                         keep_nodes=keep_nodes, 
                                         ties_enabled=FALSE,
                                         use_log=use_log)
  
  # Produces a list of all results. For i=0,1,...,n, 
  # fit_list[[1+4*i]]: pi0[[alpha=i]],
  # fit_list[[2+4*i]]: p[[alpha=i]],
  # fit_list[[3+4*i]]: theta[[alpha=i]],
  # fit_list[[4+4*i]]: num_nodes[[alpha=i]].
  
  # Create striding vectors for indexing into the results
  pi0_stride <- seq(from=1, to=4*(n+1), by=4)
  p_stride <- seq(from=2, to=4*(n+1), by=4)
  theta_stride <- seq(from=3, to=4*(n+1), by=4)
  num_nodes_stride <- seq(from=4, to=4*(n+1), by=4)
  
  # Obtain values of each type
  pi0s <- fit_list[pi0_stride]
  ps <- fit_list[p_stride]
  thetas <- fit_list[theta_stride]
  nums_nodes <- fit_list[num_nodes_stride]
  
  # Create a vector of labels R1, R2, ..., RJ
  pi0_labels <- gsub(" ", "", paste(rep("R", times=J), seq(from=1,to=J,by=1)))
  # Create a vector of labels P1, P2, ..., PJ
  p_labels <- gsub(" ", "", paste(rep("P", times=J), seq(from=1,to=J,by=1)))
  
  # Create pi0s dataframe
  pi0_df <- as.data.frame(t(as.data.frame(pi0s)))
  names(pi0_df) <- pi0_labels
  
  # Create ps dataframe
  ps_df <- as.data.frame(t(as.data.frame(ps)))
  names(ps_df) <- p_labels
  
  # Create likelihoods data frame
  M_vec <- rep(M, length(alpha_grid))
  rankings_vec <- rep(rankings, length(alpha_grid))
  ratings_vec <- rep(ratings, length(alpha_grid))
  ranking_likelihood <- componentwise_dwmb(M_vec, alpha_grid, ps, thetas, rankings_vec, pi0s, ratings_vec, log=F)#[1]
  rating_likelihood <- componentwise_dwmb(M_vec, alpha_grid, ps, thetas, rankings_vec, pi0s, ratings_vec, log=F)#[2]
  
  # Create dataframe of all results
  alpha_df <- data.frame(alpha=alpha_grid)
  theta_df <- data.frame(unlist(thetas))
  names(theta_df) <- c("theta")
  num_nodes_df <- data.frame(num_nodes=unlist(nums_nodes))
  names(num_nodes_df) <- c("num_nodes")
  
  astar_df <- cbind(alpha_df, pi0_df, ps_df, theta_df, num_nodes_df)
  
  # This is the most human-interpretable dataframe containing the data
  interpretable_df <- astar_df
  
  # Pivot the data longer to plot effectively
  astar_df <- astar_df %>% pivot_longer(cols = starts_with("P"), names_to="p_name", values_to="p_value")
  astar_df <- astar_df %>% pivot_longer(cols = starts_with("R"), names_to="r_name", values_to="r_value")
  
  ast_plot <- ggplot(data=astar_df, aes(x=alpha, y=p_value, color=p_name))# + 
  if(smooth){
    ast_plot <- ast_plot + geom_smooth(method=stats::loess, level=0.95) + scale_y_continuous(limits=c(0,1))
  } else {
    ast_plot <- ast_plot + geom_line() + scale_y_continuous(limits=c(0,1))
  }

  if(parameter=="theta") {
    ast_plot <- ggplot(data=astar_df, aes(x=alpha, y=log(theta))) + 
      geom_line()
  }
  
  plot(ast_plot)
  if(return_plot) {
    return(ast_plot)
  }
  return(interpretable_df)
}

plot_rankings <- function(rankings,palette="Greens") {
  I <- nrow(rankings)
  J <- ncol(rankings)
  
  rankings_table <- as.data.frame(rankings)
  rownames(rankings_table) <- paste0("Judge ",1:I)
  names(rankings_table) <- paste0("Rank ",1:J)
  pander(rankings_table)
  
  rankings_long <- melt(rankings)
  names(rankings_long) <- c("Judge","Rank","Object")
  rankings_plot <- ggplot(rankings_long,aes(x=Object,fill=factor(Rank)))+theme_bw(base_size=10)+
    geom_bar()+#scale_colour_steps(low = "#318354", high = "#FBF5E0") + 
    scale_fill_brewer(palette=palette)+
    labs(fill="Rank",y="Count")+ggtitle("Ranks by Object")+
    theme(panel.grid = element_blank())
    
  return(rankings_plot)
  # if(return_plot) {
  #   return(rankings_plot)
  # } else {
  #   rankings_plot
  # }
}

plot_ratings <- function(ratings) {
  I <- nrow(ratings)
  J <- ncol(ratings)
  
  ratings_table <- as.data.frame(ratings)
  rownames(ratings_table) <- paste0("Judge ",1:I)
  names(ratings_table) <- c("Object: 1",2:J)
  set.alignment("right")
  pander(ratings_table)
  
  ratings_long <- melt(ratings)
  names(ratings_long) <- c("Judge","Object","Rating")
  ratings_plot <- ggplot(ratings_long,aes(x=factor(Object),y=Rating))+theme_bw(base_size=10)+
    geom_boxplot(color="gray")+geom_jitter(height=0,width=0.4,alpha=0.75)+
    labs(x="Object",y="Rating")+ggtitle("Ratings by Object")+
    theme(panel.grid = element_blank())
  
  return(ratings_plot)
}

show_mb_estimate <- function(rankings,ratings,M,digits=3) {
  MLE_mb <- fit_mb(rankings=rankings,ratings=ratings,M=M,method="ASTAR")
  if(is.matrix(MLE_mb$pi0)){
    which_keep <- sample(1:nrow(MLE_mb$pi0),1)
    MLE_mb$pi0 <- MLE_mb$pi0[which_keep,]
    MLE_mb$p <- MLE_mb$p[which_keep,]
    MLE_mb$theta <- MLE_mb$theta[which_keep,]
  }
  # Note: I set the convention that a>b means a has a numerically lower rating
  MLE_table <- data.frame(Parameter=c("Consensus Ranking, pi_0",
                                "Object Quality Parameter, p",
                                "Consensus Scale Parameter, theta"),
                    MLE=c(paste0(MLE_mb$pi0,collapse=">"),
                          paste0("(",paste0(round(MLE_mb$p,digits),collapse=","),")"),
                          round(MLE_mb$theta,digits)))
  pander(MLE_table)
  return(MLE_mb)
}

show_wmb_estimate <- function(rankings,ratings,M,alpha,digits=3) {
  if(alpha<0 || alpha>1) {stop("alpha must be within [0, 1]")}
  MLE_wmb <- ASTAR_wmb(rankings=rankings,ratings=ratings,M=M,alpha=alpha)
  if(is.matrix(MLE_wmb$pi0)){
    which_keep <- sample(1:nrow(MLE_wmb$pi0),1)
    MLE_wmb$pi0 <- MLE_wmb$pi0[which_keep,]
    MLE_wmb$p <- MLE_wmb$p[which_keep,]
    MLE_wmb$theta <- MLE_wmb$theta[which_keep,]
  }
  # Note: I set the convention that a>b means a has a numerically lower rating
  MLE_table <- data.frame(Parameter=c("Consensus Ranking, pi_0",
                                      "Object Quality Parameter, p",
                                      "Consensus Scale Parameter, theta"),
                          MLE=c(paste0(MLE_wmb$pi0,collapse=">"), 
                                paste0("(",paste0(round(MLE_wmb$p,digits),collapse=","),")"),
                                round(MLE_wmb$theta,digits)))
  pander(MLE_table)
  return(MLE_wmb)
}

ci_mb <- function(rankings,ratings,M,interval=0.90,nsamples=50,all=FALSE,method="ASTAR"){
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  bs_pi0 <- matrix(NA,nrow=nsamples,ncol=J)
  bs_parameters <- matrix(NA,nrow=nsamples,ncol=J+1)
  
  for(sample in 1:nsamples){
    bs_sample <- sample(1:I,I,replace=T)
    bs_res <- suppressMessages(fit_mb(rankings[bs_sample,],ratings[bs_sample,],M,method))
    if(is.matrix(bs_res$pi0)){
      which_keep <- sample(1:nrow(bs_res$pi0),1)
      bs_res$pi0 <- bs_res$pi0[which_keep,]
      bs_res$p <- bs_res$p[which_keep,]
      bs_res$theta <- bs_res$theta[which_keep,]
    }
    bs_pi0[sample,] <- bs_res$pi0
    bs_parameters[sample,] <- c(bs_res$p,bs_res$theta)
  }
  ci <- as.matrix(apply(bs_parameters,2,function(parameter){quantile(parameter,probs=c((1-interval)/2,1-(1-interval)/2))}))
  colnames(ci) <- colnames(bs_parameters) <- c(paste0("p",1:J),"theta")
  ci_ranks <- matrix(unlist(lapply(1:J,function(j){
    quantile(apply(bs_pi0,1,function(pi0){which(pi0==j)}),probs=c((1-interval)/2,1-(1-interval)/2))
  })),nrow=2)
  rownames(ci_ranks) <- rownames(ci)
  colnames(ci_ranks) <- paste0("Object",1:J)
  
  if(all){return(list(ci=ci,ci_ranks=ci_ranks,
                      bootstrap_pi0=bs_pi0,
                      bootstrap_ptheta=bs_parameters))
  }else{return(list(ci=ci,ci_ranks=ci_ranks))}
}

ci_wmb <- function(rankings,ratings,M,alpha,interval=0.90,nsamples=50,all=FALSE){
  I <- nrow(rankings)
  J <- ncol(rankings)
  if(any(dim(ratings)!=c(I,J))){stop("rankings and ratings must be of the same dimension")}
  bs_pi0 <- matrix(NA,nrow=nsamples,ncol=J)
  bs_parameters <- matrix(NA,nrow=nsamples,ncol=J+1)
  
  for(sample in 1:nsamples){
    # Repeatedly try samples until one does not return null
    repeat{
      bs_sample <- sample(1:I,I,replace=T)
      bs_res <- NULL
      
      alpha_used <- alpha
      n_tries <- 0
      adjusting_failed <- FALSE
      # If an error occurs, try a slightly different alpha value
      tryCatch({bs_res <- suppressMessages(ASTAR_wmb(rankings[bs_sample,],ratings[bs_sample,],M,alpha_used))},
        error=function(e){
          warning("Caught error in ASTAR, trying a similar value.")
          epsilon <- 1e-8
          if(alpha_used<(epsilon*10)) {
            alpha_used <- alpha_used+epsilon
          } else {
            alpha_used <- alpha_used-epsilon
          }
          n_tries <- n_tries + 1
          # Put a threshold on the maximum number of attempts
          if(n_tries > 20) {
            adjusting_failed <- TRUE
            break
          }
        }
      )
      if(adjusting_failed) {
        stop("Small adjustments failed to fix ASTAR optimization inaccuracies")
      }
      if(!is.null(bs_res)) {
        break
      }
    }
    
    if(is.matrix(bs_res$pi0)){
      which_keep <- sample(1:nrow(bs_res$pi0),1)
      bs_res$pi0 <- bs_res$pi0[which_keep,]
      bs_res$p <- bs_res$p[which_keep,]
      bs_res$theta <- bs_res$theta[which_keep,]
    }
    bs_pi0[sample,] <- bs_res$pi0
    bs_parameters[sample,] <- c(bs_res$p,bs_res$theta)
  }
  ci <- as.matrix(apply(bs_parameters,2,function(parameter){quantile(parameter,probs=c((1-interval)/2,1-(1-interval)/2))}))
  colnames(ci) <- colnames(bs_parameters) <- c(paste0("p",1:J),"theta")
  ci_ranks <- matrix(unlist(lapply(1:J,function(j){
    quantile(apply(bs_pi0,1,function(pi0){which(pi0==j)}),probs=c((1-interval)/2,1-(1-interval)/2))
  })),nrow=2)
  rownames(ci_ranks) <- rownames(ci)
  colnames(ci_ranks) <- paste0("Object",1:J)
  
  if(all){return(list(ci=ci,ci_ranks=ci_ranks,
                      bootstrap_pi0=bs_pi0,
                      bootstrap_ptheta=bs_parameters))
  }else{return(list(ci=ci,ci_ranks=ci_ranks))}
}

plot_p_wmb <- function(MLE_wmb, CI_wmb, show_data=F) {
  J <- length(MLE_wmb$p)
  p_plot <- as.data.frame(cbind(1:J,MLE_wmb$p,t(CI_wmb$ci[,1:J])))
  names(p_plot) <- c("Object","PointEstimate","Lower","Upper")
  if(show_data){
    print("Data for Plot of Estimated Quality by Object:")
    print(p_plot)
  }
  ggplot(p_plot,aes(x=Object,y=PointEstimate,ymin=Lower,ymax=Upper))+
    geom_point()+geom_errorbar()+ylim(c(-0.01,1.01))+theme_bw(base_size=10)+
    labs(x="Object",y="Estimated Quality (95% CI)")+
    ggtitle("Estimated Quality by Object")+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
}

# Current unused utility function to run a full analysis on a dataset.
# Unused because it is not resilient to errors in a single step. 
analyze_rankrate <- function(rankings,ratings,M,alpha,ngrid=100,nsamples=200, name="analyzed") {
  # Visualize the ranking and rating data 
  rankings_plot <- plot_rankings(rankings)
  ratings_plot <- plot_ratings(ratings)
  
  # Calculate original MB and new WMB maximum likelihood estimators
  MLE_mb <- show_mb_estimate(rankings, ratings, M)
  MLE_wmb <- show_wmb_estimate(rankings, ratings, M, alpha=alpha)
  
  # Create plots showing how estimation results change with alpha
  p_plot <- ASTAR_plot(rankings,ratings,M,n=ngrid,parameter="p",return_plot=T)
  theta_plot <- ASTAR_plot(rankings,ratings,M,n=ngrid,parameter="theta",return_plot=T)
  
  # Create plots of confidence intervals on p
  CI_wmb <- ci_wmb(rankings=rankings,ratings=ratings,M=M,alpha=0.6,
                        interval=0.95,nsamples=nsamples)
  ci_plot <- plot_p_wmb(MLE_wmb, CI_wmb)
  
  print(paste0("Saving to location: conor/thesis/semester2/outputs/", name))
  # Save all plots to files
  png(paste0("conor/thesis/semester2/outputs/", name,"_rankings_plot.png"), width=6, height=5, units="in", res=1200)
  rankings_plot
  dev.off()
  png(paste0("conor/thesis/semester2/outputs/", name,"_ratings_plot.png"), width=6, height=5, units="in", res=1200)
  ratings_plot
  dev.off()
  png(paste0("conor/thesis/semester2/outputs/", name,"_p_plot.png"), width=6, height=5, units="in", res=1200)
  p_plot
  dev.off()
  png(paste0("conor/thesis/semester2/outputs/", name,"_theta_plot.png"), width=6, height=5, units="in", res=1200)
  theta_plot
  dev.off()
  png(paste0("conor/thesis/semester2/outputs/", name,"_ci_plot.png"), width=6, height=5, units="in", res=1200)
  ci_plot
  dev.off()
  
  print(paste0("Finished analyzing data labeled ", name))
}

