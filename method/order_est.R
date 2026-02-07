source("MDCP.R")
source("PMDCP.R")
source("SMF.R")
source("PMF.R")

order_est = function(data, method, alpha = 0.05){
  p_candidate = 1:5
  n = length(data)
  train_n = floor(n / 3) * 2
  valid_n = train_n / 2
  n = train_n + valid_n
  x <- data[1: n]
  
  cvr <- c()
  len <- c()
  if (method == "conformal"){
    for (p in p_candidate){
      cv <- logical(0)
      length <- c()
      mdcp_upper<- c()
      mdcp_lower <- c()
      
      for (t in train_n:(n - 1)) {
        x_train <- x[(t-train_n+1):t]
        interval <- MDCP(x_train, p = p, alpha)
        
        lower <- interval$lower
        upper <- interval$upper
     
        mdcp_upper <- c(mdcp_upper, upper)
        mdcp_lower <- c(mdcp_lower, lower)

        y_next <- x[t + 1]
        cv <- c(cv, (y_next >= lower && y_next <= upper))
        length <- c(length, upper - lower)
      }
      
      cvr[p] <- mean(cv)
      len[p] <- mean(length)
    }
  }  
    else if (method == "conformal_predict"){
      for (p in p_candidate){
        cv <- logical(0)
        length <- c()
        mdcp_upper<- c()
        mdcp_lower <- c()
        
        for (t in train_n:(n - 1)) {
          x_train <- x[(t-train_n+1):t]
          interval <- PDCP(x_train, p = p, alpha)
          
          lower <- interval$lower
          upper <- interval$upper
          
          mdcp_upper <- c(mdcp_upper, upper)
          mdcp_lower <- c(mdcp_lower, lower)
          
          y_next <- x[t + 1]
          cv <- c(cv, (y_next >= lower && y_next <= upper))
          length <- c(length, upper - lower)
        }
        
        cvr[p] <- mean(cv)
        len[p] <- mean(length)
      }
    }
      else if (method == "MF"){
        for (p in p_candidate){
          cv <- logical(0)
          length <- c()
          mdcp_upper<- c()
          mdcp_lower <- c()
          
          for (t in train_n:(n - 1)) {
            x_train <- x[(t-train_n+1):t]
            interval <- SMF(x = x_train, p = p, alpha = alpha)
            
            lower <- interval$lower
            upper <- interval$upper
            
            mdcp_upper <- c(mdcp_upper, upper)
            mdcp_lower <- c(mdcp_lower, lower)
            
            y_next <- x[t + 1]
            cv <- c(cv, (y_next >= lower && y_next <= upper))
            length <- c(length, upper - lower)
          }
          
          cvr[p] <- mean(cv)
          len[p] <- mean(length)
        }
      }
        else if (method == "PMF"){
          for (p in p_candidate){
            cv <- logical(0)
            length <- c()
            mdcp_upper<- c()
            mdcp_lower <- c()
            
            for (t in train_n:(n - 1)) {
              x_train <- x[(t-train_n+1):t]
              interval <- PMF(x_train, p = p, alpha)
              
              lower <- interval$lower
              upper <- interval$upper
              
              mdcp_upper <- c(mdcp_upper, upper)
              mdcp_lower <- c(mdcp_lower, lower)
              
              y_next <- x[t + 1]
              cv <- c(cv, (y_next >= lower && y_next <= upper))
              length <- c(length, upper - lower)
            }
            
            cvr[p] <- mean(cv)
            len[p] <- mean(length)
          }
  }
  return(which.min(abs(cvr - (1-alpha))*len))
}
