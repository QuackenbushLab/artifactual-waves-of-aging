library(tidyverse)
library(truncnorm)

# function from Shen et al. to optimize LOESS span
optimize_loess_span <-
  function(x, y, span_range = seq(0.2, 0.6, 0.1)) {
    span_rmse <-
      purrr::map(span_range, function(span) {
        # cat(span, " ")
        temp_data =
          data.frame(x, y)
        
        prediction <-
          purrr::map(
            2:(nrow(temp_data) - 1),
            .f = function(idx) {
              temp_result =
                loess(formula = y ~ x,
                      data = temp_data[-idx,],
                      span = span)
              prediction =
                try(predict(object = temp_result,
                            newdata = temp_data[idx, -2, drop = FALSE]))
              
              if (class(prediction) == "try-error") {
                data.frame(real = temp_data$y[idx],
                           prediction = NA)
              } else{
                data.frame(real = temp_data$y[idx],
                           prediction = as.numeric(prediction))
              }
            }
          ) %>%
          dplyr::bind_rows() %>%
          dplyr::filter(!is.na(prediction))
        
        if (all(is.na(prediction$prediction))) {
          temp_rmse = NA
        } else{
          temp_rmse = sqrt(sum((
            prediction$real - prediction$prediction
          ) ^ 2) / nrow(prediction))
        }
        
        data.frame(span = span, rmse = temp_rmse)
      }) %>%
      dplyr::bind_rows()
    
    span_rmse =
      span_rmse %>%
      dplyr::filter(!is.na(rmse))
    idx = which.min(span_rmse$rmse)
    
    list(span_rmse)
}

do_de_swan <- # response should be df of rows = genes, cols = vals at age, sorted same as age vector
  function(age,response, window_center, buckets_size) {
    p_value <-
      window_center %>%
      purrr::map(function(i) {
        cat(i, " ")
        left_range <-
          c(i - buckets_size / 2, i)
        right_range <-
          c(i, i + buckets_size / 2)
        
        left_idx <-
          which(age >= left_range[1] &
                  age < left_range[2])
        
        control_sample_id <- left_idx
        
        right_idx <-
          which(age >= right_range[1] &
                  age < right_range[2])
        
        case_sample_id <- right_idx
        expression_data = data.frame("response"=response)
        p_value <-
          seq_len(nrow(expression_data)) %>%
          purrr::map(function(i) {
            # change to wilcox
            cntrl_data = data.frame("expr"=as.numeric(expression_data[i, control_sample_id]),
                                    "type"="cntrl")
            case_data = data.frame("expr"=as.numeric(expression_data[i, case_sample_id]),
                                   "type"="case")
            both_data = rbind.data.frame(cntrl_data,case_data)
            wilcox.test(as.numeric(expression_data[i, control_sample_id]),
                        as.numeric(expression_data[i, case_sample_id]),
                        exact = T)$p.value
            
            
          }) %>%
          unlist()
        
        p_value <-
          data.frame(variable_id = row.names(expression_data),
                     p_value) %>%
          dplyr::mutate(p_value_adjust = p.adjust(p_value, "BH"),
                        center = i)
        p_value
        
      }) %>%
      dplyr::bind_rows()
    
    p_value
}


do_de_swan_first_window <- # response should be df of rows = genes, cols = vals at age, sorted same as age vector
  function(age,response, window_center = c(40)) {
    p_value <-
      window_center %>%
      purrr::map(function(i) {
        cat(i, " ")
        left_range <-
          c(i - 15, i)
        right_range <-
          c(i, i + 10)
        
        left_idx <-
          which(age >= left_range[1] &
                  age < left_range[2])
        
        control_sample_id <- left_idx
        
        right_idx <-
          which(age >= right_range[1] &
                  age < right_range[2])
        
        case_sample_id <- right_idx
        expression_data = data.frame("response"=response)
        p_value <-
          seq_len(nrow(expression_data)) %>%
          purrr::map(function(i) {
            # change to wilcox
            cntrl_data = data.frame("expr"=as.numeric(expression_data[i, control_sample_id]),
                                    "type"="cntrl")
            case_data = data.frame("expr"=as.numeric(expression_data[i, case_sample_id]),
                                   "type"="case")
            both_data = rbind.data.frame(cntrl_data,case_data)
            #pvalue(wilcox_test(both_data$expr ~ factor(both_data$type), distribution="exact"))
            
            wilcox.test(as.numeric(expression_data[i, control_sample_id]),
                        as.numeric(expression_data[i, case_sample_id]), 
                        exact = T)$p.value
          }) %>%
          unlist()
        
        p_value <-
          data.frame(variable_id = row.names(expression_data),
                     p_value) %>%
          dplyr::mutate(p_value_adjust = p.adjust(p_value, "BH"),
                        center = i)
        p_value
        
      }) %>%
      dplyr::bind_rows()
    
    p_value
  }

test_age_dist = function(my_ages,npred, bucket, 
                         optimize_loess = F, 
                         do_ipop = F,
                         do_log = F,
                         distribution = "standard_normal")
{
  windows = seq(min(my_ages) + bucket/2,
                max(my_ages) - bucket/2,
                by=1)
  
  if(distribution == "uniform")
    my_responses = matrix(runif(n=length(my_ages)*npred),nrow=npred)
  if(distribution == "standard_normal")
    my_responses = matrix(rnorm(n=length(my_ages)*npred),nrow=npred)
  
  if(!do_ipop)
    res = do_de_swan(age = my_ages, 
                        response = data.frame(my_responses), 
                        window_center = windows,
                        buckets_size = bucket)
  
  if(do_ipop)
  {
    ipop_windows = seq(41,65,by=1)
    first_window = do_de_swan_first_window(age = my_ages,
                                  response = data.frame(my_responses), 
                                  window_center = c(40))
    other_windows = do_de_swan(age = my_ages,
                                    response =data.frame(my_responses),
                                    window_center = ipop_windows,
                                    buckets_size = bucket)
    res = rbind.data.frame(first_window, other_windows)
  }
  
  # now do loess on my_responses and run de_swan on that
  loess_lattice = data.frame(x=seq(min(my_ages),
                                   max(my_ages),by=0.5))
  loessed_data = list()
  opt_spans = list()
  for(i in 1:nrow(my_responses))
  {
    if(i %% 50 == 0 & do_log) 
      log_print(paste0("processing molecule:",i))
    # get optimal loess span
    # by default, use 0.6, which is what is preferred for
    # most molecules when tested for optimal
    loess_span = 0.6
    if(optimize_loess)
    {
      optimal_span = optimize_loess_span(x=my_ages, y=my_responses[i,], span_range = seq(0.2, 0.6, 0.1))
      opt_span_val = optimal_span[[1]]$span[which.min(optimal_span[[1]]$rmse)]
      opt_spans[[i]] = opt_span_val
      loess_span = opt_span_val
    }
    
    # run loess model 
    loess_model = loess(my_responses[i,] ~ my_ages, span = loess_span)
    loessed_data[[i]] = predict(loess_model, newdata = loess_lattice$x)
  }
  
  loessed_data = Reduce(rbind.data.frame, loessed_data)
  if(!do_ipop)
    res_loess = do_de_swan(age = loess_lattice$x, response = loessed_data, window_center = windows,
                              buckets_size = bucket)
  if(do_ipop)
  {
    ipop_windows = seq(41,65,by=1)
    first_window = do_de_swan_first_window(age = loess_lattice$x,
                                           response = loessed_data,
                                           window_center = c(40))
    other_windows = do_de_swan(age = loess_lattice$x,
                                   response = loessed_data,
                                   window_center = ipop_windows,
                                   buckets_size = bucket)
    res_loess = rbind.data.frame(first_window, other_windows)
  }
  
  return(list(res,res_loess,opt_spans))
}


get_profile = function(x, id)
{
  x %>%
    group_by(center) %>%
    summarize(nhits = sum(p_value_adjust < 0.05)) 
}

get_mean_profile = function(sim_list, loess=T)
{
  idx = ifelse(loess==T,2,1)
  all_profiles = lapply(sim_list,function(x){get_profile(x[[idx]])})
  for(i in 1:length(all_profiles))
    all_profiles[[i]]$id = paste0("Trial:",i)
  
  all_profiles_df = Reduce(rbind.data.frame, all_profiles)
  mean_profile = all_profiles_df %>%
    group_by(center) %>%
    summarize(mean_nhits = mean(nhits))
  return(mean_profile)
}

plot_from_sim_list = function(sim_list, npred, loess=T){
  idx = ifelse(loess==T,2,1)
  all_profiles = lapply(sim_list,function(x){get_profile(x[[idx]])})
  for(i in 1:length(all_profiles))
    all_profiles[[i]]$id = paste0("Trial:",i)
  
  all_profiles_df = Reduce(rbind.data.frame, all_profiles)
  
  p = ggplot() +
    geom_line(data=all_profiles_df, aes(x=center,y=nhits,color=id),alpha=0.2) +
    theme(legend.position = "none")
  
  mean_profile = all_profiles_df %>%
    group_by(center) %>%
    summarize(mean_nhits = mean(nhits))
  
  if(loess)
    { p=p + 
    geom_line(data = mean_profile,aes(x=center,y=mean_nhits)) +
    geom_point(data = mean_profile,aes(x=center,y=mean_nhits)) +
    xlab("window center") +
    ylab("# molecules with q < 0.05") + #+ 
    ylim(0.5*npred,npred); return(p)}
  if(!loess)
  { p=p + 
      geom_line(data = mean_profile,aes(x=center,y=mean_nhits)) +
      geom_point(data = mean_profile,aes(x=center,y=mean_nhits)) +
      xlab("window center") +
      ylab("# molecules with q < 0.05") + #+ 
      ylim(0,4); return(p)}
  
}

loess_matrix = function(my_ages, my_responses, optimize_loess=F){
  loess_lattice = data.frame(x=seq(min(my_ages),
                                     max(my_ages),by=0.5))
  loessed_data = list()
  opt_spans = list()
  for(i in 1:nrow(my_responses))
  {
    # get optimal loess span
    # by default, use 0.6, which is what is preferred for
    # most molecules when tested for optimal
    loess_span = 0.6
    if(optimize_loess)
    {
      optimal_span = optimize_loess_span(x=my_ages, y=my_responses[i,], span_range = seq(0.2, 0.6, 0.1))
      opt_span_val = optimal_span[[1]]$span[which.min(optimal_span[[1]]$rmse)]
      opt_spans[[i]] = opt_span_val
      loess_span = opt_span_val
    }
    
    # run loess model 
    loess_model = loess(my_responses[i,] ~ my_ages, span = loess_span)
    loessed_data[[i]] = predict(loess_model, newdata = loess_lattice$x)
  }
  
  loessed_data = Reduce(rbind.data.frame, loessed_data)
  return(list("lattice_ages"=loess_lattice$x,"response"=loessed_data))
}

