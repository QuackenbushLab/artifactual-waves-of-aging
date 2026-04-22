## Copied and Adapted from https://github.com/jaspershen-lab/ipop_aging/tree/main

suppressPackageStartupMessages({
  library(plyr)
  library(purrr)
  library(dplyr)
  # library(furrr)  # optional
})

# ... optimize_loess_span ...
# ... run_loess ...

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
    
    # span_rmse
    # plot <-
    #   data.frame(x, y) %>%
    #   ggplot(aes(x, y)) +
    #   geom_point(size = 5) +
    #   # geom_line() +
    #   base_theme
    
    # span_rmse =
    #   span_rmse %>%
    #   dplyr::filter(!is.na(rmse))
    # idx = which.min(span_rmse$rmse)
    # for(i in 1:nrow(span_rmse)){
    # plot =
    #   plot +
    #   geom_smooth(
    #     se = FALSE,
    #     span = span_rmse$span[idx],
    #     color = ggsci::pal_lancet()(n = 9)[idx]
    #   )
    # # }
    
    # plot =
    #   plot +
    #   ggplot2::ggtitle(label = paste("Span: ", span_rmse$span[idx])) +
    #   theme(title = element_text(colour = ggsci::pal_lancet()(n = 9)[idx]))
    
    list(span_rmse)
  }


run_loess <- function(variable_info, sample_info, expression_data, sup_negative = TRUE) {
  
  new_expression_data <-
    vector(mode = "list", length = nrow(variable_info))
  
  for (i in seq_along(variable_info$variable_id)) {
    temp_variable_id <- variable_info$variable_id[i]
    cat(i, " ")
    temp_data <-
      data.frame(value = as.numeric(expression_data[temp_variable_id,]),
                 sample_info)
    
    optimize_span <-
      optimize_loess_span(
        x = temp_data$adjusted_age,
        y = temp_data$value,
        span_range = c(0.3, 0.4, 0.5, 0.6)
      )
    
    span <-
      optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
    
    value <- temp_data$value
    adjusted_age <- temp_data$adjusted_age
    
    ls_reg <-
      loess(value ~ adjusted_age,
            span = span)
    
     prediction_value =
        predict(ls_reg,
                newdata = data.frame(adjusted_age = seq(26, 75, by = 0.5)))
    new_expression_data[[i]] <- as.numeric(prediction_value)
  }
  
  new_expression_data <-
    new_expression_data %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
    colnames(new_expression_data) <-
    as.character(seq(26, 75, by = 0.5))
  
  rownames(new_expression_data) <- rownames(expression_data)
  
  negative_idx <-
    new_expression_data %>%
    apply(1, function(x) {
      sum(x < 0)
    }) %>%
    `>`(0) %>%
    which()
  
  if (sup_negative){
    new_expression_data[which(new_expression_data < 0, arr.ind = TRUE)] <-
    0
  }
  
  
  return(new_expression_data)
}
