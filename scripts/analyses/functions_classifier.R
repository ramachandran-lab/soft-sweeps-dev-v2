library(mclust)  # gaussian mixture
library(plotROC)  # roc curve analysis

# # # 

my_max=5
my_size=100

missingno <- function(train_x_uni, N=100) {
  for (stat in names(train_x_uni)) {
    train_x_uni[[stat]][sample(c(1:length(train_x_uni[[stat]])), N, replace=T)] <- NA
    train_x_uni[[stat]][sample(c(1:length(train_x_uni[[stat]])), N, replace=T)] <- Inf
    train_x_uni[[stat]][sample(c(1:length(train_x_uni[[stat]])), N, replace=T)] <- -Inf
  }
  return(train_x_uni)
}

fit_mix <- function(train_x, remove_outlier=NA) {
  train_x_help <- as_tibble(train_x)
  len_na <- sum(is.na(train_x_help))
  train_x_help <- na.omit(train_x_help)
  len_total <- nrow(train_x_help)
  train_x_help_fin <- train_x_help[is.finite(rowSums(train_x_help)),]
  train_x_help_inf <- train_x_help[!is.finite(rowSums(train_x_help)),]
  len_finite <- nrow(train_x_help_fin)
  fit_x <- NULL
  # add argument remove_outlier, if not NA, retain only middle 
  # that proportion of points before fitting univariate classifier
  if (!is.na(remove_outlier)) {
    train_x_help_fin <- train_x_help_fin[which(
      train_x_help_fin>=quantile(as.matrix(train_x_help_fin), (1-remove_outlier)/2) & 
      train_x_help_fin<=quantile(as.matrix(train_x_help_fin), (1+remove_outlier)/2)),]
  }
  # fitting univariate classifier
  fit_x$dens <- densityMclust(as.matrix(train_x_help_fin), G=1:my_max,
    initialization=list(sample(1:nrow(train_x_help_fin), size=my_size, replace=T)))
  if (ncol(train_x_help) > 1) {
    stop("To do: multivariate implementation of boundary densities.")
  } else {
    fit_x$spike <- NULL
    fit_x$spike$fin <- len_finite
    fit_x$spike$neg <- sum(as.matrix(train_x_help_inf)<0)
    fit_x$spike$pos <- sum(as.matrix(train_x_help_inf)>0)
    fit_x$spike$na <- len_na
  }
  return(fit_x)
}

logd_mix <- function(fit_x, pred_x) {
  pred_x_help <- as_tibble(pred_x)
  n_col <- ncol(pred_x_help)
  w_fin <- which(is.finite(rowSums(pred_x_help)))
  w_neg <- which(is.infinite(rowSums(pred_x_help)) & pred_x_help<0)
  w_pos <- which(is.infinite(rowSums(pred_x_help)) & pred_x_help>0)
  sum_total <- fit_x$spike$fin+fit_x$spike$neg+fit_x$spike$pos
  pred_temp <- log10(predict(fit_x$dens, pred_x_help[w_fin,])*fit_x$spike$fin/sum_total)
  pred_x_help$logdens <- NA
  pred_x_help$logdens[w_fin] <- pred_temp
  if (n_col > 1) {
    stop("To do: multivariate implementation of boundary densities.")
  } else {
    pred_x_help$logdens[w_neg] <- log10(fit_x$spike$neg/sum_total)
    pred_x_help$logdens[w_pos] <- log10(fit_x$spike$pos/sum_total)
  }
  logdens_x <- pred_x_help["logdens"]
  return(logdens_x)
}

fit_mix_uni <- function(train_x_uni, train_y_uni, remove_outlier=NA) {
  fit_x_uni <- NULL
  fit_x_uni$my_classes <- names(table(train_y_uni))
  fit_x_uni$my_stat_temp <- names(train_x_uni)
  for (stat in fit_x_uni$my_stat_temp) {
    fit_x_uni[[stat]] <- NULL
    for (cls in fit_x_uni$my_classes) {
      train_x <- train_x_uni[which(train_y_uni==cls), stat]
      fit_x_uni[[stat]][[cls]] <- fit_mix(train_x, remove_outlier)
    }
  }
  return(fit_x_uni)
}

logd_mix_uni <- function(fit_x_uni, pred_x_uni) {
  logd_y_uni <- NULL
  logd_y_uni$my_classes <- fit_x_uni$my_classes
  logd_y_uni$my_stat_temp <- fit_x_uni$my_stat_temp
  for (stat in logd_y_uni$my_stat_temp) {
    for (cls in logd_y_uni$my_classes) {
      if (is.null(logd_y_uni[[stat]])) {
        logd_y_uni[[stat]] <- logd_mix(fit_x_uni[[stat]][[cls]], pred_x_uni[stat])
      } else {
        logd_y_uni[[stat]] <- cbind(logd_y_uni[[stat]], logd_mix(fit_x_uni[[stat]][[cls]], pred_x_uni[stat]))
      }
    }
    names(logd_y_uni[[stat]]) <- logd_y_uni$my_classes
    logd_y_uni[[stat]] <- as_tibble(logd_y_uni[[stat]])
  }
  return(logd_y_uni)
}

pred_helper <- function(my_pred, prob_postodds) {
  my_pred_help <- my_pred
  my_pred_help <- my_pred_help %>% 
    mutate(neutral=ifelse(my_tol(neutral, 1), 1-.Machine$double.eps, neutral), 
           neutral=ifelse(my_tol(neutral, 0), 0+.Machine$double.eps, neutral)) %>%
    mutate(hardsweep=ifelse(my_tol(hardsweep, 1), 1-.Machine$double.eps, hardsweep), 
           hardsweep=ifelse(my_tol(hardsweep, 0), 0+.Machine$double.eps, hardsweep)) %>%
    mutate(softsweep=ifelse(my_tol(softsweep, 1), 1-.Machine$double.eps, softsweep), 
           softsweep=ifelse(my_tol(softsweep, 0), 0+.Machine$double.eps, softsweep))
  my_pred_help <- my_pred_help %>%
    mutate(pred_class = ifelse((neutral >= hardsweep) & (neutral >= softsweep), "neutral", NA)) %>%
    mutate(pred_class = ifelse((hardsweep >= neutral) & (hardsweep >= softsweep), "hardsweep", pred_class)) %>%
    mutate(pred_class = ifelse((softsweep >= neutral) & (softsweep >= hardsweep), "softsweep", pred_class))
  my_pred_help <- my_pred_help %>%
    mutate(pred_class_alt = pred_class) %>% 
    mutate(pred_class_alt=ifelse(pred_class_alt=="softsweep", ifelse(
      softsweep>=prob_postodds*hardsweep, "softsweep", "softsweepweak"), pred_class_alt)) %>%
    mutate(pred_class_alt=ifelse(pred_class_alt=="hardsweep", ifelse(
      hardsweep>=prob_postodds*softsweep, "hardsweep", "hardsweepweak"), pred_class_alt))
  my_pred_help <- my_pred_help %>% 
    mutate(pred_sweep_neutral = ifelse((neutral >= (softsweep+hardsweep)), "neutral", NA)) %>% 
    mutate(pred_sweep_neutral = ifelse(((softsweep+hardsweep) >= neutral), "sweep", pred_sweep_neutral))
  my_pred_help <- my_pred_help %>% 
    mutate(pred_sweep_neutral_alt = pred_sweep_neutral) %>% 
    mutate(pred_sweep_neutral_alt=ifelse(pred_sweep_neutral_alt=="sweep", ifelse(
      (softsweep+hardsweep)>=prob_postodds*neutral, "sweep", "sweepweak"), pred_sweep_neutral_alt))
  my_pred_help <- my_pred_help %>% 
    mutate(pred_soft_hard = ifelse((hardsweep >= softsweep), "hardsweep", NA)) %>% 
    mutate(pred_soft_hard = ifelse((softsweep >= hardsweep), "softsweep", pred_soft_hard))
  my_pred_help <- my_pred_help %>% 
    mutate(pred_soft_hard_alt = pred_soft_hard) %>% 
    mutate(pred_soft_hard_alt=ifelse(pred_soft_hard_alt=="softsweep", ifelse(
      softsweep>=prob_postodds*hardsweep, "softsweep", "softsweepweak"), pred_soft_hard_alt)) %>%
    mutate(pred_soft_hard_alt=ifelse(pred_soft_hard_alt=="hardsweep", ifelse(
      hardsweep>=prob_postodds*softsweep, "hardsweep", "hardsweepweak"), pred_soft_hard_alt))
  return(my_pred_help)
}

logd_mix_nb <- function(logd_y_uni, stat_sweep_final) {
  logd_y_nb <- NULL
  logd_y_nb$my_classes <- logd_y_uni$my_classes
  logd_y_nb$stat_final <- stat_sweep_final
  
  for (cls in logd_y_nb$my_classes) {
    for (stat in stat_sweep_final) {
      if (is.null(logd_y_nb$temp[[cls]])) {
        logd_y_nb$temp[[cls]] <- logd_y_uni[[stat]][cls]
      } else {
        logd_y_nb$temp[[cls]] <-
          cbind(logd_y_nb$temp[[cls]], logd_y_uni[[stat]][cls])
      }
    }
    names(logd_y_nb$temp[[cls]]) <- stat_sweep_final
    logd_y_nb$temp[[cls]] <- as_tibble(logd_y_nb$temp[[cls]])
    # ignore NAs unless all NA
    logd_y_nb$temp[[cls]] <- logd_y_nb$temp[[cls]] %>%
      mutate(NB=rowSums(logd_y_nb$temp[[cls]], na.rm=T)) %>%
      mutate(NB=ifelse(rowSums(!is.na(.))==1, NA, NB))
  }
  
  for (cls in logd_y_nb$my_classes) {
    if (is.null(logd_y_nb$main)) {
      logd_y_nb$main <- logd_y_nb$temp[[cls]]$NB
    } else {
      logd_y_nb$main <-
        cbind(logd_y_nb$main, logd_y_nb$temp[[cls]]$NB)
    }
  }
  logd_y_nb$main <- as_tibble(logd_y_nb$main)
  names(logd_y_nb$main) <- logd_y_nb$my_classes
  return(logd_y_nb)
}

pred_mix_nb <- function(logd_y_nb, prob_prior, prob_postodds) {
  pred_y_nb <- NULL
  pred_y_nb$my_classes <- logd_y_nb$my_classes
  pred_y_nb$stat_final <- logd_y_nb$stat_final
  for (cls in pred_y_nb$my_classes) {
    if (is.null(pred_y_nb$main)) {
      pred_y_nb$main <- 10**logd_y_nb$main[[cls]]*prob_prior[[cls]]
    } else {
      pred_y_nb$main <-
        cbind(pred_y_nb$main, 10**logd_y_nb$main[[cls]]*prob_prior[[cls]])
    }
  }
  pred_y_nb$main <- as_tibble(pred_y_nb$main)
  names(pred_y_nb$main) <- pred_y_nb$my_classes
  # normalize
  pred_y_nb$main <- as_tibble(pred_y_nb$main/rowSums(pred_y_nb$main))
  # argmax
  pred_y_nb$main <- pred_helper(pred_y_nb$main, prob_postodds)
  return(consolidate_nb(pred_y_nb))
}

consolidate_nb <- function(test_y_nb) {
  return(test_y_nb)
}

logit <- function(p) {
  return(log10(p/(1-p)))
}
