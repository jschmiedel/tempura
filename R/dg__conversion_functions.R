## functions to convert dG to fitness and vice versa
## for both stability and binding phenotypes

convert_dg2foldingfitness <- function(
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  fitness_scale,
  no_folded_states
) {
  if (fitness_scale == "lin" & no_folded_states == 1) {
      f_fitness <- function_folding_dg2fitness(
          f_ddg_var = f_ddg_var[[1]],
          f_dgwt = f_dgwt,
          f_fitwt = f_fitwt,
          f_fit0 = f_fit0
      )
  } else if (fitness_scale == "log" & no_folded_states == 1) {
      f_fitness <- function_folding_dg2logf(
          f_ddg_var = f_ddg_var[[1]],
          f_dgwt = f_dgwt,
          f_fitwt = f_fitwt,
          f_fit0 = f_fit0
      )
  } else if (fitness_scale == "lin" & no_folded_states == 2) {
      f_fitness <- function_folding_dg2fitness_4state(
          fA_ddg_var = f_ddg_var[[1]],
          fB_ddg_var = f_ddg_var[[2]],
          fA_dgwt = f_dgwt[1],
          fB_dgwt = f_dgwt[2],
          f_fitwt = f_fitwt,
          f_fit0 = f_fit0
      )
  } else if (fitness_scale == "log" & no_folded_states == 2) {
      f_fitness <- function_folding_dg2logf_4state(
          fA_ddg_var = f_ddg_var[[1]],
          fB_ddg_var = f_ddg_var[[2]],
          fA_dgwt = f_dgwt[1],
          fB_dgwt = f_dgwt[2],
          f_fitwt = f_fitwt,
          f_fit0 = f_fit0
      )
  }
  return(f_fitness)
}


convert_dg2foldinggradient <- function(
  f_ddg,
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  fitness,
  w,
  mutxvar,
  fitness_scale,
  no_folded_states,
  lambda,
  par
) {
  if (fitness_scale == "lin" & no_folded_states == 1) {
    gradient_f <- function_folding_dg2fitness_gradient( ##
      f_ddg = f_ddg[[1]],
      f_ddg_var = f_ddg_var[[1]],
      f_dgwt = f_dgwt,
      f_fitwt = f_fitwt,
      f_fit0 = f_fit0,
      fitness = fitness,
      w = w,
      mutxvar = mutxvar,
      lambda = lambda
    )
  } else if (fitness_scale == "log" & no_folded_states == 1) {
    gradient_f <- function_folding_dg2logf_gradient( ##
      f_ddg = f_ddg[[1]],
      f_ddg_var = f_ddg_var[[1]],
      f_dgwt = f_dgwt,
      f_fitwt = f_fitwt,
      f_fit0 = f_fit0,
      fitness = fitness,
      w = w,
      mutxvar = mutxvar,
      lambda = lambda
    )
  } else if (fitness_scale == "lin" & no_folded_states == 2) {
    gradient_f <- function_folding_dg2fitness_4state_gradient( ##
      fA_ddg = f_ddg[[1]],
      fB_ddg = f_ddg[[2]],
      fA_ddg_var = f_ddg_var[[1]],
      fB_ddg_var = f_ddg_var[[2]],
      fA_dgwt = f_dgwt[1],
      fB_dgwt = f_dgwt[2],
      f_fitwt = f_fitwt,
      f_fit0 = f_fit0,
      fitness = fitness,
      w = w,
      mutxvar = mutxvar,
      lambda = lambda
    )
  } else if (fitness_scale == "log" & no_folded_states == 2) {
    gradient_f <- function_folding_dg2logf_4state_gradient( ##
      fA_ddg = f_ddg[[1]],
      fB_ddg = f_ddg[[2]],
      fA_ddg_var = f_ddg_var[[1]],
      fB_ddg_var = f_ddg_var[[2]],
      fA_dgwt = f_dgwt[1],
      fB_dgwt = f_dgwt[2],
      f_fitwt = f_fitwt,
      f_fit0 = f_fit0,
      fitness = fitness,
      w = w,
      mutxvar = mutxvar,
      lambda = lambda
    )
  }
  return(gradient_f)
}



convert_dg2bindingfitness <- function(
    b_ddg_var,
    f_ddg_var,
    b_dgwt,
    f_dgwt,
    b_fitwt,
    b_fit0,
    fitness_scale,
    no_folded_states
) {
  if (fitness_scale == "lin" & no_folded_states == 1) {
    b_fitness <- function_binding_dg2fitness(
        b_ddg_var = b_ddg_var,
        f_ddg_var = f_ddg_var[[1]],
        b_dgwt = b_dgwt,
        f_dgwt = f_dgwt,
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0
      )
  } else if (fitness_scale == "log" & no_folded_states == 1) {
    b_fitness <- function_binding_dg2logf(
        b_ddg_var = b_ddg_var,
        f_ddg_var = f_ddg_var[[1]],
        b_dgwt = b_dgwt,
        f_dgwt = f_dgwt,
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0
      )
  } else if (fitness_scale == "lin" & no_folded_states == 2) {
    b_fitness <- function_binding_dg2fitness_4state(
        b_ddg_var = b_ddg_var,
        fA_ddg_var = f_ddg_var[[1]],
        fB_ddg_var = f_ddg_var[[2]],
        b_dgwt = b_dgwt,
        fA_dgwt = f_dgwt[1],
        fB_dgwt = f_dgwt[2],
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0
      )
  } else if (fitness_scale == "log" & no_folded_states == 2) {
    b_fitness <- function_binding_dg2logf_4state(
        b_ddg_var = b_ddg_var,
        fA_ddg_var = f_ddg_var[[1]],
        fB_ddg_var = f_ddg_var[[2]],
        b_dgwt = b_dgwt,
        fA_dgwt = f_dgwt[1],
        fB_dgwt = f_dgwt[2],
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0
      )
  }
  return(b_fitness)
}


convert_dg2bindinggradient <- function(
      b_ddg,
      f_ddg,
      b_ddg_var,
      f_ddg_var,
      b_dgwt,
      f_dgwt,
      b_fitwt,
      b_fit0,
      fitness,
      w,
      mutxvar,
      fitness_scale,
      no_folded_states,
      lambda
  ) {
  if (fitness_scale == "lin" & no_folded_states == 1) {
    gradient_b <- function_binding_dg2fitness_gradient( ##
        b_ddg = b_ddg,
        f_ddg = f_ddg[[1]],
        b_ddg_var = b_ddg_var,
        f_ddg_var = f_ddg_var[[1]],
        b_dgwt = b_dgwt,
        f_dgwt = f_dgwt,
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0,
        fitness = fitness,
        w = w,
        mutxvar = mutxvar,
        lambda = lambda
    )
  } else if (fitness_scale == "log" & no_folded_states == 1) {
    gradient_b <- function_binding_dg2logf_gradient( ##
        b_ddg = b_ddg,
        f_ddg = f_ddg[[1]],
        b_ddg_var = b_ddg_var,
        f_ddg_var = f_ddg_var[[1]],
        b_dgwt = b_dgwt,
        f_dgwt = f_dgwt,
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0,
        fitness = fitness,
        w = w,
        mutxvar = mutxvar,
        lambda = lambda
    )
  } else if (fitness_scale == "lin" & no_folded_states == 2) {
    gradient_b <- function_binding_dg2fitness_4state_gradient( ##
        b_ddg = b_ddg,
        fA_ddg = f_ddg[[1]],
        fB_ddg = f_ddg[[2]],
        b_ddg_var = b_ddg_var,
        fA_ddg_var = f_ddg_var[[1]],
        fB_ddg_var = f_ddg_var[[2]],
        b_dgwt = b_dgwt,
        fA_dgwt = f_dgwt[1],
        fB_dgwt = f_dgwt[2],
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0,
        fitness = fitness,
        w = w,
        mutxvar = mutxvar,
        lambda = lambda
    )
  } else if (fitness_scale == "log" & no_folded_states == 2) {
    gradient_b <- function_binding_dg2logf_4state_gradient( ##
        b_ddg = b_ddg,
        fA_ddg = f_ddg[[1]],
        fB_ddg = f_ddg[[2]],
        b_ddg_var = b_ddg_var,
        fA_ddg_var = f_ddg_var[[1]],
        fB_ddg_var = f_ddg_var[[2]],
        b_dgwt = b_dgwt,
        fA_dgwt = f_dgwt[1],
        fB_dgwt = f_dgwt[2],
        b_fitwt = b_fitwt,
        b_fit0 = b_fit0,
        fitness = fitness,
        w = w,
        mutxvar = mutxvar,
        lambda = lambda
    )
  }
  return(gradient_b)
}


function_folding_dg2fitness <- function(
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  rt = 1.99e-3 * 310.15
) {
  f_fitness <- f_fit0 +
    (f_fitwt - f_fit0) * (1 + exp(f_dgwt / rt)) /
      (1 + exp((f_dgwt + f_ddg_var) / rt))
}

function_folding_dg2fitness_gradient <- function(
  f_ddg,
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  par,
  rt = 1.99e-3 * 310.15
) {
  cf <- (f_fitwt - f_fit0) / rt
  ei <- exp(f_ddg_var / rt)
  ewt <- exp(f_dgwt / rt)
  eiwt <- ei * ewt
  dfitness_df_ddg <- -cf * eiwt * (1 + ewt) / (1 + eiwt)^2 #
  dfitness_df_ddg[is.na(dfitness_df_ddg)] <- 0

  dfitness_df_dgwt <- cf * (ewt - eiwt) / (1 + eiwt)^2 #
  dfitness_df_fitwt <- (1 + ewt) / (1 + eiwt) #
  dfitness_df_fit0 <- (eiwt - ewt) / (1 + eiwt) #
  fitness_pred <- function_folding_dg2fitness(
    f_ddg_var = f_ddg_var,
    f_dgwt = f_dgwt,
    f_fitwt = f_fitwt,
    f_fit0 = f_fit0
  )
  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2
  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_f_ddg <- mutxvar %*% (derr_dfitness * dfitness_df_ddg) +
                    2 * f_ddg * lambda #regularization term
  gradient_f_dgwt <- sum(derr_dfitness * dfitness_df_dgwt, na.rm = T)
  gradient_f_fitwt <- sum(derr_dfitness * dfitness_df_fitwt, na.rm = T)
  gradient_f_fit0 <- sum(derr_dfitness * dfitness_df_fit0, na.rm = T)

  return(list(
    f_ddg = gradient_f_ddg,
    f_dgwt = gradient_f_dgwt,
    f_fitwt = gradient_f_fitwt,
    f_fit0 = gradient_f_fit0
  ))
}

function_binding_dg2fitness <- function(
  b_ddg_var,
  f_ddg_var,
  b_dgwt,
  f_dgwt,
  b_fitwt,
  b_fit0,
  rt = 1.99e-3 * 310.15
) {
  b_fitness <- b_fit0 +
    (b_fitwt - b_fit0) * (1 + exp(b_dgwt / rt) * (1 + exp(f_dgwt / rt))) /
     (1 + exp((b_dgwt + b_ddg_var) / rt) * (1 + exp((f_dgwt + f_ddg_var) / rt)))
}

function_binding_dg2fitness_gradient <- function(
  f_ddg,
  b_ddg,
  f_ddg_var,
  b_ddg_var,
  f_dgwt,
  b_dgwt,
  b_fitwt,
  b_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (b_fitwt - b_fit0) / rt
  ebi <- exp(b_ddg_var / rt)
  ebwt <- exp(b_dgwt / rt)
  efi <- exp(f_ddg_var / rt)
  efwt <- exp(f_dgwt / rt)

  ebiwt <- ebi * ebwt
  efiwt <- efi * efwt
  eiwt <- ebiwt * (1 + efiwt)
  ewt <- ebwt * (1 + efwt)

  dfitness_df_ddg <- -cf * ebiwt * efiwt * (1 + ewt) /
                    (1 + eiwt) ^ 2
  dfitness_df_ddg[is.na(dfitness_df_ddg)] <- 0

  dfitness_db_ddg <- -cf * (1 + ewt) * eiwt /
                      (1 + eiwt) ^ 2
  dfitness_db_ddg[is.na(dfitness_db_ddg)] <- 0

  dfitness_df_dgwt <- cf * ebwt * efwt * (1 + eiwt - ebi * efi * (1 + ewt)) /
    (1 + eiwt) ^ 2 #
  dfitness_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!
  dfitness_db_fitwt <- (1 + ewt) / (1 + eiwt) #
  dfitness_db_fit0 <- (eiwt - ewt) / (1 + eiwt) #

  fitness_pred <- function_binding_dg2fitness(
    f_ddg_var = f_ddg_var,
    b_ddg_var = b_ddg_var,
    f_dgwt = f_dgwt,
    b_dgwt = b_dgwt,
    b_fitwt = b_fitwt,
    b_fit0 = b_fit0
  )
  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2
  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_f_ddg <- mutxvar %*% (derr_dfitness * dfitness_df_ddg) +
                    2 * f_ddg * lambda # regularization term
  gradient_b_ddg <- mutxvar %*% (derr_dfitness * dfitness_db_ddg) +
                    2 * b_ddg * lambda  # regularization term
  gradient_f_dgwt <- sum(derr_dfitness * dfitness_df_dgwt, na.rm = T)
  gradient_b_dgwt <- sum(derr_dfitness * dfitness_db_dgwt, na.rm = T)
  gradient_b_fitwt <- sum(derr_dfitness * dfitness_db_fitwt, na.rm = T)
  gradient_b_fit0 <- sum(derr_dfitness * dfitness_db_fit0, na.rm = T)

  return(list(
    f_ddg = gradient_f_ddg,
    b_ddg = gradient_b_ddg,
    f_dgwt = gradient_f_dgwt,
    b_dgwt = gradient_b_dgwt,
    b_fitwt = gradient_b_fitwt,
    b_fit0 = gradient_b_fit0
  ))
}


function_folding_dg2fitness_4state <- function(
  fA_ddg_var,
  fB_ddg_var,
  fA_dgwt,
  fB_dgwt,
  f_fitwt,
  f_fit0,
  rt = 1.99e-3 * 310.15
) {
  sf <- f_fit0 +
    (f_fitwt - f_fit0) *
    (exp(-(fA_dgwt + fA_ddg_var) / rt) + exp(-(fB_dgwt + fB_ddg_var) / rt)) *
    (1 + exp(-fA_dgwt / rt) + exp(-fB_dgwt / rt)) /
      ((1 + exp(-(fA_dgwt + fA_ddg_var) / rt) + exp(-(fB_dgwt + fB_ddg_var) / rt)) *
      (exp(-fA_dgwt / rt) + exp(-fB_dgwt / rt)))
}

function_folding_dg2fitness_4state_gradient <- function(
  fA_ddg,
  fB_ddg,
  fA_ddg_var,
  fB_ddg_var,
  fA_dgwt,
  fB_dgwt,
  f_fitwt,
  f_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (f_fitwt - f_fit0) / rt

  eAi <- exp(-fA_ddg_var / rt)
  eBi <- exp(-fB_ddg_var / rt)

  eAwt <- exp(-fA_dgwt / rt)
  eBwt <- exp(-fB_dgwt / rt)

  eAiwt <- eAi * eAwt
  eBiwt <- eBi * eBwt

  ewt <- eAwt + eBwt
  eiwt <- eAiwt + eBiwt

  dfitness_dfA_ddg <- -cf * eAiwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2
  dfitness_dfA_ddg[is.na(dfitness_dfA_ddg)] <- 0

  dfitness_dfB_ddg <- -cf * eBiwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2
  dfitness_dfB_ddg[is.na(dfitness_dfB_ddg)] <- 0

  dfitness_dfA_dgwt <- cf * (eAwt * eiwt * (1 + eiwt)  - eAiwt * ewt * (1 + ewt)) /
    ((1 + eiwt) * ewt)^2 #

  dfitness_dfB_dgwt <- cf * (eBwt * eiwt * (1 + eiwt)  - eBiwt * ewt * (1 + ewt)) /
    ((1 + eiwt) * ewt)^2 #

  dfitness_df_fitwt <- eiwt * (1 + ewt) / ((1 + eiwt) * ewt) #

  dfitness_df_fit0 <- (ewt - eiwt) / ((1 + eiwt) * ewt) #

  fitness_pred <- function_folding_dg2fitness_4state(
    fA_ddg_var = fA_ddg_var,
    fB_ddg_var = fB_ddg_var,
    fA_dgwt = fA_dgwt,
    fB_dgwt = fB_dgwt,
    f_fitwt = f_fitwt,
    f_fit0 = f_fit0
  )

  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2
  derr_dfitness[is.na(derr_dfitness)] <- 0

  gradient_fA_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfA_ddg)+
                      ((1 + sqrt(2)) * fA_ddg - (1 - sqrt(2)) * fB_ddg) * lambda   #regularization term
  gradient_fB_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfB_ddg) +
                      ((1 + sqrt(2)) * fB_ddg - (1 - sqrt(2)) * fA_ddg) * lambda  #regularization term
  gradient_fA_dgwt <- sum(derr_dfitness * dfitness_dfA_dgwt, na.rm = T)
  gradient_fB_dgwt <- sum(derr_dfitness * dfitness_dfB_dgwt, na.rm = T)
  gradient_f_fitwt <- sum(derr_dfitness * dfitness_df_fitwt, na.rm = T)
  gradient_f_fit0 <- sum(derr_dfitness * dfitness_df_fit0, na.rm = T)

  return(list(
    fA_ddg = gradient_fA_ddg,
    fB_ddg = gradient_fB_ddg,
    fA_dgwt = gradient_fA_dgwt,
    fB_dgwt = gradient_fB_dgwt,
    f_fitwt = gradient_f_fitwt,
    f_fit0 = gradient_f_fit0
  ))
}

function_binding_dg2fitness_4state <- function(
  b_ddg_var,
  fA_ddg_var,
  fB_ddg_var,
  b_dgwt,
  fA_dgwt,
  fB_dgwt,
  b_fitwt,
  b_fit0,
  rt = 1.99e-3 * 310.15
) {
  b_fitness <- b_fit0 +
    (b_fitwt - b_fit0) * (1 + exp(b_dgwt / rt) * (1 + exp(fB_dgwt / rt) * (1 + exp(-fA_dgwt / rt)))) /
     (1 + exp((b_dgwt + b_ddg_var) / rt) *
        (1 + exp((fB_dgwt + fB_ddg_var) / rt) * (1 + exp(-(fA_dgwt + fA_ddg_var) / rt))))
}

function_binding_dg2fitness_4state_gradient <- function(
  fA_ddg,
  fB_ddg,
  b_ddg,
  fA_ddg_var,
  fB_ddg_var,
  b_ddg_var,
  fA_dgwt,
  fB_dgwt,
  b_dgwt,
  b_fitwt,
  b_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (b_fitwt - b_fit0) / rt

  ebi <- exp(b_ddg_var / rt)
  ebwt <- exp(b_dgwt / rt)
  ebiwt <- ebi * ebwt

  eAfi <- exp(-fA_ddg_var / rt)
  eAfwt <- exp(-fA_dgwt / rt)
  eAfiwt <- eAfi * eAfwt

  eBfi <- exp(fB_ddg_var / rt)
  eBfwt <- exp(fB_dgwt / rt)
  eBfiwt <- eBfi * eBfwt

  eiwt <- ebiwt * (1 + eBfiwt * (1 + eAfiwt))
  ewt <- ebwt * (1 + eBfwt * (1 + eAfwt))

  dfitness_db_ddg <- -cf * (1 + ewt) * eiwt /
                      (1 + eiwt) ^ 2
  dfitness_db_ddg[is.na(dfitness_db_ddg)] <- 0

  dfitness_dfA_ddg <- cf * ebiwt * eBfiwt * eAfiwt * (1 + ewt) /
                      (1 + eiwt) ^ 2
  dfitness_dfA_ddg[is.na(dfitness_dfA_ddg)] <- 0

  dfitness_dfB_ddg <- -cf * ebiwt * eBfiwt * (1 + eAfiwt) * (1 + ewt) /
                        (1 + eiwt) ^ 2
  dfitness_dfB_ddg[is.na(dfitness_dfB_ddg)] <- 0

  dfitness_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!

  dfitness_dfA_dgwt <- cf * (-ebwt * eBfwt * eAfwt * (1 + eiwt) + (1 + ewt) * ebiwt * eBfiwt * eAfiwt) /
    (1 + eiwt) ^ 2 #

  dfitness_dfB_dgwt <- cf * (ebwt * eBfwt * (1 + eAfwt) * (1 + eiwt) -
    (1 + ewt) * ebiwt * eBfiwt * (1 + eAfiwt)) /
    (1 + eiwt) ^ 2 #

  dfitness_db_fitwt <- (1 + ewt) / (1 + eiwt) #
  dfitness_db_fit0 <- 1 - (1 + ewt) / (1 + eiwt) #

  fitness_pred <- function_binding_dg2fitness_4state(
    b_ddg_var = b_ddg_var,
    fA_ddg_var = fA_ddg_var,
    fB_ddg_var = fB_ddg_var,
    b_dgwt = b_dgwt,
    fA_dgwt = fA_dgwt,
    fB_dgwt = fB_dgwt,
    b_fitwt = b_fitwt,
    b_fit0 = b_fit0
  )

  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2
  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_b_ddg <- mutxvar %*% (derr_dfitness * dfitness_db_ddg) +
                    2 * b_ddg * lambda  #regularization term
  gradient_fA_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfA_ddg) +
                    ((1 + sqrt(2)) * fA_ddg - (1 - sqrt(2)) * fB_ddg) * lambda  #regularization term
  gradient_fB_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfB_ddg) +
                    ((1 + sqrt(2)) * fB_ddg - (1 - sqrt(2)) * fA_ddg) * lambda  #regularization term

  gradient_b_dgwt <- sum(derr_dfitness * dfitness_db_dgwt, na.rm = T)
  gradient_fA_dgwt <- sum(derr_dfitness * dfitness_dfA_dgwt, na.rm = T)
  gradient_fB_dgwt <- sum(derr_dfitness * dfitness_dfB_dgwt, na.rm = T)

  gradient_b_fitwt <- sum(derr_dfitness * dfitness_db_fitwt, na.rm = T)
  gradient_b_fit0 <- sum(derr_dfitness * dfitness_db_fit0, na.rm = T)

  return(list(
    b_ddg = gradient_b_ddg,
    fA_ddg = gradient_fA_ddg,
    fB_ddg = gradient_fB_ddg,
    b_dgwt = gradient_b_dgwt,
    fA_dgwt = gradient_fA_dgwt,
    fB_dgwt = gradient_fB_dgwt,
    b_fitwt = gradient_b_fitwt,
    b_fit0 = gradient_b_fit0
  ))
}





## functions to convert dG to log-fitness and vice versa
## for both stability and binding phenotypes

#assumes fwt and f0 values are on log scale

function_folding_dg2logf <- function(
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  rt = 1.99e-3 * 310.15
) {
  f_fitness <- log(
    exp(f_fit0) +
    (exp(f_fitwt) - exp(f_fit0)) *
    (1 + exp(f_dgwt / rt)) /
    (1 + exp((f_dgwt + f_ddg_var) / rt))
  )
}

function_folding_dg2logf_gradient <- function(
  f_ddg,
  f_ddg_var,
  f_dgwt,
  f_fitwt,
  f_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (exp(f_fitwt) - exp(f_fit0)) / rt
  ei <- exp(f_ddg_var / rt)
  ewt <- exp(f_dgwt / rt)
  eiwt <- ei * ewt

  p <- (1 + ewt) / (1 + eiwt)
  f_pred_nolog <- exp(f_fit0) + (exp(f_fitwt) - exp(f_fit0)) * p

  dfitness_df_ddg <- -cf * eiwt * (1 + ewt) / (1 + eiwt)^2 #
  dfitness_df_ddg[is.na(dfitness_df_ddg)] <- 0

  dfitness_df_dgwt <- cf * (ewt - eiwt) / (1 + eiwt)^2
  dfitness_df_fitwt <- exp(f_fitwt) * p
  dfitness_df_fit0 <- exp(f_fit0) * (1 - p)

  fitness_pred <- function_folding_dg2logf(
    f_ddg_var = f_ddg_var,
    f_dgwt = f_dgwt,
    f_fitwt = f_fitwt,
    f_fit0 = f_fit0
  )
  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2 / f_pred_nolog

  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_f_ddg <- mutxvar %*% (derr_dfitness * dfitness_df_ddg) +
                    2 * f_ddg * lambda  #regularization term
  gradient_f_dgwt <- sum(derr_dfitness * dfitness_df_dgwt, na.rm = T)
  gradient_f_fitwt <- sum(derr_dfitness * dfitness_df_fitwt, na.rm = T)
  gradient_f_fit0 <- sum(derr_dfitness * dfitness_df_fit0, na.rm = T)

  return(list(
    f_ddg = gradient_f_ddg,
    f_dgwt = gradient_f_dgwt,
    f_fitwt = gradient_f_fitwt,
    f_fit0 = gradient_f_fit0
  ))
}

function_binding_dg2logf <- function(
  b_ddg_var,
  f_ddg_var,
  b_dgwt,
  f_dgwt,
  b_fitwt,
  b_fit0,
  rt = 1.99e-3 * 310.15
) {
  b_fitness <- log(exp(b_fit0) +
    (exp(b_fitwt) - exp(b_fit0)) * (1 + exp(b_dgwt / rt) * (1 + exp(f_dgwt / rt))) /
     (1 + exp((b_dgwt + b_ddg_var) / rt) * (1 + exp((f_dgwt + f_ddg_var) / rt))))
}

function_binding_dg2logf_gradient <- function(
  f_ddg,
  b_ddg,
  f_ddg_var,
  b_ddg_var,
  f_dgwt,
  b_dgwt,
  b_fitwt,
  b_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (exp(b_fitwt) - exp(b_fit0)) / rt
  ebi <- exp(b_ddg_var / rt)
  ebwt <- exp(b_dgwt / rt)
  efi <- exp(f_ddg_var / rt)
  efwt <- exp(f_dgwt / rt)

  ebiwt <- ebi * ebwt
  efiwt <- efi * efwt
  eiwt <- ebiwt * (1 + efiwt)
  ewt <- ebwt * (1 + efwt)

  p <- (1 + ewt) / (1 + eiwt)
  f_pred_nolog <- exp(b_fit0) + (exp(b_fitwt) - exp(b_fit0)) * p

  dfitness_df_ddg <- -cf * ebiwt * efiwt * (1 + ewt) /
                    (1 + eiwt) ^ 2
  dfitness_df_ddg[is.na(dfitness_df_ddg)] <- 0

  dfitness_db_ddg <- -cf * (1 + ewt) * eiwt /
                    (1 + eiwt) ^ 2
  dfitness_db_ddg[is.na(dfitness_db_ddg)] <- 0

  dfitness_df_dgwt <- cf * ebwt * efwt * (1 + eiwt - ebi * efi * (1 + ewt)) /
    (1 + eiwt) ^ 2 #
  dfitness_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!

  dfitness_db_fitwt <- exp(b_fitwt) * p #

  dfitness_db_fit0 <- exp(b_fit0) * (1 - p)

  fitness_pred <- function_binding_dg2logf(
    f_ddg_var = f_ddg_var,
    b_ddg_var = b_ddg_var,
    f_dgwt = f_dgwt,
    b_dgwt = b_dgwt,
    b_fitwt = b_fitwt,
    b_fit0 = b_fit0
  )
  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2 / f_pred_nolog

  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_f_ddg <- mutxvar %*% (derr_dfitness * dfitness_df_ddg) +
                    2 * f_ddg * lambda  # regularization term
  gradient_b_ddg <- mutxvar %*% (derr_dfitness * dfitness_db_ddg) +
                    2 * b_ddg * lambda  # regularization term
  gradient_f_dgwt <- sum(derr_dfitness * dfitness_df_dgwt, na.rm = T)
  gradient_b_dgwt <- sum(derr_dfitness * dfitness_db_dgwt, na.rm = T)
  gradient_b_fitwt <- sum(derr_dfitness * dfitness_db_fitwt, na.rm = T)
  gradient_b_fit0 <- sum(derr_dfitness * dfitness_db_fit0, na.rm = T)

  return(list(
    f_ddg = gradient_f_ddg,
    b_ddg = gradient_b_ddg,
    f_dgwt = gradient_f_dgwt,
    b_dgwt = gradient_b_dgwt,
    b_fitwt = gradient_b_fitwt,
    b_fit0 = gradient_b_fit0
  ))
}


function_folding_dg2logf_4state <- function(
  fA_ddg_var,
  fB_ddg_var,
  fA_dgwt,
  fB_dgwt,
  f_fitwt,
  f_fit0,
  rt = 1.99e-3 * 310.15
) {
  f_fitness <- log(exp(f_fit0) +
    (exp(f_fitwt) - exp(f_fit0)) *
    (exp(-(fA_dgwt + fA_ddg_var) / rt) + exp(-(fB_dgwt + fB_ddg_var) / rt)) *
    (1 + exp(-fA_dgwt / rt) + exp(-fB_dgwt / rt)) /
      ((1 + exp(-(fA_dgwt + fA_ddg_var) / rt) + exp(-(fB_dgwt + fB_ddg_var) / rt)) *
      (exp(-fA_dgwt / rt) + exp(-fB_dgwt / rt))))
}

function_folding_dg2logf_4state_gradient <- function(
  fA_ddg,
  fB_ddg,
  fA_ddg_var,
  fB_ddg_var,
  fA_dgwt,
  fB_dgwt,
  f_fitwt,
  f_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (exp(f_fitwt) - exp(f_fit0)) / rt

  eAi <- exp(-fA_ddg_var / rt)
  eBi <- exp(-fB_ddg_var / rt)

  eAwt <- exp(-fA_dgwt / rt)
  eBwt <- exp(-fB_dgwt / rt)

  eAiwt <- eAi * eAwt
  eBiwt <- eBi * eBwt

  ewt <- eAwt + eBwt
  eiwt <- eAiwt + eBiwt

  p <- eiwt * (1 + ewt) / ((1 + eiwt) * ewt)
  f_pred_nolog <- exp(f_fit0) + (exp(f_fitwt) - exp(f_fit0)) * p


  dfitness_dfA_ddg <- -cf * eAiwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2
  dfitness_dfA_ddg[is.na(dfitness_dfA_ddg)] <- 0

  dfitness_dfB_ddg <- -cf * eBiwt * (1 + ewt) * ewt / ((1 + eiwt) * ewt)^2
  dfitness_dfB_ddg[is.na(dfitness_dfB_ddg)] <- 0

  dfitness_dfA_dgwt <- cf * (eAwt * eiwt * (1 + eiwt)  - eAiwt * ewt * (1 + ewt)) /
    ((1 + eiwt) * ewt)^2 #

  dfitness_dfB_dgwt <- cf * (eBwt * eiwt * (1 + eiwt)  - eBiwt * ewt * (1 + ewt)) /
    ((1 + eiwt) * ewt)^2 #

  dfitness_df_fitwt <-  exp(f_fitwt) * p#

  dfitness_df_fit0 <- exp(f_fit0) * (1 - p) #

  fitness_pred <- function_folding_dg2logf_4state(
    fA_ddg_var = fA_ddg_var,
    fB_ddg_var = fB_ddg_var,
    fA_dgwt = fA_dgwt,
    fB_dgwt = fB_dgwt,
    f_fitwt = f_fitwt,
    f_fit0 = f_fit0
  )

  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2 / f_pred_nolog
  derr_dfitness[is.na(derr_dfitness)] <- 0

  gradient_fA_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfA_ddg) +
                    ((1 + sqrt(2)) * fA_ddg - (1 - sqrt(2)) * fB_ddg) * lambda  #regularization term
  gradient_fB_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfB_ddg) +
                    ((1 + sqrt(2)) * fB_ddg - (1 - sqrt(2)) * fA_ddg) * lambda  #regularization term
  gradient_fA_dgwt <- sum(derr_dfitness * dfitness_dfA_dgwt, na.rm = T)
  gradient_fB_dgwt <- sum(derr_dfitness * dfitness_dfB_dgwt, na.rm = T)
  gradient_f_fitwt <- sum(derr_dfitness * dfitness_df_fitwt, na.rm = T)
  gradient_f_fit0 <- sum(derr_dfitness * dfitness_df_fit0, na.rm = T)

  return(list(
    fA_ddg = gradient_fA_ddg,
    fB_ddg = gradient_fB_ddg,
    fA_dgwt = gradient_fA_dgwt,
    fB_dgwt = gradient_fB_dgwt,
    f_fitwt = gradient_f_fitwt,
    f_fit0 = gradient_f_fit0
  ))
}

function_binding_dg2logf_4state <- function(
  b_ddg_var,
  fA_ddg_var,
  fB_ddg_var,
  b_dgwt,
  fA_dgwt,
  fB_dgwt,
  b_fitwt,
  b_fit0,
  rt = 1.99e-3 * 310.15
) {
  b_fitness <- log(exp(b_fit0) +
    (exp(b_fitwt) - exp(b_fit0)) * (1 + exp(b_dgwt / rt) * (1 + exp(fB_dgwt / rt) * (1 + exp(-fA_dgwt / rt)))) /
     (1 + exp((b_dgwt + b_ddg_var) / rt) *
        (1 + exp((fB_dgwt + fB_ddg_var) / rt) * (1 + exp(-(fA_dgwt + fA_ddg_var) / rt)))))
}

function_binding_dg2logf_4state_gradient <- function(
  fA_ddg,
  fB_ddg,
  b_ddg,
  fA_ddg_var,
  fB_ddg_var,
  b_ddg_var,
  fA_dgwt,
  fB_dgwt,
  b_dgwt,
  b_fitwt,
  b_fit0,
  fitness,
  w,
  mutxvar,
  lambda,
  rt = 1.99e-3 * 310.15
) {
  cf <- (exp(b_fitwt) - exp(b_fit0)) / rt

  ebi <- exp(b_ddg_var / rt)
  ebwt <- exp(b_dgwt / rt)
  ebiwt <- ebi * ebwt

  eAfi <- exp(-fA_ddg_var / rt)
  eAfwt <- exp(-fA_dgwt / rt)
  eAfiwt <- eAfi * eAfwt

  eBfi <- exp(fB_ddg_var / rt)
  eBfwt <- exp(fB_dgwt / rt)
  eBfiwt <- eBfi * eBfwt

  eiwt <- ebiwt * (1 + eBfiwt * (1 + eAfiwt))
  ewt <- ebwt * (1 + eBfwt * (1 + eAfwt))

  p <- (1 + ewt) / (1 + eiwt)
  f_pred_nolog <- exp(b_fit0) + (exp(b_fitwt) - exp(b_fit0)) * p


  dfitness_db_ddg <- -cf * (1 + ewt) * eiwt /
                    (1 + eiwt) ^ 2
  dfitness_db_ddg[is.na(dfitness_db_ddg)] <- 0

  dfitness_dfA_ddg <- cf * ebiwt * eBfiwt * eAfiwt * (1 + ewt) /
                      (1 + eiwt) ^ 2
  dfitness_dfA_ddg[is.na(dfitness_dfA_ddg)] <- 0

  dfitness_dfB_ddg <- -cf * ebiwt * eBfiwt * (1 + eAfiwt) * (1 + ewt) /
                      (1 + eiwt) ^ 2
  dfitness_dfB_ddg[is.na(dfitness_dfB_ddg)] <- 0

  dfitness_db_dgwt <- cf * (ewt - eiwt) /
    (1 + eiwt) ^ 2 #!

  dfitness_dfA_dgwt <- cf * (-ebwt * eBfwt * eAfwt * (1 + eiwt) + (1 + ewt) * ebiwt * eBfiwt * eAfiwt) /
    (1 + eiwt) ^ 2 #

  dfitness_dfB_dgwt <- cf * (ebwt * eBfwt * (1 + eAfwt) * (1 + eiwt) -
    (1 + ewt) * ebiwt * eBfiwt * (1 + eAfiwt)) /
    (1 + eiwt) ^ 2 #

  dfitness_db_fitwt <- exp(b_fitwt) * p #
  dfitness_db_fit0 <- exp(b_fit0) * (1 - p) #

  fitness_pred <- function_binding_dg2logf_4state(
    b_ddg_var = b_ddg_var,
    fA_ddg_var = fA_ddg_var,
    fB_ddg_var = fB_ddg_var,
    b_dgwt = b_dgwt,
    fA_dgwt = fA_dgwt,
    fB_dgwt = fB_dgwt,
    b_fitwt = b_fitwt,
    b_fit0 = b_fit0
  )

  derr_dfitness <- -2 * (fitness - fitness_pred) / w^2 / f_pred_nolog
  derr_dfitness[is.na(derr_dfitness)] <- 0
  gradient_b_ddg <- mutxvar %*% (derr_dfitness * dfitness_db_ddg) +
                    2 * b_ddg * lambda  #regularization term
  gradient_fA_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfA_ddg) +
                    ((1 + sqrt(2)) * fA_ddg - (1 - sqrt(2)) * fB_ddg) * lambda  #regularization term
  gradient_fB_ddg <- mutxvar %*% (derr_dfitness * dfitness_dfB_ddg) +
                    ((1 + sqrt(2)) * fB_ddg - (1 - sqrt(2)) * fA_ddg) * lambda  #regularization term

  gradient_b_dgwt <- sum(derr_dfitness * dfitness_db_dgwt, na.rm = T)
  gradient_fA_dgwt <- sum(derr_dfitness * dfitness_dfA_dgwt, na.rm = T)
  gradient_fB_dgwt <- sum(derr_dfitness * dfitness_dfB_dgwt, na.rm = T)

  gradient_b_fitwt <- sum(derr_dfitness * dfitness_db_fitwt, na.rm = T)
  gradient_b_fit0 <- sum(derr_dfitness * dfitness_db_fit0, na.rm = T)

  return(list(
    b_ddg = gradient_b_ddg,
    fA_ddg = gradient_fA_ddg,
    fB_ddg = gradient_fB_ddg,
    b_dgwt = gradient_b_dgwt,
    fA_dgwt = gradient_fA_dgwt,
    fB_dgwt = gradient_fB_dgwt,
    b_fitwt = gradient_b_fitwt,
    b_fit0 = gradient_b_fit0
  ))
}
