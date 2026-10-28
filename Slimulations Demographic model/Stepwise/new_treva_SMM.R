install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE,repos = "https://cloud.r-project.org/")
  }
}

required_packages <- c("phytools", "expm", "ape", "dplyr", "stringr", "GenSA")
lapply(required_packages, install_if_missing)
#####SUPER MEGA TRY 1 OF MERGE 
library(ape)
library(phytools)
library(expm)
library(stringr)
library(GenSA)


###Arguments

Model = c(0,0)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Please provide tree_dir, str_count_dir and run times.")
tree_dir <- args[1]
str_count_dir <- args[2]

samples <- args[3]
runn <- as.integer(2)
nodes<-as.integer(samples)*2-1
if (is.na(runn) || runn <= 0) stop("run times (args[4]) must be a positive integer.")

if (length(args) >= 5 && nzchar(args[5])) {
  Ex <- as.numeric(args[5])
  if (is.na(Ex)) Ex <- 2
} else {
  Ex <- 2
}

if (length(args) >= 6 && nzchar(args[6])) {
  v5 <- as.integer(args[6])
  if (!is.na(v5) && v5 %in% c(1,2)) Model[1] <- v5
}
if (length(args) >= 7 && nzchar(args[7])) {
  v6 <- as.integer(args[7])
  if (!is.na(v6) && v6 %in% c(1,2)) Model[2] <- v6
}
extract_number <- function(filename) {
  as.numeric(sub(".*STR(\\d+)_newick\\.txt", "\\1", filename))
}
tree_dir<-dirname(str_count_dir)
tree_files <- list.files(tree_dir, pattern = "*newick.txt", full.names = TRUE)
#str_count_files <- list.files(str_count_dir, pattern = "*.txt", full.names = TRUE)
tree_files <- tree_files[order(sapply(tree_files, extract_number))]
#str_count_files <- str_count_files[order(sapply(str_count_files, extract_number))]
specie=read.csv(str_count_dir)
specie[] <- lapply(specie, function(x) as.integer(gsub("\\[|\\]", "", x)))

strs<-names(specie)

str_counts <- lapply(strs, function(s) {

  str_values <- unlist(specie[s])
  # Create named factor
  factor(setNames(str_values, paste0(0:nodes)))
})

num_trees <- length(tree_files)
#tree_files
trees <- lapply(tree_files, read.tree)
print(trees)
#extract_number <- function(path) as.numeric(str_extract(basename(path), "\\d+"))
#tree_files <- list.files(tree_dir, pattern = "\\.txt$", full.names = TRUE)
#str_count_files <- list.files(str_count_dir, pattern = "\\.txt$", full.names = TRUE)
#tree_files <- tree_files[order(sapply(tree_files, extract_number))]
#str_count_files <- str_count_files[order(sapply(str_count_files, extract_number))]
trees <- lapply(tree_files, read.tree)

#str_counts <- lapply(str_count_files, function(f) {
 # x <- unlist(read.csv(f, header = FALSE))
 # factor(setNames(x, paste0(1:length(x))))
#})

s <- 0.0
u <- 0.5
v <- 0.0

mu_lower <- 1e-8
mu_upper <- .1
p_lower <- 0.0
p_upper <- 0.9
m_lower <- .1
m_upper <- 0.9
g <- 0

log_sum_exp <- function(v) {
  m <- max(v)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(v - m)))
}

###TPM Functions###

create_transition_matrix_tpm <- function(mu, k, m, s, u, v, p, lo, hi) {
  alpha_microsat <- function(i, u, v, lo) {
    max(0, min(1, u - (v * (i - lo))))
  }
  
  beta_microsat <- function(i, s, lo) {
    1 + ((i - lo) * s)
  }
  
  gamma_microsat <- function(m, i, j, lo, hi) {
    num <- m * (1 - m) ^ (abs(i - j) - 1)
    if (lo <= i && i < j && j <= hi) {
      denom <- 1 - (1 - m) ^ (hi - i)
    } else if (lo <= j && j < i && i <= hi) {
      denom <- 1 - (1 - m) ^ (i - lo)
    } else {
      stop("microsat model error: i and j must be between lo and hi")
    }
    num / denom
  }
  
  nalleles <- hi - lo + 1
  Q <- matrix(0, nrow = nalleles, ncol = nalleles)
  
  for (i in lo:hi) {
    for (j in lo:hi) {
      ii <- i - lo + 1
      jj <- j - lo + 1
      if (i == j) {
        next
      } else if (i == j - 1) {
        Q[ii, jj] <- beta_microsat(i, s, lo) * alpha_microsat(i, u, v, lo) *
          (p + (1 - p) * gamma_microsat(m, i, j, lo, hi))
      } else if (i < j - 1) {
        Q[ii, jj] <- beta_microsat(i, s, lo) * alpha_microsat(i, u, v, lo) *
          (1 - p) * gamma_microsat(m, i, j, lo, hi)
      } else if (i == j + 1) {
        Q[ii, jj] <- beta_microsat(i, s, lo) * (1 - alpha_microsat(i, u, v, lo)) *
          (p + (1 - p) * gamma_microsat(m, i, j, lo, hi))
      } else if (i > j + 1) {
        Q[ii, jj] <- beta_microsat(i, s, lo) * (1 - alpha_microsat(i, u, v, lo)) *
          (1 - p) * gamma_microsat(m, i, j, lo, hi)
      }
    }
  }
  
  for (i in 1:nalleles) {
    Q[i, i] <- -sum(Q[i, -i])
  }
  
  Q <- Q * mu
  
  return(Q)
}

pruning_tpm <- function(mu, tree, x, p, m, Ex, ...) {
  pw <- reorder(tree, "postorder")
  pw$tip.label <- gsub("n", "", pw$tip.label)
  lo <- min(as.numeric(levels(x)), na.rm = TRUE)
  hi <- max(as.numeric(levels(x)), na.rm = TRUE)
  k0 <- hi - lo + 1
  J <- floor((Ex) * k0)
  g <- 0
  
  lo <- floor(lo - (J/2))
  if ( lo < 0){
    g <- (-1)*(lo)
    lo = 0
  }
  hi <- floor(hi + (J/2)) + g
  k <- hi - lo + 1
  
  Q <- create_transition_matrix_tpm(mu, k, m, s, u, v, p, lo, hi)
  
  L <- rbind(to.matrix(x[pw$tip.label], as.character(seq(lo, hi))),
             matrix(0, tree$Nnode, k, dimnames = list(1:tree$Nnode + Ntip(tree))))
  
  L[L == 0] <- -Inf
  L[L == 1] <- 0
  
  nn <- unique(pw$edge[, 1])
  for (i in 1:length(nn)) {
    ee <- which(pw$edge[, 1] == nn[i])
    logPP <- matrix(NA, length(ee), k)
    
    for (j in 1:length(ee)) {
      P <- expm(Q * pw$edge.length[ee[j]])
      child <- pw$edge[ee[j], 2]
      
      tmp <- rep(-Inf, k)
      for (s in 1:k) {
        tmp[s] <- log_sum_exp(log(P[s, ]) + L[child, ])
      }
      logPP[j, ] <- tmp
    }
    
    L[nn[i], ] <- colSums(logPP)
  }
  
  pi <- rep(1/k, k)
  root <- nn[length(nn)]
  prob <- log_sum_exp(log(pi) + L[root, ])
  
  return(prob)
}

objective_function_tpm <- function(par) {
  mu <- 10^(par[1])
  p <- par[2]
  m <- par[3]
  total <- 0
  tree_liks <- numeric(length(trees))
  for (i in seq_along(trees)) {
    llik <- tryCatch({
      pruning_tpm(mu, trees[[i]], str_counts[[i]], p, m, Ex)
    }, error = function(e) {
      message(sprintf("Failed at mu=%.2e p=%.4f m=%.4f (tree %d): %s", mu, p, m, i, conditionMessage(e)))
      return(-1e6)
    })
    if (!is.finite(llik)) llik <- -1e6
    tree_liks[i] <- -llik
    total <- total + (-llik)
  }
  tested_params[[length(tested_params) + 1]] <<- c(mu, p, m, tree_liks, composite = total)
  return(total)
}

cycle_tpm_local <- function(runn){
  overall_best_composite_value <- Inf
  overall_best_tree_liks <- rep(Inf, length(trees))
  
  overall_best_tree_params <- matrix(NA, nrow = length(trees), ncol = 3)
  colnames(overall_best_tree_params) <- c("mu", "p", "m")
  overall_best_tree_params <- data.frame(
    mu = rep(NA, length(trees)),
    p = rep(NA, length(trees)),
    m = rep(NA, length(trees)),
    likelihood = rep(Inf, length(trees))
  )
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  
  for (rep in 1:runn) {
    mu <- 10^runif(1, log10(mu_lower), log10(mu_upper))
    p  <- runif(1, min = p_lower, max = p_upper)
    m  <- runif(1, min = m_lower, max = m_upper)
    
    x0 <- c(mu, p, m)
    print(x0)
    
    #tested_params <- list()
    tested_params <<- list()
    
    result <- optim(
      par = x0,
      fn = objective_function_tpm,
      method = "L-BFGS-B",
      lower = c(log10(mu_lower), p_lower, m_lower),
      upper = c(log10(mu_upper), p_upper, m_upper),
      control = list(trace = 1, ndeps = c(1e-2, 1e-3, 1e-3))
    )
    
    tested_df <- as.data.frame(do.call(rbind, tested_params))
    colnames(tested_df) <- c("mu", "p", "m", paste0("tree", seq_along(trees)), "composite")
    
    for (j in seq_along(trees)) {
      col <- paste0("tree", j)
      idx <- which.min(tested_df[[col]])
      if (tested_df[[col]][idx] < overall_best_tree_liks[j]) {
        overall_best_tree_liks[j] <- tested_df[[col]][idx]
        overall_best_tree_params[j, ] <- c(
          mu = tested_df$mu[idx],
          p = tested_df$p[idx],
          m = tested_df$m[idx],
          likelihood = -tested_df[[col]][idx]
        )
      }
    }
    
    composite_mu <- tested_df$mu[which.min(tested_df$composite)]
    composite_p <- tested_df$p[which.min(tested_df$composite)]
    composite_m <- tested_df$m[which.min(tested_df$composite)]
    
    if (result$value < overall_best_composite_value) {
      overall_best_composite_value <- result$value
      overall_best_composite <- c(mu = composite_mu, p = composite_p, m = composite_m)
    }
  }
  tree_df <- cbind(
    tree = rownames(overall_best_tree_params),
    overall_best_tree_params
  )
  
  composite_row <- data.frame(
    tree = "composite",
    mu = overall_best_composite["mu"],
    p = overall_best_composite["p"],
    m = overall_best_composite["m"],
    likelihood = -overall_best_composite_value
  )
  
  final_df <- rbind(tree_df, composite_row)
  
  return(final_df)
}

cycle_tpm_global <- function(runn){
  overall_best_composite_value <- Inf
  overall_best_tree_liks <- rep(Inf, length(trees))
  
  overall_best_tree_params <- matrix(NA, nrow = length(trees), ncol = 3)
  colnames(overall_best_tree_params) <- c("mu", "p", "m")
  overall_best_tree_params <- data.frame(
    mu = rep(NA, length(trees)),
    p = rep(NA, length(trees)),
    m = rep(NA, length(trees)),
    likelihood = rep(Inf, length(trees))
  )
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  
  for (rep in 1:runn) {
    mu <- 10^runif(1, log10(mu_lower), log10(mu_upper))
    p  <- runif(1, min = p_lower, max = p_upper)
    m  <- runif(1, min = m_lower, max = m_upper)
    
    x0 <- c(mu, p, m)
    print(x0)
    
    #tested_params <- list()
    tested_params <<- list()
    
    result <- GenSA(
      par = x0,
      fn = objective_function_tpm,
      lower = c(log10(mu_lower), p_lower, m_lower),
      upper = c(log10(mu_upper), p_upper, m_upper),
      control = list(
        max.call = 500,
        verbose = FALSE
      )
    )
    
    tested_df <- as.data.frame(do.call(rbind, tested_params))
    colnames(tested_df) <- c("mu", "p", "m", paste0("tree", seq_along(trees)), "composite")
    
    for (j in seq_along(trees)) {
      col <- paste0("tree", j)
      idx <- which.min(tested_df[[col]])
      if (tested_df[[col]][idx] < overall_best_tree_liks[j]) {
        overall_best_tree_liks[j] <- tested_df[[col]][idx]
        overall_best_tree_params[j, ] <- c(
          mu = tested_df$mu[idx],
          p = tested_df$p[idx],
          m = tested_df$m[idx],
          likelihood = -tested_df[[col]][idx]
        )
      }
    }
    
    composite_mu <- tested_df$mu[which.min(tested_df$composite)]
    composite_p <- tested_df$p[which.min(tested_df$composite)]
    composite_m <- tested_df$m[which.min(tested_df$composite)]
    
    if (result$value < overall_best_composite_value) {
      overall_best_composite_value <- result$value
      overall_best_composite <- c(mu = composite_mu, p = composite_p, m = composite_m)
    }
  }
  tree_df <- cbind(
    tree = rownames(overall_best_tree_params),
    overall_best_tree_params
  )
  
  composite_row <- data.frame(
    tree = "composite",
    mu = overall_best_composite["mu"],
    p = overall_best_composite["p"],
    m = overall_best_composite["m"],
    likelihood = -overall_best_composite_value
  )
  
  final_df <- rbind(tree_df, composite_row)
  
  return(final_df)
}

### SMM ###

create_transition_matrix_smm <- function(mu, k) {
  Q <- matrix(0, k, k)
  for (i in 1:k) {
    if (i == 1) {
      Q[i, i] <- -mu
      Q[i, i + 1] <- mu
    } else if (i == k) {
      Q[i, i] <- -mu
      Q[i, i - 1] <- mu
    } else {
      Q[i, i] <- -mu
      Q[i, i - 1] <- mu / 2
      Q[i, i + 1] <- mu / 2
    }
  }
  Q
}

pruning_smm <- function(mu, tree, x, Ex) {
  pw <- reorder(tree, "postorder")
  pw$tip.label <- gsub("n", "", pw$tip.label)
  lo <- min(as.numeric(levels(x)), na.rm = TRUE)
  hi <- max(as.numeric(levels(x)), na.rm = TRUE)
  k0 <- hi - lo + 1
  J <- floor((Ex-1) * k0)
  
  lo <- floor(lo - (J/2))
  if ( lo < 0){
    g <- (-1)*(lo)
    lo = 0
  }
  hi <- floor(hi + (J/2)) + g
  k <- hi - lo + 1
  
  Q <- create_transition_matrix_smm(mu, k)
  
  L <- rbind(to.matrix(x[pw$tip.label], as.character(seq(lo, hi))),
             matrix(0, tree$Nnode, k, dimnames = list(1:tree$Nnode + Ntip(tree))))
  
  L[L == 0] <- -Inf
  L[L == 1] <- 0
  
  nn <- unique(pw$edge[, 1])
  for (node_idx in seq_along(nn)) {
    node <- nn[node_idx]
    ee <- which(pw$edge[, 1] == node)
    logPP <- matrix(NA, length(ee), k)
    
    for (j in seq_along(ee)) {
      edge_id <- ee[j]
      P <- expm(Q * pw$edge.length[edge_id])
      child <- pw$edge[edge_id, 2]
      tmp <- rep(-Inf, k)
      for (s in 1:k) {
        vals <- log(P[s, ]) + L[child, ]
        tmp[s] <- log_sum_exp(vals)
      }
      logPP[j, ] <- tmp
    }
    
    L[node, ] <- colSums(logPP)
  }
  
  pi <- rep(1 / k, k)
  root <- nn[length(nn)]
  loglik_root <- log_sum_exp(log(pi) + L[root, ]) 
  
  return(loglik_root)
}

objective_function_smm <- function(par) {
  mu <- 10^(par[1])
  total_loglik <- 0
  per_tree <- numeric(length(trees))
  
  for (i in seq_along(trees)) {
    result <- tryCatch({
      llik <- pruning_smm(mu, trees[[i]], str_counts[[i]], Ex)
      if (!is.finite(llik)) stop("Non-finite likelihood")
      llik
    }, error = function(e) {
      message(sprintf("Failed at mu=%.4e (tree %d): %s", mu, i, conditionMessage(e)))
      return(-1e6)
    })
    if (!is.finite(llik)) llik <- -1e6
    per_tree[i] <- -result
    total_loglik <- total_loglik + -result
  }
  
  tested_params[[length(tested_params) + 1]] <<- c(mu = mu,
                                                   setNames(per_tree, paste0("tree", seq_along(trees))),
                                                   composite = total_loglik)
  return(total_loglik)
}

cycle_smm_local <- function(runn){
  overall_best_composite_value <- Inf
  overall_best_tree_liks <- rep(Inf, length(trees))
  
  overall_best_tree_params <- matrix(NA, nrow = length(trees), ncol = 1)
  colnames(overall_best_tree_params) <- c("mu")
  overall_best_tree_params <- data.frame(
    mu = rep(NA, length(trees)),
    likelihood = rep(Inf, length(trees))
  )
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  
  for (rep in 1:runn) {
    mu <- 10^runif(1, log10(mu_lower), log10(mu_upper))
    
    x0 <- c(mu)
    
    tested_params <<- list()
    #tested_params <- list()
    
    result <- optim(
      par = x0,
      fn = objective_function_smm,
      method = "L-BFGS-B",
      lower = c(log10(mu_lower)),
      upper = c(log10(mu_upper)),
      control = list(trace = 1, ndeps = 1e-2)
    )
    
    tested_df <- as.data.frame(do.call(rbind, tested_params))
    colnames(tested_df) <- c("mu", paste0("tree", seq_along(trees)), "composite")
    
    for (j in seq_along(trees)) {
      col <- paste0("tree", j)
      idx <- which.min(tested_df[[col]])
      if (tested_df[[col]][idx] < overall_best_tree_liks[j]) {
        overall_best_tree_liks[j] <- tested_df[[col]][idx]
        overall_best_tree_params[j, ] <- c(
          mu = tested_df$mu[idx],
          likelihood = -tested_df[[col]][idx]
        )
      }
    }
    
    composite_mu <- tested_df$mu[which.min(tested_df$composite)]
    
    if (result$value < overall_best_composite_value) {
      overall_best_composite_value <- result$value
      overall_best_composite <- c(mu = composite_mu)
    }
  }
  tree_df <- cbind(
    tree = rownames(overall_best_tree_params),
    overall_best_tree_params
  )
  
  composite_row <- data.frame(
    tree = "composite",
    mu = overall_best_composite["mu"],
    likelihood = -overall_best_composite_value
  )
  
  final_df <- rbind(tree_df, composite_row)
  
  return(final_df)
}

cycle_smm_global <- function(runn){
  overall_best_composite_value <- Inf
  overall_best_tree_liks <- rep(Inf, length(trees))
  
  overall_best_tree_params <- matrix(NA, nrow = length(trees), ncol = 1)
  colnames(overall_best_tree_params) <- c("mu")
  overall_best_tree_params <- data.frame(
    mu = rep(NA, length(trees)),
    likelihood = rep(Inf, length(trees))
  )
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  rownames(overall_best_tree_params) <- paste0("tree", seq_along(trees))
  
  for (rep in 1:runn) {
    mu <- 10^runif(1, log10(mu_lower), log10(mu_upper))
    
    x0 <- c(mu)
    
    #tested_params <- list()
    tested_params <<- list()
    
    result <- GenSA(
      par = x0,
      fn = objective_function_smm,
      lower = log10(mu_lower),
      upper = log10(mu_upper),
      control = list(
        max.call = 500,
        verbose = FALSE
      )
    )
    
    tested_df <- as.data.frame(do.call(rbind, tested_params))
    colnames(tested_df) <- c("mu", paste0("tree", seq_along(trees)), "composite")
    
    for (j in seq_along(trees)) {
      col <- paste0("tree", j)
      idx <- which.min(tested_df[[col]])
      if (tested_df[[col]][idx] < overall_best_tree_liks[j]) {
        overall_best_tree_liks[j] <- tested_df[[col]][idx]
        overall_best_tree_params[j, ] <- c(
          mu = tested_df$mu[idx],
          likelihood = -tested_df[[col]][idx]
        )
      }
    }
    
    composite_mu <- tested_df$mu[which.min(tested_df$composite)]
    
    if (result$value < overall_best_composite_value) {
      overall_best_composite_value <- result$value
      overall_best_composite <- c(mu = composite_mu)
    }
  }
  tree_df <- cbind(
    tree = rownames(overall_best_tree_params),
    overall_best_tree_params
  )
  
  composite_row <- data.frame(
    tree = "composite",
    mu = overall_best_composite["mu"],
    likelihood = -overall_best_composite_value
  )
  
  final_df <- rbind(tree_df, composite_row)
  
  return(final_df)
}

if( Model[1] == 1){
  if (Model[2] == 1){
    cj <- cycle_smm_local(runn)
  }
  if (Model[2] == 2) {
    cj <- cycle_smm_global(runn)
  }
}

if( Model[1] == 2){
  if (Model[2] < 2){
    cj <- cycle_tpm_local(runn)
  }
  if (Model[2] == 2) {
    cj <- cycle_tpm_global(runn)
  }
}

if (Model[1] == 0){
  if (Model[2] < 2){
    SMM <- cycle_smm_local(runn)
    SMM$likelihood <- 2 - 2*(SMM$likelihood) 
    TPM <- cycle_tpm_local(runn)
    TPM$likelihood <- 6 - 2*(TPM$likelihood)
    ###
    c <- vector("list", nrow(TPM))
    
    for (i in seq_len(nrow(SMM))) {
      if (SMM$likelihood[i] < TPM$likelihood[i]) {
        r <- SMM[i, ]
        r$p <- 1
        r$m <- 0
        r$source <- "SMM"
      } else {
        r <- TPM[i, ]
        r$source <- "TPM"
      }
      c[[i]] <- r
    }
    
    cj <- do.call(rbind, c)
    head(cj)
  }

  if (Model[2] == 2) {
    SMM <- cycle_smm_global(runn)
    SMM$likelihood <- 2 - 2*(SMM$likelihood) 
    TPM <- cycle_tpm_global(runn)
    TPM$likelihood <- 6 - 2*(TPM$likelihood)
    ###
    c <- vector("list", nrow(TPM))
    
    for (i in seq_len(nrow(SMM))) {
      if (SMM$likelihood[i] < TPM$likelihood[i]) {
        r <- SMM[i, ]
        r$p <- 1
        r$m <- 0
        r$source <- "SMM"
      } else {
        r <- TPM[i, ]
        r$source <- "TPM"
      }
      c[[i]] <- r
    }
    
    cj <- do.call(rbind, c)
    head(cj)
  }
}

write.table(cj,
            file = sprintf("%s/Overall_best_likelihood.tsv",tree_dir),
            sep = "\t", row.names = FALSE, quote = FALSE)