get_projections_seeded <- function(row_start, output_length, indices, u, k, x, options) {
  options$row_start <- row_start
  options$output_length <- output_length
  
  .Call(C_esrp_generate_random_projections, indices, u, k, x, options)
}

get_projections_unseeded <- function(output_length, indices, u, k, x, options) {
  options$output_length <- output_length
  
  .Call(C_esrp_generate_random_projections, indices, u, k, x, options)
}

generate_random_projections <- function(
  u,
  k,
  x = NULL,
  seed = NA_integer_,
  d = NA_integer_,
  type = c("normal", "discrete"),
  threads = NA_integer_
) {
  matched_call <- match.call()

  if (!is.list(u)) u <- list(u)
  type <- match.arg(type)
  if (!is.null(x) && !is.list(x)) x <- list(x)
  num_input_vectors <- length(u)
  
  if (!all(sapply(u, typeof) == "integer"))
    stop("all input vectors must be of integer type")
  indices <- sort(Reduce(union, u))
  if (!is.null(x)) {
    if (length(x) != length(u))
      stop("length of x must equal length of u")

    if (!all(sapply(x, typeof) == "double"))
      stop("all input x must be of double type")
    
    for (i in seq_along(u)) {
      if (length(u[[i]]) != length(x[[i]]))
        stop("x at position ", i, " does not have the same length as u at position ", i)
      o <- order(u[[i]])
      u[[i]] <- u[[i]][o]
      x[[i]] <- x[[i]][o]
    }
  } else {
    if (length(u) == 1L) {
      u <- NULL
    } else {
      u <- lapply(u, sort)
    }
  }
  
  if (type == "discrete" && (is.na(d) || d[1L] < max(lengths(u))))
    stop("for discrete type, 'd' must be a positive integer of length greater than all input vectors")
  
  seed <- as.integer(seed)
  
  threads <- as.integer(threads)
  if (is.na(threads) || threads[1L] < 1L)
    threads <- 1L
  else
    threads <- threads[1L]
  
  options <- list()

  if (!is.na(d))
    options$d <- d

  if (!is.na(seed))
    options$seed <- seed
  
  k <- as.integer(k[1L])
  if (is.na(k) || k < 0L)
    stop("k must be a positive integer")

  threads <- min(threads, k)
  
  if (threads > 1L) {
    cluster <- makeCluster(threads)

    invisible(clusterEvalQ(cluster, library(esrp)))
    
    num_projections_per_thread <- k %/% threads
    off_by_one_index <- k %% threads
    if (off_by_one_index > 0L) {
      num_projections_per_thread <- c(
        rep.int(num_projections_per_thread + 1L, off_by_one_index),
        rep.int(num_projections_per_thread, threads - off_by_one_index)
      )
    } else {
      num_projections_per_thread <- rep.int(num_projections_per_thread, threads)
    }

    if (!is.na(seed)) {
      row_start <- c(0L, cumsum(num_projections_per_thread[-length(num_projections_per_thread)]))
      
      clusterExport(cluster, "get_projections_seeded", asNamespace("esrp"))

      cluster_results <- clusterMap(cluster, get_projections_seeded, row_start = row_start, output_length = num_projections_per_thread, MoreArgs = list(indices = indices, u = u, k = k, x = x, options = options))
    } else {
      clusterExport(cluster, "get_projections_unseeded", asNamespace("esrp"))

      cluster_results <- clusterMap(cluster, get_projections_unseeded, output_length = num_projections_per_thread, MoreArgs = list(indices = indices, u = u, k = k, x = x, options = options))
    }
    stopCluster(cluster)
    
    lapply(seq_len(num_input_vectors), function(i) unlist(lapply(seq_len(threads), function(j) cluster_results[[j]][[i]])))
  } else {
    .Call(C_esrp_generate_random_projections,
      indices,
      u,
      k,
      x,
      options
    )
  }
}

