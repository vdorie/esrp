\name{random_projections}

\alias{generate_random_projections}

\title{Efficient Sparse Random Projections}
\description{
  Generate random projections of sparse input vectors using
  multiple threads, allowing for millions of dimensions of projected
  coordinates on consumer-grade hardware.
}
\usage{
generate_random_projections(
  u,
  k,
  x = NULL,
  seed = NA_integer_,
  d = NA_integer_,
  type = c("normal", "discrete"),
  threads = NA_integer_
)
}
\arguments{
  \item{u}{a list of input indices vectors}
  \item{k}{integer giving the desired number of output dimensions}
  \item{x}{optional list of real vectors of length equal to \code{u} where
  each element must have the same length as \code{u[[i]]}. These are
  the values that are projected, and if absent are assumed to be 1s}
  \item{seed}{an optional integer giving a seed used when creating pseudo
  random number generators for each input dimension}
  \item{d}{an optional integer giving the size of the input space. Since the
  inputs are coded sparsely, it is possible to supply \code{u} without giving
  the size of the space in which those vectors are embedded. However,
  \code{d} is required when using the \code{"discrete"} \code{type}, as it is
  used to assign probabilities to the possible values of each discrete random
  variable.}
  \item{type}{one of \code{"normal"} or \code{"discrete"}, which determines
  the distribution of the random variables. Discrete refers to the distribution
  specified by Li, Hastie, and Church in the references.}
  \item{threads}{an optional integer specifying how many threads to run in
  parallel}
}
\value{
  Returns a list of vectors, one for each input. Each output vector is the
  random projection of the corresponding input into \code{k} dimensions.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\references{
Li, P., Hastie, T.J., and Church, K.W. (2006)
   Very Sparse Random Projections.
   Proceedings of the 12th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining,
   Association for Computing Machinery, New York, NY, 287-296.
}
\examples{
# Example with fixed SSNs
u1 <- c(
  "000-00-0000", "123-56-6789", "555-55-5555", "424-24-2424", "636-55-3226"
)
u2 <- c(
   "000-00-0000", "555-55-5555", "716-76-2323"
)

# turn SSNs to integers by first removing dashes
# and then parsing strings
print(u1 <- as.integer(gsub("-", "", u1)))
print(u2 <- as.integer(gsub("-", "", u2)))

v_12 <- esrp::generate_random_projections(
  u = list(u1, u2),
  k = 7,
  seed = 0
)

u3 <- c(
  "000-00-0000", "555-55-5555", "646-92-6614"
)
u3 <- as.integer(gsub("-", "", u3))

v_13 <- esrp::generate_random_projections(
  u = list(u1, u3),
  k = 7,
  seed = 0
)

# Assert that method is consistent with a fixed seed
stopifnot(all(v_12[[1]] == v_13[[1]]))


# Multithreaded example with random data
set.seed(0)
# Subtract 1 so that numbers start at 0
u1 <- sample(100000L, 1000L, replace = FALSE) - 1L
u2 <- sample(u1, 200L, replace = FALSE)
u2 <- c(u2, sample(100000L, 30L, replace = FALSE))

x <- list(
  rnorm(length(u1)),
  rnorm(length(u2))
)

v_mt <- esrp::generate_random_projections(
  u = list(u1, u2),
  k = 50001,
  x = x,
  threads = 2,
  seed = 0
)

# repeat single threaded
v_st <- esrp::generate_random_projections(
  u = list(u1, u2),
  k = 50001,
  x = x,
  threads = 1,
  seed = 0
)

# Assert that method is consistent with a fixed seed
for (i in 1:2)
  stopifnot(v_mt[[i]] == v_st[[i]])
}

