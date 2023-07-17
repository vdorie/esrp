Efficient Sparse Random Projections
-----------------------------------

## Installation

1. Install your platform's R toolchain.
   * For Mac, install [Xcode](https://apps.apple.com/us/app/xcode/id497799835?mt=12).
   * For Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html).
   * For Linux, use your package manager to install gcc or clang (e.g. on [Debian](https://linuxize.com/post/how-to-install-gcc-compiler-on-debian-10/)). It is extremely likely that this is installed already.
2. From within R, install the `remotes` package, i.e. `install.packages("remotes")`.
3. Execute `install_github("vdorie/esrp")`

## Example

Most of the documentation is in the [manual](man/generate_random_projections.Rd) file. An example from that is recreated below:

```R
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
```
