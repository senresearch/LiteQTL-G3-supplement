
# What are the required packages. 
packages = c("tidyverse","qtl","mice","tictoc","qtl2", "parallel")

# Specify CRAN mirror
local({r <- getOption("repos")
   r["CRAN"] <- "https://mirrors.nics.utk.edu/cran"
   options(repos=r)
})

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
