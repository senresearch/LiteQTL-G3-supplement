local({r <- getOption("repos")
   r["CRAN"] <- "https://mirrors.nics.utk.edu/cran"
   options(repos=r)
})

install.packages("tidyverse")
install.packages("qtl")
install.packages("tictoc")
install.packages("qtl2")
