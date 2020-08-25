#!/bin/bash

# If you have error: guix locale can't be found. Use the following line.
# export GUIX_LOCPATH=$HOME/.guix-profile/lib/locale

# Install R 
guix package -i r
GUIX_PROFILE="/home/$USER/.guix-profile"

# Install R packages. You have to install R packages uisng GUIX, install.packages does not work under this environment. 
# The convention is r-[packagename]. You can search for the package before installing by `guix search r-[packagename]`
# reference1: https://guix.gnu.org/manual/en/html_node/Invoking-guix-environment.html
# reference2: https://guix.gnu.org/manual/en/html_node/Invoking-guix-package.html
# reference3: https://elephly.net/posts/2017-03-24-r-with-guix.html
guix install r-tidyverse
guix install r-mice
guix install r-tictoc 
guix install r-qtl 
guix install r-qtl2 

R_LIB="/home/$USER/.guix-profile/site-library"
