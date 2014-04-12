## 1. R requirement and installation

R (http://www.r-project.org) is a language and environment for statistical computing and graphics. We assume R (version 3.0.3) has been installed in your local machine. It can be installed following quick instructions below for different platforms (Windows, Mac, and Linux).

* Quick link for `Windows` and `Mac`: [Download R for Windows](http://www.stats.bris.ac.uk/R/bin/windows/base/R-3.0.3-win.exe) and [Download R for Mac](http://www.stats.bris.ac.uk/R/bin/macosx/R-latest.pkg).

* Below are `shell command lines in Terminal` (for `Linux`):

Assume you have a `ROOT (sudo)` privilege:
    
    sudo su
    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.0.3.tar.gz
    tar xvfz R-3.0.3.tar.gz
    cd R-3.0.3
    ./configure
    make
    make check
    make install
    R # start R

Assume you do not have a ROOT privilege and want R installation under your home directory (below `/home/hfang` should be replaced with yours):

    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.0.3.tar.gz
    tar xvfz R-3.0.3.tar.gz
    cd R-3.0.3
    ./configure --prefix=/home/hfang/R-3.0.3
    make
    make check
    make install
    /home/hfang/R-3.0.3/bin/R # start R

## 2. Install the packages from CRAN (Bioconductor & R-Forge)

Notes: below are `R command lines (NOT shell command lines in Terminal)`.

To install [stable release version](http://cran.r-project.org/package=dnet), run:

    source("http://bioconductor.org/biocLite.R")
    biocLite(pkgs=c("supraHex","graph","Rgraphviz"))
    install.packages("dnet",repos="http://cran.r-project.org",type="source")

To install [latest development version](http://r-forge.r-project.org/projects/dnet) (`highly recommended` for benefits of latest improvements), run:

    ## First, install from CRAN
    list.pkg <- c("ape","igraph")
    for(pkg in list.pkg){
        if(!require(pkg, character.only=T)){
            install.packages(pkg,repos="http://www.stats.bris.ac.uk/R",dependencies=T)
        }
    }
    
    ## Second, install from Bioconductor
    list.pkg <- c("graph","Rgraphviz","affy","limma")
    source("http://bioconductor.org/biocLite.R")
    for(pkg in list.pkg){
        if(!require(pkg, character.only=T)) biocLite(pkg)
    }
    
    ## Last, install from R-Forge
    list.pkg <- c("supraHex","dnet")
    for(pkg in list.pkg){
        install.packages(pkg,repos="http://R-Forge.R-project.org", type="source")
    }