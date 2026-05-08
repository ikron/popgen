# Population genetics
This repository contains materials for population genetics course I teach at the University of Jyväskylä

Load the scripts into R by running
```
urlfile <- "https://raw.github.com/ikron/popgen/master/scripts/simpopgen.R"
source(urlfile)
```
# Polymorphism data for Neurospora
You can load the Neurospora polymorphism data by running
```
datafile <- "https://raw.github.com/ikron/popgen/master/data/ncrassa.RData"
load(url(datafile))
```
# Data for Darwin's finches
```
finchurl <- "https://raw.github.com/ikron/popgen/master/data/finch.csv"
finch <- read.csv(finchurl, header = T, sep = ",", dec = ".")
```

# Data for human heights
```
galtonurl <- "https://raw.github.com/ikron/popgen/master/data/galton1889.csv"
galton <- read.csv(galtonurl, header = T, sep = ",", dec = ".")
```

# Data for association mapping
```
myY  <- read.table("https://raw.github.com/ikron/popgen/master/data/mdp_traits.txt", head = TRUE)
myG <- read.table("https://raw.github.com/ikron/popgen/master/data/mdp_genotype_test.hmp.txt", head = FALSE)
#A simulated trait
simY  <- read.table("https://raw.github.com/ikron/popgen/master/data/mdp_traits_validation.txt", head = TRUE)[,-c(2:4)]
```
## Setting up GAPIT
Load the GAPIT functions by running
```
source("http://zzlab.net/GAPIT/gapit_functions.txt")
```
I tested this in the computer class, and at least for me installing GAPIT worked. GAPIT requires a lot of dependencies, so installing all of them takes a bit of time. 
