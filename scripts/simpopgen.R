#Ilkka Kronholm
#Course about population genetics

#Functions used in population genetic simulations with R


#This function simulates genetic drift
genetic.drift <- function(N, nsim, generations, p0 = 0.5, size.change = F, size.vec = NULL) {
    #Make the results matrix
    res.mat <- matrix(rep(0, (nsim+1)*(generations+1)), ncol = nsim + 1)
    colnames(res.mat) <- c("Generation", paste(rep("sim", nsim), 1:nsim, sep = ""))
    res.mat[1,2:ncol(res.mat)] <- p0 #Initialize starting frequency
    res.mat[,1] <- 0:generations #Number of generations

    if(size.change == FALSE) {
    #Wright-Fisher model of genetic drfit
    for(i in 2:nrow(res.mat)) {
        for(j in 2:ncol(res.mat)) {
            pop <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat[(i-1),j], res.mat[(i-1),j] ))
            res.mat[i,j] <- sum(pop)/(2*N)
        }
    }
    }

    #Simulation with changes in population size
    if(size.change == TRUE) {
        #Check that simulated generations and population size vector match
        if(length(size.vec) != generations) { stop(cat("Length of population size vector and number of generations don't match!\n")) }
        #Wright-Fisher model of genetic drfit
        for(i in 2:nrow(res.mat)) {
            #Update population size
            N <- size.vec[i-1]
            for(j in 2:ncol(res.mat)) {
                pop <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat[(i-1),j], res.mat[(i-1),j] ))
                res.mat[i,j] <- sum(pop)/(2*N)
            }
        }
    } #Done
        
    return(res.mat)
}

#This function simulates selection and drfit
selection.drift <- function(N = 100, nsim, s = 0.1, generations, p0 = 0.2, h = 0.5, infinite.size = F, custom.geno.w = F, geno.w.vec = NULL) {
    #Make the results matrix
    res.mat <- matrix(rep(0, (nsim+1)*(generations+1)), ncol = nsim + 1)
    colnames(res.mat) <- c("Generation", paste(rep("sim", nsim), 1:nsim, sep = ""))
    res.mat[1,2:ncol(res.mat)] <- p0 #Initialize starting frequency
    res.mat[,1] <- 0:generations #Number of generations
    #Genotype fitnesses
    if(custom.geno.w == FALSE) { geno.w <- c(1, 1-h*s, 1-s) }
    if(custom.geno.w == TRUE) { geno.w <- geno.w.vec }

    #Wright-Fisher model of genetic drfit
    for(i in 2:nrow(res.mat)) {
        for(j in 2:ncol(res.mat)) {
            #Effect of selection
            p <- res.mat[(i-1),j] #Alle freq in the previous generation
            w.bar <- (p^2)*geno.w[1] + 2*p*(1-p)*geno.w[2] + ((1-p)^2)*geno.w[3] #Average pop W
            #Change in allele frequency due to selection
            delta.p.sel <- (p*(1-p)*(p*(geno.w[1] - geno.w[2]) + (1-p)*(geno.w[2] - geno.w[3])) ) / w.bar
            if(infinite.size == FALSE) {
            #Effect of drift
            popsample <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat[(i-1),j], res.mat[(i-1),j] ))
            p.tplus1 <- sum(popsample)/(2*N)
            delta.p.drift <- p.tplus1 - res.mat[(i-1),j]

            #Change in allele frequency due to all forces
            res.mat[i,j] <- res.mat[(i-1),j] + delta.p.sel + delta.p.drift
            if(res.mat[i,j] < 0) { res.mat[i,j] <- 0 } #Allele freq cannot be < 0
            if(res.mat[i,j] > 1) { res.mat[i,j] <- 1 } #Allele freq cannot be > 1
            }
            #If population size is infinite, only effect of selection
            if(infinite.size == TRUE) { res.mat[i,j] <- res.mat[(i-1),j] + delta.p.sel }
        }
    }

    #return results
    return(res.mat)
}


#Function to plot the results
plot.allele.freq <- function(data) {
    #Reformatting the data frame for plotting
    data <- data.frame(data)
    data <- pivot_longer(data, cols = colnames(data)[-1], names_to = "Simulation", values_to = "Frequency")

    #Making the plot
    ggplot(data, aes(x = Generation, y = Frequency, group = Simulation)) +
        geom_line(col = "red") +
        xlab("Generation") +
        ylab("p") +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        scale_x_continuous(limits = c(0, max(data$Generation)), expand = c(0,0)) +
        theme(plot.margin = unit(c(15.5, 15.5, 5.5, 5.5), "point"), text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

}

#This function simulates selection and drfit
population.structure <- function(N = 100, m = 0.1, npop, generations, p0 = 0.5) {
    #Make the results matrix
    res.mat <- matrix(rep(0, (npop+1)*(generations+1)), ncol = npop + 1)
    colnames(res.mat) <- c("Generation", paste(rep("pop", npop), 1:npop, sep = ""))
    res.mat[1,2:ncol(res.mat)] <- p0 #Initialize starting frequency
    res.mat[,1] <- 0:generations #Number of generations
    #Genotype fitnesses
    #if(custom.geno.w == FALSE) { geno.w <- c(1, 1-h*s, 1-s) }
    #if(custom.geno.w == TRUE) { geno.w <- geno.w.vec }

    #Wright-Fisher model of genetic drfit
    for(i in 2:nrow(res.mat)) {
        #for(j in 2:ncol(res.mat)) {
            #Effect of selection
            #p <- res.mat[(i-1),j] #Alle freq in the previous generation
            #w.bar <- (p^2)*geno.w[1] + 2*p*(1-p)*geno.w[2] + ((1-p)^2)*geno.w[3] #Average pop W
            #Change in allele frequency due to selection
            #delta.p.sel <- (p*(1-p)*(p*(geno.w[1] - geno.w[2]) + (1-p)*(geno.w[2] - geno.w[3])) ) / w.bar
            #if(infinite.size == FALSE) {
            #Effect of drift
            p.tplus1 <- rep(0, npop) #Initialize effect of drift
            for(j in 1:npop) {
                popsample <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat[(i-1),(j+1)], res.mat[(i-1),(j+1)]))
                p.tplus1[j] <- sum(popsample)/(2*N)
            }
            delta.p.drift <- p.tplus1 - res.mat[(i-1),2:ncol(res.mat)]
            #Effect of migration
            mean.p <- mean(res.mat[(i-1),2:ncol(res.mat)])
            delta.p.mig <- m*(mean.p-res.mat[(i-1),2:ncol(res.mat)])

            #Change in allele frequency due to all forces
            res.mat[i,2:ncol(res.mat)] <- res.mat[(i-1),2:ncol(res.mat)] + delta.p.mig + delta.p.drift
            res.mat[i,2:ncol(res.mat)][res.mat[i,2:ncol(res.mat)] < 0] <- 0
            res.mat[i,2:ncol(res.mat)][res.mat[i,2:ncol(res.mat)] > 1] <- 1
            #if(res.mat[i,j] < 0) { res.mat[i,j] <- 0 } #Allele freq cannot be < 0
            #if(res.mat[i,j] > 1) { res.mat[i,j] <- 1 } #Allele freq cannot be > 1
            
            #If population size is infinite, only effect of selection
            #if(infinite.size == TRUE) { res.mat[i,j] <- res.mat[(i-1),j] + delta.p.sel }
        
    }

    #return results
    return(res.mat)
}

#plot.alleles.pop <- 

#Function to plot the results
plot.allele.pop <- function(data) {
    #Reformatting the data frame for plotting
    data <- data.frame(data)
    data <- pivot_longer(data, cols = colnames(data)[-1], names_to = "Population", values_to = "Frequency")

    #Making the plot
     ggplot(data, aes(x = Generation, y = Frequency, group = Population, color = Population)) +
        geom_line() +
        xlab("Generation") +
        ylab("p") +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        scale_x_continuous(limits = c(0, max(data$Generation)), expand = c(0,0)) +
        theme(plot.margin = unit(c(15.5, 15.5, 5.5, 5.5), "point"), text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

}

#Function to calculate summary statistics from structures populations
pop.summary.stat <- function(data) {
    npop <- ncol(data)-1 #Number of populations
    ngen <- nrow(data) #Number of generations
    res.mat <- cbind(data[,1], matrix(rep(0, ngen*3), ncol = 3))
    colnames(res.mat) <- c("Generation", "HS", "HT", "FST") #column names

    het <- function(p) {2*p*(1-p)} #Calculate H-W heterozygosity
    pophet <- apply(data[,-1], MARGIN = 2, het) #Heterozygosity for each pop
    Hs <- apply(pophet, MARGIN = 1, mean) #Hs
    pbar <- apply(data[,-1], MARGIN = 1, mean) #average allele frequency
    Ht <- het(pbar) #Ht

    Fst <- (Ht - Hs) / Ht #Calculate Fst
    Fst[Fst < 0] <- 0
    
    res.mat[,2] <- Hs
    res.mat[,3] <- Ht
    res.mat[,4] <- Fst

    return(res.mat)
}
        
plot.summary.stat <- function(summary) {
    summary <- data.frame(summary) #Changing to a data.frame
    #Making the plot
    ggplot(summary, aes(x = Generation, y = FST)) +
        geom_line() +
        xlab("Generation") +
        ylab("F_ST") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(limits = c(0, max(summary$Generation)), expand = c(0,0)) +
        theme(plot.margin = unit(c(15.5, 15.5, 5.5, 5.5), "point"), text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
}


#Function to simulate local adaptation
#This function simulates selection and drfit when there are two populations
local.adaptation <- function(N = 100, m = 0.1, s, generations, mu = 1E-6, n.p0 = 0.5, s.p0 = 0.999) {
    npop <- 2 #Two populations
    
    #Make the results matrix
    #For the selected locus
    res.mat.s <- matrix(rep(0, (npop+3)*(generations+1)), ncol = npop + 3)
    colnames(res.mat.s) <- c("Generation", paste(rep("s.pop", npop), 1:npop, sep = ""), "w.bar1", "w.bar2")
    res.mat.s[1,2:(npop+1)] <- s.p0 #Initialize starting frequency
    res.mat.s[,1] <- 0:generations #Number of generations
    #Genotype fitnesses
    geno.w.pop1 <- c(1, 1-0.5*s, 1-s)
    geno.w.pop2 <- c(1-s, 1-0.5*s, 1)

    #For the neutral locus
    res.mat.n <- matrix(rep(0, (npop+1)*(generations+1)), ncol = npop + 1)
    colnames(res.mat.n) <- c("Generation", paste(rep("n.pop", npop), 1:npop, sep = ""))
    res.mat.n[1,2:ncol(res.mat.n)] <- n.p0 #Initialize starting frequency
    res.mat.n[,1] <- 0:generations #Number of generations
    

    #Wright-Fisher model of genetic drfit
    for(i in 2:(generations+1)) {
            #Effect of selection
            #Population 1
            p1 <- res.mat.s[(i-1),2] #Allele freq in the previous generation
            w.bar1 <- (p1^2)*geno.w.pop1[1] + 2*p1*(1-p1)*geno.w.pop1[2] + ((1-p1)^2)*geno.w.pop1[3] #Average pop W
            res.mat.s[(i-1),4] <- w.bar1
            #Change in allele frequency due to selection
            delta.p.sel1 <- (p1*(1-p1)*(p1*(geno.w.pop1[1] - geno.w.pop1[2]) + (1-p1)*(geno.w.pop1[2] - geno.w.pop1[3])) ) / w.bar1

            #Population 2
            p2 <- res.mat.s[(i-1),3]
            w.bar2 <- (p2^2)*geno.w.pop2[1] + 2*p2*(1-p2)*geno.w.pop2[2] + ((1-p2)^2)*geno.w.pop2[3] #Average pop W
            res.mat.s[(i-1),5] <- w.bar2
            #Change in allele frequency due to selection
            delta.p.sel2 <- (p2*(1-p2)*(p2*(geno.w.pop2[1] - geno.w.pop2[2]) + (1-p2)*(geno.w.pop2[2] - geno.w.pop2[3])) ) / w.bar2

            
            #Effect of drift
            p.tplus1 <- rep(0, npop) #Initialize effect of drift (selected locus)
            p.tplus1.n <- rep(0, npop) #Initialize the effect of drift (neutral locus)
            for(j in 1:npop) {
                popsample <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat.s[(i-1),(j+1)], res.mat.s[(i-1),(j+1)]))
                p.tplus1[j] <- sum(popsample)/(2*N)
                popsample <- sample(c(0,1), size = 2*N, replace = TRUE, prob = c(1-res.mat.n[(i-1),(j+1)], res.mat.n[(i-1),(j+1)]))
                p.tplus1.n[j] <- sum(popsample)/(2*N)
            }
            delta.p.drift.s <- p.tplus1 - res.mat.s[(i-1),2:(npop+1)]
            delta.p.drift.n <- p.tplus1.n - res.mat.n[(i-1),2:(npop+1)]
            #Effect of migration
            mean.p.s <- mean(res.mat.s[(i-1),2:(npop+1)])
            mean.p.n <- mean(res.mat.n[(i-1),2:(npop+1)])
            delta.p.mig.s <- m*(mean.p.s-res.mat.s[(i-1),2:(npop+1)])
            delta.p.mig.n <- m*(mean.p.n-res.mat.n[(i-1),2:(npop+1)])
            #Effect of mutations
            #−μp t−1 + ν(1 − p t−1 )
            delta.p.mut.s <- -mu*res.mat.s[(i-1),2:(npop+1)] + mu*(1-res.mat.s[(i-1),2:(npop+1)])
            delta.p.mut.n <- -mu*res.mat.n[(i-1),2:(npop+1)] + mu*(1-res.mat.n[(i-1),2:(npop+1)])

            #Change in allele frequency due to all forces
            #Selected locus
            res.mat.s[i,2] <- res.mat.s[(i-1),2] + delta.p.mig.s[1] + delta.p.drift.s[1] + delta.p.mut.s[1] + delta.p.sel1
            res.mat.s[i,3] <- res.mat.s[(i-1),3] + delta.p.mig.s[2] + delta.p.drift.s[2] + delta.p.mut.s[2] + delta.p.sel2
            res.mat.s[i,2:(npop+1)][res.mat.s[i,2:(npop+1)] < 0] <- 0
            res.mat.s[i,2:(npop+1)][res.mat.s[i,2:(npop+1)] > 1] <- 1
            #Neutral locus
            res.mat.n[i,2:(npop+1)] <- res.mat.n[(i-1),2:(npop+1)] + delta.p.mig.n + delta.p.drift.n + delta.p.mut.n
            res.mat.n[i,2:(npop+1)][res.mat.n[i,2:(npop+1)] < 0] <- 0
            res.mat.n[i,2:(npop+1)][res.mat.n[i,2:(npop+1)] > 1] <- 1
            #res.mat.s[i,2] <- 
            
            #if(res.mat[i,j] < 0) { res.mat[i,j] <- 0 } #Allele freq cannot be < 0
            #if(res.mat[i,j] > 1) { res.mat[i,j] <- 1 } #Allele freq cannot be > 1
            
            #If population size is infinite, only effect of selection
            #if(infinite.size == TRUE) { res.mat[i,j] <- res.mat[(i-1),j] + delta.p.sel }
        
    }

#Average fitness for the final generation
p1 <- res.mat.s[(generations-1),2]
res.mat.s[(generations+1),4] <- (p1^2)*geno.w.pop1[1] + 2*p1*(1-p1)*geno.w.pop1[2] + ((1-p1)^2)*geno.w.pop1[3] #Average pop W
#Population 2
p2 <- res.mat.s[(generations-1),3]
res.mat.s[(generations+1),5] <- (p2^2)*geno.w.pop2[1] + 2*p2*(1-p2)*geno.w.pop2[2] + ((1-p2)^2)*geno.w.pop2[3]

    res.mat <- cbind(res.mat.s, res.mat.n[,2:3]) 

    #return results
    return(res.mat)
}

#Function to plot the results
plot.allele.local <- function(data) {
    #Reformatting the data frame for plotting
    data <- data.frame(data)
    data <- pivot_longer(data, cols = colnames(data)[-c(1,4,5)], names_to = "Population", values_to = "Frequency")
    data$type <- ifelse(grepl("s", data$Population), "selected", "neutral")
    data$pop <- ifelse(grepl("1", data$Population), "Population 1", "Population 2")

    #Making the plot
     ggplot(data, aes(x = Generation, y = Frequency, group = Population, color = type)) +
        geom_line() +
        xlab("Generation") +
        ylab("p") +
        scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
        scale_x_continuous(limits = c(0, max(data$Generation)), expand = c(0,0)) +
        facet_grid(pop ~ .) +
        theme(plot.margin = unit(c(15.5, 15.5, 5.5, 5.5), "point"), panel.spacing = unit(12, "point"), text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

}

#This function calculates summary statistics 
local.adapt.summary.stat <- function(data) {
    ngen <- nrow(data)
    res.mat <- cbind(data[,1], matrix(rep(0, ngen*2), ncol = 2))
    colnames(res.mat) <- c("Generation", "FST_sel", "FST_neut")

    #Calculate FST for the selected locus
    het <- function(p) {2*p*(1-p)} #Calculate H-W heterozygosity
    pophet.s <- apply(data[,2:3], MARGIN = 2, het) #Heterozygosity for each pop
    Hs.s <- apply(pophet.s, MARGIN = 1, mean) #Hs
    pbar.s <- apply(data[,2:3], MARGIN = 1, mean) #average allele frequency
    Ht.s <- het(pbar.s) #Ht

    Fst.s <- (Ht.s - Hs.s) / Ht.s #Calculate Fst
    Fst.s[Fst.s < 0] <- 0

    res.mat[,2] <- Fst.s

    #Calculate FST for the neutral locus
    pophet.n <- apply(data[,6:7], MARGIN = 2, het) #Heterozygosity for each pop
    Hs.n <- apply(pophet.n, MARGIN = 1, mean) #Hs
    pbar.n <- apply(data[,6:7], MARGIN = 1, mean) #average allele frequency
    Ht.n <- het(pbar.n) #Ht

    Fst.n <- (Ht.n - Hs.n) / Ht.n #Calculate Fst
    Fst.n[Fst.n < 0] <- 0

    res.mat[,3] <- Fst.n
    return(res.mat)

}

plot.local.adapt.fst <- function(summary) {
    summary <- data.frame(summary) #Changing to a data.frame
    summary <- pivot_longer(summary, cols = colnames(summary)[-1], names_to = "type", values_to = "FST")
    #Making the plot
    ggplot(summary, aes(x = Generation, y = FST, group = type, color = type)) +
        geom_line() +
        xlab("Generation") +
        ylab("F_ST") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(limits = c(0, max(summary$Generation)), expand = c(0,0)) +
        theme(plot.margin = unit(c(15.5, 15.5, 5.5, 5.5), "point"), text = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
}

