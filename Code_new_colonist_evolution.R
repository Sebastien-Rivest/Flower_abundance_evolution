# R Code for the simulation of pollination system evolution in plant and pollinator communities, as reported in:
# "Interspecific variation in plant abundance can contribute to the diversity and evolutionary assembly of floral communities".

# Packages used in the model
library(sads)
require(bipartite)
library(progress)
library(ggplot2)
library(ggthemes)
library(magrittr)
library(ggpubr)
library(extrafont)
library(maxnodf)
library(truncnorm)

# Settings for the set of simulations
plantsp.final <- seq(1,12,by = 1) # List of plant species
pollisp <- seq(1,10,by = 1) # List of pollinator species
num.sims = 100 # Number of simulations
num.gen = 1000 # Number of generations
num.gen.new <-200

# Input parameters:
colonization.rate = 25 # Rate of colonization by new plant species
av.plant.abund = 250 # Average plant abundance in the communities
visits.per.reward = 5 # Number of pollinator visits that one unit of reward sustains
K.pollinators = 1000 # Average value of Ki for the pollinators (limitation by factors other than reward availability)
av.Pj = 50000 # Average number of pollen grains produced per flower
av.Oj = 10 # Average number of ovules produced per flower
av.tj = 0.01 # Average proportion of conspecific pollen on pollinator bodies transported to conspecific stigmas
av.bj = 0.5 # Average proportion of received pollen grains resulting in fertilized seeds
av.qrew = 1 # Average reward production rate

# Creating empty objects in which the data will go.
plantpopsizes.det <- array(0,dim=c(length(plantsp.final),num.sims)) # Maximal population sizes of the plant species
pollpopsizes <- array(0,dim=c(length(pollisp),num.gen+1,num.sims)) # Population sizes of the pollinators
pollpopsizes.max <- array(0,dim=c(length(pollisp),num.sims)) # Values of Ki for the pollinators
pollicarry <- array(0,dim=c(length(pollisp),num.sims)) # Carryover capacity of the pollinators
polliremo <- array(0,dim=c(length(pollisp),num.sims)) # Pollen removal rate of the pollinators
poll.dep <- array(0,dim=c(length(plantsp.final),length(pollisp),num.gen+1,num.sims)) # Male fitness
femalfit <- array(0,dim=c(length(plantsp.final),num.gen+1,num.sims)) # Female fitness
selection.gradient <- array(0,dim=c(length(plantsp.final),length(pollisp))) # Selection gradient resulting from male fitness
selection.gradient.f <- array(0,dim=c(length(plantsp.final),length(pollisp))) # Selection gradient resulting from female fitness
fg <- array(0,dim=c(length(plantsp.final),length(pollisp))) # Total selection gradient
genome <- array(numeric(0),dim=c(length(plantsp.final),length(pollisp),num.gen+1,num.sims)) # Attractiveness values
network <- lapply(1:num.sims, matrix, data= 0, nrow=length(plantsp.final), ncol=length(pollisp)) # Plant-pollinator networks

adapt.for <- array(1,dim=c(length(plantsp.final),num.sims)) # Values of adaptive foraging

# Simulations of the evolutionary assembly of flower communities
for (j in 1:num.sims){ # Loop for each simulation
  #Sampling parameter values for this simulation
  plantpopsizes.det[,j] <- rpoilog(length(plantsp.final), mu = 5, sig =1) %>% +1 %>% # Randomly sampling the maximal plant abundances (kj).
    "/" (mean(.)) %>% "*"(av.plant.abund) %>% ceiling(.)
  plantpopsizes <-1 # initial plant abundance = 1
  plantsp <- c(1) # List of plant species in the community, starting with one species (the list expands as other species colonize the community).
  pollpopsizes[,1,j] <- rpoilog(length(pollisp), mu = 5, sig = 0.7) %>% # Randomly sampling the initial number of visits by each pollinator species.
    "/" (mean(.)) %>% "*"(visits.per.reward) %>% ceiling(.) 
  pollpopsizes.max[,j] <- rpoilog(length(pollisp), mu = 5, sig = 0.7) %>% # Randomly sampling the values of Ki (limitation by factors other than reward availability).
    "/" (mean(.)) %>% "*"(K.pollinators) %>% ceiling(.) 
  pollicarry[,j] <- runif(length(pollisp), 0.55, 0.9) # Sampling carryover values of the different pollinators from a uniform distribution.
  polliremo[,j] <- runif(length(pollisp), 0.3, 0.7) # Sampling removal values from a uniform distribution.
  samply<-runif(100,0.5,0.5) # Initial values of attractiveness are set at 0.5, but could also be sampled from an uniform distribution.
  genome[plantsp,,1,j] <- sample(samply,length(genome[plantsp,,1,j]),replace = TRUE)
  # Sampling values of Pj, Oj, tj, bj, and qrew from normal distribution truncated at zero.
  # The distributions of tj and bj are also truncated at 1 because this is the maximal value that they can take.
  Pj <- rtruncnorm(length(plantsp.final),a=0, b=Inf, av.Pj, av.Pj/4)
  Oj <- rtruncnorm(length(plantsp.final),a=0, b=Inf, av.Oj, av.Oj/4)
  tj <- rtruncnorm(length(plantsp.final),a=0, b=1, av.tj, av.tj/4)
  bj <- rtruncnorm(length(plantsp.final),a=0, b=1, av.bj, av.bj/4)
  qrew <- rtruncnorm(length(plantsp.final),a=0, b=Inf, av.qrew, av.qrew/4)
  
  reward <- array(0,dim=c(length(plantsp.final))) # Values of total reward per plant species
  reward.f <- array(0,dim=c(length(plantsp.final))) # Values of reward per flower for each species
  
  for (g in 1:num.gen){ # Loop for each generation
    # New colonization events
    if (length(plantsp) < length(plantsp.final)){ # Stop potential colonization when all plant species have colonized.
      # New colonization events occur with a probability given by the colonization rate.
      if (sample(1000,1) > 1000-1000/colonization.rate){ # Sampling probability of a new colonization.
        plantsp <- c(plantsp,tail(plantsp,1)+1) # Adding the new colonist.
        plantpopsizes <- c(plantpopsizes,plantpopsizes.det[tail(plantsp,1)]) # Adding abundance of the new colonist.
        genome[tail(plantsp,1),,g,j] <- sample(samply,length(genome[tail(plantsp,1),,g,j]),replace = TRUE) # Adding the initial values of attractiveness.
      }
    }
    # Calculating the number of visits by each pollinator to each plant species to build the plant-pollinator-network.
    for (i in 1:length(plantsp)){
      for (m in 1:length(pollisp)){
        network[[j]][i,m] <- pollpopsizes[m,g,j]*((plantpopsizes[i]*genome[i,m,g,j]*adapt.for[i,j])/sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp,j]))
      }
    }
    # Measuring reward availability for each plant species.
    for (i in 1:length(plantsp)){
      reward.f[i] <- 1/(1+(sum(network[[j]][i,])/plantpopsizes[i])/qrew[i]) # Reward content per flower for each species
      reward[i] <- (1/(1+(sum(network[[j]][i,])/plantpopsizes[i])/qrew[i])) * plantpopsizes[i] # Reward content for the whole plant populations
    }
    # Measuring the relative reward content of each plant species which determines adaptive foraging.
    # The values of adaptive foraging are updated at each generation to adjust to changes in reward availability between species.
    for (i in 1:length(plantsp)){
      adapt.for[i,j] <- adapt.for[i,j] * (reward.f[i]/mean(reward.f[plantsp]))
    }
    
    # Measuring for each plant species male and female fitness and the selection gradients for attractiveness to each pollinator species.
    for(i in seq_along(plantsp)){ # Loop for each plant species
      # Measuring the number of pollen grains removed by each pollinator alone ("info.removed") and the proportion of pollen removed by each pollinator ("H").
      info.removed <- matrix(rep(0,length(pollisp)),nrow=1)
      H <- matrix(rep(0,length(pollisp)),nrow=1)
      for (m in 1:length(pollisp)) {
        info.removed[1,m] <- (1-polliremo[m,j]) ^ ((pollpopsizes[m,g,j] * genome[i,m,g,j] *adapt.for[i,j]) / sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp,j]))
        H[1,m] <-  polliremo[m,j] * pollpopsizes[m,g,j] * genome[i,m,g,j] * adapt.for[i,j] / sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp,j])
      }
      # Measuring, for each pollinator, its contribution to male fitness and the selection gradient for attractiveness to this pollinator.
      for(t in seq_along(pollisp)){ # Loop for each pollinator species
        sel.g <- numeric(0) # vector of the effect of variation in attractiveness to a given pollinator on pollination by each pollinator
        # Calculating the contribution of each pollinator to the focal plant's male fitness.
        poll.dep[i,t,g,j] <- (Pj[i]*tj[i]*(1-prod(info.removed)) * (H[1,t]/sum(H)))  * pollicarry[t,j]^
          (sum(plantpopsizes*genome[plantsp,t,g,j]*adapt.for[plantsp,j])/(plantpopsizes[i]*genome[i,t,g,j]*adapt.for[i,j]))
        # Measuring the selection gradient resulting from male fitness using the partial derivative of male fitness as a function of attractiveness to a given pollinator.
        # The partial derivative of equation (2) (see article) in respect to attractiveness to a given pollinator is the sum of... 
        # the effect of attractiveness to a given pollinator on pollination by this pollinator (first component),
        # and its effect on pollination by the other pollinators (second component).
        for(m in seq_along(pollisp)){
          if (m == t){ # The effect of attractiveness to a given pollinator on pollination by this pollinator is calculated first. 
            # In the loop, this corresponds to when pollinator "m" = pollinator "t" (m == t).
            # Input parameters for Equation (2) first component, see Input parameters above.
            Z <- prod(info.removed)
            r <- polliremo[m,j]
            v <- pollpopsizes[m,g,j]
            x <- genome[i,m,g,j]
            y <- genome[i,m,g,j] # Average attractiveness, which is equal to the individual attractiveness
            a <- adapt.for[i,j]
            A <- sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp,j])
            G <- sum(H)
            B <- plantpopsizes[i]
            C <- pollicarry[m,j]
            # Partial derivative of Equation (2) first component
            sel.g <- c(sel.g, ((Pj[i] * tj[i] * (1 - Z)) * (v * r * a/A/G - (v * x * r * a/A) * (v * r * a/(A))/G^2)
                               - Pj[i] * tj[i] * (Z * (log(1 - r) * (v * a/A))) * ((v * x * r * a/A)/G)) * C^(A /(B * y * a)))
          }
          else { # The effect of attractiveness to a given pollinator on pollination by the other pollinators
            # Input parameters for Equation (2) second component, see Input parameters above.
            Z <- prod(info.removed)
            r <- polliremo[m,j]
            v <- pollpopsizes[m,g,j]
            y <- genome[i,m,g,j]
            A <- sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp,j])
            G <- sum(H)
            C <- pollicarry[m,j]
            B <- plantpopsizes[i]
            D <- sum(plantpopsizes*genome[plantsp,t,g,j]*adapt.for[plantsp,j])
            x <- genome[i,t,g,j]
            a <- adapt.for[i,j]
            o <- polliremo[t,j]
            q <- pollpopsizes[t,g,j]
            # Partial derivative of Equation (2) second component
            sel.g <- c(sel.g, -((((Pj[i] * tj[i] * (1 - Z)) * ((v * y * a/(A)) * r) * (o * q * a/D)/(G^2) +
                                    Pj[i] * tj[i] * (Z * (log((1 - o)) * (q * a/D))) * (((v * y * a/(A)) * r)/G)) * C^(A/(B * y * a)))))
          }
        }
        selection.gradient[i,t] <- sum(sel.g) # The selection gradient of attractiveness to a given pollinator is the sum of the effect of variation
        # of this trait on pollination by each pollinator.
      }
      # Calculating female fitness.
      
      femalfit[i,g,j] <- Oj[i] * (1-exp(-bj[i]*sum(poll.dep[i,,g,j])/Oj[i]))
      # Measuring the selection gradient resulting from female fitness.
      for(t in seq_along(pollisp)){
        # Input parameters for Equation (3), see Input parameters above.
        z <- sum(poll.dep[i,,g,j])-poll.dep[i,t,g,j]
        M <- poll.dep[i,t,g,j]
        x <- genome[i,t,g,j]
        a <- genome[i,t,g,j] # Average attractiveness, which is equal to the individual attractiveness
        # Partial derivative of Equation (3)
        selection.gradient.f[i,t] <- Oj[i] * (exp(-bj[i] * (z + (M * x/a))/Oj[i]) * (bj[i] * (M/a)/Oj[i]))
      }
    }
    # Measuring changes in attractiveness to each pollinator for each plant species in this generation.
    for (i in 1:length(plantsp)){
      for (m in 1:length(pollisp)){
        if (femalfit[i,g,j] == 0) { # If female fitness is zero (i.e no pollination occurred) then selection on attractiveness is zero (avoid division by zero).
          fg[i,m] <- 0
        }
        else { # Otherwise measure selection gradient
          # The selection gradient is divided by the average total fitness (2*femalfit[i,g,j]) to obtain relative rather than absolute fitness.
          fg[i,m] <- selection.gradient[i,m]/sum(poll.dep[i,,g,j]) + selection.gradient.f[i,m]/femalfit[i,g,j]
        }
        # Change in attractiveness
        # The additive genetic variance is a non-negative function of attractiveness that approaches zero when the mean attractiveness value is close to its minimum or maximum.
        genome[i,m,g+1,j] <- genome[i,m,g,j] + (genome[i,m,g,j]*(1-genome[i,m,g,j]) * fg[i,m]/2) 
        # Preventing attractiveness from reaching fixation (i.e. reach 0 or 1) (simulating the input of genetic variability).
        if (genome[i,m,g+1,j] > 0.99) { 
          genome[i,m,g+1,j] <- 0.99
        }
        if (genome[i,m,g+1,j] < 0.01) {
          genome[i,m,g+1,j] <- 0.01
        }
      } 
    }

    # Measuring pollinator populations growth.
    for (m in 1:length(pollisp)){
      visprop <- rep(0,length(plantsp))
      for (i in 1:length(plantsp)){
        visprop[i] <- network[[j]][i,m]/sum(network[[j]][i,]) # Proportion of visits made by each pollinator to each plant species
        # This measure is used to determine the quantity of reward that each pollinator can acquire.
      }
      # Pollinator population sizes
      pollpopsizes[m,g+1,j] <- pollpopsizes.max[m,j]*(1-exp(-((visits.per.reward+visits.per.reward^2)*sum(reward[plantsp]*visprop)/pollpopsizes.max[m,j])))
      if (pollpopsizes[m,g+1,j] < 1){ # Preventing pollinator populations from going completely extinct.
        pollpopsizes[m,g+1,j] <- 1 # Pollinator populations cannot be lower than 1, simulating low but constant migration.
      }
    }
    # Plant population sizes (constant in the basic model, but can be modeled to change)
    for (i in plantsp){
      plantpopsizes[i] <- plantpopsizes.det[i,j]
    }
    
  }
} # End of simulation

# Creating empty objects in which the data will go for the new colonist's evolution
ar.competition <- array(0,dim=c(300,num.sims)) # Average level of competition on the pollinators on which the new colonist evolved.
ar.abundance <- array(0,dim=c(300,num.sims)) # Average total number of visits made by each pollinator on which the new colonist evolved
ar.carryover <- array(0,dim=c(300,num.sims)) # Average pollen carryover of the pollinators on which the new colonist evolved
ar.removal <- array(0,dim=c(300,num.sims)) # Average removal rate of the pollinators on which the new colonist evolved
ar.generalization <- array(0,dim=c(300,num.sims)) # Degree of specialization of the new colonist
# New objects used during the simulations including the new colonist
plantpopsizes.new<-rbind(plantpopsizes.det,rep(0,num.sims))
adapt.for.new<-rbind(adapt.for,rep(1,num.sims))
genome.new <- array(numeric(0),dim=c(length(plantsp.final)+1,length(pollisp),num.gen+1,num.sims))
pollpopsizes.new <- array(0,dim=c(length(pollisp),num.gen.new+1,num.sims))
reward.new <- array(0,dim=c(length(plantsp)+1))
reward.f.new <- array(0,dim=c(length(plantsp)+1))
network.new <-  lapply(1:num.sims, matrix, data= 0, nrow=length(plantsp.final)+1, ncol=length(pollisp)) # Plant-pollinator networks
poll.dep <- array(0,dim=c(length(pollisp),num.gen+1,num.sims)) # Male fitness
femalfit <- array(0,dim=c(num.gen+1,num.sims)) 
selection.gradient <- array(0,dim=c(length(pollisp))) 
selection.gradient.f <- array(0,dim=c(length(pollisp))) 
fg <- array(0,dim=c(length(pollisp))) 
qrew <-c(qrew,1)

Pj = 50000 # Number of pollen grains produced per flower
Oj = 10 # Number of ovules produced per flower
tj = 0.01 # Proportion of conspecific pollen on pollinator bodies transported to conspecific stigmas
bj = 0.5 # Proportion of received pollen grains resulting in fertilized seeds

# Simulations of the evolution of a new plant colonist in function of its abundance
for (j in 1:num.sims){ # Loop of simulation
  for (incr in 1:300){ # The simulation is repeated from abundances of 10 to 2000 with increments of 10 between the sets of 100 simulations.
    pollpopsizes.new[,1,j] <- pollpopsizes[,num.gen,j]
    adapt.for.new[,j] <- c(adapt.for[,j],1) # Adding the adaptive foraging for the new colonist which start at one (pollinators have not visited the flower yet).
    i <- length(plantsp)+1 # identity of the new colonist
    genome.new[,,1,j] <- rbind(genome[,,num.gen,j],rep(0.5,length(pollisp))) # Adding the starting values of attractiveness of the new colonist.
    plantpopsizes.det.new <- incr*10 # Increments of 10
    plantpopsizes.new[i,j] <- plantpopsizes.det.new
    
    for (g in 1:num.gen.new){ # Loop for each generation
      # Flower evolution is measured only for the new colonist, using the same process as for the evolutionary assembly of the flower community above.
      for (ij in 1:length(plantpopsizes.new[,j])){
        for (m in 1:length(pollisp)){
          network.new[[j]][ij,m] <- pollpopsizes.new[m,g,j]* ((plantpopsizes.new[ij,j]*genome.new[ij,m,g,j]*adapt.for.new[ij,j])/sum(plantpopsizes.new[,j]*genome.new[,m,g,j]*adapt.for.new[,j]))
        }
      }
      for (ij in 1:length(plantpopsizes.new[,j])){
        reward.f.new[ij] <- 1/(1+(sum(network.new[[j]][ij,])/plantpopsizes.new[ij,j])/qrew[ij])
        reward.new[ij] <- (1/(1+(sum(network.new[[j]][ij,])/plantpopsizes.new[ij,j])/qrew[ij]))*plantpopsizes.new[ij,j]
        
      }
      for (ij in 1:length(plantpopsizes.new[,j])){
        adapt.for.new[ij,j] <- adapt.for.new[ij,j] * reward.f.new[ij]/mean(reward.f.new)
      }
      
      info.removed<- matrix(rep(0,length(pollisp)),nrow=1)
      H <- matrix(rep(0,length(pollisp)),nrow=1)
      for (m in 1:length(pollisp)) {
        info.removed[1,m] <- (1-polliremo[m,j]) ^ ((pollpopsizes.new[m,g,j]*genome.new[i,m,g,j]*adapt.for.new[i,j]) / sum(plantpopsizes.new[,j]*genome.new[,m,g,j]*adapt.for.new[,j]))
        H[1,m] <-  polliremo[m,j] * pollpopsizes.new[m,g,j] * genome.new[i,m,g,j] * adapt.for.new[i,j] / sum(plantpopsizes.new[,j]*genome.new[,m,g,j]*adapt.for.new[,j])
      }
      for(t in seq_along(pollisp)){ # Loop for each pollinator species
        sel.g <- numeric(0) # vector of the effect of variation in attractiveness to a given pollinator on pollination by each pollinator
        poll.dep[t,g,j] <- (Pj*tj*(1-prod(info.removed)) * (H[1,t]/sum(H)))  * pollicarry[t,j]^
          (sum(plantpopsizes.new[,j]*genome.new[,t,g,j]*adapt.for.new[,j])/(plantpopsizes.new[i,j]*genome.new[i,t,g,j]*adapt.for.new[i,j]))
        for(m in seq_along(pollisp)){
          if (m == t){ 
            Z <- prod(info.removed)
            r <- polliremo[m,j]
            v <- pollpopsizes.new[m,g,j]
            x <- genome.new[i,m,g,j]
            y <- genome.new[i,m,g,j] # Average attractiveness, which is equal to the individual attractiveness
            a <- adapt.for.new[i,j]
            A <- sum(plantpopsizes.new[,j]*genome.new[,m,g,j]*adapt.for.new[,j])
            G <- sum(H)
            B <- plantpopsizes.new[i,j]
            C <- pollicarry[m,j]
            sel.g <- c(sel.g, ((Pj * tj * (1 - Z)) * (v * r * a/A/G - (v * x * r * a/A) * (v * r * a/(A))/G^2)
                               - Pj * tj * (Z * (log(1 - r) * (v * a/A))) * ((v * x * r * a/A)/G)) * C^(A/(B * y * a)))
          }
          else { 
            Z <- prod(info.removed)
            r <- polliremo[m,j]
            v <- pollpopsizes.new[m,g,j]
            y <- genome.new[i,m,g,j]
            A <- sum(plantpopsizes.new[,j]*genome.new[,m,g,j]*adapt.for.new[,j])
            G <- sum(H)
            C <- pollicarry[m,j]
            B <- plantpopsizes.new[i,j]
            D <- sum(plantpopsizes.new[,j]*genome.new[,t,g,j]*adapt.for.new[,j])
            x <- genome.new[i,t,g,j]
            a <- adapt.for.new[i,j]
            o <- polliremo[t,j]
            q <- pollpopsizes.new[t,g,j]
            sel.g <- c(sel.g, -((((Pj * tj * (1 - Z)) * ((v * y * a/(A)) * r) * (o * q * a/D)/(G^2) +
                                    Pj * tj * (Z * (log((1 - o)) * (q * a/D))) * (((v * y * a/(A)) * r)/G)) * C^(A/(B * y * a)))))
          }
        }
        selection.gradient[t] <- sum(sel.g) 
      }
      femalfit[g,j] <- Oj * (1-exp(-bj*sum(poll.dep[,g,j])/Oj))
      for(t in seq_along(pollisp)){
        z <- sum(poll.dep[,g,j])-poll.dep[t,g,j]
        M <- poll.dep[t,g,j]
        x <- genome.new[i,t,g,j]
        a <- genome.new[i,t,g,j]
        selection.gradient.f[t] <- Oj * (exp(-bj * (z + (M * x/a))/Oj) * (bj * (M/a)/Oj))
      }
      
      genome.new[,,g+1,j] <- genome.new[,,1,j] # Expect for the new colonist, the attractiveness of the other plant species remains the same.
      
      for (m in 1:length(pollisp)){
        if (femalfit[g,j] == 0) {
          fg[m] <- 0
        }
        else {
          fg[m] <- selection.gradient[m]/sum(poll.dep[,g,j]) + selection.gradient.f[m]/femalfit[g,j]
        }
        genome.new[i,m,g+1,j] <-genome.new[i,m,g,j] + (genome.new[i,m,g,j]*(1-genome.new[i,m,g,j]) * fg[m]/2)
        if (genome.new[i,m,g+1,j] > 0.99) {
          genome.new[i,m,g+1,j] <- 0.99
        }
        if (genome.new[i,m,g+1,j] < 0.01) {
          genome.new[i,m,g+1,j] <- 0.01
        }
      } 

      pollpopsizes.new[,g+1,j] <- pollpopsizes.new[,g,j] # Pollinator populations are assumed to have reached their maximum to avoid confunding efffects.
      plantpopsizes.new[i,j] <- plantpopsizes.det.new
    }
    # Measuring the strength of competition for the exploitation of each pollinator (similar to pollinator specialization).
    # Each pollinator for which a plant's attractiveness was more than 0.2 was considered a pollinator of this species.
    comp <- numeric(0) # Empty object
    for (k in 1:length(which(genome.new[i,,g,j] > 0.2))){ # Length = number of pollinators of the new colonist
      pol <- which(genome.new[i,,g,j] > 0.2)[k] # Object of the pollinators of the new colonist
      # Measuring the total abundance of the other plant species exploiting this pollinator (the measure of competition).
      comp <- c(comp,((sum(plantpopsizes.new[which(genome.new[,pol,g,j] > 0.2),j])) - plantpopsizes.new[i,j]))
    }
    # Stocking the attributes of the pollinators toward which the new colonist evolved.
    ar.competition[incr,j] <- mean(comp) # Average level of competition on the pollinators on which the new colonist evolved.
    ar.abundance[incr,j] <- mean(pollpopsizes.new[which(genome.new[i,,g,j] > 0.2),g,j]) # Average total number of visits made by each pollinator on which the new colonist evolved
    ar.carryover[incr,j] <- mean(pollicarry[which(genome.new[i,,g,j] > 0.2),j]) # Average pollen carryover of the pollinators on which the new colonist evolved
    ar.generalization[incr,j] <- length(which(genome.new[i,,g,j] > 0.2)) # Degree of specialization of the new colonist
    ar.removal[incr,j] <- mean(polliremo[which(genome.new[i,,g,j] > 0.2),j]) # Average removal rate of the pollinators on which the new colonist evolved
  }
}

# Visualization of the pollination systems evolved for the new colonist in function of its floral abundance.
# Create empty objects in which the mean and standard error (SE) between simulations at each value of the new colonist abundance will go.
mean.competition <- array(0,dim=c(300))
se.competition <- array(0,dim=c(300))
mean.abundance <- array(0,dim=c(300))
se.abundance <- array(0,dim=c(300))
mean.carryover <- array(0,dim=c(300))
se.carryover <- array(0,dim=c(300))
mean.removal <- array(0,dim=c(300))
se.removal <- array(0,dim=c(300))
mean.generalization <- array(0,dim=c(300))
se.generalization <- array(0,dim=c(300))

for (a in 1:300) {
  mean.competition[a] <- mean(ar.competition[a,]) # Mean value of pollinator specialization.
  se.competition[a] <- sd(ar.competition[a,])/sqrt(num.sims) # SE of pollinator specialization.
  mean.abundance[a] <- mean(ar.abundance[a,]) # Mean value of pollinator abundance.
  se.abundance[a] <- sd(ar.abundance[a,])/sqrt(num.sims) # SE of pollinator abundance.
  mean.carryover[a] <- mean(ar.carryover[a,]) # Mean value of pollinator carryover.
  se.carryover[a] <- sd(ar.carryover[a,])/sqrt(num.sims) # SE of pollinator carryover.
  mean.removal[a] <- mean(ar.removal[a,]) # Mean value of pollinator carryover.
  se.removal[a] <- sd(ar.removal[a,])/sqrt(num.sims) # SE of pollinator carryover.
  mean.generalization[a] <- mean(ar.generalization[a,]) # Mean value of plant generalization in pollination system.
  se.generalization[a] <- sd(ar.generalization[a,])/sqrt(num.sims) # SE of plant generalization in pollination system.
}

# Standardizing values with maximum and minimum values of 1 and 0 respectively.
mcomp <- 1-(mean.competition-(min(mean.competition)))/(max(mean.competition)-min(mean.competition))
scomp <- se.competition/(max(mean.competition)-min(mean.competition))
mabun <- (mean.abundance-min(mean.abundance))/(max(mean.abundance)-min(mean.abundance))
sabun <- se.abundance/(max(mean.abundance)-min(mean.abundance))
mcarr <- (mean.carryover- min(mean.carryover))/(max(mean.carryover)- min(mean.carryover))
scarr <- se.carryover/(max(mean.carryover)- min(mean.carryover))
mremo <- (mean.removal- min(mean.removal))/(max(mean.removal)- min(mean.removal))
sremo <- se.removal/(max(mean.removal)- min(mean.removal))
mgene <- (mean.generalization- min(mean.generalization))/(max(mean.generalization)- min(mean.generalization))
sgene <- se.generalization/(max(mean.generalization)- min(mean.generalization))

# Creating data frame.
x<- seq(from = 10, to = 3000, by = 10)
parameter <- c(mcomp,mabun,mcarr,mgene,mremo)
abund <- rep(x/(x+3000),5)
spec <- c(rep("Pollinator specialization", 300),rep("Pollinator abundance", 300),rep("Carryover", 300),rep("Generalization", 300),rep("Removal", 300))
datafr<- data.frame(parameter,abund,spec)

# Plot of pollinator attributes and degree of generalization in function of the new colonist's abundance
p5 <- ggplot(datafr,aes(y = parameter, x = abund, colour = spec)) + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + geom_line(size = 1.8)
colour <- c("#009E73","black","#edae49","dodgerblue4","grey")
p5 <- p5 + coord_cartesian(xlim =c(0,0.5), ylim = c(0,1)) +
  labs(x="Relative abundance", y="Parameter value") +
  scale_colour_manual(values=colour) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_text(size=15), axis.title.y = element_text(size=15)) +
  theme(legend.text=element_text(size=13)) +
  theme(legend.position = "top") +
  theme(panel.border = element_rect(colour = "gray20", fill=NA, size=1.3)) +
  geom_ribbon(data = datafr[datafr$spec=="Pollinator specialization", ],aes(ymin = mcomp-scomp, ymax = mcomp+scomp),alpha = 0.3, colour = NA,fill = "dodgerblue4") +
  geom_ribbon(data = datafr[datafr$spec=="Pollinator abundance", ],aes(ymin = mabun-sabun, ymax = mabun+sabun),alpha = 0.3, colour = NA,fill = "#edae49") +
  geom_ribbon(data = datafr[datafr$spec=="Carryover", ],aes(ymin = mcarr-scarr, ymax = mcarr+scarr),alpha = 0.3, colour = NA,fill = "#009E73") +
  geom_ribbon(data = datafr[datafr$spec=="Generalization", ],aes(ymin = mgene-sgene, ymax = mgene+sgene),alpha = 0.3, colour = NA,fill = "black") +
  geom_ribbon(data = datafr[datafr$spec=="Removal", ],aes(ymin = mremo-sremo, ymax = mremo+sremo),alpha = 0.3, colour = NA,fill = "grey")
p5
