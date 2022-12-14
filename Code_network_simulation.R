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

# Simulations of the evolutionary assembly of flower communities
for (j in 1:num.sims){ # Loop for each simulation
  #Sampling parameter values for this simulation
  plantpopsizes.det[,j] <- rpoilog(length(plantsp.final), mu = 5, sig =1) %>% +1 %>% # Randomly sampling the maximal plant abundances (kj).
    "/" (mean(.)) %>% "*"(av.plant.abund) %>% ceiling(.)
  plantpopsizes <- plantpopsizes.det[1]
  plantsp <- c(1) # List of plant species in the community, starting with one species (the list expands as other species colonize the community).
  pollpopsizes[,1,j] <- rpoilog(length(pollisp), mu = 5, sig = 0.7) %>% # Randomly sampling the initial number of visits by each pollinator species.
    "/" (mean(.)) %>% "*"(visits.per.reward) %>% ceiling(.) 
  pollpopsizes.max[,j] <- rpoilog(length(pollisp), mu = 5, sig = 0.7) %>% # Sampling the values of Ki (limitation by factors other than reward availability).
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
  
  adapt.for <- array(1,dim=c(length(plantsp.final))) # Values of adaptive foraging
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
        network[[j]][i,m] <- pollpopsizes[m,g,j]*((plantpopsizes[i]*genome[i,m,g,j]*adapt.for[i])/sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp]))
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
      adapt.for[i] <- adapt.for[i] * (reward.f[i]/mean(reward.f[plantsp]))
    }
    # Measuring for each plant species male and female fitness and the selection gradients for attractiveness to each pollinator species.
    for(i in seq_along(plantsp)){ # Loop for each plant species
      # Measuring the number of pollen grains removed by each pollinator alone ("info.removed") and the proportion of pollen removed by each pollinator ("H").
      info.removed <- matrix(rep(0,length(pollisp)),nrow=1)
      H <- matrix(rep(0,length(pollisp)),nrow=1)
      for (m in 1:length(pollisp)) {
        info.removed[1,m] <- (1-polliremo[m,j]) ^ ((pollpopsizes[m,g,j] * genome[i,m,g,j] *adapt.for[i]) /
                                                     sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp]))
        H[1,m] <-  polliremo[m,j] * pollpopsizes[m,g,j] * genome[i,m,g,j] * adapt.for[i] / sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp])
      }
      # Measuring, for each pollinator, its contribution to male fitness and the selection gradient for attractiveness to this pollinator.
      for(t in seq_along(pollisp)){ # Loop for each pollinator species
        sel.g <- numeric(0) # vector of the effect of variation in attractiveness to a given pollinator on pollination by each pollinator
        # Calculating the contribution of each pollinator to the focal plant's male fitness.
        poll.dep[i,t,g,j] <- (Pj[i]*tj[i]*(1-prod(info.removed)) * (H[1,t]/sum(H)))  * pollicarry[t,j]^
          (sum(plantpopsizes*genome[plantsp,t,g,j]*adapt.for[plantsp])/(plantpopsizes[i]*genome[i,t,g,j]*adapt.for[i]))
        # Measuring selection gradient from male fitness using the partial derivative of male fitness as a function of attractiveness to a given pollinator.
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
            a <- adapt.for[i]
            A <- sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp])
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
            A <- sum(plantpopsizes*genome[plantsp,m,g,j]*adapt.for[plantsp])
            G <- sum(H)
            C <- pollicarry[m,j]
            B <- plantpopsizes[i]
            D <- sum(plantpopsizes*genome[plantsp,t,g,j]*adapt.for[plantsp])
            x <- genome[i,t,g,j]
            a <- adapt.for[i]
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
        if (femalfit[i,g,j] == 0) { # If female fitness is zero (i.e no pollination occurred) selection on attractiveness is zero (avoid division by zero).
          fg[i,m] <- 0
        }
        else { # Otherwise measure selection gradient
          # Measure selection gradient 
          fg[i,m] <- selection.gradient[i,m]/sum(poll.dep[i,,g,j]) + selection.gradient.f[i,m]/femalfit[i,g,j]
        }
        # Change in attractiveness
        # The additive genetic variance is a non-negative function of attractiveness,
        #that approaches zero when the mean attractiveness value is close to its minimum or maximum.
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

# Because the values of attractiveness are prevented from being lower than 0.01 (simulating the input of genetic variability),
# the values of attractiveness of 0.01 are now changed back to zero to prevent the networks from containing every possible link.
for (j in 1:num.sims){
  for (m in 1:length(pollisp)){
    for (i in 1:length(plantsp.final)){
      if (genome[i,m,g,j] == 0.01){
        genome[i,m,g,j] = 0
        network[[j]][i,m] = 0
      }
    }
  }
}

# Removing pollinators that are effectively extinct in the community (have no link with any plant species, population size = 1).
# It is necessary to remove those pollinators not participating in the network to measure NODFc.
# Pollinator extinctions were rare under most parameter values, but could occur when the pollinator-to-plant ratio was high.
for (j in 1:num.sims){
  if (length(which(colSums(network[[j]]) == 0)) > 0){
    network[[j]] <- network[[j]][,-which(colSums(network[[j]]) == 0)]
  }
}

# Visualization of the resulting plant-pollinator networks
plotweb(network[[j]],col.high ="#990000",bor.col.high= "#990000",
        col.low = "#003366",bor.col.low="#003366",
        col.interaction = "#999999",bor.col.interaction = "#999999")
visweb(network[[j]])

# Plot of degree of generalization in function of plant abundance
# The pollinators of a plant species were determined as the ones for which the plant's attractiveness was more than 0.2.
x <- numeric(0)
y <- numeric(0)
for (j in 1:num.sims){
  for (i in 1:length(plantsp.final)){
    y <- c(y,length(which(genome[i,,num.gen,j] > 0.2)))
    x <- c(x,plantpopsizes.det[i,j]/sum(plantpopsizes.det[,j]))
  }
}
gen <- data.frame(x,y)
ggplot(gen, aes(x=x, y=y)) +
  geom_point(alpha=0.3, color="#c0392b") +
  geom_smooth(alpha=0.25, color="black", fill="black") +
  labs(x="Relative abundance", y="Generalization") +
  coord_cartesian(xlim =c(0,0.5), ylim = c(0,length(pollisp))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# Measuring average and variability in degree of generalization in the plant communities.
specialization<-array(0,dim=c(length(plantsp),num.sims))
for (j in 1:num.sims) {
  for (i in 1:length(plantsp.final)) {
    specialization[i,j] <- length(which(genome[i,,num.gen,j] > 0.2)) # degree of generalization of each plant species
  }
}
specialization.mean <- array(0,dim=c(num.sims)) # Average degree of generalization in each community
specialization.sd <- array(0,dim=c(num.sims)) # Standard deviation in degrees of generalization in each community
for (j in 1:num.sims) {
  specialization.mean[j] <- mean(specialization[,j])
  specialization.sd[j] <- sd(specialization[,j])
}
mean(specialization.mean) # Average degree of generalization within communities across simulations
sd(specialization.mean) # Standard deviation in average degree of generalization within communities
mean(specialization.sd) # Standard deviation in degree of generalization within communities across simulations
sd(specialization.sd) # Standard deviation of this measure of variation

# Measuring network structure
NODFc <- array(0,dim=c(num.sims)) # NODFc
nestedness <- array(0,dim=c(num.sims)) # NODF
connectance <- array(0,dim=c(num.sims)) # Connectance
nsharedpartners <- array(0,dim=c(num.sims)) # Number of shared partners
for (j in 1:num.sims) {
  NODFc[j] <- NODFc(web = network[[j]])
  nestedness[j] <- networklevel(network[[j]],index = 'NODF')
  connectance[j] <- networklevel(network[[j]],index = 'connectance')
  nsharedpartners[j] <- networklevel(network[[j]],index = 'mean number of shared partners')[2]
}
mean(NODFc) # Average NODFc across simulations
sd(NODFc) # Standard deviation of NODFc
mean(nestedness) # Average NODF across simulations
sd(nestedness) # Standard deviation of NODF
mean(connectance) # Average connectance across simulations
sd(connectance) # Standard deviation of connectance
mean(nsharedpartners) # Average of average number of shared partners across simulations
sd(nsharedpartners) # Standard deviation of of average number of shared partners