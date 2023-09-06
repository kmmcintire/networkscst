# example of hosts with overlapping epidemic

rm(list=ls())

setwd("~/Dropbox (University of Michigan)/Manuscripts/McIntireetal_Networks/Data&Code/daphnia-parasite-networks/2023/networks")

setwd("c:/Users/krism/Desktop/networks/networks")

# NOTE: requires all of the functions in "centrality-lake-functions.R"
source("centrality-lake-functions.R")

# packages
require(tidyverse)
require(zoo)
require(igraph)
require(bipartite)
require(dunn.test)
require(conover.test)
require(cowplot)


# looping through and plotting everything, then picking ones to use as examples ================

# general code ---------------------------------------------------------------------
clean_auc <- read.csv("Clean-Data-2014-2016.csv", header = TRUE)

# rename columns; easier to "gather" this way
colnames(clean_auc)[which(names(clean_auc) == "pasteuria.prev")] <- "pasteuria"
colnames(clean_auc)[which(names(clean_auc) == "metsch.prev")] <- "metsch"
colnames(clean_auc)[which(names(clean_auc) == "scarlet.prev")] <- "scarlet"
colnames(clean_auc)[which(names(clean_auc) == "brood.prev")] <- "brood"
colnames(clean_auc)[which(names(clean_auc) == "larssonia.prev")] <- "larssonia"
colnames(clean_auc)[which(names(clean_auc) == "spider.prev")] <- "spider"
colnames(clean_auc)[which(names(clean_auc) == "gurleya.prev")] <- "gurleya"

# split data depending on total number of hosts counted in parasite sample
total_less <- clean_auc %>%
  filter(Total < 20)

# if less than 20 total, set prevalence to 0
total_less$pasteuria <- 0
total_less$metsch <- 0
total_less$scarlet <- 0
total_less$brood <- 0
total_less$larssonia <- 0
total_less$spider <- 0
total_less$gurleya <- 0

# don't do anything if total is greater/equal to 20
total_greater <- clean_auc %>%
  filter(Total >= 20)

# join all data back together
loop_data <- full_join(total_less, total_greater)

# gather into tall format; select appropriate columns
loop_data <- loop_data %>%
  gather(Parasite.Species, Prevalence, pasteuria:gurleya) %>%
  select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Prevalence)

# arrange data 
loop_data <- loop_data %>%
  group_by(Host.Species, Year, Lake, Parasite.Species, Julian.Day) %>%
  arrange()

# set year and julian day as correct class
loop_data$Year <- as.factor(loop_data$Year)
loop_data$Julian.Day <- as.integer(loop_data$Julian.Day)

# order everything by julian day
loop_data <- loop_data[order(loop_data$Julian.Day), ] #make sure julian days are in order

# make sure it's in data frame format
loop_data <- as.data.frame(loop_data)

loop_data$Parasite.Species <- as.factor(loop_data$Parasite.Species)

# MAD is adding these next bits to convert from character to factor because it's not working with them as characters
loop_data$Parasite.Species <- factor(loop_data$Parasite.Species)
loop_data$Lake <- factor(loop_data$Lake)
loop_data$Focal.Host <- factor(loop_data$Focal.Host)###<-KMM- error!!!!!!!!

parasite_vec <- c("pasteuria", "brood", "scarlet") # set up vectors of names for loop
lake_vec <- c(levels(loop_data$Lake))
year_vec <- c(2014, 2015, 2016)



plots <- vector("list", 135) 

q <- 1

for (i in 1:length(parasite_vec)) {
  
  paradat <- loop_data %>% filter(Parasite.Species == paste0(parasite_vec[i]))
  
  for (j in 1:length(lake_vec)) {
    lakedat <- paradat %>% filter(Lake == paste0(lake_vec[j]))
    
    for (k in 1:length(year_vec)) {
      plotdat <- lakedat %>% filter(Year == paste0(year_vec[k]))
      
      p1 <- ggplot(data = plotdat, aes(x = Julian.Day, y = Prevalence, color = Host.Species)) +
        geom_point() +
        geom_line() +
        ggtitle(paste0(parasite_vec[i], lake_vec[j], year_vec[k]))
      
      plots[[q]] <- p1
      
      q <- q + 1
    }
  }
}



pdf(file = "exampleplots.pdf", width = 6, height = 5, useDingbats=F) 
plots
dev.off()

i <- 1
j <- 1
k <- 1


# ones to use for prevalence

bruin2014scar <- loop_data %>%
  filter(Lake == "Bruin" & Year == 2014 & Parasite.Species == "scarlet") %>%
  filter(Host.Species == "retrocurva"|Host.Species == "pulicaria")

crookedp2016brood <- loop_data %>%
  filter(Lake == "CrookedP" & Year == 2016 & Parasite.Species == "brood") %>%
  filter(Host.Species == "dentifera"|Host.Species == "ceriodaphnia")



# join all data back together
loop_data <- full_join(total_less, total_greater)

str(loop_data)




# gather into tall format; select appropriate columns
loop_data <- loop_data %>%
  gather(Parasite.Species, Prevalence, pasteuria:gurleya) %>%
  select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Prevalence, Density)

inf_host_density <- loop_data %>%
  mutate(infected_density = Prevalence*Density)


# arrange data 
inf_host_density <- inf_host_density %>%
  group_by(Host.Species, Year, Lake, Parasite.Species, Julian.Day) %>%
  arrange()

# set year and julian day as correct class
inf_host_density$Year <- as.factor(inf_host_density$Year)
inf_host_density$Julian.Day <- as.integer(inf_host_density$Julian.Day)

# order everything by julian day
inf_host_density <- inf_host_density[order(inf_host_density$Julian.Day), ] #make sure julian days are in order

# make sure it's in data frame format
inf_host_density <- as.data.frame(inf_host_density)

inf_host_density$Parasite.Species <- as.factor(inf_host_density$Parasite.Species)

# MAD is adding these next bits to convert from character to factor because it's not working with them as characters
inf_host_density$Lake <- factor(inf_host_density$Lake)


parasite_vec <- c("pasteuria", "brood", "scarlet") # set up vectors of names for loop
lake_vec <- c(levels(inf_host_density$Lake))
year_vec <- c(2014, 2015, 2016)



dens_plots <- vector("list", 135) 

q <- 1

for (i in 1:length(parasite_vec)) {
  
  paradat <- inf_host_density %>% filter(Parasite.Species == paste0(parasite_vec[i]))
  
  for (j in 1:length(lake_vec)) {
    lakedat <- paradat %>% filter(Lake == paste0(lake_vec[j]))
    
    for (k in 1:length(year_vec)) {
      plotdat <- lakedat %>% filter(Year == paste0(year_vec[k]))
      
      p1 <- ggplot(data = plotdat, aes(x = Julian.Day, y = infected_density, color = Host.Species)) +
        geom_point() +
        geom_line() +
        ggtitle(paste0(parasite_vec[i], lake_vec[j], year_vec[k]))
      
      dens_plots[[q]] <- p1
      
      q <- q + 1
    }
  }
}



pdf(file = "example_dens_plots.pdf", width = 6, height = 5, useDingbats=F) 
dens_plots
dev.off()



# SCARLET IN BRUIN 2014

inf_host_density

bruin2014scardens <- inf_host_density %>%
  filter(Lake == "Bruin" & Year == 2014 & Parasite.Species == "scarlet") %>%
  filter(Host.Species == "retrocurva"|Host.Species == "pulicaria") %>%
  mutate(Days = Julian.Day - 2014204)



min(bruin2014scardens$Julian.Day)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73")


retro <- bruin2014scardens %>% filter(Host.Species == "retrocurva")
y1 <- retro$infected_density
pulic <- bruin2014scardens %>% filter(Host.Species == "pulicaria")
y2 <- pulic$infected_density


y1[y1 <= 0] <- NA
y2[y2 <= 0] <- NA

ttt <- cbind(y1,y2)

#avg_epi <- rowMeans(ttt, na.rm = FALSE)

avg_epi <- apply(ttt[, c("y1", "y2")], 1, function(x) (prod(x[x!=0]))^(1/sum(x!=0)))
avg_epi[is.na(avg_epi)] <- 0

dday <-retro$Days

bruin14plot <- ggplot() +
  #geom_area(aes(x = retro$Days, y = pmin(y1, y2)), fill = 'black', alpha = 0.1) +
  geom_area(aes(x = dday, y = avg_epi), fill = 'black', alpha = 0.1) +
  geom_point(data = bruin2014scardens, 
             aes(x = Days, y = infected_density, color = Host.Species, shape = Host.Species),
             size = 3) +
  geom_line(data = bruin2014scardens, 
            aes(x = Days, y = infected_density, color = Host.Species),
            size = 1) +
  #ylab(bquote("Infected host density (animals/"~m^2,")")) +
  ylab(expression(paste("Infected host density (animals/",m^2,")", sep = ""))) +
  xlab("Days") +
  scale_fill_manual(values = c("#d95f02", "#7570b3")) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100, 125),
                     labels = c("0", "25", "50", "75", "100", "125")) +
  guides(color = guide_legend(title = "Host Species"),
         shape = guide_legend(title = "Host Species")) +
  theme(legend.position = "top") +
  theme_bw()
bruin14plot  


pdf(file = "Bruin_2014_epidemic_overlap.pdf", width = 6, height = 4, useDingbats=F) 
bruin14plot
dev.off()



