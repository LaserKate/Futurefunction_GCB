# LIBRARIES
library("vegan")
library("spatstat")
library("psych")
library("raster")
library("plotrix")
library("dplyr")

#Data source and function source: Madin et al. 2023 Journal of Applied Ecology
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2664.14447

# FUNCTIONS
source("functions.R")

# DATA PREPARATION
source("data_prep.R")
#inside

#data prep----
##########

###***Analysis for every Reef ----
#Calculate rate of change per reef region with DHW from Gonzalez-Barrios----
RateChange_Gonzalez #download raw data directly from Gonzalez-Barrios article
library(ggpubr)
library(dplyr)

RateofChange_Aus<-RateChange_Gonzalez %>% filter(Ecoregion == c("Exmouth_to_Broome", "Ningaloo","Torres_Strait_Northern_Great_Barrier_Reef", "Central_and_Southern_Great_Barrier_Reef")) %>% ggplot(aes(x = DHW, y = absolute_rate_t)) +
  geom_point(aes( color = Ecoregion)) + geom_smooth(method = "lm", aes(DHW, absolute_rate_t), se=T,lty=3, color="black")+
  labs(title="Australian reefs", x="Degree Heating Weeks", y= "Mean annual rate of change")+theme_bw()+
  stat_regline_equation(aes(DHW, absolute_rate_t))

###***END Analysis for every Reef ----
#****-----
#*****----

#LHI Analysis ***----
#A. ---
#Original data McWilliam et al. 2018 but data csv available within Madin et al. R workspace in the data folder
#QAQC as per Methods outlined in Quigley & Baird 2024
ctd <- read.csv("data", as.is=TRUE) #

dat <- unique(ctd[c("specie_name", "family_molecules")])
names(dat) <- c("species", "family_molecules")
dat$genus <- sapply(strsplit(dat$species, " "), "[", 1)

# range size (Hughes et al.)
rs <- ctd[ctd$trait_name=="Range size" & ctd$resource_id==605, c("specie_name", "value")]
rs$value <- as.numeric(rs$value)
names(rs) <- c("species", "range")
rs <- aggregate(range ~ species, rs, mean)
dat <- merge(dat, rs, all.x=TRUE)
dat$range <- dat$range / max(dat$range, na.rm=TRUE) # normalize

# Local abundance
ab <- ctd[ctd$trait_name=="Abundance GBR", c("specie_name", "value")]
names(ab) <- c("species", "abund")
dat <- merge(dat, ab, all.y=TRUE)
dat$abund[dat$abund=="common"] <- 1  # normalize
dat$abund[dat$abund=="uncommon"] <- 0.5 # normalize
dat$abund[dat$abund=="rare"] <- 0.25 # normalize
dat$abund <- as.numeric(dat$abund)

# Bleaching susceptibility
names(bri) <- c("species", "BRI")
BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$bleach <- bri$BRI[match(dat$species, bri$species)]
dat$bleach[is.na(dat$bleach)] <- BRI_genus[is.na(dat$bleach)]
dat$bleach <- 1 - (dat$bleach / 100) # normalize

# Growth form typical
gf <- ctd[ctd$trait_name=="Growth form typical", c("specie_name", "value")]
names(gf) <- c("species", "growth_form")
dat <- merge(dat, gf, all.x=TRUE)

# Life history
lh <- ctd[ctd$trait_name=="Life history strategy", c("specie_name", "value")]
names(lh) <- c("species", "life_history")
lh$life_history[lh$life_history=="weedy"] <- "Weedy"
lh$life_history[lh$life_history=="generalist"] <- "Generalist"
lh$life_history[lh$life_history=="competitive"] <- "Competitive"
lh$life_history[lh$life_history=="stress-tolerant"] <- "Stress-tolerant"
dat <- merge(dat, lh, all.x=TRUE) #400

# Traits - McWilliam et al. 2018 PNAS
dim(tra[tra$domain=="pacific",])
tra <- tra[,c("species", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat <- merge(dat, tra) #391
dim(dat)

# Goreau Madin ellipses
lh <- ctd[ctd$location_name %in% c("Great Barrier Reef (GBR)", "Indo-Pacific (unspecified)"), c("specie_name", "goreau_Madin")]
names(lh) <- c("species", "goreau_Madin")
lh$goreau_Madin[lh$goreau_Madin=="builders"] <- "Builders"
lh$goreau_Madin[lh$goreau_Madin=="cementers"] <- "Cementers"
lh$goreau_Madin[lh$goreau_Madin=="fillers"] <- "Fillers"
lh$goreau_Madin[lh$goreau_Madin=="cementers&fillers"] <- "Cementers & Fillers"
dat <- merge(dat, lh, all.x=TRUE) #496
dim(dat)

dat$goreau <- NA
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_spacesize %in% c(4, 5)] <- "Fillers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(1, 2)] <- "Cementers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(3, 4, 5) & dat$cat_SA_vol  %in% c(1, 2)] <- "Builders"

# Restoration potential from PLoS paper
dat$restore <- 0
dat$restore[dat$growth_form %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$growth_form %in% c("massive", "submassive")] <- 5
dat$restore[dat$growth_form %in% c("laminar")] <- 4
dat$restore[dat$growth_form %in% c("columnar", "encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$growth_form %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$growth_form %in% c("hispidose")] <- 1
dat$restore <- (dat$restore) / 6 # normalize
dim(dat)

dat <- dat[!duplicated(dat$species), ]
dim(dat)

#data prep----

#need to make sure that dat has the correct PCs

# TRAIT SPACE
#start at 358
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) #

## all species Lord Howe----
lhi <- read.csv("lordhowe_updateAB23.csv", as.is=TRUE)
dim(lhi) 
lhi$species[!(lhi$species %in% dat$species)]
dat$lhi <- dat$species %in% lhi$species
sum(dat$lhi) 

# Identify rows with duplicate names in the first column
duplicate_names <- duplicated(dat$species)
# Remove rows with duplicate names
dat.test <- dat[!duplicate_names, ]
sum(dat.test$lhi) #

# export
dat_lhi<- dat[dat$lhi,] 
dim(dat_lhi)

#B. Lord Howe Abundances (Keith et al.)----
lhi_Keith_Abun ##Transect data from Keith et al. 2014 transect data. Contact A.Baird for this or find in publication
head(lhi_Keith_Abun)
lhi_Keith_Abun<-lhi_Keith_Abun %>% filter(Reef == c("Lord Howe Island"))
library(reshape)
melted_data <- melt(lhi_Keith_Abun, id.vars = c("date", "observer", "Reef", "Habitat", "species", "site", "transect"))
print(melted_data)

means_per_species_persite <- melted_data %>%
  group_by(species, site, transect) %>% 
  summarize(
    Mean_ = mean(value,na.rm = TRUE))

# View the resulting dataframe with means
print(means_per_species_persite)

means_per_species_persite %>% #filter(dat_lhi.genus == "Acropora") %>% 
  ggplot(aes(y = Mean_, x = site, fill=species))+ 
  geom_bar(stat = "identity", position="fill")

means_per_species_persite %>% #filter(dat_lhi.genus == "Acropora") %>% 
  ggplot(aes(y = Mean_, x = species, fill=site))+ 
  geom_bar(stat = "identity", position="fill")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

means_per_species <- melted_data %>%
  group_by(species) %>% 
  summarize(
    Mean_2 = mean(value,na.rm = TRUE),
  )
print(means_per_species)

means_per_species %>% #filter(dat_lhi.genus == "Acropora") %>% 
  ggplot(aes(y = Mean_2, x = species, fill=Mean_2))+ 
  geom_bar(stat = "identity")

####

#C. Adding in new relative abundances into the Choice function----
#Ecological persistence Lord Howe score-----

dat_lhi<- dat[dat$lhi,] 
dim(dat_lhi)

merged_df_lhi <- merge(means_per_species, dat_lhi, by = "species", all = TRUE) #

# Add conditional logic
merged_df_lhi$CombinedRelAbud <- ifelse(!is.na(merged_df_lhi$Mean_2), merged_df_lhi$Mean_2, merged_df_lhi$abund)

merged_df_lhi$abund_1.5C <- merged_df_lhi$CombinedRelAbud + (-5.78/100)
merged_df_lhi$abund_2C <- merged_df_lhi$CombinedRelAbud + (-15.86/100)

Abundance_lineplot_LH<-merged_df_lhi
dat_lhi2<-merged_df_lhi #just to save the dataset
dim(merged_df_lhi) #

#change all negatives to zeros
dat_lhi2$CombinedRelAbud <- ifelse(dat_lhi2$CombinedRelAbud < 0, 0, dat_lhi2$CombinedRelAbud)
dat_lhi2$abund_1.5C <- ifelse(dat_lhi2$abund_1.5C < 0, 0, dat_lhi2$abund_1.5C)
dat_lhi2$abund_2C <- ifelse(dat_lhi2$abund_2C < 0, 0, dat_lhi2$abund_2C)

#normalize Common, Uncommon, Rare
dat_lhi2$CombinedRelAbud[dat_lhi2$CombinedRelAbud > 0.5 & dat_lhi2$CombinedRelAbud <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_1.5C[dat_lhi2$abund_1.5C > 0.5 & dat_lhi2$abund_1.5C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_2C[dat_lhi2$abund_2C > 0.5 & dat_lhi2$abund_2C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare

dat_lhi2$CombinedRelAbud[dat_lhi2$CombinedRelAbud > 0.25 & dat_lhi2$CombinedRelAbud <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_1.5C[dat_lhi2$abund_1.5C > 0.25 & dat_lhi2$abund_1.5C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_2C[dat_lhi2$abund_2C > 0.25 & dat_lhi2$abund_2C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare

dat_lhi2$CombinedRelAbud[dat_lhi2$CombinedRelAbud > 0 & dat_lhi2$CombinedRelAbud <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_1.5C[dat_lhi2$abund_1.5C > 0.0 & dat_lhi2$abund_1.5C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_lhi2$abund_2C[dat_lhi2$abund_2C > 0.0 & dat_lhi2$abund_2C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare

####
library(tidyverse)
df_long <- pivot_longer(
  dat_lhi2,
  cols = c("CombinedRelAbud", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

df_long <- df_long %>%
  filter(!(Scenario == "CombinedRelAbud" & Abund == 0))

df_long$Scenario <- factor(df_long$Scenario, levels = c("CombinedRelAbud", "abund_1C", "abund_1.5C", "abund_2C"))

RelativeAbund_Allscenarios<- ggplot(df_long, aes(Scenario, species, fill= Abund)) + 
  geom_tile()+scale_fill_viridis_c()


#D.Estimate Ecological persistance from normalized abundance, bleaching susceptibility, and range size----
s = 82
dim(dat_lhi2)
dat_lhi2$abund_2C <- ifelse(dat_lhi2$abund_2C > 0, -1, dat_lhi2$abund_2C)

#calculate ecological pers Madin et al. 2023
#input to equation must be called, needs to be abund called for input to work
dat_lhi2_Norm<-dat_lhi2
dat_lhi2_Norm$abund<-dat_lhi2$CombinedRelAbud

dat_lhi2_1.5<-dat_lhi2
dat_lhi2_1.5$abund<-dat_lhi2$abund_1.5

dat_lhi2_2<-dat_lhi2
dat_lhi2_2$abund<-dat_lhi2$abund_2

choice_1 <- choice(dat_lhi2_Norm, s, vars=c("abund","range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_1) <- c("Species", "Ecological persistence Present")
choice_2 <- choice(dat_lhi2_1.5, s, vars=c("abund","range","bleach"), trait=FALSE)[c("species", "value")]
names(choice_2) <- c("Species", "Ecological persistence 1.5C")
choice_3 <- choice(dat_lhi2_2, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_3) <- c("Species", "Ecological persistence 2C")

matx <- merge(choice_1, choice_2, all=TRUE)
matx <- merge(matx, choice_3, all=TRUE)
matx$`Ecological persistence Present`[is.na(matx$`Ecological persistence Present`)] <- 0
matx$`Ecological persistence 1.5C`[is.na(matx$`Ecological persistence 1.5C`)] <- 0
matx$`Ecological persistence 2C`[is.na(matx$`Ecological persistence 2C`)] <- 0
matx <- matx[order(matx["Ecological persistence Present"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 1.5C"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 2C"], decreasing=FALSE),]

image(t(as.matrix(matx[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(matx)-1))), labels=matx$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"), cex.axis=0.8, las=2)


#E. Plot Ecological Persistence with ggplot----
mat_lhi <- merge(choice_1, choice_2, all=TRUE)
mat_lhi <- merge(mat_lhi, choice_3, all=TRUE)
mat_lhi$Ecol.Persis<-mat_lhi$`Ecological persistence`

#Now determine which species are lost with different scenarios and re-run the trait space

df_long_heatmap <- pivot_longer(
  mat_lhi,
  cols = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"),
  names_to = "Scenario",
  values_to = "Ecol.Persistance"
)

df_long_heatmap$Scenario <- factor(df_long_heatmap$Scenario, levels = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"))

(Ecol.Pers_Allscenarios<- ggplot(df_long_heatmap, aes(x=Scenario, y = Species, fill = Ecol.Persistance)) +
  geom_tile() + scale_fill_viridis_c())

#MODEL output saved -----
#save(mat_lhi, file = "LHI_model_pers.RData")

#F. Trait space----

#trait space made on only the corals we have data for for the persistence metric
#run from very top from Madin et al. data
dim(dat) 
duplicate_names <- duplicated(dat$species)
dat.test <- dat[!duplicate_names, ]
dim(dat.test) #

## all spp

lhi_all #use model output from above to create species list for present that have data needed for this step
lhi_all$species[!(lhi_all$species %in% dat$species)]
dat$lhi_all <- dat$species %in% lhi_all$species
sum(dat$lhi_all)
duplicate_names <- duplicated(dat$species)
# Remove rows with duplicate names
dat.test <- dat[!duplicate_names, ]
sum(dat.test$lhi_all) #

lhi_all_1.5 #use model output from above to create species list for present that have data needed for this step
dim(lhi_all_1.5) #
lhi_all_1.5$species[!(lhi_all_1.5$species %in% dat.test$species)]
dat.test$lhi_all_1.5 <- dat.test$species %in% lhi_all_1.5$species
sum(dat.test$lhi_all_1.5) #
#

#G. brought over ----
dat<-dat.test 
dat <- dat[!duplicated(dat$species), ] 
dim(dat)

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) #

#H. create convex hull for all data----
TABS<-dat[, 25:26] 
CHS<- chull(TABS)

#dat <- dat[, -c(19:20)] #remove extra if needed, check your dataframe


##I. create convex hulls for each goreau class
lhi_h <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & lhi_all == 'TRUE',
         select=c('species','goreau_Madin', 'lhi_all','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#present
tab7 <- lhi_h[1:20, 4:5] #builders
tab8 <- lhi_h[21:34, 4:5] #Cementer
tab9 <- lhi_h[35:38, 4:5] #Cementers & Fillers
tab901 <- lhi_h[39:47, 4:5] #fillers

lhi_h.15 <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & lhi_all_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'lhi_all_1.5','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#1.5C
tab10 <- lhi_h.15[1:13, 4:5] #builders
tab11 <- lhi_h.15[14:24, 4:5] #Cementer
tab12 <- lhi_h.15[25:25, 4:5] #Cementers & Fillers
tab1201 <- lhi_h.15[26:31, 4:5] #fillers

#---

lhi_ <- dat%>% 
  subset(lhi_all == 'TRUE',
         select=c('species','goreau_Madin', 'lhi_all','PC1','PC2','PC3'))
tab000 <- lhi_[, 4:5]

ch000 <- chull(tab000)
ch7 <- chull(tab7)
ch8 <- chull(tab8)
ch9 <- chull(tab9)
ch901 <- chull(tab901)

lhi_15 <- dat%>% 
  subset(lhi_all_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'lhi_all_1.5','PC1','PC2','PC3'))
tab000_1.5 <- lhi_15[, 4:5]

ch000_1.5 <- chull(tab000_1.5)
ch7_1.5 <- chull(tab10)
ch8_1.5 <- chull(tab11)
ch9_1.5 <- chull(tab12)
ch901_1.5 <- chull(tab1201)


#I.Plot PCA for lord howe island climate scenarios----
plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
lhi_pres <- dat[dat$lhi_all,] #
lhi_1.5 <- dat[dat$lhi_all_1.5,] #
mtext("Lord Howe Island under different climate scenarios", line=-1, adj=0, cex=1.5)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("PC1", 1, 2, cex=0.8)
mtext("PC2", 2, 2, cex=0.8)
text(2, -2.5, "Present: n=62, 1.5C: n=45", cex=1.5)
polygon(tab000[ch000, ], border = "#73c391", col = rgb(115, 195, 145, 128, maxColorValue = 255), lwd = 2) #present day
polygon(tab000_1.5[ch000_1.5, ], border="#f36e45", col = rgb(243, 110, 69, 100, maxColorValue = 255), lwd=2) #1.5C
polygon(TABS[CHS, ], border="black", lwd=2)
points(PC2 ~ PC1, lhi_pres, col="#73c391", pch=20, cex=1)
points(PC2 ~ PC1, lhi_1.5, col="#f36e45", pch=20, cex=1.5)

#save(dat, lhi_pres,lhi_1.5, TABS,file="LHI_traitspaceGCB.RData")
#Save MODEL data ----

dat %>% group_by(lhi_all)  %>% group_by(goreau_Madin) %>% summarise(n())
dat %>% group_by(lhi_all_1.5)  %>% group_by(goreau_Madin) %>% summarise(n())
lhi_pres %>% group_by(lhi_all, genus)  %>% group_by(goreau_Madin) %>% summarise(n())
lhi_pres %>% group_by(lhi_all, genus,goreau_Madin) %>% summarise(n())

#J. Loss in function ----
lhi_pres$genus 
class(lhi_pres$genus)
lhi_pres$genus <- as.factor(lhi_pres$genus)

result_df_pres <- lhi_pres %>%
  group_by(lhi_all, genus,goreau_Madin) %>% summarise(proportionPres = n())

result_df_1.5 <- lhi_1.5 %>%
  group_by(lhi_all_1.5, genus, goreau_Madin) %>% summarise(proportion1.5 = n())


#-----
combined_df <- merge(result_df_pres, result_df_1.5, by = c("genus", "goreau_Madin"), all = TRUE)

combined_df <- combined_df %>%mutate(proportion1.5 = ifelse(is.na(proportion1.5), 0, proportion1.5))

melted_df <- melt(combined_df,
                  id.vars = c("genus", "goreau_Madin"),
                  measure.vars = c("proportionPres", "proportion1.5"),
                  variable.name = c("time"),
                  value.name = "value")

melted_df$goreau_Madin <- ifelse(is.na(melted_df$goreau_Madin) | melted_df$goreau_Madin == "NA", "NAs", melted_df$goreau_Madin)

custom_order_LHI <- c("NAs", "Cementers", "Fillers", "Builders", "Cementers & Fillers")

# Change the order of the factor levels in goreau_Madin
melted_df$goreau_Madin <- factor(melted_df$goreau_Madin, levels = custom_order_LHI)

# Calculate percentage change using dplyr
percentage_change_df <- melted_df %>%
  group_by(goreau_Madin, variable) %>% summarise(proportion = sum(value))

(percentage_change_result <- percentage_change_df %>%
  group_by(goreau_Madin) %>% arrange(variable) %>%
  mutate(percentage_change = (proportion / lag(proportion) - 1) * 100))

#K. Change in Volume LHI ----
# Assuming 'tab000' and 'tab000_1.5' are your data frames containing PC1 and PC2 scores,
# and 'ch000' and 'ch000_1.5' contain the indices

# Function to calculate the area of a convex hull in 2D space
calculate_convex_hull_area <- function(data, indices) {
  subset_points <- data[indices, ]
  n <- nrow(subset_points)
  
  # Duplicate the first point to close the polygon
  subset_points <- rbind(subset_points, subset_points[1, ])
  
  # Calculate the area using the Shoelace formula
  area <- 0.5 * sum(subset_points[1:(n-1), 1] * subset_points[2:n, 2] -
                      subset_points[2:n, 1] * subset_points[1:(n-1), 2])
  
  return(abs(area))
}

# Calculate areas for present day and 1.5C scenarios
area_present_day <- calculate_convex_hull_area(tab000, ch000)
area_1.5C <- calculate_convex_hull_area(tab000_1.5, ch000_1.5)

# Assuming the convex hull is in 2D, volume is then area in 2D
volume_present_day <- area_present_day #
volume_1.5C <- area_1.5C #

# Given values
value1 <- volume_present_day
value2 <- volume_1.5C

# Calculate percent difference
(percent_difference <- abs((value1 - value2) / ((value1 + value2) / 2)) * 100)

# Calculate the difference in volumes
volume_difference <- abs(volume_present_day - volume_1.5C)

# Print the result
cat("Difference in Convex Hull Volumes:", volume_difference, "\n")


#*****
#*
#*
#Ningaloo Analysis ****----
#A. ---
#Original data McWilliam et al. 2018 but data.csv available within Madin et al. R workspace
#QAQC as per Methods in Quigley & Baird 2024
ctd <- read.csv("data", as.is=TRUE) #

dat <- unique(ctd[c("specie_name", "family_molecules")])
names(dat) <- c("species", "family_molecules")
dat$genus <- sapply(strsplit(dat$species, " "), "[", 1)

# range size (hughes)
rs <- ctd[ctd$trait_name=="Range size" & ctd$resource_id==605, c("specie_name", "value")]
rs$value <- as.numeric(rs$value)
names(rs) <- c("species", "range")
rs <- aggregate(range ~ species, rs, mean)
dat <- merge(dat, rs, all.x=TRUE)
dat$range <- dat$range / max(dat$range, na.rm=TRUE) # normalize

# Local abundance
ab <- ctd[ctd$trait_name=="Abundance GBR", c("specie_name", "value")]
names(ab) <- c("species", "abund")
dat <- merge(dat, ab, all.y=TRUE)
dat$abund[dat$abund=="common"] <- 1  # normalize
dat$abund[dat$abund=="uncommon"] <- 0.5 # normalize
dat$abund[dat$abund=="rare"] <- 0.25 # normalize
dat$abund <- as.numeric(dat$abund)

# Bleaching susceptibility
names(bri) <- c("species", "BRI")
BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$bleach <- bri$BRI[match(dat$species, bri$species)]
dat$bleach[is.na(dat$bleach)] <- BRI_genus[is.na(dat$bleach)]
dat$bleach <- 1 - (dat$bleach / 100) # normalize

# Growth form typical
gf <- ctd[ctd$trait_name=="Growth form typical", c("specie_name", "value")]
names(gf) <- c("species", "growth_form")
dat <- merge(dat, gf, all.x=TRUE)

# Life history
lh <- ctd[ctd$trait_name=="Life history strategy", c("specie_name", "value")]
names(lh) <- c("species", "life_history")
lh$life_history[lh$life_history=="weedy"] <- "Weedy"
lh$life_history[lh$life_history=="generalist"] <- "Generalist"
lh$life_history[lh$life_history=="competitive"] <- "Competitive"
lh$life_history[lh$life_history=="stress-tolerant"] <- "Stress-tolerant"
dat <- merge(dat, lh, all.x=TRUE) 

dim(tra[tra$domain=="pacific",])
tra <- tra[,c("species", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat <- merge(dat, tra) 

lh <- ctd[ctd$location_name %in% c("Great Barrier Reef (GBR)", "Indo-Pacific (unspecified)"), c("specie_name", "goreau_Madin")]
names(lh) <- c("species", "goreau_Madin")
lh$goreau_Madin[lh$goreau_Madin=="builders"] <- "Builders"
lh$goreau_Madin[lh$goreau_Madin=="cementers"] <- "Cementers"
lh$goreau_Madin[lh$goreau_Madin=="fillers"] <- "Fillers"
lh$goreau_Madin[lh$goreau_Madin=="cementers&fillers"] <- "Cementers & Fillers"
dat <- merge(dat, lh, all.x=TRUE) 

dat$goreau <- NA
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_spacesize %in% c(4, 5)] <- "Fillers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(1, 2)] <- "Cementers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(3, 4, 5) & dat$cat_SA_vol  %in% c(1, 2)] <- "Builders"

# Restoration potential from PLoS paper
dat$restore <- 0
dat$restore[dat$growth_form %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$growth_form %in% c("massive", "submassive")] <- 5
dat$restore[dat$growth_form %in% c("laminar")] <- 4
dat$restore[dat$growth_form %in% c("columnar", "encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$growth_form %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$growth_form %in% c("hispidose")] <- 1
dat$restore <- (dat$restore) / 6 # normalize

dat <- dat[!duplicated(dat$species), ] 
dim(dat )

#data prep----

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat)

## all species ningaloo---

ning_all<- read.csv("ningaloo_all.csv", as.is=TRUE)

dim(ning_all)
ning_all$species[!(ning_all$species %in% dat$species)]
dat$ning_all <- dat$species %in% ning_all$species 
sum(dat$ning_all) 
# Identify rows with duplicate names in the first column
duplicate_names <- duplicated(dat$species)
# Remove rows with duplicate names
dat.test <- dat[!duplicate_names, ]
sum(dat.test$ning_all) 

dat_ning<- dat[dat$ning_all,] 
dim(dat_ning)
head(dat_ning)

#B. Ning abundances----

DBCA_Ning #<- ###Contact DBCA for waiver of data usage as per information in the Data Availability Statement 

head(DBCA_Ning)
Ning_Abun<-DBCA_Ning %>% filter(Level2Class == c("Hard coral")) 
library(reshape)
melted_data.ning <- melt(Ning_Abun, id.vars = c("Site", "Replicate", "Date", "Level3Class", "Level4Class", "LVL4_Count", "Percent_cover"))

#data clean up
melted_data.ning2 <- melted_data.ning %>%
  mutate(Level4Class = gsub("Acroporidae", "Acropora spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Lobophylliidae", "Lobophyllia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Pectiniidae", "Pectinia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Merulinidae", "Merulina spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Pocilloporidae", "Pocillopora spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Poritidae", "Porites spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Euphyllidae", "Euphyllia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Fungiidae", "Fungia spp.", Level4Class)) %>%
  filter(!(Level4Class %in% c("Hard coral")))

(means_per_species_persite.ning <- melted_data.ning2 %>%
  group_by(Level4Class, Site, Replicate) %>% 
  summarize(
    Mean_ = mean(Percent_cover,na.rm = TRUE)))

means_per_species.ning <- melted_data.ning2 %>%
  group_by(Level4Class) %>% 
  summarize(
    Mean_2 = mean(Percent_cover,na.rm = TRUE),
  )
print(means_per_species.ning)

means_per_species.ning.again <- means_per_species_persite.ning %>%
  group_by(Level4Class) %>% 
  summarize(
    Mean_2 = mean(Mean_,na.rm = TRUE),
  )

#C.Adding in the field collected relative abundances to the species list----

# Create a new column with NA values
new_column <- rep(NA, nrow(dat_ning))

# Insert the new column next to the 5th column
dat_ning <- cbind(dat_ning[,1:4], new_column, dat_ning[,5:ncol(dat_ning)])

# Rename the new column to abund_DBCA
colnames(dat_ning)[5] <- "abund_DBCA"

dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Acropora", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Alveopora", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Acanthastrea", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Astrea", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Astreopora", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Catalaphyllia", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Caulastraea", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Coelastrea", 1.124228e-0, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Coeloseris", 5.216049e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Coscinaraea", 5.710403e-04, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Ctenactis", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Cycloseris", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Cyphastrea", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Diploastrea", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Dipsastraea", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Echinophyllia", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Echinopora", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Euphyllia", 8.385310e-03, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Favites", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Fungia", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Galaxea", 8.385310e-03, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Gardineroseris", 5.216049e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Goniastrea", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Goniopora", 3.003522e-0, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Herpolitha", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Hydnophora", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Isopora", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Leptastrea", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Leptoria", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Leptoseris", 5.216049e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Lithophyllon", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Lobactis", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Lobophyllia", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Merulina", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Montipora", 15.8127, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Moseleya", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Mycedium", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Oulophyllia", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Oxypora", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Pachyseris", 3.774939e-03, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Palauastrea", 3.359307e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Pavona", 5.216049e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Pectinia", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Physogyra", 5.019929e-04, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Platygyra", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Plerogyra", 5.120328e-04, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Plesiastrea", 5.120328e-04, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Pleuractis", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Pocillopora", 3.359307e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Podabacia", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Porites", 3.003522e-0, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Sandalolitha", 1.151355e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Scapophyllia", 1.124228e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Psammocora", 2.022482e-03, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Seriatopora", 3.359307e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Stylocoeniella", 3.359307e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Stylophora", 3.359307e-01, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Symphyllia", 7.215818e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Turbinaria", 1.069773e-02, dat_ning$abund_DBCA)
dat_ning$abund_DBCA <- ifelse(dat_ning$genus == "Acanthastrea", 7.215818e-02, dat_ning$abund_DBCA)

dat_ning.test<-dat_ning[,c(1,5)]

as.numeric(dat_ning$abund_DBCA) 

#these values are already percentages
dat_ning$abund_1.5C <- (dat_ning$abund_DBCA -7.22) 
dat_ning$abund_2C <- (dat_ning$abund_DBCA -18.02) 

Abundance_lineplot_Ning<-dat_ning #save in % abundances from 0 to 100

df_long.ning <- pivot_longer(
  dat_ning,
  cols = c("abund_DBCA", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)
df_long.ning.Abun <- df_long.ning %>%
  filter(!(Scenario == "abund_DBCA" & Abund == 0))
df_long.ning.Abun$Scenario <- factor(df_long.ning.Abun$Scenario, levels = c("abund_DBCA", "abund_1C", "abund_1.5C", "abund_2C"))
df_long.ning.Abun$Abund <- ifelse(df_long.ning.Abun$Abund < 0, 0, df_long.ning.Abun$Abund)

#now convert to 1.0 to 0 scale
dat_ning$abund_DBCA.R <- (dat_ning$abund_DBCA) /100
dat_ning$abund_1.5C.R <- (dat_ning$abund_1.5C) /100
dat_ning$abund_2C.R <- (dat_ning$abund_2C) /100

##check the plot now that it is coverted to decmils
df_long.ning <- pivot_longer(
  dat_ning,
  cols = c("abund_DBCA.R", "abund_1.5C.R", "abund_2C.R"),
  names_to = "Scenario",
  values_to = "Abund"
)
df_long.ning.RelAb <- df_long.ning %>%
  filter(!(Scenario == "abund_DBCA.R" & Abund == 0))
df_long.ning.RelAb$Scenario <- factor(df_long.ning.RelAb$Scenario, levels = c("abund_DBCA.R", "abund_1C.R", "abund_1.5C.R", "abund_2C.R"))
df_long.ning.RelAb$Abund <- ifelse(df_long.ning.RelAb$Abund < 0, 0, df_long.ning.RelAb$Abund)

dat_ning2<-dat_ning #just to save the dataset

#change all negatives to zeros
dat_ning2$abund_DBCA.R <- ifelse(dat_ning2$abund_DBCA.R < 0, 0, dat_ning2$abund_DBCA.R)
dat_ning2$abund_1.5C.R <- ifelse(dat_ning2$abund_1.5C.R < 0, 0, dat_ning2$abund_1.5C.R)
dat_ning2$abund_2C.R <- ifelse(dat_ning2$abund_2C.R < 0, 0, dat_ning2$abund_2C.R)

#normalize 
dat_ning2$abund_DBCA.R[dat_ning2$abund_DBCA.R > 0.5 & dat_ning2$abund_DBCA.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_1.5C.R[dat_ning2$abund_1.5C.R > 0.5 & dat_ning2$abund_1.5C.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_2C.R[dat_ning2$abund_2C.R > 0.5 & dat_ning2$abund_2C.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare

dat_ning2$abund_DBCA.R[dat_ning2$abund_DBCA.R > 0.25 & dat_ning2$abund_DBCA.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_1.5C.R[dat_ning2$abund_1.5C.R > 0.25 & dat_ning2$abund_1.5C.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_2C.R[dat_ning2$abund_2C.R > 0.25 & dat_ning2$abund_2C.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare

dat_ning2$abund_DBCA.R[dat_ning2$abund_DBCA.R > 0 & dat_ning2$abund_DBCA.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_1.5C.R[dat_ning2$abund_1.5C.R > 0.0 & dat_ning2$abund_1.5C.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_ning2$abund_2C.R[dat_ning2$abund_2C.R > 0.0 & dat_ning2$abund_2C.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare

####
df_long.ning <- pivot_longer(
  dat_ning2,
  cols = c("abund_DBCA.R", "abund_1.5C.R", "abund_2C.R"),
  names_to = "Scenario",
  values_to = "Abund"
)

#filter out any that were "zero" in the Present day
df_long.ning <- df_long.ning %>%
  filter(!(Scenario == "abund_DBCA.R" & Abund == 0))

df_long.ning$Scenario <- factor(df_long.ning$Scenario, levels = c("abund_DBCA.R", "abund_1C", "abund_1.5C.R", "abund_2C.R"))


#D.Estimate Ecolgoical persistance from normalized abundance, bleaching susceptibility, and range size----

s = 283
dim(dat_ning2)
dat_ning2$abund_2C.R <- ifelse(dat_ning2$abund_2C.R > 0, -1, dat_ning2$abund_2C.R)

#calculate ecological pers
dat_ning_Norm<-dat_ning2
dat_ning_Norm$abund<-dat_ning2$abund_DBCA.R

dat_ning_1.5<-dat_ning2
dat_ning_1.5$abund<-dat_ning2$abund_1.5C.R

dat_ning_2<-dat_ning2
dat_ning_2$abund<-dat_ning2$abund_2C.R

choice_1 <- choice(dat_ning_Norm, s, vars=c("abund","range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_1) <- c("Species", "Ecological persistence Present")
choice_2 <- choice(dat_ning_1.5, s, vars=c("abund","range","bleach"), trait=FALSE)[c("species", "value")]
names(choice_2) <- c("Species", "Ecological persistence 1.5C")
choice_3 <- choice(dat_ning_2, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_3) <- c("Species", "Ecological persistence 2C")

matx <- merge(choice_1, choice_2, all=TRUE)
matx <- merge(matx, choice_3, all=TRUE)
matx$`Ecological persistence Present`[is.na(matx$`Ecological persistence Present`)] <- 0
matx$`Ecological persistence 1.5C`[is.na(matx$`Ecological persistence 1.5C`)] <- 0
matx$`Ecological persistence 2C`[is.na(matx$`Ecological persistence 2C`)] <- 0
matx <- matx[order(matx["Ecological persistence Present"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 1.5C"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 2C"], decreasing=FALSE),]

image(t(as.matrix(matx[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(matx)-1))), labels=matx$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"), cex.axis=0.8, las=2)


#E. Plot Ecological Persistence with ggplot----
mat_ning <- merge(choice_1, choice_2, all=TRUE)
mat_ning <- merge(mat_ning, choice_3, all=TRUE)
mat_ning$Ecol.Persis<-mat_ning$`Ecological persistence`

##Now determine which species are lost with different scenarios and re-run the trait space

df_long_heatmap <- pivot_longer(
  mat_ning,
  cols = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"),
  names_to = "Scenario",
  values_to = "Ecol.Persistance"
)

df_long_heatmap$Scenario <- factor(df_long_heatmap$Scenario, levels = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"))

(Ecol.Pers_Allscenarios<- ggplot(df_long_heatmap, aes(x=Scenario, y = Species, fill = Ecol.Persistance)) +
    geom_tile() + scale_fill_viridis_c())

#MODEL output saved -----
save(mat_ning, file = "Ning_model_pers.RData")

#F. Trait space----

#trait space made on only the corals we have data for for the persistence metric
#run from very top from Madin et al. data in "all species Ninglaoo line"
dim(dat) 
duplicate_names <- duplicated(dat$species)
dat.test <- dat[!duplicate_names, ]
dim(dat.test)

#Present day
Ningaloo_Present #use model output from above to create species list for present that have data needed for this step
dim(Ningaloo_Present) 
Ningaloo_Present$species[!(Ningaloo_Present$species %in% dat.test$species)]
dat$Ningaloo_Present <- dat$species %in% Ningaloo_Present$species
duplicate_names1 <- duplicated(dat$species)
dat.test <- dat[!duplicate_names1, ]
sum(dat.test$Ningaloo_Present) 

Ning_all_1.5 #<- #use model output from above to create species list for present that have data needed for this step
Ning_all_1.5$species[!(Ning_all_1.5$species %in% dat$species)]
dat$Ning_all_1.5 <- dat$species %in% Ning_all_1.5$species
duplicate_names2 <- duplicated(dat$species)
dat.test <- dat[!duplicate_names2, ]
sum(dat.test$Ning_all_1.5) 

#G. brought over f----
#then create dat.test object there

dat<-dat.test 
dat <- dat[!duplicated(dat$species), ] 
dim(dat)

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) #

#H. create convex hull for all data----

TABS<-dat[, 22:23] #make sure correct PC1 and PC2 for Ning.
CHS<- chull(TABS)

#dat <- dat[, -c(22:23)]

Ning_h <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & Ningaloo_Present == 'TRUE',
         select=c('species','goreau_Madin', 'Ningaloo_Present','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#present
tab7.n <- Ning_h[1:89, 4:5] #builders
tab8.n <- Ning_h[90:132, 4:5] #Cementer
tab9.n <- Ning_h[133:147, 4:5] #Cementers & Fillers
tab901.n <- Ning_h[148:193, 4:5] #fillers

Ning_h.15 <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & Ning_all_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'Ning_all_1.5','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#1.5C
tab10.n <- Ning_h.15[1:10, 4:5] #builders
tab11.n <- Ning_h.15[11:24, 4:5] #Cementer
tab12.n <- Ning_h.15[25:32, 4:5] #Cementers & Fillers
tab1201.n <- Ning_h.15[33:73, 4:5] #fillers

#---
Ning_ <- dat%>% 
  subset(Ningaloo_Present == 'TRUE',
         select=c('species','goreau_Madin', 'Ningaloo_Present','PC1','PC2','PC3'))
tab000.n <- Ning_[, 4:5]

ch000.n <- chull(tab000.n)
ch7.n <- chull(tab7.n)
ch8.n <- chull(tab8.n)
ch9.n <- chull(tab9.n)
ch901.n <- chull(tab901.n)

Ning_15 <- dat%>% 
  subset(Ning_all_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'Ning_all_1.5','PC1','PC2','PC3'))
tab000_1.5.n <- Ning_15[, 4:5]

ch000_1.5.n <- chull(tab000_1.5.n)
ch7_1.5.n <- chull(tab10.n)
ch8_1.5.n <- chull(tab11.n)
ch9_1.5.n <- chull(tab12.n)
ch901_1.5.n <- chull(tab1201.n)

#I. Plot PCA for Ningaloo climate scenarios----
plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
Nings_pres <- dat[dat$Ningaloo_Present,] 
Nings_1.5 <- dat[dat$Ning_all_1.5,] 
mtext("Ningaloo under different climate scenarios", line=-1, adj=0, cex=1.5)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("PC1", 1, 2, cex=0.8)
mtext("PC2", 2, 2, cex=0.8)
text(2, -2.5, cex=1.5)
polygon(tab000.n[ch000.n, ], border = "#73c391", col = rgb(115, 195, 145, 128, maxColorValue = 255), lwd = 2) #present day
polygon(tab000_1.5.n[ch000_1.5.n, ], border="#f36e45", col = rgb(243, 110, 69, 100, maxColorValue = 255), lwd=2) #1.5C
polygon(TABS[CHS, ], border="black", lwd=2)
points(PC2 ~ PC1, Nings_pres, col="#73c391", pch=20, cex=1)
points(PC2 ~ PC1, Nings_1.5, col="#f36e45", pch=20, cex=1.5)

#Save MODEL data ----
save(dat, Nings_pres,Nings_1.5, TABS,file="Ningaloo_traitspaceGCB.RData")

#J. Loss in function ----
#Nings_pres
#Nings_1.5
Nings_pres$genus <- as.factor(Nings_pres$genus)

result_df_pres.Ning <- Nings_pres %>%
  group_by(Ningaloo_Present, genus,goreau_Madin) %>% summarise(proportionPres = n())

result_df_1.5.Ning <- Nings_1.5 %>%
  group_by(Ning_all_1.5, genus, goreau_Madin) %>% summarise(proportion1.5 = n())

#-----
combined_df.Ning <- merge(result_df_pres.Ning, result_df_1.5.Ning, by = c("genus", "goreau_Madin"), all = TRUE)

combined_df.Ning <- combined_df.Ning %>%mutate(proportion1.5 = ifelse(is.na(proportion1.5), 0, proportion1.5))

melted_df.Ning <- melt(combined_df.Ning,
                  id.vars = c("genus", "goreau_Madin"),
                  measure.vars = c("proportionPres", "proportion1.5"),
                  variable.name = c("time"),
                  value.name = "value")

#custom_order_LHI <- c("Builders", "Cementers", "Cementers & Fillers", "Fillers", "NA")
melted_df.Ning$goreau_Madin <- ifelse(is.na(melted_df.Ning$goreau_Madin) | melted_df.Ning$goreau_Madin == "NA", "NAs", melted_df.Ning$goreau_Madin)

custom_order_Ning <- c("Fillers", "Cementers & Fillers", "NAs", "Cementers", "Builders")

# Change the order of the factor levels in goreau_Madin
melted_df.Ning$goreau_Madin <- factor(melted_df.Ning$goreau_Madin, levels = custom_order_Ning)

# Calculate percentage change using dplyr
percentage_change_df.Ning <- melted_df.Ning %>%
  group_by(goreau_Madin, variable) %>% summarise(proportion = sum(value))

(percentage_change_result.Ning <- percentage_change_df.Ning %>%
    group_by(goreau_Madin) %>% arrange(variable) %>%
    mutate(percentage_change = (proportion / lag(proportion) - 1) * 100))

#K. Change in Volume LHI ----
# Assuming 'tab000' and 'tab000_1.5' are your data frames containing PC1 and PC2 scores,
# and 'ch000' and 'ch000_1.5' contain the indices

# Function to calculate the area of a convex hull in 2D space
calculate_convex_hull_area <- function(data, indices) {
  subset_points <- data[indices, ]
  n <- nrow(subset_points)
  
  # Duplicate the first point to close the polygon
  subset_points <- rbind(subset_points, subset_points[1, ])
  
  # Calculate the area using the Shoelace formula
  area <- 0.5 * sum(subset_points[1:(n-1), 1] * subset_points[2:n, 2] -
                      subset_points[2:n, 1] * subset_points[1:(n-1), 2])
  
  return(abs(area))
}

# Calculate areas for present day and 1.5C scenarios
area_present_day.Ning <- calculate_convex_hull_area(tab000.n, ch000.n)
area_1.5C.Ning <- calculate_convex_hull_area(tab000_1.5.n, ch000_1.5.n)

# Assuming the convex hull is in 2D, volume is then area in 2D
(volume_present_day.Ning <- area_present_day.Ning) 
(volume_1.5C.Ning <- area_1.5C.Ning) 

# Given values
value1 <- volume_present_day.Ning
value2 <- volume_1.5C.Ning

# Calculate percent difference
(percent_difference <- abs((value1 - value2) / ((value1 + value2) / 2)) * 100)

# Print the result
cat("Percent Difference:", percent_difference, "%\n")

# Calculate the difference in volumes
(volume_difference.Ning <- abs(volume_present_day.Ning - volume_1.5C.Ning))


#***----
#***----
#***----

#Shark Bay Analysis analysis ***----
#A. ---
#Original data McWilliam et al. 2018 but data.csv available within Madin et al. R workspace
#QAQC as per Methods in Quigley & Baird 2024
ctd <- read.csv("data", as.is=TRUE) #

dat <- unique(ctd[c("specie_name", "family_molecules")])
names(dat) <- c("species", "family_molecules")
dat$genus <- sapply(strsplit(dat$species, " "), "[", 1)

# range size (hughes)
rs <- ctd[ctd$trait_name=="Range size" & ctd$resource_id==605, c("specie_name", "value")]
rs$value <- as.numeric(rs$value)
names(rs) <- c("species", "range")
rs <- aggregate(range ~ species, rs, mean)
dat <- merge(dat, rs, all.x=TRUE)
dat$range <- dat$range / max(dat$range, na.rm=TRUE) # normalize

# Local abundance
ab <- ctd[ctd$trait_name=="Abundance GBR", c("specie_name", "value")]
names(ab) <- c("species", "abund")
dat <- merge(dat, ab, all.y=TRUE)
dat$abund[dat$abund=="common"] <- 1  # normalize
dat$abund[dat$abund=="uncommon"] <- 0.5 # normalize
dat$abund[dat$abund=="rare"] <- 0.25 # normalize
dat$abund <- as.numeric(dat$abund)

# Bleaching susceptibility
names(bri) <- c("species", "BRI")
BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$bleach <- bri$BRI[match(dat$species, bri$species)]
dat$bleach[is.na(dat$bleach)] <- BRI_genus[is.na(dat$bleach)]
dat$bleach <- 1 - (dat$bleach / 100) # normalize

# Growth form typical
gf <- ctd[ctd$trait_name=="Growth form typical", c("specie_name", "value")]
names(gf) <- c("species", "growth_form")
dat <- merge(dat, gf, all.x=TRUE)

# Life history
lh <- ctd[ctd$trait_name=="Life history strategy", c("specie_name", "value")]
names(lh) <- c("species", "life_history")
lh$life_history[lh$life_history=="weedy"] <- "Weedy"
lh$life_history[lh$life_history=="generalist"] <- "Generalist"
lh$life_history[lh$life_history=="competitive"] <- "Competitive"
lh$life_history[lh$life_history=="stress-tolerant"] <- "Stress-tolerant"
dat <- merge(dat, lh, all.x=TRUE) #

# Traits - McWilliam et al. 2018 PNAS
dim(tra[tra$domain=="pacific",])
tra <- tra[,c("species", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat <- merge(dat, tra) 

# Goreau Madin ellipses
lh <- ctd[ctd$location_name %in% c("Great Barrier Reef (GBR)", "Indo-Pacific (unspecified)"), c("specie_name", "goreau_Madin")]
names(lh) <- c("species", "goreau_Madin")
lh$goreau_Madin[lh$goreau_Madin=="builders"] <- "Builders"
lh$goreau_Madin[lh$goreau_Madin=="cementers"] <- "Cementers"
lh$goreau_Madin[lh$goreau_Madin=="fillers"] <- "Fillers"
lh$goreau_Madin[lh$goreau_Madin=="cementers&fillers"] <- "Cementers & Fillers"
dat <- merge(dat, lh, all.x=TRUE) #

dat$goreau <- NA
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_spacesize %in% c(4, 5)] <- "Fillers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(1, 2)] <- "Cementers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(3, 4, 5) & dat$cat_SA_vol  %in% c(1, 2)] <- "Builders"

# Restoration potential from PLoS paper
dat$restore <- 0
dat$restore[dat$growth_form %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$growth_form %in% c("massive", "submassive")] <- 5
dat$restore[dat$growth_form %in% c("laminar")] <- 4
dat$restore[dat$growth_form %in% c("columnar", "encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$growth_form %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$growth_form %in% c("hispidose")] <- 1
dat$restore <- (dat$restore) / 6 # normalize
dim(dat)

dat <- dat[!duplicated(dat$species), ] 
dim(dat )

#data prep----

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) 

## all species Shark Bay----
shark_all<- read.csv("sharkbay_all.csv", as.is=TRUE)
shark_all$species[!(shark_all$species %in% dat$species)]
dat$shark_all <- dat$species %in% shark_all$species
sum(dat$shark_all) 

dat_shark<- dat[dat$shark_all,] 
dim(dat_shark) 

#need to make sure that dat has the PCs 

#B. Shark Bay abundance ----

DBCA_Shark #<- ###Contact DBCA for waiver of data usage as per information in the Data Availability Statement 
head(DBCA_Shark)
Shark_Abun<-DBCA_Shark %>% filter(Level2Class == c("Hard coral")) #don't need "CORAL"
#library(reshape)
melted_data.shark <- melt(Shark_Abun, id.vars = c("Site", "Replicate", "Date", "Level3Class", "Level4Class", "LVL4_Count", "Percent_cover"))
print(melted_data.shark)

#data clean up
melted_data.shark2 <- melted_data.shark %>%
  mutate(Level4Class = gsub("Acroporidae", "Acropora spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Lobophylliidae", "Lobophyllia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Pectiniidae", "Pectinia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Merulinidae", "Merulina spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Pocilloporidae", "Pocillopora spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Poritidae", "Porites spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Euphyllidae", "Euphyllia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Fungiidae", "Fungia spp.", Level4Class)) %>%
  mutate(Level4Class = gsub("Faviidae", "Favia spp.", Level4Class)) %>%
  filter(!(Level4Class %in% c("Hard coral")))

(means_per_species_persite.shark <- melted_data.shark2 %>%
    group_by(Level4Class, Site, Replicate) %>% 
    summarize(
      Mean_ = mean(Percent_cover,na.rm = TRUE)))

# View the resulting dataframe with means
print(means_per_species_persite.shark)

means_per_species.shark <- melted_data.shark2 %>%
  group_by(Level4Class) %>% 
  summarize(
    Mean_2 = mean(Percent_cover,na.rm = TRUE),
  )
print(means_per_species.shark)

#C.Adding in the field collected relative abundances to the species list----

# Create a new column with NA values
new_column <- rep(NA, nrow(dat_shark))

# Insert the new column next to the 5th column
dat_shark2 <- cbind(dat_shark[,1:4], new_column, dat_shark[,5:ncol(dat_shark)])

# Rename the new column to abund_DBCA
colnames(dat_shark2)[5] <- "abund_DBCA"

dat_shark3 <- dat_shark2 %>%
  mutate(abund_DBCA = case_when(
    genus == "Acropora" ~  + 4.457079631,
    family_molecules == "Fungiidae" ~ + 0.00412,
    genus == "Cyphastrea" ~  + 0.0323,
    genus == "Dipsastraea" ~  + 0.0267,
    genus == "Favites" ~  + 0.163,
    genus == "Goniastrea" ~  + 0.0512,
    genus == "Hydnophora" ~  + 0.00669,
    genus == "Montipora" ~  + 0.773,
    genus == "Oulophyllia" ~  + 0.0184,
    genus == "Pavona" ~  + 0.0503,
    genus == "Platygyra" ~  + 0.128,
    genus == "Pocillopora" ~  + 0.234,
    genus == "Porites" ~  + 0.318,
    genus == "Seriatopora" ~  + 0.00412,
    genus == "Stylophora" ~  + 0.0199,
    genus == "Turbinaria" ~  + 12.0,
    TRUE ~ NA_real_
  ))

dat_shark4 <- dat_shark3 %>%
  mutate(abund_DBCA = case_when(
    family_molecules == "Acroporidae" ~  + 4.46, #iso, Alveo, Astreo
    family_molecules == "Lobophyllidae" ~ + 0.0323, #realted to Mussidae  and Merulinidae (so Cyphastrea value)
    family_molecules == "Merulinidae" ~  + 0.0323, #in Merulinidae (so Cyphastrea value)
    family_molecules == "Fungiidae" ~  + 0.00412, 
    family_molecules == "Agariciidae" ~  + 0.0503, #gets Pavona value bc same family
    family_molecules == "Poritidae" ~  + 0.318,
    family_molecules == "Pachyseridae" ~  + 0.0503, #should say Agariciidae as family, gets Pavona value bc same family
    family_molecules == "Pocilloporidae" ~  + 0.234,
    genus == "Turbinaria" ~ + 12.0,
    #family_molecules == "Plesiastreidae" ~  + #genus and monophylotic family
    #family_molecules == "Psammocoridae" ~  + #genus and monophylotic family
    #family_molecules == "Siderastreidae" ~  + #genus and monophylotic family
    #family_molecules == "Coscinaraeidae" ~  + #genus and monophylotic family
    TRUE ~ NA_real_
  ))
##some species are still NA

Abundance_lineplot.Shark<-dat_shark4 #this is in % from 100

#as.numeric(dat_shark4$abund_DBCA) 

dat_shark4$abund_1.5C <- (dat_shark4$abund_DBCA -7.7) 
dat_shark4$abund_2C <- (dat_shark4$abund_DBCA -18.02) 

df_long.sh.perct.abun <- pivot_longer(
  dat_shark4,
  cols = c("abund_DBCA", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)
df_long.sh.perct.abun.Abun <- df_long.sh.perct.abun %>%
  filter(!(Scenario == "abund_DBCA" & Abund == 0))
df_long.sh.perct.abun.Abun$Scenario <- factor(df_long.sh.perct.abun.Abun$Scenario, levels = c("abund_DBCA", "abund_1C", "abund_1.5C", "abund_2C"))
df_long.sh.perct.abun.Abun$Abund <- ifelse(df_long.sh.perct.abun.Abun$Abund < 0, 0, df_long.sh.perct.abun.Abun$Abund)

#now convert to 1.0 to 0 scale
dat_shark4$abund_DBCA.R <- (dat_shark4$abund_DBCA) /100
dat_shark4$abund_1.5C.R <- (dat_shark4$abund_1.5C) /100
dat_shark4$abund_2C.R <- (dat_shark4$abund_2C) /100

##check the plot now that it is coverted to decmils
df_long.shark <- pivot_longer(
  dat_shark4,
  cols = c("abund_DBCA.R", "abund_1.5C.R", "abund_2C.R"),
  names_to = "Scenario",
  values_to = "Abund"
)
df_long.Shark.RelAb <- df_long.shark %>%
  filter(!(Scenario == "abund_DBCA.R" & Abund == 0))
df_long.Shark.RelAb$Scenario <- factor(df_long.Shark.RelAb$Scenario, levels = c("abund_DBCA.R", "abund_1C.R", "abund_1.5C.R", "abund_2C.R"))
df_long.Shark.RelAb$Abund <- ifelse(df_long.Shark.RelAb$Abund < 0, 0, df_long.Shark.RelAb$Abund)

#need to covert these to relative abundances----

#change all negatives to zeros
dat_shark4$abund_DBCA.R <- ifelse(dat_shark4$abund_DBCA.R < 0, 0, dat_shark4$abund_DBCA.R)
dat_shark4$abund_1.5C.R <- ifelse(dat_shark4$abund_1.5C.R < 0, 0, dat_shark4$abund_1.5C.R)
dat_shark4$abund_2C.R <- ifelse(dat_shark4$abund_2C.R < 0, 0, dat_shark4$abund_2C.R)

#normalize 
dat_shark4$abund_DBCA.R[dat_shark4$abund_DBCA.R > 0.5 & dat_shark4$abund_DBCA.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_1.5C.R[dat_shark4$abund_1.5C.R > 0.5 & dat_shark4$abund_1.5C.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_2C.R[dat_shark4$abund_2C.R > 0.5 & dat_shark4$abund_2C.R <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare

dat_shark4$abund_DBCA.R[dat_shark4$abund_DBCA.R > 0.25 & dat_shark4$abund_DBCA.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_1.5C.R[dat_shark4$abund_1.5C.R > 0.25 & dat_shark4$abund_1.5C.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_2C.R[dat_shark4$abund_2C.R > 0.25 & dat_shark4$abund_2C.R <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare

dat_shark4$abund_DBCA.R[dat_shark4$abund_DBCA.R > 0 & dat_shark4$abund_DBCA.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_1.5C.R[dat_shark4$abund_1.5C.R > 0.0 & dat_shark4$abund_1.5C.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
dat_shark4$abund_2C.R[dat_shark4$abund_2C.R > 0.0 & dat_shark4$abund_2C.R <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare

####

df_long.shark <- pivot_longer(
  dat_shark4,
  cols = c("abund_DBCA.R", "abund_1.5C.R", "abund_2C.R"),
  names_to = "Scenario",
  values_to = "Abund"
)

#filter out any that were "zero" in the Present day
#
df_long.shark <- df_long.shark %>%
  filter(!(Scenario == "abund_DBCA.R" & Abund == 0))

df_long.shark$Scenario <- factor(df_long.shark$Scenario, levels = c("abund_DBCA.R", "abund_1.5C.R", "abund_2C.R"))

#D.Estimate Ecolgoical persistance from normalized abundance, bleaching susceptibility, and range size----
s = 251
dim(dat_shark4)
dat_shark4$abund_2C.R <- ifelse(dat_shark4$abund_2C.R > 0, -1, dat_shark4$abund_2C.R)

#calculate ecological persistance Madin et al. 2023
#input to equation must be called, needs to be abund called for input to work
dat_shark4_Norm<-dat_shark4
dat_shark4_Norm$abund<-dat_shark4$abund_DBCA.R

dat_shark4_1.5<-dat_shark4
dat_shark4_1.5$abund<-dat_shark4$abund_1.5C.R

dat_shark4_2<-dat_shark4
dat_shark4_2$abund<-dat_shark4$abund_2C.R

choice_1 <- choice(dat_shark4_Norm, s, vars=c("abund","range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_1) <- c("Species", "Ecological persistence Present")
choice_2 <- choice(dat_shark4_1.5, s, vars=c("abund","range","bleach"), trait=FALSE)[c("species", "value")]
names(choice_2) <- c("Species", "Ecological persistence 1.5C")
choice_3 <- choice(dat_shark4_2, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_3) <- c("Species", "Ecological persistence 2C")

matx <- merge(choice_1, choice_2, all=TRUE)
matx <- merge(matx, choice_3, all=TRUE)
matx$`Ecological persistence Present`[is.na(matx$`Ecological persistence Present`)] <- 0
matx$`Ecological persistence 1.5C`[is.na(matx$`Ecological persistence 1.5C`)] <- 0
matx$`Ecological persistence 2C`[is.na(matx$`Ecological persistence 2C`)] <- 0
matx <- matx[order(matx["Ecological persistence Present"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 1.5C"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 2C"], decreasing=FALSE),]

image(t(as.matrix(matx[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(matx)-1))), labels=matx$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"), cex.axis=0.8, las=2)


#E. Plot Ecological Persistence with ggplot----
mat_shark <- merge(choice_1, choice_2, all=TRUE)
mat_shark <- merge(mat_shark, choice_3, all=TRUE)
mat_shark$Ecol.Persis<-mat_shark$`Ecological persistence`

#Now determine which species are lost with different scenarios and re-run the trait space

df_long_heatmap.shark <- pivot_longer(
  mat_shark,
  cols = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"),
  names_to = "Scenario",
  values_to = "Ecol.Persistance"
)

df_long_heatmap.shark$Scenario <- factor(df_long_heatmap.shark$Scenario, levels = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"))

(Ecol.Pers_Allscenarios.shark<- ggplot(df_long_heatmap.shark, aes(x=Scenario, y = Species, fill = Ecol.Persistance)) +
    geom_tile() + scale_fill_viridis_c())

#MODEL output saved -----
save(mat_shark, file = "SB_model_pers.RData")

#F. Trait space----

#trait space made on only the corals we have data for for the persistence metric
#run from very top from Madin et al. data
dim(dat) 
duplicate_names <- duplicated(dat$species)
dat.test <- dat[!duplicate_names, ]
dim(dat.test) 

#Present day
shark_Present #use model output from above to create species list for present that have data needed for this step
dim(shark_Present) 
shark_Present$species[!(shark_Present$species %in% dat.test$species)]
dat.test$shark_Present <- dat.test$species %in% shark_Present$species
sum(dat.test$shark_Present) 

#add 1.5C here
shark_1.5 #use model output from above to create species list for present that have data needed for this step
dim(shark_1.5) 
shark_1.5$species[!(shark_1.5$species %in% dat.test$species)]
dat.test$shark_1.5 <- dat.test$species %in% shark_1.5$species
sum(dat.test$shark_1.5) 

#G. brought over f----
#go back up to top of data_prep.R and run from the very top
#then create dat.test object there

dat<-dat.test #
dat <- dat[!duplicated(dat$species), ] #

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) 

#H. create convex hull for all data----
TABS<-dat[, 21:22] 
CHS<- chull(TABS)

Sh_h <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & shark_Present == 'TRUE',
         select=c('species','goreau_Madin', 'shark_Present','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#present
tab7.s <- Sh_h[1:75, 4:5] #builders
tab8.s <- Sh_h[76:113, 4:5] #Cementer
tab9.s <- Sh_h[114:128, 4:5] #Cementers & Fillers
tab901.s <- Sh_h[129:173, 4:5] #fillers

Sh_h.15 <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & shark_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'shark_1.5','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#1.5C
tab10.s <- Sh_h.15[1, 4:5] #builders
tab11.s <- Sh_h.15[2:7, 4:5] #Cementer

#---

Shark_ <- dat%>% 
  subset(shark_Present == 'TRUE',
         select=c('species','goreau_Madin', 'shark_Present','PC1','PC2','PC3'))
tab000.s <- Shark_[, 4:5]

ch000.s <- chull(tab000.s)
ch7.s <- chull(tab7.s)
ch8.s <- chull(tab8.s)
ch9.s <- chull(tab9.s)
ch901.s <- chull(tab901.s)

Shark_15 <- dat%>% 
  subset(shark_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'shark_1.5','PC1','PC2','PC3'))
tab000_1.5.s <- Shark_15[, 4:5]

ch000_1.5.s <- chull(tab000_1.5.s)
ch7_1.5.s <- chull(tab10.s)
ch8_1.5.s <- chull(tab11.s)

#I. Plot PCA for Shark Bay climate scenarios----
plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
Shark_pres <- dat[dat$shark_Present,] 
Shark_1.5 <- dat[dat$shark_1.5,] 
mtext("Shark Bay under different climate scenarios", line=-1, adj=0, cex=1.5)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("PC1", 1, 2, cex=0.8)
mtext("PC2", 2, 2, cex=0.8)
text(2, -2.5, "Present: n=222, 1.5C: n=8", cex=1.5)
points(PC2 ~ PC1, Shark_pres, col="#73c391", pch=20, cex=1)
points(PC2 ~ PC1, Shark_1.5, col="#f36e45", pch=20, cex=1.5)
polygon(tab000.s[ch000.s, ], border = "#73c391", col = rgb(115, 195, 145, 128, maxColorValue = 255), lwd = 2) #present day
polygon(tab000_1.5.s[ch000_1.5.s, ], border="#f36e45", col = rgb(243, 110, 69, 100, maxColorValue = 255), lwd=2) #1.5C
polygon(TABS[CHS, ], border="black", lwd=2)

#save MODEL data----
#save(dat, Shark_pres,Shark_1.5, TABS,file="Shark_traitspaceGCB.RData")

#J. Loss in function ----
#Shark_pres
#Shark_1.5
Shark_pres$genus <- as.factor(Shark_pres$genus)

result_df_pres.SB <- Shark_pres %>%
  group_by(shark_Present, genus,goreau_Madin) %>% summarise(proportionPres = n())

result_df_1.5.SB <- Shark_1.5 %>%
  group_by(shark_1.5, genus, goreau_Madin) %>% summarise(proportion1.5 = n())

#-----
combined_df.SB <- merge(result_df_pres.SB, result_df_1.5.SB, by = c("genus", "goreau_Madin"), all = TRUE)

combined_df.SB <- combined_df.SB %>%mutate(proportion1.5 = ifelse(is.na(proportion1.5), 0, proportion1.5))

melted_df.SB <- melt(combined_df.SB,
                       id.vars = c("genus", "goreau_Madin"),
                       measure.vars = c("proportionPres", "proportion1.5"),
                       variable.name = c("time"),
                       value.name = "value")


melted_df.SB$goreau_Madin <- ifelse(is.na(melted_df.SB$goreau_Madin) | melted_df.SB$goreau_Madin == "NA", "NAs", melted_df.SB$goreau_Madin)

custom_order_SB <- c("Cementers", "NAs","Builders", "Cementers & Fillers","Fillers")

# Change the order of the factor levels in goreau_Madin
melted_df.SB$goreau_Madin <- factor(melted_df.SB$goreau_Madin, levels = custom_order_SB)

# Calculate percentage change using dplyr
percentage_change_df.SB <- melted_df.SB %>%
  group_by(goreau_Madin, variable) %>% summarise(proportion = sum(value))

(percentage_change_result.SB <- percentage_change_df.SB %>%
    group_by(goreau_Madin) %>% arrange(variable) %>%
    mutate(percentage_change = (proportion / lag(proportion) - 1) * 100))

#K. Change in Volume SB ----
# Assuming 'tab000' and 'tab000_1.5' are your data frames containing PC1 and PC2 scores,
# and 'ch000' and 'ch000_1.5' contain the indices

# Function to calculate the area of a convex hull in 2D space
calculate_convex_hull_area <- function(data, indices) {
  subset_points <- data[indices, ]
  n <- nrow(subset_points)
  
  # Duplicate the first point to close the polygon
  subset_points <- rbind(subset_points, subset_points[1, ])
  
  # Calculate the area using the Shoelace formula
  area <- 0.5 * sum(subset_points[1:(n-1), 1] * subset_points[2:n, 2] -
                      subset_points[2:n, 1] * subset_points[1:(n-1), 2])
  
  return(abs(area))
}

# Calculate areas for present day and 1.5C scenarios
area_present_day.SB <- calculate_convex_hull_area(tab000.s, ch000.s)
area_1.5C.SB <- calculate_convex_hull_area(tab000_1.5.s, ch000_1.5.s)

# Assuming the convex hull is in 2D, volume is then area in 2D
(volume_present_day.SB <- area_present_day.SB) 
(volume_1.5C.SB <- area_1.5C.SB) 

# Given values
value1 <- volume_present_day.SB
value2 <- volume_1.5C.SB

# Calculate percent difference
(percent_difference <- abs((value1 - value2) / ((value1 + value2) / 2)) * 100)

# Print the result
cat("Percent Difference:", percent_difference, "%\n")

# Calculate the difference in volumes
(volume_difference.SB <- abs(volume_present_day.SB - volume_1.5C.SB))


#***----
#***----
#***----

#Great Barrier Reef ----
#A. ---
#Original data McWilliam et al. 2018 in data.csv available within Madin et al. R workspace
#QAQC as per Methods in Quigley & Baird 2024
ctd <- read.csv("data", as.is=TRUE) 

dat <- unique(ctd[c("specie_name", "family_molecules")])
names(dat) <- c("species", "family_molecules")
dat$genus <- sapply(strsplit(dat$species, " "), "[", 1)

# range size (Hughes et al. originally as outlined in Madin et al. )
rs <- ctd[ctd$trait_name=="Range size" & ctd$resource_id==605, c("specie_name", "value")]
rs$value <- as.numeric(rs$value)
names(rs) <- c("species", "range")
rs <- aggregate(range ~ species, rs, mean)
dat <- merge(dat, rs, all.x=TRUE)
dat$range <- dat$range / max(dat$range, na.rm=TRUE) # normalize

# Local abundance
ab <- ctd[ctd$trait_name=="Abundance GBR", c("specie_name", "value")]
names(ab) <- c("species", "abund")
dat <- merge(dat, ab, all.y=TRUE)
dat$abund[dat$abund=="common"] <- 1  # normalize
dat$abund[dat$abund=="uncommon"] <- 0.5 # normalize
dat$abund[dat$abund=="rare"] <- 0.25 # normalize
dat$abund <- as.numeric(dat$abund)

# Bleaching susceptibility
names(bri) <- c("species", "BRI")
BRI_genus <- bri$BRI[match(dat$genus, bri$species)]
dat$bleach <- bri$BRI[match(dat$species, bri$species)]
dat$bleach[is.na(dat$bleach)] <- BRI_genus[is.na(dat$bleach)]
dat$bleach <- 1 - (dat$bleach / 100) # normalize

# Growth form typical
gf <- ctd[ctd$trait_name=="Growth form typical", c("specie_name", "value")]
names(gf) <- c("species", "growth_form")
dat <- merge(dat, gf, all.x=TRUE)

# Life history
lh <- ctd[ctd$trait_name=="Life history strategy", c("specie_name", "value")]
names(lh) <- c("species", "life_history")
lh$life_history[lh$life_history=="weedy"] <- "Weedy"
lh$life_history[lh$life_history=="generalist"] <- "Generalist"
lh$life_history[lh$life_history=="competitive"] <- "Competitive"
lh$life_history[lh$life_history=="stress-tolerant"] <- "Stress-tolerant"
dat <- merge(dat, lh, all.x=TRUE) 

# Traits - McWilliam et al. 2018 PNAS as outlined in Madin et al. 
dim(tra[tra$domain=="pacific",])
tra <- tra[,c("species", "cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")]
dat <- merge(dat, tra) 

# Goreau Madin ellipses
lh <- ctd[ctd$location_name %in% c("Great Barrier Reef (GBR)", "Indo-Pacific (unspecified)"), c("specie_name", "goreau_Madin")]
names(lh) <- c("species", "goreau_Madin")
lh$goreau_Madin[lh$goreau_Madin=="builders"] <- "Builders"
lh$goreau_Madin[lh$goreau_Madin=="cementers"] <- "Cementers"
lh$goreau_Madin[lh$goreau_Madin=="fillers"] <- "Fillers"
lh$goreau_Madin[lh$goreau_Madin=="cementers&fillers"] <- "Cementers & Fillers"
dat <- merge(dat, lh, all.x=TRUE) 
dim(dat)
dat$goreau <- NA
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_spacesize %in% c(4, 5)] <- "Fillers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(1, 2)] <- "Cementers"
dat$goreau[dat$cat_colonydiameter %in% c(4, 5) & dat$cat_colonyheight %in% c(3, 4, 5) & dat$cat_SA_vol  %in% c(1, 2)] <- "Builders"

# Restoration potential from PLoS paper as outlined in Madin et al. 
dat$restore <- 0
dat$restore[dat$growth_form %in% c("branching_open", "branching_closed")] <- 6
dat$restore[dat$growth_form %in% c("massive", "submassive")] <- 5
dat$restore[dat$growth_form %in% c("laminar")] <- 4
dat$restore[dat$growth_form %in% c("columnar", "encrusting", "encrusting_long_uprights")] <- 3
dat$restore[dat$growth_form %in% c("digitate", "corymbose", "tables_or_plates")] <- 2
dat$restore[dat$growth_form %in% c("hispidose")] <- 1
dat$restore <- (dat$restore) / 6 # normalize
dim(dat)

dat <- dat[!duplicated(dat$species), ] 
dim(dat )

#stop here if running for second part at G and F

#data prep----

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) 


#B. GBR abundances for whole GBR----
GBR_Baird_Abun_all <- read.csv("LIT_GBR_2013-16_shallow_02042020_wLatitude.csv", as.is=TRUE)

GBR_Baird_Abun_all <- GBR_Baird_Abun_all %>%
  mutate_at(vars(starts_with("intercept_")), as.numeric)

class(GBR_Baird_Abun_all$intercept_8)

result <- GBR_Baird_Abun_all %>%
  group_by(site, transect, reef) %>%
  summarise(across(starts_with("intercept_"), 
                   list(sum_column_name = ~sum(., na.rm = TRUE)), 
                   .names = "{col}_sum"))

merged_data <- merge(GBR_Baird_Abun_all, result, by = c("site", "transect", "reef"), all.x = TRUE)

relative_percentages <- merged_data %>%
  group_by(reef, transect) %>%
  rowwise() %>%
  mutate(relative_percentage_intercept_1 = intercept_1 / intercept_1_sum * 100) %>%
  mutate(relative_percentage_intercept_2 = intercept_2 / intercept_2_sum * 100) %>%
  mutate(relative_percentage_intercept_3 = intercept_3 / intercept_3_sum * 100) %>%
  mutate(relative_percentage_intercept_4 = intercept_4 / intercept_4_sum * 100) %>%
  mutate(relative_percentage_intercept_5 = intercept_5 / intercept_5_sum * 100) %>%
  mutate(relative_percentage_intercept_6 = intercept_6 / intercept_6_sum * 100) %>%
  mutate(relative_percentage_intercept_7 = intercept_7 / intercept_7_sum * 100) %>%
  mutate(relative_percentage_intercept_8 = intercept_8 / intercept_8_sum * 100) %>%
  mutate(relative_percentage_intercept_9 = intercept_9 / intercept_9_sum * 100) %>%
  mutate(relative_percentage_intercept_10 = intercept_10 / intercept_10_sum * 100) %>%
  mutate(relative_percentage_intercept_11 = intercept_11 / intercept_11_sum * 100) %>%
  ungroup() 

subset_merged_data <- relative_percentages[, c(1:3,9, 35:45) ]
subset_merged_data[, 5:15] <- lapply(subset_merged_data[, 5:15], as.numeric)
mean_values <- rowMeans(subset_merged_data[, 5:15], na.rm = TRUE)

subset_merged_data$mean_values <- mean_values

cleanGBR<-subset_merged_data[, c(1:4,16) ]

#library(reshape)

means_per_species_persite <- cleanGBR %>%
  group_by(species, site) %>% 
  summarize(
    Mean_ = mean(mean_values,na.rm = TRUE))

means_per_species_persite %>% #filter(dat_lhi.genus == "Acropora") %>% 
  ggplot(aes(y = Mean_, x = site, fill=species))+ 
  geom_bar(stat = "identity", position="fill")

means_per_species <- cleanGBR %>%
  group_by(species) %>% 
  summarize(
    Mean_ = mean(mean_values,na.rm = TRUE))

##Mean abundances now

#from data_prep.R. Need to run the dat from the top and gbr_all to get below:
# Great Barrier Reef----
## all spp
gbr_all ##Download data from Kuo 2017 and Kuo et al. 2023 from references
dim(gbr_all)
gbr_all$species[!(gbr_all$species %in% dat$species)]
dat$gbr_all <- dat$species %in% gbr_all$species
sum(dat$gbr_all) 
duplicate_names <- duplicated(dat$species)
# Remove rows with duplicate names
dat.test <- dat[!duplicate_names, ]
sum(dat.test$gbr_all) 

#trait space made on only the corals we have data for for the persistence metric
#run from very top from Madin et al. data
dim(dat) 
duplicate_names <- duplicated(dat$species)
dat.test <- dat[!duplicate_names, ]
dim(dat.test) 
##end data prep ----

dim(gbr_all) 
dim(means_per_species) 
merged_df_gbr <- merge(gbr_all, means_per_species, by = "species", all = TRUE) #
dim(merged_df_gbr) 
# Identify species present in means_per_species but not in gbr_all
species_to_drop <- setdiff(means_per_species$species, gbr_all$species)

# Filter out rows with species present in means_per_species but not in gbr_all
merged_df_gbr <- merged_df_gbr[!merged_df_gbr$species %in% species_to_drop, ]
dim(merged_df_gbr) #

# Split the 'species_genus' column by space
split_names <- strsplit((merged_df_gbr$species), " ")

# Extract genus and species from the split_names list
genus <- sapply(split_names, `[`, 1)
species <- sapply(split_names, `[`, 2)

# Create a new data frame with species and genus columns
new_data <- data.frame(species = species, genus = genus)

gbr_all_final <- cbind(merged_df_gbr, new_data)

colnames(gbr_all_final)[3] <- "species1"
head(gbr_all_final)

Genus_means <- gbr_all_final %>%
  group_by(genus) %>% 
  summarize(
    Mean_ = mean(Mean_,na.rm = TRUE))

merged_df_gbr2 <- merge(gbr_all_final, Genus_means, by = "genus", all = TRUE) #
merged_df_gbr2$mergeColumn <- ifelse(is.na(merged_df_gbr2$Mean_.x), merged_df_gbr2$Mean_.y, merged_df_gbr2$Mean_.x)

dat_small<-dat[,c(1:3)]
merged_df_gbr3 <- merge(merged_df_gbr2, dat_small, by = "species", all = TRUE) #

species_to_drop2 <- setdiff(merged_df_gbr3$species, merged_df_gbr2$species)
# Filter out rows with species present in means_per_species but not in gbr_all
merged_df_gbr_final <- merged_df_gbr3[!merged_df_gbr3$species %in% species_to_drop2, ]
dim(merged_df_gbr_final) 

as.data.frame(family_means <- merged_df_gbr_final %>%
  group_by(family_molecules) %>% 
  summarize(
    Mean_fam = mean(Mean_.x,na.rm = TRUE)))

merged_df_gbr4 <- merge(merged_df_gbr_final, family_means, by = "family_molecules", all = TRUE) #
dim(merged_df_gbr4) 

#mean of family
merged_df_gbr4$mergeColumn2 <- ifelse(is.na(merged_df_gbr4$mergeColumn), merged_df_gbr4$Mean_fam, merged_df_gbr4$mergeColumn)

merged_df_gbr4 %>%
  #filter(gbr_all == "TRUE") %>%
  ggplot(aes(y = mergeColumn2, x = reorder(species, -mergeColumn2), fill = mergeColumn2)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

na_count <- sum(is.na(merged_df_gbr4$Mean_.x))
sum(is.na(merged_df_gbr4$mergeColumn)) 
sum(is.na(merged_df_gbr4$mergeColumn2)) 

#C.change in relative abundances under diff climate scenarios----

# Add conditional logic to use "abund" column if "Mean_2" is missing or within [0, 0.25]
#merged_df_gbr.single$ <- ifelse(!is.na(merged_df_gbr.single$Mean_2), merged_df_gbr.single$Mean_2, merged_df_gbr.single$abund)

#min abundance of family-----
#as.data.frame(family_means <- merged_df_gbr_final %>%
#                group_by(family_molecules) %>% 
#                summarize(
#                  Min_fam = min(Mean_.x,na.rm = TRUE)))

#merged_df_gbr4 <- merge(merged_df_gbr_final, family_means, by = "family_molecules", all = TRUE) #
#dim(merged_df_gbr4) #402

#merged_df_gbr4$mergeColumn2 <- ifelse(is.na(merged_df_gbr4$mergeColumn), merged_df_gbr4$Min_fam, merged_df_gbr4$mergeColumn)

#merged_df_gbr4 %>%
  #filter(gbr_all == "TRUE") %>%
#  ggplot(aes(y = mergeColumn2, x = reorder(species, -mergeColumn2), fill = mergeColumn2)) +
#  geom_bar(stat = "identity") +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#sum(is.na(merged_df_gbr4$mergeColumn2))
#sum(is.infinite(merged_df_gbr4$mergeColumn2)) #6
#merged_df_gbr4$mergeColumn2 <- replace(merged_df_gbr4$mergeColumn2, is.infinite(merged_df_gbr4$mergeColumn2), NA)
#end min abundance of genus end-----

###restart below

#mean of the genus
merged_df_gbr4$Present <- (merged_df_gbr4$mergeColumn2) / 100
merged_df_gbr4$abund_1.5C <- (merged_df_gbr4$mergeColumn2 -6.74) /100
merged_df_gbr4$abund_2C <- (merged_df_gbr4$mergeColumn2 -20.18) /100

Abundance_lineplot<-merged_df_gbr4 #for lineplot below

#change all negatives to zeros
merged_df_gbr4$Present <- ifelse(merged_df_gbr4$Present < 0, 0, merged_df_gbr4$Present)
merged_df_gbr4$abund_1.5C <- ifelse(merged_df_gbr4$abund_1.5C < 0, 0, merged_df_gbr4$abund_1.5C)
merged_df_gbr4$abund_2C <- ifelse(merged_df_gbr4$abund_2C < 0, 0, merged_df_gbr4$abund_2C)

#normalize 
merged_df_gbr4$Present[merged_df_gbr4$Present > 0.5 & merged_df_gbr4$Present <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_1.5C[merged_df_gbr4$abund_1.5C > 0.5 & merged_df_gbr4$abund_1.5C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_2C[merged_df_gbr4$abund_2C > 0.5 & merged_df_gbr4$abund_2C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare

merged_df_gbr4$Present[merged_df_gbr4$Present > 0.25 & merged_df_gbr4$Present <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_1.5C[merged_df_gbr4$abund_1.5C > 0.25 & merged_df_gbr4$abund_1.5C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_2C[merged_df_gbr4$abund_2C > 0.25 & merged_df_gbr4$abund_2C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare

merged_df_gbr4$Present[merged_df_gbr4$Present > 0 & merged_df_gbr4$Present <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_1.5C[merged_df_gbr4$abund_1.5C > 0.0 & merged_df_gbr4$abund_1.5C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
merged_df_gbr4$abund_2C[merged_df_gbr4$abund_2C > 0.0 & merged_df_gbr4$abund_2C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare

####

df_long.GBR_all <- pivot_longer(
  merged_df_gbr4,
  cols = c("Present", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

#filter out any that were "zero" in the Present day
#n =
df_long.GBR_all <- df_long.GBR_all %>%
  filter(!(Scenario == "Present" & Abund == 0))

df_long.GBR_all$Scenario <- factor(df_long.GBR_all$Scenario, levels = c("Present", "abund_1C", "abund_1.5C", "abund_2C"))

RelativeAbund_Allscenarios.GBR_All<- ggplot(df_long.GBR_all, aes(Scenario, species, fill= Abund)) + 
  geom_tile()+scale_fill_viridis_c()


#D. Estimate Ecolgoical persistance from normalized abundance, bleaching susceptibility, and range size----
s = 402
dim(merged_df_gbr4)#
merged_df_gbr4$abund_2C <- ifelse(merged_df_gbr4$abund_2C > 0, -1, merged_df_gbr4$abund_2C)
View(merged_df_gbr4)

#calculate ecological persistance Madin et al. 2023

duplicate_names1 <- duplicated(merged_df_gbr4$species)#none
duplicate_names2 <- duplicated(dat$species) #
dat <- dat[!duplicated(dat$species), ] 
dim(dat) #

dat_merged_df_gbr4 <- merge(merged_df_gbr4, dat,  by = "species", all = TRUE) #
dim(dat_merged_df_gbr4) #
species_to_drop3 <- setdiff(dat$species, merged_df_gbr4$species)
# Filter out rows with species present in means_per_species but not in gbr_all
merged_df_gbr_final_final <- dat_merged_df_gbr4[!dat_merged_df_gbr4$species %in% species_to_drop3, ]
dim(merged_df_gbr_final_final) #

#input to equation must be called, needs to be abund called for input to work

merged_df_gbr4_Norm<-merged_df_gbr_final_final
merged_df_gbr4_Norm$abund<-merged_df_gbr_final_final$mergeColumn2

merged_df_gbr4_1.5<-merged_df_gbr_final_final
merged_df_gbr4_1.5$abund<-merged_df_gbr_final_final$abund_1.5

merged_df_gbr4_2<-merged_df_gbr_final_final
merged_df_gbr4_2$abund<-merged_df_gbr_final_final$abund_2C

#make to run sure PCA run at the top----

choice_1 <- choice(merged_df_gbr4_Norm, s, vars=c("abund","range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_1) <- c("Species", "Ecological persistence Present")
choice_2 <- choice(merged_df_gbr4_1.5, s, vars=c("abund","range","bleach"), trait=FALSE)[c("species", "value")]
names(choice_2) <- c("Species", "Ecological persistence 1.5C")
choice_3 <- choice(merged_df_gbr4_2, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_3) <- c("Species", "Ecological persistence 2C")

matx <- merge(choice_1, choice_2, all=TRUE)
matx <- merge(matx, choice_3, all=TRUE)
matx$`Ecological persistence Present`[is.na(matx$`Ecological persistence Present`)] <- 0
matx$`Ecological persistence 1.5C`[is.na(matx$`Ecological persistence 1.5C`)] <- 0
matx$`Ecological persistence 2C`[is.na(matx$`Ecological persistence 2C`)] <- 0
matx <- matx[order(matx["Ecological persistence Present"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 1.5C"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 2C"], decreasing=FALSE),]

image(t(as.matrix(matx[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(matx)-1))), labels=matx$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"), cex.axis=0.8, las=2)

#E. Plot Ecological Persistence with ggplot----
mat_gbr <- merge(choice_1, choice_2, all=TRUE)
mat_gbr <- merge(mat_gbr, choice_3, all=TRUE)
mat_gbr$Ecol.Persis<-mat_gbr$`Ecological persistence`

#Now determine which species are lost with different scenarios and re-run the trait space

df_long_heatmap.gbr <- pivot_longer(
  mat_gbr,
  cols = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"),
  names_to = "Scenario",
  values_to = "Ecol.Persistance"
)

df_long_heatmap.gbr$Scenario <- factor(df_long_heatmap.gbr$Scenario, levels = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"))

(Ecol.Pers_Allscenarios.gbr<- ggplot(df_long_heatmap.gbr, aes(x=Scenario, y = Species, fill = Ecol.Persistance)) +
    geom_tile() + scale_fill_viridis_c())

#MODEL output saved -----
save(mat_gbr, file = "GBR_model_pers.RData")


#end-----

#F. Trait space----

#trait space made on only the corals we have data for for the persistence metric
#run from very top from Madin et al. data
dim(dat) 
duplicate_names <- duplicated(dat$species)
dat.test <- dat[!duplicate_names, ]
dim(dat.test) 

#Present day
gbr_Present#use model output from above to create species list for present that have data needed for this step
dim(gbr_Present) #
gbr_Present$species[!(gbr_Present$species %in% dat.test$species)]
dat.test$gbr_Present <- dat.test$species %in% gbr_Present$species
sum(dat.test$gbr_Present) #

#add 1.5C here
gbr_1.5#use model output from above to create species list for present that have data needed for this step
dim(gbr_1.5) #
gbr_1.5$species[!(gbr_1.5$species %in% dat.test$species)]
dat.test$gbr_1.5 <- dat.test$species %in% gbr_1.5$species
sum(dat.test$gbr_1.5) #


#G. brought over ----
#go back up to top of data_prep.R and run from the very top
#then create dat.test object there

dat<-dat.test #from Ninglaoo data_prep #derep should be 391, lord howe 391, GBR 391
dat <- dat[!duplicated(dat$species), ] #391 Ninglaoo, 391 lord howe, GBR 391
dim(dat)

# TRAIT SPACE
cats <- c("cat_growthrate", "cat_corallitesize", "cat_colonydiameter" , "cat_skeletaldensity", "cat_colonyheight","cat_SA_vol", "cat_spacesize")
pca <- princomp(dat[cats], cor=TRUE, scores=TRUE)
space <- data.frame(PC1=jitter(pca$scores[,1], amount=0.1), PC2=jitter(pca$scores[,2], amount=0.1), PC3=jitter(pca$scores[,3], amount=0.1))
dat <- cbind(dat, space)
dim(dat) #

#H. create convex hull for all data----
TABS<-dat[, 21:22] #PC1 and PC2 for GBR
#TABS<-dat[, 23:24] #new with min family value
CHS<- chull(TABS)

#dat <- dat[, -c(19:20)]
#dat <- dat[, -c(23:25)]#new with min family value

GBR_h <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & gbr_Present == 'TRUE',
         select=c('species','goreau_Madin', 'gbr_Present','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#present
tab7.g <- GBR_h[1:105, 4:5] #builders
tab8.g <- GBR_h[106:164, 4:5] #Cementer
tab9.g <- GBR_h[165:186, 4:5] #Cementers & Fillers
tab901.g <- GBR_h[187:247, 4:5] #fillers

GBR_h.15 <- dat%>% 
  subset(goreau_Madin %in% c('Fillers','Builders', 'Cementers', 'Cementers & Fillers') & gbr_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'gbr_1.5','PC1','PC2','PC3')) %>% 
  arrange(goreau_Madin)

#1.5C
tab10.g <- GBR_h.15[1:77, 4:5] #builders
tab11.g <- GBR_h.15[78:125, 4:5] #Cementer
tab12.g <- GBR_h.15[126:146, 4:5] #Cementers & Fillers
tab1201.g <- GBR_h.15[147:202, 4:5] #fillers

#minimum family
#tab10.g <- GBR_h.15[1:69, 4:5] #builders
#tab11.g <- GBR_h.15[70:109, 4:5] #Cementer
#tab12.g <- GBR_h.15[110:130, 4:5] #Cementers & Fillers
#tab1201.g <- GBR_h.15[131:183, 4:5] #fillers


#---

GBR_ <- dat%>% 
  subset(gbr_Present == 'TRUE',
         select=c('species','goreau_Madin', 'gbr_Present','PC1','PC2','PC3'))
tab000.g <- GBR_[, 4:5]

ch000.g <- chull(tab000.g)
ch7.g <- chull(tab7.g)
ch8.g <- chull(tab8.g)
ch9.g <- chull(tab9.g)
ch901.g <- chull(tab901.g)

GBR_15 <- dat%>% 
  subset(gbr_1.5 == 'TRUE',
         select=c('species','goreau_Madin', 'gbr_1.5','PC1','PC2','PC3'))
tab000_1.5.g <- GBR_15[, 4:5]

ch000_1.5.g <- chull(tab000_1.5.g)
ch7_1.5.g <- chull(tab10.g)
ch8_1.5.g <- chull(tab11.g)
ch9_1.5.g <- chull(tab12.g)
ch901_1.5.g <- chull(tab1201.g)

#I. Plot PCA for GBR climate scenarios----

plot(PC2 ~ PC1, dat, pch=20, cex=0.5, col="grey", ylab="", xlab="", cex.axis=0.75, ylim=c(-3, 3.5), xlim=c(-4, 5), axes=FALSE)
GBR_pres <- dat[dat$gbr_Present,] #340
GBR_1.5 <- dat[dat$gbr_1.5,] #251
mtext("GBR under different climate scenarios", line=-1, adj=0, cex=1.5)
axis(1, cex.axis=0.8)
axis(2, las=2, cex.axis=0.8)
mtext("PC1", 1, 2, cex=0.8)
mtext("PC2", 2, 2, cex=0.8)
text(2, -2.5, "Present: n=340, 1.5C: n=251", cex=1.5)
polygon(tab000.g[ch000.g, ], border = "#73c391", col = rgb(115, 195, 145, 128, maxColorValue = 255), lwd = 2) #present day
polygon(tab000_1.5.g[ch000_1.5.g, ], border="#f36e45", col = rgb(243, 110, 69, 100, maxColorValue = 255), lwd=2) #1.5C
polygon(TABS[CHS, ], border="black", lwd=2)
points(PC2 ~ PC1, GBR_pres, col="#73c391", pch=20, cex=1)
points(PC2 ~ PC1, GBR_1.5, col="#f36e45", pch=20, cex=1.5)


#save MODEL data----
#save(dat, GBR_pres, GBR_1.5, TABS,file="GBR_traitspaceGCB.RData")


#J. Loss in function ----
#GBR_pres
#GBR_1.5

result_df_pres.GBR <- GBR_pres %>%
  group_by(gbr_Present, genus,goreau_Madin) %>% summarise(proportionPres = n())

result_df_1.5.GBR <- GBR_1.5 %>%
  group_by(gbr_1.5, genus, goreau_Madin) %>% summarise(proportion1.5 = n())

#-----
combined_df.GBR <- merge(result_df_pres.GBR, result_df_1.5.GBR, by = c("genus", "goreau_Madin"), all = TRUE)

combined_df.GBR <- combined_df.GBR %>%mutate(proportion1.5 = ifelse(is.na(proportion1.5), 0, proportion1.5))

melted_df.GBR <- melt(combined_df.GBR,
                     id.vars = c("genus", "goreau_Madin"),
                     measure.vars = c("proportionPres", "proportion1.5"),
                     variable.name = c("time"),
                     value.name = "value")

melted_df.GBR$goreau_Madin <- ifelse(is.na(melted_df.GBR$goreau_Madin) | melted_df.GBR$goreau_Madin == "NA", "NAs", melted_df.GBR$goreau_Madin)

custom_order_GBR <- c("Cementers & Fillers","Fillers","Cementers", "Builders", "NAs")

# Change the order of the factor levels in goreau_Madin
melted_df.GBR$goreau_Madin <- factor(melted_df.GBR$goreau_Madin, levels = custom_order_GBR)

# Calculate percentage change using dplyr
percentage_change_df.GBR <- melted_df.GBR %>%
  group_by(goreau_Madin, variable) %>% summarise(proportion = sum(value))

(percentage_change_result.GBR <- percentage_change_df.GBR %>%
    group_by(goreau_Madin) %>% arrange(variable) %>%
    mutate(percentage_change = (proportion / lag(proportion) - 1) * 100))

#K. Change in Volume GBR ----
# Assuming 'tab000' and 'tab000_1.5' are your data frames containing PC1 and PC2 scores,
# and 'ch000' and 'ch000_1.5' contain the indices

# Function to calculate the area of a convex hull in 2D space
calculate_convex_hull_area <- function(data, indices) {
  subset_points <- data[indices, ]
  n <- nrow(subset_points)
  
  # Duplicate the first point to close the polygon
  subset_points <- rbind(subset_points, subset_points[1, ])
  
  # Calculate the area using the Shoelace formula
  area <- 0.5 * sum(subset_points[1:(n-1), 1] * subset_points[2:n, 2] -
                      subset_points[2:n, 1] * subset_points[1:(n-1), 2])
  
  return(abs(area))
}

# Calculate areas for present day and 1.5C scenarios
area_present_day.GBR <- calculate_convex_hull_area(tab000.g, ch000.g)
area_1.5C.GBR <- calculate_convex_hull_area(tab000_1.5.g, ch000_1.5.g)

# Assuming the convex hull is in 2D, volume is then area in 2D
(volume_present_day.GBR <- area_present_day.GBR) 
(volume_1.5C.GBR <- area_1.5C.GBR) 

# Given values
value1 <- volume_present_day.GBR
value2 <- volume_1.5C.GBR

# Calculate percent difference
(percent_difference <- abs((value1 - value2) / ((value1 + value2) / 2)) * 100)

# Print the result
cat("Percent Difference:", percent_difference, "%\n")

# Calculate the difference in volumes
(volume_difference.GBR <- abs(volume_present_day.GBR - volume_1.5C.GBR))

# Print the result
cat("Difference in Convex Hull Volumes:", volume_difference, "\n")

######end-----

#Calculating differences in trait space at present day


# Given values
values <- c(
  GBR = 20.53593,
  LHI = 16.19585,
  Ning = 20.23883,
  SB = 19.90475
)

# Create a matrix of all combinations
combinations <- expand.grid(value1 = names(values), value2 = names(values))

# Calculate percent difference
combinations$PercentDifference <- apply(combinations, 1, function(row) {
  value1 <- values[row["value1"]]
  value2 <- values[row["value2"]]
  abs((value1 - value2) / ((value1 + value2) / 2)) * 100
})

# Plot the percent differences
melted_df <- melt(combinations, id.vars = c("value1", "value2"))
#SI Figure 12----
ggplot(melted_df, aes(x = value2, y = value1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Symmetric Matrix Heatmap",
       x = "Percent (%) overall in present day trait space", y = "Percent (%) overall in present day trait space") +
  theme_minimal()

#   A. Main figure for functional loss----
melted_df$Reef<-"Lord Howe Island"
melted_df.Ning$Reef <-"Ningaloo"
melted_df.SB$Reef <-"Shark Bay"
melted_df.GBR$Reef <-"Great Barrier Reef"

combined_df.Allreefs1 <- merge(melted_df,melted_df.Ning, by = c("genus", "goreau_Madin", "variable", "Reef", "value"), all = TRUE)
combined_df.Allreefs2 <- merge(melted_df.SB,melted_df.GBR, by = c("genus", "goreau_Madin", "variable", "Reef", "value"), all = TRUE)
combined_df.Allreefs3 <- merge(combined_df.Allreefs1,combined_df.Allreefs2, by = c("genus", "goreau_Madin", "variable",  "Reef", "value"), all = TRUE)
head(combined_df.Allreefs3)

#----
#----

#Bubble plot----
my_colors <- c("#73c391", "#f36e45")
exclude_categories <- c("NAs", "Cementers & Fillers")
filtered_data <- combined_df.Allreefs3 %>%
  filter(!goreau_Madin %in% exclude_categories)
custom_order_Goreau <- c( "Builders", "Cementers", "Fillers")
filtered_data$goreau_Madin <- factor(filtered_data$goreau_Madin, levels = custom_order_Goreau)

final_bubble<-ggplot(filtered_data, aes(x = reorder(genus, -value), y = value, size = value, color = variable)) +
  geom_point(alpha = 0.7)+
  facet_grid(goreau_Madin~Reef, scales = "free_y") +
  labs(title = "Bubble Graph Faceted by Reef",
       x = "Genus", y = "Proportion of species",
       size = "Value", color = "Variable") +
  theme_minimal() +coord_flip()+scale_color_manual(values = my_colors)

#-----

#Total species lost?-----
GBR_1.5.TRUE<-GBR_1.5 %>% filter(gbr_1.5 == c("TRUE"))
lhi_1.5.TRUE<-lhi_1.5%>% filter(lhi_all_1.5 == c("TRUE"))
Nings_1.5.TRUE<-Nings_1.5%>% filter(Ning_all_1.5== c("TRUE"))
Shark_1.5.TRUE<-Shark_1.5%>% filter(shark_1.5== c("TRUE"))

GBR_1.5.TRUE <- GBR_1.5.TRUE[, "species", drop = FALSE]
lhi_1.5.TRUE <- lhi_1.5.TRUE[, "species", drop = FALSE]
Nings_1.5.TRUE <- Nings_1.5.TRUE[, "species", drop = FALSE]
Shark_1.5.TRUE <- Shark_1.5.TRUE[, "species", drop = FALSE]

TRUE.bind<-rbind(GBR_1.5.TRUE,lhi_1.5.TRUE,Nings_1.5.TRUE,Shark_1.5.TRUE)

unique_species_Aus <- TRUE.bind %>%
  distinct(species, .keep_all = TRUE) 

GBR_P.TRUE<-GBR_pres %>% filter(gbr_Present == c("TRUE"))
lhi_P.TRUE<-lhi_pres%>% filter(lhi_all == c("TRUE"))
Nings_P.TRUE<-Nings_pres%>% filter(Ningaloo_Present== c("TRUE"))
Shark_P.TRUE<-Shark_pres%>% filter(shark_Present== c("TRUE"))

GBR_P.TRUE <- GBR_P.TRUE[, "species", drop = FALSE]
lhi_P.TRUE <- lhi_P.TRUE[, "species", drop = FALSE]
Nings_P.TRUE <- Nings_P.TRUE[, "species", drop = FALSE]
Shark_P.TRUE <- Shark_P.TRUE[, "species", drop = FALSE]

TRUE.bind.Present<-rbind(GBR_P.TRUE,lhi_P.TRUE,Nings_P.TRUE,Shark_P.TRUE)#

unique_species_Aus_present <- TRUE.bind.Present %>%
  distinct(species, .keep_all = TRUE) 


#Below is only for Central GBR abundances----

C.GBR_Baird_Abun ###Baird et al. central GBR Abundances, Baird et al. 2013 (https://zenodo.org/records/4310394)

library(reshape)

#data clean up
C.GBR_Baird_Abun <- C.GBR_Baird_Abun %>%
  mutate(genus = gsub("Lobophytum", "Lobophyllia", genus))

head(total_sums <- C.GBR_Baird_Abun %>%
  group_by(Site, transect, Year_month) %>%
  summarize(Total_Sum = sum(Total)))

head(GBR.means_per_species_persite <- C.GBR_Baird_Abun  %>%
  # Join with the total_sums data frame based on transect and Year_month
  left_join(total_sums, by = c("Site", "transect", "Year_month")) %>%
  # Calculate the percentage for each genus within each transect and Year_month combination
  mutate(Percentage = (Total / Total_Sum) * 100))
  
head(GBR.means_per_species_persite2 <- GBR.means_per_species_persite %>%
  group_by(Site, genus) %>% 
  summarize(
    Mean_ = mean(Percentage,na.rm = TRUE)))

head(GBR.means_per_species_persite3 <- GBR.means_per_species_persite %>%
       group_by(genus) %>% 
       summarize(
         Mean_ = mean(Percentage,na.rm = TRUE)))

##Mean abundances now----

#follow same instructions as above to run as outlined in full GBR section

# Create a new column with NA values
new_column <- rep(NA, nrow(gbr_all)) #

# Insert the new column next to the column
gbr_all <- cbind(gbr_all[,1], new_column)

# Rename the new column to abund_DBCA
colnames(gbr_all)[2] <- "mergeColumn2"
colnames(gbr_all)[1] <- "species"

gbr_all<-as.data.frame(gbr_all)

# Split the 'species_genus' column by space
split_names <- strsplit((gbr_all$species), " ")

# Extract genus and species from the split_names list
genus <- sapply(split_names, `[`, 1)
species <- sapply(split_names, `[`, 2)

# Create a new data frame with species and genus columns
new_data <- data.frame(species = species, genus = genus)

gbr_all_final <- cbind(gbr_all[,1], new_data, new_column)
colnames(gbr_all_final)[4] <- "abund_C.GBR.Baird"
head(gbr_all_final)


gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Acropora", 17.414866, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Alveopora", 2.158104, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Acanthastrea", 3.566270, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Astreopora", 3.447649, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Coeloseris", 1.709619, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus =="Coscinaraea", 1.674697, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Cyphastrea", 3.790106, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Diploastrea", 1.968635, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Echinophyllia", 3.2740572, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Echinopora", 2.897029, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Favites", 7.429746, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Fungia", 2.741290, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Galaxea", 3.878136, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Goniastrea", 4.196856, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Goniopora", 3.003522e-0, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Hydnophora", 2.038686, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Isopora", 2.019209, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Leptastrea", 2.727695, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Leptoria", 2.042803, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Leptoseris", 2.200977, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Lobophyllia", 3.474478, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Merulina", 2.900909, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Montipora", 23.024598, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Moseleya", 2.562834, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Mycedium", 2.286701, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Oulophyllia", 1.330439, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Oxypora", 2.465484, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Pachyseris", 4.649683, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Pavona", 2.499159, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Pectinia", 3.212615, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Platygyra", 3.988627, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Plesiastrea", 2.229877, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Pocillopora", 2.959614, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Porites", 8.735734, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Psammocora", 3.006104, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Seriatopora", 2.400261, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Stylocoeniella", 3.448276, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Stylophora", 1.937902, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Symphyllia", 2.962917, gbr_all_final$abund_C.GBR.Baird)
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Turbinaria", 12.463169, gbr_all_final$abund_C.GBR.Baird)

unique_genus_with_na <- unique(subset(gbr_all_final, is.na(abund_C.GBR.Baird))$genus)

# Print the unique species names with NA values in 'abund_C.GBR.Baird'
print(unique_genus_with_na)

####
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Anacropora", 17.414866, gbr_all_final$abund_C.GBR.Baird) #used Acroporidae
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Astrea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Australogyra", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Bernardpora", 8.735734, gbr_all_final$abund_C.GBR.Baird) #Porites 8.735734
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Cantharellus", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Caulastraea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Coelastrea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Ctenactis", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Cycloseris", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Cynarina", 3.47447, gbr_all_final$abund_C.GBR.Baird) #lobo
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Danafungia", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Dipsastraea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Echinomorpha", 3.212615, gbr_all_final$abund_C.GBR.Baird) #Pectinea
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Gardineroseris", 2.741290, gbr_all_final$abund_C.GBR.Baird) #agarace, suborder fungiina
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Halomitra", 2.741290, gbr_all_final$abund_C.GBR.Baird) ##fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Heliofungia", 2.741290, gbr_all_final$abund_C.GBR.Baird) ##fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Herpolitha", 2.741290, gbr_all_final$abund_C.GBR.Baird) ##fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Homophyllia", 3.47447, gbr_all_final$abund_C.GBR.Baird) #lobo
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Lithophyllon", 2.741290, gbr_all_final$abund_C.GBR.Baird)#fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Lobactis", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Madracis", 2.959614, gbr_all_final$abund_C.GBR.Baird) #Pocoillo
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Micromussa", 3.474478, gbr_all_final$abund_C.GBR.Baird) #Lobophyllidae
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Oulastrea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #Merulinidea
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Palauastrea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae 
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Paragoniastrea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Paramontastraea", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Pleuractis", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Podabacia", 2.741290, gbr_all_final$abund_C.GBR.Baird) #fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Polyphyllia", 2.741290, gbr_all_final$abund_C.GBR.Baird)#fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Sandalolitha", 2.741290, gbr_all_final$abund_C.GBR.Baird)#fungia
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Scapophyllia", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Stylaraea", 8.735734, gbr_all_final$abund_C.GBR.Baird) #Porites 8.735734
gbr_all_final$abund_C.GBR.Baird <- ifelse(gbr_all_final$genus == "Trachyphyllia", 2.900909, gbr_all_final$abund_C.GBR.Baird) #in family Merulinidae , used related Merulina 2.900909

print(unique_genus_with_na2 <- unique(subset(gbr_all_final, is.na(abund_C.GBR.Baird))$genus)) #402 but some with NA
#1] "Blastomussa"       "Catalaphyllia"     "Euphyllia"         "Physogyra"         "Plerogyra"         "Pseudosiderastrea" "Siderastrea"      

dat3<-dat #created from PCA after creating dat.test in data_prep.R
sum(dat3$gbr_all) #

dat_gbr<- dat3[dat3$gbr_all,] #

colnames(gbr_all_final)[1] <- "Name"
dat_gbr$Name<-dat_gbr$species
#
merged_df_gbr <- merge(gbr_all_final, dat_gbr, by = "Name", all = TRUE) #

duplicate_names.gbr <- duplicated(merged_df_gbr$Name)
# Remove rows with duplicate names
merged_df_gbr.single <- merged_df_gbr[!duplicate_names.gbr, ]
dim(merged_df_gbr.single) #

##check the relative abundances again 

merged_df_gbr %>%
  filter(gbr_all == "TRUE") %>%
  ggplot(aes(y = abund_C.GBR.Baird, x = reorder(Name, -abund_C.GBR.Baird), fill = abund_C.GBR.Baird)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

merged_df_gbr %>%
  filter(gbr_all == "TRUE") %>%
  ggplot(aes(y = abund_C.GBR.Baird, x = reorder(Name, -abund_C.GBR.Baird), fill = abund_C.GBR.Baird)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap(.~genus.x) +
  scale_x_discrete(breaks = function(x) unique(merged_df_gbr$Name[merged_df_gbr$genus.x == x]))

merged_df_gbr.single$abund_1.5C <- merged_df_gbr.single$abund_C.GBR.Baird + (-6.74/100)
merged_df_gbr.single$abund_2C <- merged_df_gbr.single$abund_C.GBR.Baird + (-20.18/100)

merged_df_gbr.single2<-merged_df_gbr.single #just to save the dataset

#change all negatives to zeros
merged_df_gbr.single$abund_C.GBR.Baird <- ifelse(merged_df_gbr.single$abund_C.GBR.Baird < 0, 0, merged_df_gbr.single$abund_C.GBR.Baird)
merged_df_gbr.single$abund_1.5C <- ifelse(merged_df_gbr.single$abund_1.5C < 0, 0, merged_df_gbr.single$abund_1.5C)
merged_df_gbr.single$abund_2C <- ifelse(merged_df_gbr.single$abund_2C < 0, 0, merged_df_gbr.single$abund_2C)

#normalize 
merged_df_gbr.single$abund_C.GBR.Baird[merged_df_gbr.single$abund_C.GBR.Baird > 0.5 & merged_df_gbr.single$abund_C.GBR.Baird <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_1.5C[merged_df_gbr.single$abund_1.5C > 0.5 & merged_df_gbr.single$abund_1.5C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_2C[merged_df_gbr.single$abund_2C > 0.5 & merged_df_gbr.single$abund_2C <= 1] <- 1.0 #normalize 0.25 to 0.5 is rare

merged_df_gbr.single$abund_C.GBR.Baird[merged_df_gbr.single$abund_C.GBR.Baird > 0.25 & merged_df_gbr.single$abund_C.GBR.Baird <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_1.5C[merged_df_gbr.single$abund_1.5C > 0.25 & merged_df_gbr.single$abund_1.5C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_2C[merged_df_gbr.single$abund_2C > 0.25 & merged_df_gbr.single$abund_2C <= 0.5] <- 0.5 #normalize 0.25 to 0.5 is rare

merged_df_gbr.single$abund_C.GBR.Baird[merged_df_gbr.single$abund_C.GBR.Baird > 0 & merged_df_gbr.single$abund_C.GBR.Baird <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_1.5C[merged_df_gbr.single$abund_1.5C > 0.0 & merged_df_gbr.single$abund_1.5C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare
merged_df_gbr.single$abund_2C[merged_df_gbr.single$abund_2C > 0.0 & merged_df_gbr.single$abund_2C <= 0.25] <- 0.25 #normalize 0.25 to 0.5 is rare

####

df_long.GBR <- pivot_longer(
  merged_df_gbr.single,
  cols = c("abund_C.GBR.Baird", "abund_1C", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

#filter out any that were "zero" in the Present day
#n =
df_long.GBR <- df_long.GBR %>%
  filter(!(Scenario == "abund_C.GBR.Baird" & Abund == 0))

df_long.GBR$Scenario <- factor(df_long.GBR$Scenario, levels = c("abund_C.GBR.Baird", "abund_1C", "abund_1.5C", "abund_2C"))

RelativeAbund_Allscenarios.GBR<- ggplot(df_long.GBR, aes(Scenario, Name, fill= Abund)) + 
  geom_tile()+scale_fill_viridis_c()

s = 402
dim(merged_df_gbr.single)
merged_df_gbr.single$abund_2C <- ifelse(merged_df_gbr.single$abund_2C > 0, -1, merged_df_gbr.single$abund_2C)
View(merged_df_gbr.single)

colnames(merged_df_gbr.single)[5] <- "species"

merged_df_gbr.single_Norm<-merged_df_gbr.single
merged_df_gbr.single_Norm$abund<-merged_df_gbr.single$abund_C.GBR.Baird

merged_df_gbr.single_1.5<-merged_df_gbr.single
merged_df_gbr.single_1.5$abund<-merged_df_gbr.single$abund_1.5

merged_df_gbr.single_2<-merged_df_gbr.single
merged_df_gbr.single_2$abund<-merged_df_gbr.single$abund_2C

#make sure PCA run
choice_1 <- choice(merged_df_gbr.single_Norm, s, vars=c("abund","range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_1) <- c("Species", "Ecological persistence Present")
choice_2 <- choice(merged_df_gbr.single_1.5, s, vars=c("abund","range","bleach"), trait=FALSE)[c("species", "value")]
names(choice_2) <- c("Species", "Ecological persistence 1.5C")
choice_3 <- choice(merged_df_gbr.single_2, s, vars=c("abund", "range", "bleach"), trait=FALSE)[c("species", "value")]
names(choice_3) <- c("Species", "Ecological persistence 2C")

matx <- merge(choice_1, choice_2, all=TRUE)
matx <- merge(matx, choice_3, all=TRUE)
matx$`Ecological persistence Present`[is.na(matx$`Ecological persistence Present`)] <- 0
matx$`Ecological persistence 1.5C`[is.na(matx$`Ecological persistence 1.5C`)] <- 0
matx$`Ecological persistence 2C`[is.na(matx$`Ecological persistence 2C`)] <- 0
matx <- matx[order(matx["Ecological persistence Present"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 1.5C"], decreasing=FALSE),]
matx <- matx[order(matx["Ecological persistence 2C"], decreasing=FALSE),]

image(t(as.matrix(matx[,2:4])), axes=FALSE, col=hcl.colors(20, "Heat", rev = TRUE))
mtext("A", line=4, cex=1.5, at=-1)

text(rep(0, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,1]), 2)), cex=0.6)
text(rep(0.5, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,2]), 2)), cex=0.6)
text(rep(1, nrow(matx[,2:4])), seq(0, 1, (1/(nrow(matx[,2:4])-1))), (round(as.matrix(matx[,2:4][,3]), 2)), cex=0.6)

axis(2, at=seq(0, 1, (1/(nrow(matx)-1))), labels=matx$Species, las=2, cex.axis=0.7, font=3)
axis(3, at=c(0, 0.5, 1), labels=c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"), cex.axis=0.8, las=2)

#Plot Ecological Persistence with ggplot
mat_gbr <- merge(choice_1, choice_2, all=TRUE)
mat_gbr <- merge(mat_gbr, choice_3, all=TRUE)
mat_gbr$Ecol.Persis<-mat_gbr$`Ecological persistence`

#Now determine which species are lost with different scenarios and re-run the trait space

df_long_heatmap.gbr <- pivot_longer(
  mat_gbr,
  cols = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"),
  names_to = "Scenario",
  values_to = "Ecol.Persistance"
)

df_long_heatmap.gbr$Scenario <- factor(df_long_heatmap.gbr$Scenario, levels = c("Ecological persistence Present", "Ecological persistence 1.5C", "Ecological persistence 2C"))

(Ecol.Pers_Allscenarios.gbr<- ggplot(df_long_heatmap.gbr, aes(x=Scenario, y = Species, fill = Ecol.Persistance)) +
    geom_tile() + scale_fill_viridis_c())

#end of Central GBR-----


#Modelled relative abundances line chart----

#GBR
Abundance_lineplot_long <- pivot_longer(
  Abundance_lineplot,
  cols = c("Present", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

Abundance_lineplot_Ning$Present<-Abundance_lineplot_Ning$abund_DBCA
Abundance_lineplot_Ning_long <- pivot_longer(
  Abundance_lineplot_Ning,
  cols = c("Present", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

Abundance_lineplot.Shark$Present<-Abundance_lineplot.Shark$abund_DBCA
Abundance_lineplot.Shark.long <- pivot_longer(
  Abundance_lineplot.Shark,
  cols = c("Present", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

Abundance_lineplot_LH$Present<-Abundance_lineplot_LH$CombinedRelAbud
Abundance_lineplot_LH.long <- pivot_longer(
  Abundance_lineplot_LH,
  cols = c("Present", "abund_1.5C", "abund_2C"),
  names_to = "Scenario",
  values_to = "Abund"
)

Abundance_lineplot_LH.long$Region<-"LordHowe"
Abundance_lineplot.Shark.long$Region<-"Shark"
Abundance_lineplot_long$Region<-"GBR"
Abundance_lineplot_Ning_long$Region<-"Ning"
Abundance_lineplot_long$Abund<-Abundance_lineplot_long$Abund*100
Abundance_lineplot_LH.long$Abund<-Abundance_lineplot_LH.long$Abund*100

merged_df <- full_join(Abundance_lineplot.Shark.long, Abundance_lineplot_long, by = c("species", "Region", "Scenario", "Abund"))
merged_df <- full_join(merged_df,Abundance_lineplot_Ning_long,  by = c("species", "Region", "Scenario", "Abund"))
merged_df <- full_join(merged_df,Abundance_lineplot_LH.long,  by = c("species", "Region", "Scenario", "Abund"))

merged_df$Scenario <- factor(merged_df$Scenario, levels = c("Present", "abund_1.5C", "abund_2C"))

merged_df.small<-merged_df[1:3039,]
merged_df.small<-merged_df.small[,c(1,24:26)]
merged_df.small <- merged_df.small %>%
  mutate(Scenario = recode(Scenario, "Present" = "1", "abund_1.5C" = "2", "abund_2C" = "3"))
merged_df.small$Scenario<-as.numeric(merged_df.small$Scenario)

#All jittered----
add_jitter <- function(vector, magnitude = 0.4) {
  jittered_vector <- vector + runif(length(vector), min = -magnitude, max = magnitude)
  return(jittered_vector)
}
# Jitter the Abund values for each species
merged_df.small$jittered_Abund <- with(merged_df.small, ave(Abund, species, FUN = add_jitter))

#*Final jittered plot for modelled abundance models----
modelled.Abundance<-ggplot(merged_df.small, aes(x = Scenario, y = jittered_Abund, color= Region, alpha = 0.1)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, aes(group = species), color= "darkgrey", linetype = "solid", size=0.5, alpha = 0.1) +
  labs(x = "Scenarios (°C)", y = "Modelled relative abundance (%)") +
  scale_fill_manual(values = c("GBR" = "#0e3a5e", "Ning" = "#53b2a8", "Shark" ="#dd9e5f", "LordHowe" = "#d1c8c1"))+
  scale_color_manual(values = c("GBR" = "#0e3a5e", "Ning" = "#53b2a8", "Shark" ="#dd9e5f", "LordHowe" = "#d1c8c1")) +
  geom_smooth(method = "loess", aes(group = Region, color=Region, fill=Region), se = TRUE, alpha = 0.4)+
  theme_bw()+geom_hline(yintercept = 0 , color = "black", linetype = "dashed", size = 0.5)+
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Present", "+1.5 (°C)", "+2 (°C)"))+geom_point(size=0.1, alpha = 0.7)

#Sensitivity analysis and non-linear functions-----

#filtering sensitivity analysis----
RateChange_Gonzalez <- RateChange_Gonzalez %>% filter(disturbance_type == c("bleaching"))
aus_bleaching.filt<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Exmouth_to_Broome", "Ningaloo",
                          "Torres_Strait_Northern_Great_Barrier_Reef", 
                          "Central_and_Southern_Great_Barrier_Reef"))
Zero_to_ten_bleaching.filt<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Western_and_Northern_Madagascar",
                          "Western and Northern Madagascar", 
                          "Mascarene_Islands", "Mascarene Islands", 
                          "Exmouth_to_Broome", "Ningaloo",
                          "Torres_Strait_Northern_Great_Barrier_Reef", 
                          "Central_and_Southern_Great_Barrier_Reef"))
neg10_to_neg20_bleaching.filt<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Panama_Bight", "Nicoya", "Cortezian", 
                          "East_African_Coral_Coast", "Maldives",
                          "Palawan/North_Borneo", "Phoenix/Tokelau/Northern_Cook_Islands", 
                          "Western_Caribbean", "Tropical_Eastern_Pacific", 
                          "Western_and_Northern_Madagascar", "Western and Northern Madagascar", 
                          "Mascarene_Islands", "Mascarene Islands", 
                          "Exmouth_to_Broome", "Ningaloo",
                          "Torres_Strait_Northern_Great_Barrier_Reef", 
                          "Central_and_Southern_Great_Barrier_Reef"))
neg30_to_neg40_bleaching.filt<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Eastern_Philippines", "South_China_Sea_Oceanic_Islands",
                          "Panama_Bight", "Nicoya", "Cortezian", 
                          "East_African_Coral_Coast", "Maldives",
                          "Palawan/North_Borneo", "Phoenix/Tokelau/Northern_Cook_Islands", 
                          "Western_Caribbean", "Tropical_Eastern_Pacific", 
                          "Western_and_Northern_Madagascar", "Western and Northern Madagascar", "Mascarene_Islands", "Mascarene Islands",
                          "Exmouth_to_Broome", "Ningaloo",
                          "Torres_Strait_Northern_Great_Barrier_Reef", 
                          "Central_and_Southern_Great_Barrier_Reef"))

aus_bleaching.filt <- aus_bleaching.filt %>%
  mutate(Sensitivity = "aus_bleaching.filt")

Zero_to_ten_bleaching.filt <- Zero_to_ten_bleaching.filt %>%
  mutate(Sensitivity = "Zero_to_ten_bleaching.filt")

neg10_to_neg20_bleaching.filt <- neg10_to_neg20_bleaching.filt %>%
  mutate(Sensitivity = "neg10_to_neg20_bleaching.filt")

neg30_to_neg40_bleaching.filt <- neg30_to_neg40_bleaching.filt %>%
  mutate(Sensitivity = "neg30_to_neg40_bleaching.filt")

All_sensitivity<-rbind(aus_bleaching.filt,Zero_to_ten_bleaching.filt,neg10_to_neg20_bleaching.filt,neg30_to_neg40_bleaching.filt)

#end of filtering sensitivity----

#figure sensitivity ----
Zero_to_ten_bleaching<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Western_and_Northern_Madagascar",
                        "Western and Northern Madagascar", 
                        "Mascarene_Islands", "Mascarene Islands", 
                        "Exmouth_to_Broome", "Ningaloo",
                        "Torres_Strait_Northern_Great_Barrier_Reef", 
                        "Central_and_Southern_Great_Barrier_Reef")) %>% ggplot(aes(x = DHW, y = absolute_rate_t)) +
  geom_point(aes( color = Ecoregion)) + geom_smooth(method = "lm", aes(DHW, absolute_rate_t), se=T,lty=3, color="black")+
  labs(title="Zero_to_ten_bleaching", x="Degree Heating Weeks", y= "Mean annual rate of change")+theme_bw()+
  stat_regline_equation(aes(DHW, absolute_rate_t))+stat_cor(label.y = max(RateChange_Gonzalez$absolute_rate_t))+scale_x_continuous(limits = c(0, 20))

neg10_to_neg20_bleaching<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Panama_Bight", "Nicoya", "Cortezian", 
                        "East_African_Coral_Coast", "Maldives",
                        "Palawan/North_Borneo", "Phoenix/Tokelau/Northern_Cook_Islands", 
                        "Western_Caribbean", "Tropical_Eastern_Pacific", 
                        "Western_and_Northern_Madagascar", "Western and Northern Madagascar", 
                        "Mascarene_Islands", "Mascarene Islands", 
                        "Exmouth_to_Broome", "Ningaloo",
                        "Torres_Strait_Northern_Great_Barrier_Reef", 
                        "Central_and_Southern_Great_Barrier_Reef")) %>% 
  ggplot(aes(x = DHW, y = absolute_rate_t)) +
  geom_point(aes( color = Ecoregion)) + geom_smooth(method = "lm", aes(DHW, absolute_rate_t), se=T,lty=3, color="black")+
  labs(title="neg10_to_neg20_bleaching", x="Degree Heating Weeks", y= "Mean annual rate of change")+theme_bw()+
  stat_regline_equation(aes(DHW, absolute_rate_t))+stat_cor(label.y = max(RateChange_Gonzalez$absolute_rate_t))+scale_x_continuous(limits = c(0, 20))

neg30_to_neg40_bleaching<-RateChange_Gonzalez %>% 
  filter(Ecoregion %in% c("Eastern_Philippines", "South_China_Sea_Oceanic_Islands",
                                                  "Panama_Bight", "Nicoya", "Cortezian", 
                                                  "East_African_Coral_Coast", "Maldives",
                                                  "Palawan/North_Borneo", "Phoenix/Tokelau/Northern_Cook_Islands", 
                                                  "Western_Caribbean", "Tropical_Eastern_Pacific", 
                                                  "Western_and_Northern_Madagascar", "Western and Northern Madagascar", "Mascarene_Islands", "Mascarene Islands",
                                                  "Exmouth_to_Broome", "Ningaloo",
                                                  "Torres_Strait_Northern_Great_Barrier_Reef", 
                                                  "Central_and_Southern_Great_Barrier_Reef")) %>% 
  ggplot(aes(x = DHW, y = absolute_rate_t)) +
  geom_point(aes(color = Ecoregion)) + geom_smooth(method = "lm", aes(DHW, absolute_rate_t), se=T,lty=3, color="black")+
  labs(title="neg30_to_neg40_bleaching", x="Degree Heating Weeks", y= "Mean annual rate of change")+theme_bw()+
  stat_regline_equation(aes(DHW, absolute_rate_t))+stat_cor(label.y = max(RateChange_Gonzalez$absolute_rate_t))+scale_x_continuous(limits = c(0, 20))

allreefs_bleaching<-RateChange_Gonzalez %>% 
  ggplot(aes(x = DHW, y = absolute_rate_t)) +
  geom_point(aes(color = Ecoregion)) + geom_smooth(method = "lm", aes(DHW, absolute_rate_t), se=T,lty=3, color="black")+
  labs(title="allreefs_bleaching", x="Degree Heating Weeks", y= "Mean annual rate of change")+theme_bw()+
  stat_regline_equation(aes(DHW, absolute_rate_t))+stat_cor(label.y = max(RateChange_Gonzalez$absolute_rate_t))+scale_x_continuous(limits = c(0, 20))

#library(gridExtra)
spatial_test_sensitivity<-grid.arrange(Zero_to_ten_bleaching, neg10_to_neg20_bleaching, neg30_to_neg40_bleaching, ncol = 1)

grid.arrange(allreefs_bleaching,allreefs_bleaching_loess)

new_order <- c("aus_bleaching.filt", "Zero_to_ten_bleaching.filt", "neg10_to_neg20_bleaching.filt", "neg30_to_neg40_bleaching.filt")

# Change the factor order
cbPalette <- c("#009E73", "#56B4E9", "#CC79A7", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#99999")

All_sensitivity <- All_sensitivity %>%
  mutate(Sensitivity = factor(Sensitivity, levels = new_order))

allsens_bleaching_loess<-ggplot(All_sensitivity) + 
  geom_point(aes(DHW, absolute_rate_t), alpha = 0.4, size = 3) +
  geom_smooth(method = "loess", aes(DHW, absolute_rate_t, color = Sensitivity), se = TRUE) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(y = "RateofChange%") +geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values=cbPalette)+
  theme_classic()

allsens_bleaching_lm<-ggplot(All_sensitivity) + 
  geom_point(aes(DHW, absolute_rate_t), alpha = 0.4, size = 3) +
  geom_smooth(method = "lm", aes(DHW, absolute_rate_t, color = Sensitivity), se = TRUE) +
  scale_x_continuous(limits = c(0, 20)) +
  labs(y = "RateofChange%") +geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()+scale_colour_manual(values=cbPalette)

grid.arrange(allsens_bleaching_loess,allsens_bleaching_lm, nrow =1)

#CMIP3 vs. CMIP6----

RateChange_sensitivity <- read.csv("rateChange_sensitivity_species.csv", as.is=TRUE) #provided in Github by author
RateChange_sensitivity <- RateChange_sensitivity %>% filter(IPCC %in% c("CMIP6"))
#library(ggpubr)

LHI_diff<-RateChange_sensitivity %>% filter(Reef == "LHI", Model != "Aus-only") %>%
ggplot(aes(x = Scenario, y = RateChange)) +  
  geom_point(size = 5, alpha = 0.5, aes(color=Model)) + 
  geom_hline(yintercept = -5.78, linetype = "dashed")+
  geom_hline(yintercept = -15.86, col= "red", linetype = "dashed")+
  labs(x = "LHI") + scale_y_continuous(limits = c(-100, 5)) +
  theme_classic()

Ning_diff<-RateChange_sensitivity %>% filter(Reef == "Ning", Model != "Aus-only") %>%
  ggplot(aes(x = Scenario, y = RateChange)) +  
  geom_point(size = 5, alpha = 0.5, aes(color=Model)) + 
  geom_hline(yintercept = -7.22, linetype = "dashed")+
  geom_hline(yintercept = -18.02, col= "red", linetype = "dashed")+
  labs(x = "Ning") + scale_y_continuous(limits = c(-100, 5)) +
  theme_classic()

Shark_diff<-RateChange_sensitivity %>% filter(Reef == "Shark", Model != "Aus-only") %>%
  ggplot(aes(x = Scenario, y = RateChange)) +  
  geom_point(size = 5, alpha = 0.5, aes(color=Model)) + 
  geom_hline(yintercept = -7.7, linetype = "dashed")+
  geom_hline(yintercept = -18.02, col= "red", linetype = "dashed")+
  labs(x = "Shark") + scale_y_continuous(limits = c(-100, 5)) + 
  theme_classic()

GBR_diff<-RateChange_sensitivity %>% filter(Reef == "GBR", Model != "Aus-only") %>%
  ggplot(aes(x = Scenario, y = RateChange)) +  
  geom_point(size = 5, alpha = 0.5, aes(color=Model)) + 
  geom_hline(yintercept = -6.74, linetype = "dashed")+
  geom_hline(yintercept = -20.18, col= "red", linetype = "dashed")+
  labs(x = "GBR") + scale_y_continuous(limits = c(-100, 5)) + 
  theme_classic()

ggarrange(LHI_diff, Ning_diff, Shark_diff,GBR_diff, nrow = 1, common.legend = TRUE, legend="bottom")
grid.arrange(LHI_diff, Ning_diff, Shark_diff,GBR_diff, nrow =1, common.legend = TRUE)

model_order <- c("Aus-only", unique(RateChange_sensitivity$Model[RateChange_sensitivity$Model != "Aus-only"]))
# Plot the data

RateChange_sensitivity <- na.omit(RateChange_sensitivity)

new_order2 <- c("neg30_to_neg40_bleaching.filt", "neg10_to_neg20_bleaching.filt", "Zero_to_ten_bleaching.filt", "aus_bleaching.filt", "AllReefs", "All Central Indo-Pacific")
RateChange_sensitivity$Difference_plot <- ifelse(RateChange_sensitivity$Difference >= 0, -1 * RateChange_sensitivity$Difference, abs(RateChange_sensitivity$Difference))

RateChange_sensitivity$Value_Color <- ifelse(RateChange_sensitivity$Difference_plot >= 0, "Positive", "Negative")

Absolute_difference_sensitivity<-ggplot(RateChange_sensitivity, aes(x = factor(Model, levels = new_order2), y = Difference_plot, color = Value_Color)) +
  geom_segment(aes(x = factor(Model, levels = new_order2), xend = factor(Model, levels = new_order2), y = 0, yend = Difference_plot), position = position_dodge(width = 0.5), color = "grey50")+
  geom_point(position = position_dodge(width = 0.5), size = 4) + # Dots
  facet_grid(Scenario ~ Reef) +
  labs(x = "Model", y = "Absolute difference in estimate rate of change") +
  theme_minimal() +
  geom_text(aes(label = round(Difference_plot)), hjust = -0.8, color="black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(2, "lines")) +
  geom_hline(yintercept = 0, color = "black")+
  scale_color_manual(values = c("Positive" = "blue", "Negative" = "red")) + # Manually set colors
  geom_vline(xintercept = which(new_order2 == "AllReefs") - 0.5, color = "black", linetype = "dashed")

RateChange_sensitivity <- read.csv("rateChange_sensitivity_species.csv", as.is=TRUE)
RateChange_sensitivity_IPCC <- RateChange_sensitivity %>% filter(Model %in% c("Aus-only"))

# Create a plot with points and segments
CMIP3_CMIP6_rateofchange<-ggplot(RateChange_sensitivity_IPCC, aes(x = IPCC, y = RateChange, color = Reef)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +facet_grid(. ~ Scenario)+
  labs(x = "IPCC", y = "Rate of change in coral cover (%)", color = "Model") +
  theme_minimal()

#Map ----
#--- Loading data ---#
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(sdmpredictors)
library(raster)
library(sp)
library(dismo)
library(ggspatial)

sf_use_s2(FALSE) 

world <- ne_countries(scale = "medium", returnclass = "sf")

#Bio-ORACLE - RCCP Pathways----
datasets = list_datasets(terrestrial = FALSE, marine = TRUE)
rcp = c("RCP45")
# Variables of interest
variables = c("temp")
# Extract future data sets
future = list_layers_future(datasets) %>%
  
dplyr::filter(grepl(paste(rcp, collapse = "|"), scenario)) %>% 
  
dplyr::filter(year == 2024 | year == 2050 | year == 2100) %>% 
# keep variables of interest using a regular expression
dplyr::filter(grepl(paste(variables, collapse = "|"), layer_code))

temp.present = c("BO2_tempmean_bdmean","BO2_temprange_bdmean","BO2_tempmean_ss","BO2_temprange_ss")
temp.present = gsub("BO2", "BO21", temp.present)
temp.future = c("BO21_RCP45_2050_tempmean_bdmean",
                "BO21_RCP45_2100_tempmean_bdmean",
                "BO21_RCP45_2050_temprange_bdmean",
                "BO21_RCP45_2100_temprange_bdmean",
                "BO21_RCP45_2050_tempmean_ss",
                "BO21_RCP45_2100_tempmean_ss",
                "BO21_RCP45_2050_temprange_ss",
                "BO21_RCP45_2100_temprange_ss")

temp = c(temp.present, temp.future)

temp.rasters = load_layers(temp)
temp.rasters$BO21_RCP45_2100_temprange_ss

# Define a boundary
boundary.l = extent(c(xmin = 145, xmax = 160, ymin = -40, ymax = -10))
temp.rasters2 = crop(temp.rasters, boundary.l)

cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Plot mean bottom mean temperature
rast.l.p <- raster(temp.rasters2, layer=3)
rast.l.50 <- raster(temp.rasters2, layer=9)
rast.l.100 <- raster(temp.rasters2, layer=10)

crs(rast.l.p) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.p.l <- rasterToPoints(rast.l.p)

crs(rast.l.50) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.50.l <- rasterToPoints(rast.l.50)

crs(rast.l.100) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.100.l <- rasterToPoints(rast.l.100)

df.present <-data.frame(map.p.l)

df.50 <- data.frame(map.50.l) 

df.100 <- data.frame(map.100.l)

#make appropriate column headings
colnames(df.present) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.50) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.100) <- c('Longitude', 'Latitude', 'Mean_temp')

#---# Lord Howe----
df.present$Variable<-"Present"
df.50$Variable<-"2050"
df.100$Variable<-"2100"
df.plot.LordHowe<-rbind(df.present,df.50, df.100)

df.plot.LordHowe$Variable <- factor(df.plot.LordHowe$Variable, levels = c("Present", "2050", "2100"))

LordHowe.Temp.plot<-ggplot(data = world) +
  geom_tile(data=df.plot.LordHowe, aes(x=Longitude, y=Latitude, fill=Mean_temp))+
  geom_sf(color = '#414754', fill = '#414754') + xlab("Longitude") + ylab("Latitude") +
  scale_fill_gradientn('Mean_temp', guide="colorbar", #trans = 'sqrt',
                       colors = c("#4D7DB5","#76C4C9","#CBE7C8","#F6E1A2","#F6AC6D","#DF4D51")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
  #                       style = north_arrow_nautical) +
  coord_sf(xlim = c(145, 160), ylim = c(-40, -20), expand = FALSE) + #edge south is mid way between Lack Mcleod
  theme_bw() + facet_wrap(.~Variable)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

LordHowe.Temp.plot


#Ningaloo Plot ----
boundary.N = extent(c(xmin = 110, xmax = 120, ymin = -25, ymax = -20))

temp.rasters3 = crop(temp.rasters, boundary.N)

cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Plot mean bottom mean temperature

rast.p.n <- raster(temp.rasters3, layer=3) 
rast.50.n <- raster(temp.rasters3, layer=9)
rast.100.n <- raster(temp.rasters3, layer=10) 

crs(rast.p.n) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.p.n <- rasterToPoints(rast.p.n)

crs(rast.50.n) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.50.n <- rasterToPoints(rast.50.n)

crs(rast.100.n) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.100.n <- rasterToPoints(rast.100.n)

df.present.n <- data.frame(map.p.n)

df.50.n <- data.frame(map.50.n) 

df.100.n <- data.frame(map.100.n)

#make appropriate column headings
colnames(df.present.n) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.50.n) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.100.n) <- c('Longitude', 'Latitude', 'Mean_temp')

df.present.n$Variable<-"Present"
df.50.n$Variable<-"2050"
df.100.n$Variable<-"2100"
df.plot.Ning<-rbind(df.present.n,df.50.n, df.100.n)

df.plot.Ning$Variable <- factor(df.plot.Ning$Variable, levels = c("Present", "2050", "2100"))

Ninglaoo.Temp.plot<-ggplot(data = world) +
  geom_tile(data=df.plot.Ning, aes(x=Longitude, y=Latitude, fill=Mean_temp))+
  geom_sf(color = '#414754', fill = '#414754') + xlab("Longitude") + ylab("Latitude") +
  scale_fill_gradientn('Mean_temp', guide="colorbar", #trans = 'sqrt',
                       colors = c("#4D7DB5","#76C4C9","#CBE7C8","#F6E1A2","#F6AC6D","#DF4D51")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
  #                       style = north_arrow_nautical) +
  coord_sf(xlim = c(113, 115), ylim = c(-24.5, -21.5), expand = FALSE) + #edge south is mid way between Lack Mcleod
  theme_bw() + facet_wrap(.~Variable)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

Ninglaoo.Temp.plot


#Shark Bay Plot ----

#This is the full plotting area
boundary.S = extent(c(xmin = 111, xmax = 115, ymin = -28, ymax = -24.5))

temp.rasters4 = crop(temp.rasters, boundary.S)

cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Plot mean bottom mean temperature
rast.s.p <- raster(temp.rasters4, layer=3)
rast.s.50 <- raster(temp.rasters4, layer=9)
rast.s.100 <- raster(temp.rasters4, layer=10) 

crs(rast.s.p) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.p.s <- rasterToPoints(rast.s.p)

crs(rast.s.50) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.50.s <- rasterToPoints(rast.s.50)

crs(rast.s.100) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.100.s <- rasterToPoints(rast.s.100)

#make a dataframe of points for ggplot
df.p.s <- data.frame(map.p.s)

df.50.s <- data.frame(map.50.s)

df.100.s <- data.frame(map.100.s)

#make appropriate column headings
colnames(df.p.s) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.50.s) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.100.s) <- c('Longitude', 'Latitude', 'Mean_temp')

df.p.s$Variable<-"Present"
df.50.s$Variable<-"2050"
df.100.s$Variable<-"2100"
df.plot.Shark<-rbind(df.p.s,df.50.s, df.100.s)

df.plot.Shark$Variable <- factor(df.plot.Shark$Variable, levels = c("Present", "2050", "2100"))

Shark.Temp.plot<-ggplot(data = world) +
  geom_tile(data=df.plot.Shark, aes(x=Longitude, y=Latitude, fill=Mean_temp))+
  geom_sf(color = '#414754', fill = '#414754') + xlab("Longitude") + ylab("Latitude") +
  scale_fill_gradientn('Mean_temp', guide="colorbar", #trans = 'sqrt',
                       colors = c("#4D7DB5","#76C4C9","#CBE7C8","#F6E1A2","#F6AC6D","#DF4D51")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #                       pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
  #                       style = north_arrow_nautical) +
  coord_sf(xlim = c(111, 115), ylim = c(-28, -24.5), expand = FALSE) + #edge south is mid way between Lack Mcleod
  theme_bw() + facet_wrap(.~Variable)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

Shark.Temp.plot

#Great Barrier Reef ----

boundary.G = extent(c(xmin = 142, xmax = 152, ymin = -25, ymax = -10.5))

# Crop rasters to boundary extent

temp.rasters5 = crop(temp.rasters, boundary.G)

cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Plot mean bottom mean temperature
rast.g.p <- raster(temp.rasters5, layer=3) 
rast.g.50 <- raster(temp.rasters5, layer=9)
rast.g.100 <- raster(temp.rasters5, layer=10) 

crs(rast.g.p) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.p.g <- rasterToPoints(rast.g.p)

crs(rast.g.50) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.50.g <- rasterToPoints(rast.g.50)

crs(rast.g.100) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
map.100.g <- rasterToPoints(rast.g.100)

df.p.g <- data.frame(map.p.g)

df.50.g <- data.frame(map.50.g) 

df.100.g <- data.frame(map.100.g)

#make appropriate column headings
colnames(df.p.g) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.50.g) <- c('Longitude', 'Latitude', 'Mean_temp')
colnames(df.100.g) <- c('Longitude', 'Latitude', 'Mean_temp')

df.p.g$Variable<-"Present"
df.50.g$Variable<-"2050"
df.100.g$Variable<-"2100"
df.plot.GBR<-rbind(df.p.g,df.50.g, df.100.g)

df.plot.GBR$Variable <- factor(df.plot.GBR$Variable, levels = c("Present", "2050", "2100"))

GBR.Temp.plot<-ggplot(data = world) +
  geom_tile(data=df.plot.GBR, aes(x=Longitude, y=Latitude, fill=Mean_temp))+
  geom_sf(color = '#414754', fill = '#414754') + xlab("Longitude") + ylab("Latitude") +
  scale_fill_gradientn('Mean_temp', guide="colorbar", #trans = 'sqrt',
                       colors = c("#4D7DB5","#76C4C9","#CBE7C8","#F6E1A2","#F6AC6D","#DF4D51")) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #pad_x = unit(0.1, "in"), pad_y = unit(0.25, "in"),
  #style = north_arrow_nautical) +
  coord_sf(xlim = c(142, 152), ylim = c(-25, -10.5), expand = FALSE) + #edge south is mid way between Lack Mcleod
  theme_bw() + facet_wrap(.~Variable)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


GBR.Temp.plot

#Main plot together----
#base plot for Figure 1 -----
library(cowplot)
WHS_plot <- plot_grid(
  LordHowe.Temp.plot +
    theme(legend.position="top", #theme(legend.position="none"
          plot.margin = unit(rep(0.5, 4), "cm")),
  Ninglaoo.Temp.plot +
    theme(legend.position="top",
          plot.margin = unit(rep(0.5, 4), "cm")),
  Shark.Temp.plot+
    theme(legend.position="top",
          plot.margin = unit(rep(0.5, 4), "cm")),
  GBR.Temp.plot +
    theme(legend.position="top",
          plot.margin = unit(rep(0.5, 4), "cm")),
  align = 'hv',
  labels = c("A.", "B.", 'C.', "D."),
  #hjust = -1,
  #vjust = 0,
  nrow = 2,
  label_size = 18
)


#Filter by Park boundaries----
LordHowe.Filt.boundries <- df.plot.LordHowe  %>%
  # Filter the data based on specific Longitude and Latitude values
  filter(Longitude >= 158.9583 & Longitude <= 159.1250 &
           Latitude >= -31.62500 & Latitude <= -31.45833) 

LordHowe.Filt.boundries %>%
  # Group by Region and Variable, then summarize the data
  group_by(Variable) %>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

####
Ningaloo.Filt.boundries <- df.plot.Ning  %>%
  # Filter the data based on specific Longitude and Latitude values
  filter(Longitude >= 113.39 & Longitude <= 114.412 &
           Latitude >= -24.5 & Latitude <= -21.0) 

Ningaloo.Filt.boundries %>%
  # Group by Region and Variable, then summarize the data
  group_by(Variable) %>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

Shark.Bay.Filt.boundries <- df.plot.Shark  %>%
  # Filter the data based on specific Longitude and Latitude values
  filter(Longitude >= 112.5 & Longitude <= 114.5 &
           Latitude >= -27.00 & Latitude <= -25.5) 

Shark.Bay.Filt.boundries %>%
  # Group by Region and Variable, then summarize the data
  group_by(Variable) %>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

GBR.Filt.boundries <- df.plot.GBR  %>%
  # Filter the data based on specific Longitude and Latitude values
  filter(Longitude >= 142 & Longitude <= 152.89 &
           Latitude >= -24.2 & Latitude <= -9) 

GBR.Filt.boundries %>%
  # Group by Region and Variable, then summarize the data
  group_by(Variable) %>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

######

LordHowe.Filt.boundries$Region<-"Lord.Howe"
Ningaloo.Filt.boundries$Region<-"Ningaloo"
GBR.Filt.boundries$Region<-"GBR"
Shark.Bay.Filt.boundries$Region<-"Shark.Bay"
df.All<-rbind(LordHowe.Filt.boundries, Ningaloo.Filt.boundries, GBR.Filt.boundries, Shark.Bay.Filt.boundries)

df.All %>% group_by(Region, Variable)%>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

df.All$Region <- factor(df.All$Region, levels = c("GBR", "Ningaloo", "Shark.Bay", "Lord.Howe"))

df.All.Present<-df.All %>%filter(Variable == "Present") 
df.All.2050<-df.All %>%filter(Variable == "2050")
df.All.2100<-df.All %>%filter(Variable == "2100")

my_comparisions.region<-list( c("GBR", "Ningaloo"), c("GBR", "Shark.Bay"), 
                              c("Ningaloo", "Shark.Bay"),c("Lord.Howe", "Shark.Bay"),
                              c("Lord.Howe", "Ningaloo"),c("GBR", "Lord.Howe"))
#library(ggpubr)

df.All.Present %>% group_by(Region, Variable)%>%
  summarise(lowest = min(Mean_temp), highest = max(Mean_temp))

df.All.2050%>% group_by(Region, Variable)%>%
  summarise(lowest = min(Mean_temp), mean= mean(Mean_temp), highest = max(Mean_temp))

df.All.2100%>% group_by(Region, Variable)%>%
  summarise(lowest = min(Mean_temp), mean= mean(Mean_temp), highest = max(Mean_temp))

# Subset data for "Present" and "2050"
present_data <- subset(df.All, Variable == "Present")
future_2050_data <- subset(df.All, Variable == "2050")

# Merge data for "Present" and "2050" on Longitude, Latitude, and Region
merged_2050_data <- merge(present_data, future_2050_data, by = c("Longitude", "Latitude", "Region"), suffixes = c("_present", "_2050"))

# Calculate relative change in Mean_temp for 2050
merged_2050_data$Relative_change_2050 <- ((merged_2050_data$Mean_temp_2050 - merged_2050_data$Mean_temp_present))

# Subset data for "2100"
future_2100_data <- subset(df.All, Variable == "2100")

# Merge data for "Present" and "2100" on Longitude, Latitude, and Region
merged_2100_data <- merge(present_data, future_2100_data, by = c("Longitude", "Latitude", "Region"), suffixes = c("_present", "_2100"))

# Calculate relative change in Mean_temp for 2100
merged_2100_data$Relative_change_2100 <- ((merged_2100_data$Mean_temp_2100 - merged_2100_data$Mean_temp_present))

# Print the results
print("Relative change in Mean_temp for 2050:")
print(merged_2050_data)

print("Relative change in Mean_temp for 2100:")
print(merged_2100_data)

GBR_MMM_2050<-merged_2050_data %>%
  filter(Mean_temp_2050 >= 28.4 & Region == "GBR") 
SB_MMM_2050<-merged_2050_data %>%
  filter(Mean_temp_2050 >= 24.5 & Region == "Shark.Bay") 
LH_MMM_2050<-merged_2050_data %>%
  filter(Mean_temp_2050 >= 24.5 & Region == "Lord.Howe") 
Ning_MMM_2050<-merged_2050_data %>%
  filter(Mean_temp_2050 >= 27.5 & Region == "Ningaloo") 

#Relative_change_2050
rel.ch.50.plot<-ggplot(merged_2050_data, aes(y = Relative_change_2050, x = Mean_temp_2050, fill = Region, color = Region)) +
  geom_point(size = 3, shape = 21, alpha = 0.2) +
  labs(
    y = "Relative Change (%)",
    x = "Mean Temperature (°C) 2050") +
  #theme_minimal(base_size = 11) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 360, hjust = 1, size = 9)) +
  scale_fill_manual(values = c("#0e3a5e", "#53b2a8", "#dd9e5f", "#d1c8c1")) +
  scale_color_manual(values = c("#000000", "#000000", "#000000", "#000000")) +
  guides(color = guide_legend(title = "Region")) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 28.4, color = "#0e3a5e", linetype = "dashed", size = 0.5) +  
  geom_vline(xintercept = 24.5, color = "#d1c8c1", linetype = "dashed", size = 0.5) +  
  geom_vline(xintercept = 24.6, color = "#dd9e5f", linetype = "dashed", size = 0.5) +  #this value is technically 24.5 but I have jittered it a bit so it is not completely overlapping w the other line
  geom_vline(xintercept = 27.5, color = "#53b2a8", linetype = "dashed", size = 0.5)+
  ylim(0.7, 1.5)+ xlim(21, 29.5)+
  geom_point(aes(x = 27.54725, y = 0.8016706), size = 3, shape = 23, fill = "#2774b3", color = "black")+
  geom_point(aes(x = 22.66526, y = 1.0635547), size = 3, shape = 23, fill = "#d1c8c1", color = "black")+
  geom_point(aes(x = 24.33001, y = 1.1702215), size = 3, shape = 23, fill = "#dd9e5f", color = "black")+
  geom_point(aes(x = 26.61903, y = 1.0626179), size = 3, shape = 23, fill = "#53b2a8", color = "black")+
  
  geom_errorbar(aes(x = 27.54725,y = 0.8016706, xmin =26.67359, xmax =28.42091), width = 0, color = "black")+ #size = 2.5
  geom_errorbar(aes(x = 22.66526, y = 1.0635547, xmin =26.10442, xmax =27.13364), width = 0, color = "black")+
  geom_errorbar(aes(x = 24.33001, y = 1.1702215, xmin =23.95488, xmax =24.70514), width = 0, color = "black")+
  geom_errorbar(aes(x = 26.61903, y = 1.0626179, xmin =22.63598, xmax =22.69454), width = 0, color = "black")+
  
  geom_errorbar(aes(x = 27.54725, y = 0.8016706, ymin =0.7596897, ymax =0.8436515), width = 0, color = "black")+
  geom_errorbar(aes(x = 22.66526, y = 1.0635547, ymin =1.047932, ymax =1.077304), width = 0, color = "black")+
  geom_errorbar(aes(x = 24.33001, y = 1.1702215, ymin =1.160906, ymax =1.179537), width = 0, color = "black")+
  geom_errorbar(aes(x = 26.61903, y = 1.0626179, ymin =1.051374, ymax =1.075735), width = 0, color = "black")

###
#% MMMs 2100----
GBR_MMM_2100<-merged_2100_data %>%
  filter(Mean_temp_2100 >= 28.4 & Region == "GBR") #2969 of GBR original 11,216 = 26.47111
SB_MMM_2100<-merged_2100_data %>%
  filter(Mean_temp_2100 >= 24.5 & Region == "Shark.Bay") #141 of original 309 = 45.6%
dim(SB_MMM_2100) #141
LH_MMM_2100<-merged_2100_data %>%
  filter(Mean_temp_2100 >= 24.5 & Region == "Lord.Howe") # 0 of original
Ning_MMM_2100<-merged_2100_data %>%
  filter(Mean_temp_2100 >= 27.5 & Region == "Ningaloo") # 11 of original 263 at #4.18%

rel.ch.100.plot<-ggplot(merged_2100_data, aes(y = Relative_change_2100, x = Mean_temp_2100, fill = Region, color = Region)) +
  geom_point(size = 3, shape = 21, alpha = 0.2) +
  labs(
    y = "Relative Change (%)",
    x = "Mean Temperature (°C) 2100") +
  #theme_minimal(base_size = 11) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 360, hjust = 1, size = 9)) +
  scale_fill_manual(values = c("#0e3a5e", "#53b2a8", "#dd9e5f", "#d1c8c1")) +
  scale_color_manual(values = c("#000000", "#000000", "#000000", "#000000")) +
  guides(color = guide_legend(title = "Region")) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 28.4, color = "#0e3a5e", linetype = "dashed", size = 0.5) +  
  geom_vline(xintercept = 24.6, color = "#d1c8c1", linetype = "dashed", size = 0.5) + #jittered, see comment above 
  geom_vline(xintercept = 24.5, color = "#dd9e5f", linetype = "dashed", size = 0.5) +  
  geom_vline(xintercept = 27.5, color = "#53b2a8", linetype = "dashed", size = 0.5) +
  ylim(0.7, 1.5)+ xlim(21, 29.5)+
  geom_point(aes(x = 27.77073, y = 1.025159), size = 3, shape = 23, fill = "#2774b3", color = "black")+
  geom_point(aes(x = 26.79884, y = 1.242432), size = 3, shape = 23, fill = "#53b2a8", color = "black")+
  geom_point(aes(x = 24.45298, y = 1.293195), size = 3, shape = 23, fill = "#dd9e5f", color = "black")+
  geom_point(aes(x = 22.72493, y = 1.123225), size = 3, shape = 23, fill = "#d1c8c1", color = "black")+
  
  geom_errorbar(aes(x = 27.77073,y = 1.025159, xmin =26.89476, xmax =28.6467), width = 0, color = "black")+ #size = 2.5
  geom_errorbar(aes(x = 26.79884, y = 1.242432, xmin =26.27837, xmax =27.31931), width = 0, color = "black")+
  geom_errorbar(aes(x = 24.45298, y = 1.293195, xmin =24.07634, xmax =24.82962), width = 0, color = "black")+
  geom_errorbar(aes(x = 22.72493, y = 1.123225, xmin =22.70151, xmax =22.74835), width = 0, color = "black")+
  
  geom_errorbar(aes(x = 27.77073, y = 1.025159, ymin =0.9864271, ymax =1.063891), width = 0, color = "black")+
  geom_errorbar(aes(x = 26.79884, y = 1.242432, ymin =1.236326, ymax =1.248538), width = 0, color = "black")+
  geom_errorbar(aes(x = 24.45298, y = 1.293195, ymin =1.284712, ymax =1.301678), width = 0, color = "black")+
  geom_errorbar(aes(x = 22.72493, y = 1.123225, ymin =1.103935, ymax =1.142515), width = 0, color = "black")

Mean_byRelChange<-plot_grid(
  rel.ch.50.plot +
    theme(legend.position="top", #theme(legend.position="none"
          plot.margin = unit(rep(0.5, 4), "cm")),
  rel.ch.100.plot +
    theme(legend.position="top",
          plot.margin = unit(rep(0.5, 4), "cm")),
  align = 'hv',
  labels = c("A.", "B."),
  #hjust = -1,
  #vjust = 0,
  nrow = 1,
  label_size = 18
)

####
# Combine merged_2050_data and merged_2100_data into one data frame
combined_data <- merge(merged_2050_data, merged_2100_data, by = c("Longitude", "Latitude", "Region"), suffixes = c("_2050", "_2100"))
combined_data <- merge(combined_data,present_data, by = c("Longitude", "Latitude", "Region"))

# Create a new column "relative.change" by combining relative change values from both columns

result <- combined_data %>%
  gather(key = "Variable", value = "Mean.Temp.Total", Mean_temp, Mean_temp_present_2050, Mean_temp_present_2100)

result2 <- result %>%
  gather(key = "Variable2", value = "Rel.Change.Total", Relative_change_2050, Relative_change_2100) 

#library(ggridges)

Result2_Rel.Change<-ggplot(result2, aes(x = Rel.Change.Total, y = Variable2, fill = Region)) +
  geom_density_ridges(alpha = 0.5) +
  labs(title = "Ridgeline Plot of Rel.Change.Total by Geographic and Variable2",
       x = "Rel.Change.Total",
       y = "Variable2")+theme_bw()+
  theme(axis.text.x = element_text(angle = 360, hjust = 1, size = 9)) +
  scale_fill_manual(values = c("#0e3a5e", "#53b2a8", "#dd9e5f", "#d1c8c1"))

# Assuming 'result2' is your data frame
# Perform ANOVA for each level of Variable2
anova_results <- lapply(levels(result2$Variable2), function(level) {
  subset_data <- result2[result2$Variable2 == level, ]
  aov_result <- aov(Rel.Change.Total ~ Region, data = subset_data)
  return(summary(aov_result))
})

# Print ANOVA results
print(anova_results)

