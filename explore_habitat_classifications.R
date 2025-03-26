library(phytools)
library(geiger)

classification=read.table("final_classification_k_7.csv", header = FALSE, sep = ",")
rownames(classification) <- classification[,1]
classification$V1<- NULL

environment=read.table("medians_data.csv", header = FALSE, sep = ",")
rownames(environment) <- environment[,1]
environment$V1<- NULL

tribe=read.table("tribe_classification.csv", header = FALSE, sep = ",")
rownames(tribe) <- tribe[,1]
tribe$V1<- NULL

combined <- merge(classification, environment, by="row.names", all=TRUE)
combined <- na.omit(combined)


combined[,2] <- as.factor(as.numeric(combined[,2]))
names(combined) <- c("species", "habitat", "BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9", "BIO10", "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19", "aspect", "elevation", "slope", "bulk_density", "clay_percent", "coarse_fragment_percent", "phx10", "sand_percent", "siltpercent", "soil_organic_content", "LandCover_1_needle_leaf", "LandCover_2_evergreen", "LandCover_3_deciduous", "LandCover_4_mixed", "LandCover_5_shrub", "LandCover_6_herb")
rownames(combined) <- combined[,1]
combined$species <- NULL

combined <- merge(combined, tribe, by="row.names", all=TRUE)
combined <- na.omit(combined)

ag <- aggregate(. ~ habitat, combined, function(x) c(mean = mean(x), sd = sd(x)))
write.csv(ag, file = "mean_stdev.csv")

# Tribe plots
ggplot(aes(y = as.numeric(BIO1/10), x = as.factor(V2)), data = combined) + geom_violin() + ylab("Mean annual temperature")


library(ggplot2)
# Mean annual temperature
p1 <- ggplot(aes(y = as.numeric(BIO1/10), x = as.factor(V2)), data = combined) + geom_violin() + ylab("Mean annual temperature") + xlab("") + geom_boxplot(width=0.1, fill = "gray") + coord_flip()
# Temperature seasonality
p2 <- ggplot(aes(y = as.numeric(BIO4/10), x = as.factor(V2)), data = combined) + geom_violin() + ylab("Temperature seasonality") + xlab("") + geom_boxplot(width=0.1, fill = "gray") + coord_flip()
# Precipitation of wettest month
p3 <- ggplot(aes(y = as.numeric(BIO12), x = as.factor(V2)), data = combined) + geom_violin() + ylab("Mean annual precipitation") + xlab("") + geom_boxplot(width=0.1, fill = "gray") + coord_flip()
# Elevation
p4 <- ggplot(aes(y = as.numeric(elevation), x = as.factor(V2)), data = combined) + geom_violin() + ylab("Elevation") + xlab("") + geom_boxplot(width=0.1, fill = "gray") + coord_flip()



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(p1, p2, p3, p4, cols=2)