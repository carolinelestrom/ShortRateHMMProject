####################################################################################################################
####################################################################################################################
#------------------------------------------ Plotting SMM -----------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Dataframe
Rate3M <- as.numeric(yields_df$"3M")
M3data <- as.data.frame(Rate3M)
M3data$DateCont <- yields_df$DateCont
M3data$Date <- yields_df$Date
head(M3data)


M3data$ViterbiCont <- exp(Viterbi_sigma_theta_kappa)*ssmmod_sigma_theta_kappa$estimate[4]

### Scatter plot with Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = ViterbiCont)) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("Continuous SSM w/ 3M Rates") + 
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_gradientn(colours = brewer.pal(n = 8, name = "Pastel2")) +
  #scale_color_gradientn(colours = brewer.pal(n = 11, name = "PRGn")) +
  #scale_color_gradientn(colours = brewer.pal(n = 11, name = "BrBG")) +
  #scale_color_gradientn(colours = brewer.pal(n = 11, name = "PuOr")) +
  #scale_color_gradientn(colours = brewer.pal(n = 11, name = "RdBu")) +
  #scale_color_gradientn(colours = brewer.pal(n = 11, name = "Spectral")) +
  #scale_color_gradient(low = "#997A8D", high = "#87A96B", na.value = NA) +
  labs(color = "States") + # Change legend label
  guides(color = guide_colourbar(barwidth = unit(1, "inches"),  # bar width
                                 barheight = unit(1.7, "inches"),  # bar height
                                 draw.ulim = T,  # show the tick at the upper limit
                                 draw.llim = T)) +  # show the tick at the lower limit
  theme(legend.position = c(0.95, 0.95),  # Adjust the position
        legend.justification = c("right", "top"), # Align the legend box
        legend.title = element_text(size = 30), # Change the size of the legend title
        legend.text = element_text(size = 20)) # Change the size of the legend text








### Plotting Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", alpha = 0.3, size = 3) +
  geom_point(aes(x = DateCont, y = ViterbiCont), alpha = 1, color = "#404080", size = 0.7) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("Continuous SSM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  annotate(geom="text", x=2014.9, y=0.11, label=paste("Observed 3M Rates"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=2009, y=0.10, label=paste("Viterbi Decoded State Sequence"),
           color="#404080", size = 17)









### State-dependent distributions
point <- seq(0.00, 0.115, length.out = 1000)
grid <- seq(from = min(Viterbi_sigma_theta_kappa), to = max(Viterbi_sigma_theta_kappa), length = 200)
delta_plot <- as.numeric(table(factor(Viterbi_sigma_theta_kappa, levels = bstar)) / length(viterbi_sigma_theta_kappa))
prob <- rowSums(matrix(delta_plot, ncol = 1, byrow = TRUE))


theta <- numeric(length(grid))
kappa <- numeric(length(grid))
sigma <- numeric(length(grid))
for (i in (1:length(grid))){
  theta[i] <- exp(grid[i])*ssmmod_sigma_theta_kappa$estimate[4]
  kappa[i] <- exp(grid[i])*exp(ssmmod_sigma_theta_kappa$estimate[3])
  sigma[i] <- exp(grid[i])*exp(ssmmod_sigma_theta_kappa$estimate[5])
}


ContDens <- matrix(0, ncol = 2 + length(grid), nrow = length(point))
for (i in (2:(length(grid) + 1))){
  ContDens[,i] <- prob[i-1] * dnorm(point, theta[i - 1], sqrt(sigma[i - 1]^2 / (2 * kappa[i - 1])))
}
ContDens[,(length(grid) + 2)] <- rowSums(ContDens)
ContDens[,1] <- point
ContDensMix <- as.data.frame(ContDens)



ContDensPlot <- data.frame(Point = point,
                           Dens = prob[1] * dnorm(point, theta[1], sqrt(sigma[1]^2 / (2 * kappa[1]))),
                           State = 1)

for (i in 2:(length(grid))){
  ContDensPlot <- rbind(ContDensPlot, data.frame(Point = point,
                                                 Dens = prob[i] * dnorm(point, theta[i], sqrt(sigma[i]^2 / (2 * kappa[i]))),
                                                 State = i))
}



#mycolors <- colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 8, name = "Pastel2"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 8, name = "Pastel1"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 8, name = "Accent"))(length(grid))
mycolors <- colorRampPalette(brewer.pal(n = 8, name = "Set3"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 11, name = "PRGn"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 11, name = "BrBG"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 11, name = "PuOr"))(length(grid))
#mycolors <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(length(grid))

labs <- c("1", rep("", 48), "50", rep("",49), "100", rep("", 49), "150", rep("", 48), "200") # c("1", rep("", 98), "100", rep("", 98), "199")



ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6", alpha=0.7, position = 'identity') +
  theme_bw() +
  #theme(legend.position = "none") +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("Continuous SSM w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = ContDensPlot, aes(x = Point, y = Dens, color = factor(State)), size = 1, inherit.aes = FALSE, alpha = 0.7) +
  scale_color_manual(values = mycolors, guide = guide_colorsteps(ticks = FALSE,
                                                                 barwidth = unit(1, "inches"),  # bar width
                                                                 barheight = unit(1.7, "inches")),  # bar height)
                     labels = labs) +
  geom_line(data = ContDensMix, aes(x = V1, y = V202), color = "black", size = 1.3, linetype = "dotdash", inherit.aes = FALSE) +
  labs(color = "States") + # Change legend label
  #guides(color = guide_colourbar(barwidth = unit(1, "inches"),  # bar width
  #                               barheight = unit(1.7, "inches"),  # bar height
  #                               draw.ulim = T,  # show the tick at the upper limit
  #                               draw.llim = T)) +  # show the tick at the lower limit
  theme(legend.position = c(0.95, 0.95),  # Adjust the position
        legend.justification = c("right", "top"), # Align the legend box
        legend.title = element_text(size = 30), # Change the size of the legend title
        legend.text = element_text(size = 20)) # Change the size of the legend text









