####################################################################################################################
####################################################################################################################
#---------------------------------------- Plotting Vasicek ---------------------------------------------------------
####################################################################################################################
####################################################################################################################




### Dataframe
Rate3M <- as.numeric(yields_df$"3M")
M3data <- as.data.frame(Rate3M)
M3data$DateCont <- yields_df$DateCont
M3data$Date <- yields_df$Date
head(M3data)



### Scatter plot
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("Standard Vasicek Model w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  annotate(geom="text", x=2016.7, y=0.107, label=paste("Vasicek : 100%"),
           color="#87A96B", size = 17)



### State-dependent distributions
point <- seq(0.00, 0.115, length.out = 1000)
densVasi <- dnorm(point, b_hat, sqrt(sigma_hat^2 / (2 * a_hat)))
plot_Vasi <- data.frame(x = point, densVasi = densVasi)


ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6",alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#FF8C00", "#404080", "#87A96B", "#FFBCD9")) +
  theme_bw() +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("Standard Vasicek Model w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = plot_Vasi, aes(x = x, y = densVasi), color = "black", linetype = "dotdash", size = 1.7, inherit.aes = FALSE) +
  annotate(geom="text", x=0.103, y=100, label=paste("Fitted Model"),
           color="black", size = 17)


