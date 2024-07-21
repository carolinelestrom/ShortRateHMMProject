####################################################################################################################
####################################################################################################################
#------------------------------------------ Plotting HMM -----------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Dataframe
Rate3M <- as.numeric(yields_df$"3M")
M3data <- as.data.frame(Rate3M)
M3data$DateCont <- yields_df$DateCont
M3data$Date <- yields_df$Date
head(M3data)




####################################################################################################################
####################################################################################################################
#----------------------------------------- 2-state Plots -----------------------------------------------------------
####################################################################################################################
####################################################################################################################


M3data$Viterbi2 <- as.factor(viterbi2)

### Scatter plot with Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = factor(Viterbi2)), show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("2-state HMM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_manual(values=c("#CC5500", "#404080")) +
  annotate(geom="text", x=2016.7, y=0.11, label=paste("State 1 : 42.69%"),
           color="#CC5500", size = 17) +
  annotate(geom="text", x=2016.7, y=0.10, label=paste("State 2 : 57.31%"),
           color="#404080", size = 17)


M3data$Viterbi2Plot <- b2_hat[M3data$Viterbi2]

### Plotting Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", alpha = 0.3, size = 3) +
  #geom_point(aes(x = DateCont, y = Viterbi2Plot), alpha = 1, color = "#404080", size = 1.7) +
  geom_line(aes(x = DateCont, y = Viterbi2Plot), alpha = 1, color = "#404080", size = 1.7, linetype = "longdash") +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("2-state HMM w/ 3M Rates") +
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
dens2_1 <- delta2_hat[1] * dnorm(point, b2_hat[1], sqrt(sigma2_hat[1]^2 / (2 * a2_hat[1])))
dens2_2 <- delta2_hat[2] * dnorm(point, b2_hat[2], sqrt(sigma2_hat[2]^2 / (2 * a2_hat[2])))

dens2_mix <- dens2_1 + dens2_2
plot_dens2 <- data.frame(x = point, dens2_1 = dens2_1, dens2_2 = dens2_2,
                         dens2_mix = dens2_mix)


ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6",alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#FF8C00", "#404080", "#87A96B", "#FFBCD9")) +
  theme_bw() +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("2-state HMM w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = plot_dens2, aes(x = x, y = dens2_1), color = "#FF8C00", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens2, aes(x = x, y = dens2_2), color = "#404080", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens2, aes(x = x, y = dens2_mix), color = "black", linetype = "dotdash", size = 1.7, inherit.aes = FALSE) +
  annotate(geom="text", x=0.0979, y=100, label=paste("State 1 -"),
           color="#FF8C00", size = 17) +
  annotate(geom="text", x=0.0979, y=90, label=paste("State 2 -"),
           color="#404080", size = 17) +
  annotate(geom="text", x=0.103, y=80, label=paste("Fitted Model"),
           color="black", size = 17)








####################################################################################################################
####################################################################################################################
#----------------------------------------- 3-state Plots -----------------------------------------------------------
####################################################################################################################
####################################################################################################################




### Dataframe
M3data$Viterbi3 <- as.factor(viterbi3)
head(M3data)


### Scatter plot with Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = factor(Viterbi3)), show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("3-state HMM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B")) +
  annotate(geom="text", x=2016.7, y=0.11, label=paste("State 1 : 22.26%"),
           color="#CC5500", size = 17) +
  annotate(geom="text", x=2016.7, y=0.10, label=paste("State 2 : 24.22%"),
           color="#404080", size = 17) +
  annotate(geom="text", x=2016.7, y=0.09, label=paste("State 3 : 53.52%"),
           color="#87A96B", size = 17)


M3data$Viterbi3Plot <- b3_hat[M3data$Viterbi3]

### Plotting Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", alpha = 0.3, size = 3) +
  #geom_point(aes(x = DateCont, y = Viterbi3Plot), alpha = 1, color = "#404080", size = 1.7) +
  geom_line(aes(x = DateCont, y = Viterbi3Plot), alpha = 1, color = "#404080", size = 1.7, linetype = "longdash") +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("3-state HMM w/ 3M Rates") +
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
dens3_1 <- delta3_hat[1] * dnorm(point, b3_hat[1], sqrt(sigma3_hat[1]^2 / (2 * a3_hat[1])))
dens3_2 <- delta3_hat[2] * dnorm(point, b3_hat[2], sqrt(sigma3_hat[2]^2 / (2 * a3_hat[2])))
dens3_3 <- delta3_hat[3] * dnorm(point, b3_hat[3], sqrt(sigma3_hat[3]^2 / (2 * a3_hat[3])))

dens3_mix <- dens3_1 + dens3_2 + dens3_3
plot_dens3 <- data.frame(x = point, dens3_1 = dens3_1, dens3_2 = dens3_2, dens3_3 = dens3_3,
                         dens3_mix = dens3_mix)


ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6",alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#FF8C00", "#404080", "#87A96B", "#FFBCD9")) +
  theme_bw() +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("3-state HMM w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = plot_dens3, aes(x = x, y = dens3_1), color = "#FF8C00", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens3, aes(x = x, y = dens3_2), color = "#404080", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens3, aes(x = x, y = dens3_3), color = "#87A96B", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens3, aes(x = x, y = dens3_mix), color = "black", linetype = "dotdash", size = 1.7, inherit.aes = FALSE) +
  annotate(geom="text", x=0.0979, y=100, label=paste("State 1 -"),
           color="#FF8C00", size = 17) +
  annotate(geom="text", x=0.0979, y=90, label=paste("State 2 -"),
           color="#404080", size = 17) +
  annotate(geom="text", x=0.0979, y=80, label=paste("State 3 -"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=0.103, y=70, label=paste("Fitted Model"),
           color="black", size = 17)






####################################################################################################################
####################################################################################################################
#----------------------------------------- 4-state Plots -----------------------------------------------------------
####################################################################################################################
####################################################################################################################




### Dataframe
M3data$Viterbi4 <- as.factor(viterbi4)
head(M3data)


### Scatter plot with Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = factor(Viterbi4)), show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("4-state HMM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_manual(values=c("#CC5500", "#404080", "#FFBCD9", "#87A96B")) +
  annotate(geom="text", x=2016.7, y=0.11, label=paste("State 1 : 20.87%"),
           color="#CC5500", size = 17) +
  annotate(geom="text", x=2016.7, y=0.10, label=paste("State 2 : 30.78%"),
           color="#404080", size = 17) +
  annotate(geom="text", x=2016.7, y=0.09, label=paste("State 3 : 31.51%"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=2016.7, y=0.08, label=paste("State 4 : 16.84%"),
           color="#FFBCD9", size = 17)




M3data$Viterbi4Plot <- b4_hat[M3data$Viterbi4]

### Plotting Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", alpha = 0.3, size = 3) +
  geom_point(aes(x = DateCont, y = Viterbi4Plot), alpha = 1, color = "#404080", size = 1.7) +
  #geom_line(aes(x = DateCont, y = Viterbi4Plot), alpha = 1, color = "#404080", size = 1.7, linetype = "longdash") +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("4-state HMM w/ 3M Rates") +
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
dens4_1 <- delta4_hat[1] * dnorm(point, b4_hat[1], sqrt(sigma4_hat[1]^2 / (2 * a4_hat[1])))
dens4_2 <- delta4_hat[2] * dnorm(point, b4_hat[2], sqrt(sigma4_hat[2]^2 / (2 * a4_hat[2])))
dens4_3 <- delta4_hat[3] * dnorm(point, b4_hat[3], sqrt(sigma4_hat[3]^2 / (2 * a4_hat[3])))
dens4_4 <- delta4_hat[4] * dnorm(point, b4_hat[4], sqrt(sigma4_hat[4]^2 / (2 * a4_hat[4])))

dens4_mix <- dens4_1 + dens4_2 + dens4_3 + dens4_4
plot_dens4 <- data.frame(x = point, dens4_1 = dens4_1, dens4_2 = dens4_2, dens4_3 = dens4_3, dens4_4 = dens4_4, 
                         dens4_mix = dens4_mix)


ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6",alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#FF8C00", "#404080", "#87A96B", "#FFBCD9")) +
  theme_bw() +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("4-state HMM w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = plot_dens4, aes(x = x, y = dens4_1), color = "#FF8C00", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens4, aes(x = x, y = dens4_2), color = "#404080", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens4, aes(x = x, y = dens4_3), color = "#FFBCD9", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens4, aes(x = x, y = dens4_4), color = "#87A96B", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens4, aes(x = x, y = dens4_mix), color = "black", linetype = "dotdash", size = 1.7, inherit.aes = FALSE) +
  annotate(geom="text", x=0.0979, y=100, label=paste("State 1 -"),
           color="#FF8C00", size = 17) +
  annotate(geom="text", x=0.0979, y=90, label=paste("State 2 -"),
           color="#404080", size = 17) +
  annotate(geom="text", x=0.0979, y=80, label=paste("State 3 -"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=0.0979, y=70, label=paste("State 4 -"),
           color="#FFBCD9", size = 17) +
  annotate(geom="text", x=0.103, y=60, label=paste("Fitted Model"),
           color="black", size = 17)





####################################################################################################################
####################################################################################################################
#----------------------------------------- 5-state Plots -----------------------------------------------------------
####################################################################################################################
####################################################################################################################




### Dataframe
M3data$Viterbi5 <- as.factor(viterbi5)
head(M3data)


### Scatter plot with Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = factor(Viterbi5)), show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("5-state HMM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_manual(values=c("#CC5500", "#404080", "#996666", "#FFBCD9", "#87A96B")) +
  annotate(geom="text", x=2016.7, y=0.11, label=paste("State 1 : 20.80%"),
           color="#CC5500", size = 17) +
  annotate(geom="text", x=2016.7, y=0.10, label=paste("State 2 : 19.40%"),
           color="#404080", size = 17) +
  annotate(geom="text", x=2016.7, y=0.09, label=paste("State 3 : 13.78%"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=2016.7, y=0.08, label=paste("State 4 : 28.80%"),
           color="#FFBCD9", size = 17) +
  annotate(geom="text", x=2016.7, y=0.07, label=paste("State 5 : 17.22%"),
           color="#996666", size = 17)






M3data$Viterbi5Plot <- b5_hat[M3data$Viterbi5]

### Plotting Viterbi
M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(color = "#87A96B", alpha = 0.3, size = 3) +
  #geom_point(aes(x = DateCont, y = Viterbi5Plot), alpha = 1, color = "#404080", size = 1.7) +
  geom_line(aes(x = DateCont, y = Viterbi5Plot), alpha = 1, color = "#404080", size = 1.7, linetype = "longdash") +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("5-state HMM w/ 3M Rates") +
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
dens5_1 <- delta5_hat[1] * dnorm(point, b5_hat[1], sqrt(sigma5_hat[1]^2 / (2 * a5_hat[1])))
dens5_2 <- delta5_hat[2] * dnorm(point, b5_hat[2], sqrt(sigma5_hat[2]^2 / (2 * a5_hat[2])))
dens5_3 <- delta5_hat[3] * dnorm(point, b5_hat[3], sqrt(sigma5_hat[3]^2 / (2 * a5_hat[3])))
dens5_4 <- delta5_hat[4] * dnorm(point, b5_hat[4], sqrt(sigma5_hat[4]^2 / (2 * a5_hat[4])))
dens5_5 <- delta5_hat[5] * dnorm(point, b5_hat[5], sqrt(sigma5_hat[5]^2 / (2 * a5_hat[5])))

dens5_mix <- dens5_1 + dens5_2 + dens5_3 + dens5_4 + dens5_5
plot_dens5 <- data.frame(x = point, dens5_1 = dens5_1, dens5_2 = dens5_2, 
                         dens5_3 = dens5_3, dens5_4 = dens5_4,
                         dens5_5 = dens5_5,
                         dens5_mix = dens5_mix)


ggplot(M3data, aes(x=Rate3M, y = ..density..)) +
  geom_histogram(bins = 70, color="#4D5D53", fill = "#F0EAD6",alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#FF8C00", "#404080", "#87A96B", "#FFBCD9")) +
  theme_bw() +
  xlab("3M Rates") +
  ylab("Density") +
  ggtitle("5-state HMM w/ 3M Rates") +
  scale_y_continuous(breaks=NULL) +
  coord_cartesian(ylim = c(0, 100)) + # Set y-axis limits without removing data
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_1), color = "#FF8C00", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_2), color = "#404080", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_3), color = "#996666", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_4), color = "#FFBCD9", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_5), color = "#87A96B", size = 1.7, inherit.aes = FALSE) +
  geom_line(data = plot_dens5, aes(x = x, y = dens5_mix), color = "black", linetype = "dotdash", size = 1.7, inherit.aes = FALSE) +
  annotate(geom="text", x=0.0979, y=100, label=paste("State 1 -"),
           color="#FF8C00", size = 17) +
  annotate(geom="text", x=0.0979, y=90, label=paste("State 2 -"),
           color="#404080", size = 17) +
  annotate(geom="text", x=0.0979, y=80, label=paste("State 3 -"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=0.0979, y=70, label=paste("State 4 -"),
           color="#FFBCD9", size = 17) +
  annotate(geom="text", x=0.0979, y=60, label=paste("State 5 -"),
           color="#996666", size = 17) +
  annotate(geom="text", x=0.103, y=50, label=paste("Fitted Model"),
           color="black", size = 17)







### Scatter plot with Viterbi and states for prediction
nervous <- M3data$DateCont[9535]
calm <- M3data$DateCont[9200]
moderate <- M3data$DateCont[8760]

M3data %>%
  ggplot(aes(x = DateCont, y = Rate3M)) +
  geom_point(aes(color = factor(Viterbi5)), show.legend = FALSE) +
  theme_bw() +
  xlab("Time") +
  ylab("3M Rates") +
  ggtitle("5-state HMM w/ 3M Rates") +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  scale_x_continuous(breaks = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  scale_y_continuous(breaks=c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10), limits = c(0, 0.11)) +
  scale_color_manual(values=c("#CC5500", "#404080", "#996666", "#FFBCD9", "#87A96B")) + 
  geom_vline(xintercept = nervous, col = "#CB4154", linewidth = 1.7, linetype = "longdash") + 
  geom_vline(xintercept = moderate, col = "#004225", linewidth = 1.7, linetype = "longdash") + 
  geom_vline(xintercept = calm, col = "#5F9EA0", linewidth = 1.7, linetype = "longdash") +
  annotate(geom="text", x=(2015), y=0.09, label=paste("~ 2022"),
           color="#CB4154", size = 17) +
  annotate(geom="text", x=(2015), y=0.11, label=paste("~ 2021"),
           color="#5F9EA0", size = 17) +
  annotate(geom="text", x=2015, y=0.10, label=paste("~ 2019"),
           color="#004225", size = 17) +
  annotate(geom="text", x=1997, y=0.11, label=paste("State 1 -"),
           color="#CC5500", size = 17) +
  annotate(geom="text", x=1997, y=0.10, label=paste("State 2 -"),
           color="#404080", size = 17) +
  annotate(geom="text", x=1997, y=0.09, label=paste("State 3 -"),
           color="#87A96B", size = 17) +
  annotate(geom="text", x=1997, y=0.08, label=paste("State 4 -"),
           color="#FFBCD9", size = 17) +
  annotate(geom="text", x=1997, y=0.07, label=paste("State 5 -"),
           color="#996666", size = 17)


