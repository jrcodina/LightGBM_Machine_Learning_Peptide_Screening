best <- read.csv("./20230717_Results_ADCP_Best.csv")
worse <- read.csv("./20230717_Results_ADCP_Worse.csv")
initial <- read.csv("./20230706_Results_ADCP_Random.csv")

# Load the ggplot2 package
library(ggplot2)

# Create a new data frame, marking the source of each row
best$Group <- 'Most Probable (8,000)'
worse$Group <- 'Less Probable (8,000)'
initial$Group <- 'Random Sample (1,600)'

# Combine the data frames
combined <- rbind(best, worse, initial)
combined <- na.omit(combined[,-4])

# Plot the histograms
p <- ggplot(combined, aes(x=affinity, fill=Group)) +
  geom_histogram(position="identity", alpha=0.7, bins=30) +
  labs(x = "Docking Score", y = "Count") +
  scale_fill_manual(values = c('red', 'steelblue', 'black'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        text = element_text(size = 36),
        axis.title = element_text(size = 36),
        axis.text = element_text(size = 36, color = "black"),
        panel.border = element_rect(linewidth = 3),
        panel.grid = element_line(color = "gray"))


ggsave(plot = p, device = "png", "./20230717_Figure_Hist_ADCP.png",width = 15, height = 8.75, dpi = 600)


# Is our data normally distributed? - YES
ggplot(combined, aes(sample = affinity)) +
  facet_wrap(~Group) +
  stat_qq() +
  stat_qq_line()

# Have our groups the same variance? - YES
var(na.rm = T, worse$affinity)/var(na.rm = T, best$affinity) # Less than 4

library(stats)

# Calculate the significance
t_test_result <- t.test(na.omit(worse$affinity), na.omit(best$affinity), var.equal = TRUE)
t_test_result


reference <- initial[order(initial$affinity),]

reference_best <- max(reference[1:(1600*.2),3]) #-9.2

