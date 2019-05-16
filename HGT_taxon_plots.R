library(tidyverse)

go_rsoc = read_csv("rsoc_GO.csv")

plot_data = go_rsoc %>%
  count(`best.HGT.taxon`)

ggplot(plot_data,
       aes(x=reorder(`best.HGT.taxon`, n), y = n, fill = n)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme_classic() +
  labs(x = "Most likely origin of horizontally transferred gene",
       y = "Count of genes",
       fill = "Count of genes")


go_rmag = read_csv("rmag_GO.csv")

plot_data = go_rmag %>%
  count(`best.HGT.taxon`)

ggplot(plot_data,
       aes(x=reorder(`best.HGT.taxon`, n), y = n, fill = n)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "yellow", high = "red") +
  theme_classic() +
  labs(x = "Most likely origin of horizontally transferred gene",
       y = "Count of genes",
       fill = "Count of genes")

