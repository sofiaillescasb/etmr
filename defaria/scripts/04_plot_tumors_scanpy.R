





etmr_cell_plot <- ggplot(all_umap %>% filter(condition == "etmr"), aes(x = UMAP_1, y = UMAP_2, color = CellClass)) +
  geom_point(alpha = 0.5, size = 0.5) +
  theme_classic() 

all_meta %>%
  filter(condition== "etmr") %>%
  mutate(Age = ordered(as.character(Age), levels = c("5", "8", "11.5", "14"))) %>%
  arrange(Age) %>%
  ggplot(aes(x = leiden, fill = Age)) +
  geom_bar(position = "fill", alpha= .7,color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme_classic() +
  ggtitle("ETMR data: scanpy integration")

all_meta %>%
  mutate(Age = ordered(as.character(Age), levels = c("5", "8", "11.5", "14"))) %>%
  arrange(Age) %>%
  ggplot(aes(x = leiden, fill = Age)) +
  geom_bar(position = "fill", color = "darkgrey")  +
  scale_fill_manual(values = c("#c9cba3","#ffe1a8", "#e26d5c", "#723d46")) +
  theme_classic() +
  ggtitle("Reference and ETMR: scanpy integration")


all_meta %>%
  filter(condition== "etmr") %>%
  ggplot(aes(x = leiden,  fill = CellClass)) +
  geom_bar(position = "fill",  color = "gray") +
  theme_classic() +
  scale_fill_manual(values = c("#ff9aa2", "#ffb7b2", "#ffdac1", "#e2f0cb", "#b5ead7", "#c7ceea"))  +
  ggtitle("ETMR data: scanpy integration")
