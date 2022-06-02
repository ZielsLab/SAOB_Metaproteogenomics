library(tidyverse)
library(ggpubr)

# ideel test checks for assemblies
# read in all files into one dataframe where a column is the name of the file and hence the name of the assembly

path <- "results/ideel/"
files <- dir(path, pattern="*.txt")
diamond_tables <- data_frame(filename = files) %>% 
  mutate(file_contents = map(filename, ~ read.table(file.path(path, .)))
    ) %>% 
  unnest(cols = c(file_contents))

names(diamond_tables) <- c("assembly", "qlen", "slen")

ideel_plot <- ggplot(diamond_tables, aes(x=qlen/slen)) + 
  geom_histogram(fill='white', color='grey25', bins=20) +
  xlab('query length / hit length') +
  ylab('frequency') +
  # scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  facet_wrap(~ assembly) +
  theme_bw()

ideel_plot

ggsave("figures/ideel-assembly-comparisons.png", ideel_plot, width=25, height=25, units=c("cm"))

x# raconx3 medakax3 ilm
diamond_tables %>% 
  filter(assembly == "raconx3_medakax3_raconilmx1_contigs.txt") %>% 
  count()

diamond_tables %>% 
  filter(assembly == "raconx3_medakax3_raconilmx1_contigs.txt") %>%
  mutate(proportion = qlen / slen) %>% 
  filter(proportion > .90) %>% 
  count()

49972 / 65029

# raconx3 medakax3 ilm polypolish
diamond_tables %>% 
  filter(assembly == "raconx3_medakax3_raconilmx1_polypolish_contigs.txt") %>% 
  count()

diamond_tables %>% 
  filter(assembly == 'raconx3_medakax3_raconilmx1_polypolish_contigs.txt') %>% 
  mutate(proportion = qlen / slen) %>% 
  filter(proportion > .90) %>% 
  count()

51987 / 63187

# raconx3 medakax3 

racon3medaka3 <- read.table("results/ideel/raconx3_medakax3_contigs.txt", header=FALSE, sep='\t')
names(racon3medaka3) <- c('qlen', 'slen')

racon3medaka3 %>% 
  mutate(proportion = qlen / slen) %>% 
  count()

racon3medaka3 %>% 
  mutate(proportion = qlen / slen)
filter(proportion > .90) %>% 
  count()

45161 / 70342 

pseudogenes <- sum(racon3medaka3$qlen / racon3medaka3$slen)
print(paste0('Encountered genes < 0.9 reference length: ', pseudogenes))
str(racon3medaka3)

p1 <- ggplot(racon3medaka3, aes(x=qlen/slen)) + 
  geom_histogram(fill='white', color='grey25', bins=20) +
  xlab('query length / hit length') +
  ylab('frequency') +
  # scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  theme_minimal() + 
  ggtitle("Raconx3 + Medakax3 (IDEEL = 64%)")

# raconx3 medakax3 polypolish 
racon3medaka3poly <- read.table("results/ideel/raconx3_medakax3_polypolish_contigs.txt", header=FALSE, sep='\t')
names(racon3medaka3poly) <- c("qlen", "slen")

racon3medaka3poly %>% 
  mutate(proportion = qlen / slen) %>% 
  filter(proportion > 0.90) %>% 
  count()

49867 / 65421

p2 <- ggplot(racon3medaka3poly, aes(x=qlen/slen)) + 
  geom_histogram(fill='white', color='grey25', bins=20) +
  xlab('query length / hit length') +
  ylab('frequency') +
  # scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  theme_minimal() +
  ggtitle("Raconx3 + Medakax3 + Polypolish (IDEEL = 76%)")


