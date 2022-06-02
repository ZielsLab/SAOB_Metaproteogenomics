library(tidyverse)

# ideel test checks for assemblies
# read in all files into one dataframe where a column is the name of the file and hence the name of the assembly



# theme function
theme_min = function (
    size=10, font=NA, face='plain', 
    panelColor=backgroundColor, axisColor='#999999', 
    gridColor=gridLinesColor, textColor='black') 
{   
  theme_text = function(...)
    ggplot2::theme_text(family=font, face=face, colour=textColor, size=size, ...)
  
  opts(
    axis.text.x = theme_text(),
    axis.text.y = theme_text(),
    #axis.line = theme_blank(),
    axis.ticks = theme_segment(colour=axisColor, size=0.25),
    
    panel.border = theme_rect(colour=backgroundColor),
    # panel.border = theme_blank(),
    
    legend.background = theme_blank(),
    legend.key = theme_blank(),
    legend.key.size = unit(1.5, 'lines'),
    legend.text = theme_text(hjust=0),
    legend.title = theme_text(hjust=0),
    
    # panel.background = theme_rect(fill=panelColor, colour=NA),
    panel.background = element_blank(),
    
    # panel.grid.major = theme_line(colour=gridColor, size=0.33),
    panel.grid.major = element_blank(),
    
    # panel.grid.minor = theme_blank(),
    panel.grid.minor = element_blank(),
    
    strip.background = theme_rect(fill=NA, colour=NA),
    strip.text.x = theme_text(hjust=0),
    strip.text.y = theme_text(angle=-90),
    plot.title = theme_text(hjust=0),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'lines'))
}

# plot all together faceting by the input assembly


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

ggplot(racon3medaka3, aes(x=qlen/slen)) + 
  geom_histogram(fill='white', color='grey25', bins=20) +
  xlab('query length / hit length') +
  ylab('frequency') +
  # scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  theme_minimal()

