library(tidyverse)
set.seed(20200308)

# blank placeholder plot for now
ggplot() + theme_void() + annotate('text', x=1, y=1, 
                                   label='insert algorithm diagram here', 
                                   size=8, angle=45)
dims <- eval(parse(text=snakemake@params[['dim']]))
ggsave(snakemake@output[['tiff']],
       device = 'tiff', dpi=300,
       width=dims[1], height=dims[2], units='in')
