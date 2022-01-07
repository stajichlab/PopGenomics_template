library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(stringr)
library(dplyr)

chrlist = read_csv("genome/chrom_nums.csv",col_names=c("Chr","CHR"),col_types="ci")
mosdepthdir = "coverage/mosdepth"

windows = c(5000,10000, 50000)
alternating_colors = rep( c("red", "black"), times = length(chrlist$Chr))

plot_strain <- function(strain, data) {
    l = data %>% filter(Strain == strain)
    Title = sprintf("Chr coverage plot for %s", strain)
    p <- ggplot(l, aes(x = pos, y = NormDepth, color = factor(CHR))) +
        scale_colour_manual(values = alternating_colors) +
        geom_point(alpha = 0.9, size = 0.8, shape = 16, show.legend = FALSE) +
        labs(title = Title,xlab = "Position", y = "Normalized Read Depth") +
        scale_x_continuous(name = "Chromosome",	expand = c(0, 0), breaks = ticks, labels = unique(l$CHR) ) +
        scale_y_continuous(name = "Normalized Read Depth",expand = c(0, 0), limits = c(0, 3)) +
        theme_classic()
                                        #+ guides(fill = guide_legend(keywidth = 3,keyheight = 1))
}

plot_chrs <- function(chrom, data) {
    Title = sprintf("Chr%s depth of coverage", chrom)
    l = data %>% filter(CHR == chrom)
    l$bp <- l$Start
    p <- ggplot(l, aes(x = bp, y = NormDepth, color = factor(Strain))) +
	geom_point(alpha = 0.7,	size = 0.8, shape = 16, show.legend = FALSE) + # scale_color_brewer(palette='RdYlBu',type='seq') +
	labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
	scale_x_continuous(expand = c(0,	0), name = "Position") +
	scale_y_continuous(name = "Normalized Read Depth", expand = c(0, 0), limits = c(0, 5)) +
	theme_classic()
                                        #guides(fill = guide_legend(keywidth = 3,keyheight = 1))
}

ChromDepths = tibble(strain=c(),ChromDepth=c())
for (window in windows) {
    inpattern = sprintf(".%sbp.regions.bed.gz$", window)
    file_list <- list.files(path = mosdepthdir, pattern = inpattern)
    bedwindows <- data.frame()
    for (i in 1:length(file_list)) {
        infile = sprintf("%s/%s", mosdepthdir, file_list[i])
        strain = str_replace(file_list[i], inpattern, "")
        t = read_tsv(infile, col_names = c("Chr", "Start", "End", "Depth"), col_types ="ciid")
        t$Strain = c(strain)
        medianDepth = median(t$Depth)
        t$NormDepth = t$Depth / medianDepth
        t = t %>% inner_join(.,chrlist)
        ChromDepths <- bind_rows(ChromDepths,t %>% group_by(Strain,CHR) %>% summarize(meddpth=median(Depth),.groups="keep"))
        bedwindows <- bind_rows(bedwindows, t)
    }
    length(bedwindows$Strain)
    unique(bedwindows$Strain)
    length(unique(bedwindows$Strain))
    unique(bedwindows$CHR)
    d = bedwindows %>% filter(CHR %in% chrlist$CHR) %>% arrange(CHR,Start,Strain)

    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$Start, d$CHR, length))

    d$pos = NA

    nchr = length(unique(chrlist$CHR))
    lastbase = 0
    ticks = NULL
    minor = vector(, 8)

    for (i in 1:nchr) {
        if (i == 1) {
            d[d$index == i, ]$pos = d[d$index == i, ]$Start
        } else {
            ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
            lastbase = lastbase + max(d[d$index == (i - 1), "Start"])
            minor[i] = lastbase
            d[d$index == i, "Start"] = d[d$index == i, "Start"] - min(d[d$index == i, "Start"]) + 1
            d[d$index == i, "End"] = lastbase
            d[d$index == i, "pos"] = d[d$index == i, "Start"] + lastbase
        }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
                                        # ticks
    minorB <- tapply(d$End, d$index, max, probs = 0.5)
                                        # minorB minor
                                        #xmax = ceiling(max(d$pos) * 1.03)
                                        #xmin = floor(max(d$pos) * -0.03)

                                        #pdf(pdffile, width = 16, height = 4)
    Title = "Depth of sequence coverage"

                                        # scale_color_brewer(palette='RdYlBu',type='seq') +
    p <- ggplot(d, aes(x = pos, y = Depth, color = Strain)) + geom_vline(mapping = NULL,
                                                                         xintercept = minorB, alpha = 0.5, size = 0.1, colour = "grey15") +
        geom_point(alpha = 0.8, size = 0.4, shape = 16, show.legend = FALSE) +
  	labs(title = Title, xlab = "Position", y = "Normalized Read Depth") +
        scale_x_continuous(name = "Chromosome",  expand = c(0, 0), breaks = ticks, labels = unique(d$CHR)) +
        scale_y_continuous(name = "Normalized Read Depth",expand = c(0, 0), limits = c(0, 5)) +
        theme_classic()
#guides(fill = guide_legend(keywidth = 3,keyheight = 1))

    pdffile = sprintf("plots/genomewide_coverage_%dkb.pdf", window/1000)
    ggsave(p,file=pdffile,width=16,height=4)

    plts <- lapply(unique(d$Strain), plot_strain, data = d)

    strains = unique(d$Strain)
    for(i in 1:length(strains ) ) {
      pdffile=sprintf("plots/StrainPlot_%dkb.%s.pdf", window/1000,strains[[i]])
  	  ggsave(plot = plts[[i]], file = pdffile,width=16)
    }

    plts <- lapply(1:nchr, plot_chrs, data = d)
    for(i in 1:nchr ) {
        pdffile=sprintf("plots/ChrPlot_%dkb.Chr%s.pdf", window/1000,i)
        ggsave(plot = plts[[i]], file = pdffile)
    }
}
cd<-ChromDepths %>% pivot_wider(names_from=CHR, values_from=meddpth,names_prefix="Chr")
write_tsv(cd,"coverage/Chromosome_Depths.tsv")
