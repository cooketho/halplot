# halplot.r
# Tom Cooke
# 2017-04-03
# Read in a table of HAL graph edges from haltraverse and plot them.

library('ggplot2')
library('plyr')
library('reshape2')
library('RColorBrewer')
library('IRanges')

calcOffset <- function(df, padding) {
    # Calculates offset for a sequence in a genome based on the order and lengths of sequences.
    
    df2 <- df[order(df$order),]
    df2$seqLen <- df2$max - df2$min + 1
    df2$seqOffset <- c(0, cumsum(df2$seqLen + padding)[1:nrow(df2) - 1])
    return(df2[,c('genome', 'seq', 'seqLen', 'seqOffset', 'seqRev')])
}

calcSegmentOffset <- function(df) {
    # Compresses sequences for legibility by deleting unused segments.
    
    # Don't want duplicated segments during gap calculations.
    df2 <- unique(df[,1:4])
    df2 <- df2[order(df2$start),]
    # Make IRanges interval for region spanned by segments (1-based coords).
    range <- IRanges(1, end=max(df2$end))
    # Make Iranges intervals for segments.
    seg <- IRanges(start=df2$start, end=df2$end)
    # Calculate gaps.
    gaps <- as.data.frame(setdiff(range, union(seg, seg)))
    # ?
    gaps$end <- gaps$end + 1
    # Rename columns, as gaps will be matched to segments by gap end to segment start.
    gaps2 <- gaps[,c('end', 'width')]
    names(gaps2) <- c('start', 'width')
    df2$dummy <- 1:nrow(df2)
    df3 <- ddply(df2, .(start), head, n=1)
    # Join preceding gap widths to segments.
    df4 <- join(df3, gaps2, by='start', type='left', match='all')
    df5 <- join(df2, df4[,c('dummy', 'width')], by='dummy', type='left', match='all')
    df5$width[is.na(df5$width)] <- 0
    # Calculate offset for removing gaps.
    df5$segOffset <- -cumsum(df5$width)
    # Join offsets with original data frame.
    df6 <- join(df, df5[,c('genome', 'seq', 'start', 'segOffset')], by=c('genome', 'seq', 'start'), type='left', match='first')
    # Calculate gap-corrected coords.
    df6$start <- df6$start + df6$segOffset
    df6$end <- df6$end + df6$segOffset
    return(df6[,1:10])
}

# Read in the data from haltraverse.
data <- read.table('wga6.edges', header=F)
names(data) <- c('parentGenome', 'parentSeq', 'parentStart', 'parentEnd', 'childGenome', 'childSeq', 'childStart', 'childEnd', 'reversed', 'feature')
data$pseq <- data$parentSeq
data$cseq <- data$childSeq

# Convert the wide-format data to long-format.
topSegments <- data[,c(1:4, 9:12)]
names(topSegments) <- c('genome', 'seq', 'start', 'end', 'rev', 'feature', 'pseq', 'cseq')
topSegments$id <- 1:nrow(topSegments)
topSegments$bot <- rep(0, nrow(topSegments))    # '0' denotes top segment.
botSegments <- data[,5:12]
names(botSegments) <- c('genome', 'seq', 'start', 'end', 'rev', 'feature', 'pseq', 'cseq')
botSegments$id <- 1:nrow(botSegments)
botSegments$bot <- rep(1, nrow(topSegments))    # '1' denotes bottom segment.
segments <- rbind(topSegments, botSegments)

# Calculate x offset for removing gaps between segments.
segments2 <- ddply(segments, .(genome, seq), .fun=function(df) calcSegmentOffset(df))
# Summarize min and max sequence positions for each genome & sequence.
sequences <- ddply(segments2, .(genome, seq), summarize, min=min(start), max=max(end))

# Read in optimized sequence orders.
config <- read.table('wga6.config')
names(config) <- c('genome', 'seq', 'order', 'seqRev')
# Manually flip one of the sequences.
config[config$seq=='Anc2refChr1825', 'seqRev'] <- 0

# Join sequence data with sequence order and orientation data.
sequences3 <- join(sequences, config, by=c('genome', 'seq'), type='left', match='all')

# Join the segment data with the sequence orientation data.                    
sequences4 <- sequences3[,c('seq', 'seqRev')]
names(sequences4) <- c('cseq', 'cSeqRev')
sequences5 <- sequences3[,c('seq', 'seqRev')]
names(sequences5) <- c('pseq', 'pSeqRev')
segments2b <- join(segments2, sequences4, by=c('cseq'), type='left', match='all')
segments2c <- join(segments2b, sequences5, by=c('pseq'), type='left', match='all')

# Calculate x offset for plotting sequences in the given order.
offsets <- ddply(sequences3, .(genome), .fun=function(df) calcOffset(df, 0))

# Join the segment data with the sequence order data.
segments3 <- join(segments2c, offsets, by=c('genome', 'seq'), type='left', match='all')

# Calculate the segment starts and ends (still need to apply genome offsets).
segments3$startCorrected <- segments3$seqOffset + (-1)^segments3$seqRev * segments3$start + segments3$seqLen * segments3$seqRev
segments3$endCorrected <- segments3$seqOffset + (-1)^segments3$seqRev * segments3$end + segments3$seqLen * segments3$seqRev
# Make column denoting whether child and parent seq have same orientation.
segments3$seqRev <- as.integer(segments3$cSeqRev==segments3$pSeqRev)

# Manually set the x and y offsets for each genome in R plot coordinates.
yOffset <- data.frame(genome=c('human', 'opossum', 'anole', 'budgie', 'medaka', 'Anc0', 'Anc1', 'Anc2', 'Anc3'),
                      y=c(0, 0, 0, 0, 0, 2.5, 1.8, 1, 1))

xOffset <- data.frame(genome=c('human', 'opossum', 'anole', 'budgie', 'medaka', 'Anc0', 'Anc1', 'Anc2', 'Anc3'),
                      xOffset=c(0,40000,80000,160000,200000,100000,70000,25000,120000))

# Join segment data with offsets.
columns <- c('genome', 'seq', 'startCorrected', 'endCorrected', 'rev', 'seqRev', 'id', 'bot', 'pseq', 'cseq', 'feature')
data5 <- join(segments3[,columns], yOffset, by='genome', type='left', match='all')
# Convert data to long-form (start and end are now on separate rows).
data6 <- melt(data5, id=c('genome', 'seq', 'pseq', 'cseq', 'id', 'rev', 'seqRev', 'bot', 'y', 'feature'), value.name='x')
# Join with x offsets.
data7 <- join(data6, xOffset, by='genome', type='left', match='all')

# Polygon vertexes are plotted in order they appear in data frame.
# Hence, need to re-order data frame to ensure correct plotting of polygons.
data8 <- data7[order(data7$id, data7$y, data7$x * ((-1)^data7$rev * (-1)^data7$seqRev)^data7$bot),]
data8$id <- factor(data8$id)

# Make plotOrder column so that polygons get plotted in order of child genome.
botGenome <- ddply(data8, .(id), summarize, botGenome=genome[1])
data9 <- join(data8, botGenome, by='id', type='left', match='all')
data10 <- data9[order(data9$botGenome),]
plotOrder <- data.frame(id=unique(data10$id), plotOrder=1:length(unique(data10$id)))
data11 <- join(data10, plotOrder, by='id', type='left', match='all')

# Re-order factor levels for genes.
lev <- c('PDSS1', 'ABI1', 'ANKRD26', 'ACBD5', 'MASTL', 'YME1L1', 'PKS', 'PTCHD3', 'RAB18', '.')
data11$feature <- factor(data11$feature, levels=lev)

# Set colors and themes.
col <- brewer.pal(n=12, 'Paired')[c(1, 12, 10, 3:4, 2, 6:9)]
o <- theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
o <- o + theme(panel.border=element_rect(fill=NA, color='black'), panel.background=element_rect(fill=NA, color='black'))

# Make HAL graph plot.
p <- ggplot(data11, aes(x=x + xOffset, y=y, group=plotOrder, fill=feature))
p <- p + geom_polygon(alpha=0.2)
p <- p + scale_fill_manual(values=col)
p

# Save to pdf.
pdf('wga6.hal.pdf', height=4, width=6)
p + o
dev.off()

