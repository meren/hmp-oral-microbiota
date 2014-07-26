#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridBase))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(Biostrings))



PLOT_OLIGO_PER_SITE_STACKBAR <- function(df, oligo, sites = c('BM', 'HP', 'KG', 'PT', 'ST', 'SUBP', 'SUPP', 'SV', 'TD', 'TH'), title = 'no title', labels = NA, colors = NA){
    df_sub_genus <- df[df$oligo %in% oligo, ]
    df_sub_genus <- df_sub_genus[df_sub_genus$site %in% sites, ]
    
    dfx <- data.frame(
            site=character(), 
            oligo=character(), 
            abundance=numeric(), 
            stringsAsFactors=FALSE)
    names(dfx) <- c('site', 'oligo', 'abundance')
    
    N <- 0
    for (site in sites){
        site_reduced <- df_sub_genus[df_sub_genus$site == site, ]
        for(oligo_in_site in oligo){
            N <- N + 1
            oligo_reduced <- site_reduced[site_reduced$oligo == oligo_in_site, ]
            oligo_in_site_abundance <- sum(oligo_reduced$abundance)
            dfx[N, ] <- c(site=site, oligo = oligo_in_site, abundance = oligo_in_site_abundance)
        }
    }
    
    
    # replace oligos with labels in legend
    if(!is.na(labels)){
        dfx$oligo <- factor(dfx$oligo, levels = oligo, labels = labels)
    }
    
    dfx$abundance <- as.numeric(as.character(dfx$abundance))
    dfx$oligo <- factor(dfx$oligo)
    
    p <- ggplot(dfx, aes(x=factor(site), y=abundance, group= oligo, fill = oligo))
    p <- p + geom_bar(position="fill", stat = "identity", width=0.90, colour = 'black')
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks.y = element_blank(), legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + labs(title=title)
    # color manually:
    if(is.na(colors)){
        colors <- sample(rainbow(length(levels(dfx$oligo)), s = 0.7, v = 0.9))
    }
    p <- p + scale_fill_manual(limits = levels(dfx$oligo), values = colors)
    print(p)
    return(colors)
}


PLOT_OLIGO_EVERYWHERE <- function(df, oligo, cluster_site='TH', sites = c('BM', 'HP', 'KG', 'PT', 'ST', 'SUBP', 'SUPP', 'SV', 'TD', 'TH'), metric='horn', method='average', title = 'no title', labels = NA, colors = NA){
    # get genus
    df_sub_genus <- df[df$oligo %in% oligo, ]
    df_sub_genus <- df_sub_genus[df_sub_genus$site %in% sites, ]
    
    #†get clustering based on SITE
    df_sub_genus_sub_site <- df_sub_genus[df_sub_genus$site == cluster_site, ]
    sub_matrix <- acast(df_sub_genus_sub_site, sample~oligo, value.var="abundance", fill = 0)
    distance <- metric #"manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "morisita", "horn", "mountford", "raup" , "binomial" or "chao"
    d <- vegdist(sub_matrix, method=distance) 
    fit <- hclust(d, method=method) # "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
    samples_order <- fit$labels[fit$order]
    
    # set the order of samples based on clustering results:
    df_sub_genus$sample <- factor(df_sub_genus$sample, levels=samples_order)
    df_sub_genus <- df_sub_genus[!df_sub_genus$sample %in% c(NA), ]
    
    # this is very cool:
    df_sub_genus$oligo <- reorder(df_sub_genus$oligo, df_sub_genus$abundance, FUN=sum)
    df_sub_genus$oligo <- factor(df_sub_genus$oligo, levels=levels(df_sub_genus$oligo))
    df_sub_genus <- df_sub_genus[order(df_sub_genus$oligo), ]
  
    # replace oligos with labels in legend
    if(!is.na(labels)){
        df_sub_genus$oligo <- factor(df_sub_genus$oligo, levels = oligo, labels = labels)
    }
    
    df_sub_genus <- df_sub_genus[with(df_sub_genus, order(site, oligo)), ]
    df_sub_genus$site <- factor(df_sub_genus$site)
    
    p <- ggplot(df_sub_genus, aes(x=factor(sample), y=abundance, labels, color, fill = oligo))
    p <- p + geom_bar(position="fill", stat = "identity", width=0.90)
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks.y = element_blank(), legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + labs(title=title)
    p <- p + labs(x='', y='')
    p <- p + scale_y_continuous(breaks = NULL)
    p <- p + facet_grid(site ~ .)
    # color manually:
    if(is.na(colors)){
        colors <- sample(rainbow(length(levels(factor(df_sub_genus$oligo))), s = 0.5, v = 0.7))
    }
    p <- p + scale_fill_manual(limits = levels(factor(df_sub_genus$oligo)), values = colors)
    print(p)
    
    return(samples_order)
}


PLOT_OLIGO_PERCENT <- function(df, oligos, samples_order, site = 'KG', title = 'no title'){
    # this function plots the percentage of all reads combined in 'oligos' list to all reads
    # observed in the 'site'
    
    dfx <- data.frame(sample=character(),
            site=character(), 
            abundance_of_oligos=numeric(), 
            abundance_of_all=numeric(), 
            p_abundance_of_oligos=numeric(), 
            stringsAsFactors=FALSE)
    names(dfx) <- c('sample', 'site', 'abundance_of_oligos', 'abundance_of_all', 'p_abundance_of_oligos')
    
    N <- 0
    dfx_sub_site <- df[df$site == site, ]
    for(sample in samples_order){
        dfx_sample_reduced <- dfx_sub_site[dfx_sub_site$sample == sample, ]
        dfx_oligo_reduced <- dfx_sample_reduced[dfx_sample_reduced$oligo %in% oligos, ]
        
        abundance_of_all <- sum(dfx_sample_reduced$abundance)
        abundance_of_oligos <- sum(dfx_oligo_reduced$abundance)
        p_abundance_of_oligos <- abundance_of_oligos * 100.0 / abundance_of_all
        
        N <- N + 1
        dfx[N, ] <- c(sample = sample, site=site, abundance_of_oligos = abundance_of_oligos, abundance_of_all = abundance_of_all, p_abundance_of_oligos = p_abundance_of_oligos)
    }
    
    dfx$p_abundance_of_oligos <- as.numeric(as.character(dfx$p_abundance_of_oligos))
    
    dfx$sample <- ordered(dfx$sample,levels=samples_order)
    
    #†prepare ggplot object
    q <- ggplot(dfx, aes(x=factor(sample), y=p_abundance_of_oligos))
    
    q <- q + stat_summary(fun.y = sum, geom = "bar")
    q <- q + theme(axis.text.x = element_text(angle = 90, size=10, vjust=0, face = 'bold'), legend.position = 'none', axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    q <- q + labs(title=title)
    q <- q + labs(x='', y=paste("Percentage in ", site, sep=''))
    q <- q + scale_y_sqrt()
    print(q)
    
}


DIST_BETWEEN_TWO_SEQS <- function(s1, s2){
    s1 <- gsub("-", "", s1)
    s2 <- gsub("-", "", s2)
    
    alignment <- pairwiseAlignment(s1, s2, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
    s1a <- as.character(pattern(alignment))
    s2a <- as.character(subject(alignment))
    mismatch_map <- strsplit(c(s1a, s2a), split= '')
    num_mismatches <- length(which(mismatch_map[[1]] != mismatch_map[[2]]))
    return(100 - (num_mismatches * 100 / nchar(s1a)))
}


OLIGO_DIST <- function(df, oligos, labels=c(), otu_limits=TRUE){
    df_subset <- df[df$OLIGO %in% oligos, ]

    dfx <- data.frame(OLIGO=character(),
            REP_SEQ=character(), 
            stringsAsFactors=FALSE)
    names(dfx) <- c('OLIGO', 'REP_SEQ')
    
    N <- 0
    for(i in seq(1, length(oligos))){
        oligo <- oligos[i]
        if(length(labels) > 0){
            label <- labels[i]
        } else {
            label <- oligos[i]
        }

        N <- N + 1
        dfx[N, ] <- c(OLIGO = label, REP_SEQ = as.character(df_subset[df_subset$OLIGO == oligo, ]$REP_SEQ))
    }
    
    df_subset <- dfx

    if(length(labels) > 0){
        df_subset$OLIGO <- ordered(df_subset$OLIGO,levels=labels)
    } else {
        df_subset$OLIGO <- ordered(df_subset$OLIGO,levels=oligos)
    }
    
    
    df_subset$OLIGO <- factor(df_subset$OLIGO)
    
    num_oligos <- nrow(df_subset)
    dist_mat <- matrix(nrow=num_oligos, ncol=num_oligos)
    colnames(dist_mat) <- df_subset$OLIGO
    rownames(dist_mat) <- df_subset$OLIGO
    
    for(i in 1:num_oligos) {
        for(j in 1:num_oligos){
            o1 <- df_subset$OLIGO[i]
            o2 <- df_subset$OLIGO[j]
            dist_mat[i,j] <- DIST_BETWEEN_TWO_SEQS(df_subset[o1, ]$REP_SEQ, df_subset[o2, ]$REP_SEQ)
        }
    }
    
    updated_df <- melt(dist_mat)
    names(updated_df) <- c('OLIGO1', 'OLIGO2', 'DIST')

    # find the best order
    d <- vegdist(dist_mat, method="canberra") 
    fit <- hclust(d, method="complete") #, "single", "complete", "average", "mcquitty", "median" or "centroid"
    oligos_order <- fit$labels[fit$order]
    
    # set the order of samples based on clustering results:
    updated_df$OLIGO1 <- factor(updated_df$OLIGO1, levels=rev(oligos_order))
    updated_df$OLIGO2 <- factor(updated_df$OLIGO2, levels=oligos_order)
    
    p <- ggplot(updated_df, aes(OLIGO1, OLIGO2, fill = DIST))
    p <- p + geom_tile()
    if(otu_limits){
        p <- p + scale_fill_gradient2(low = "steelblue", high = "red", mid="white", midpoint=98.5, limits=c(97,100), na.value='#bbfdc6')
    } else {
        p <- p + scale_fill_gradient2(low = "steelblue", high = "red", mid="white")
    }
    p <- p + theme(axis.text.x = element_blank(), legend.position = 'right', axis.ticks.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(p)
    return(levels(updated_df$OLIGO1))
}


STANDARD_FIG <- function(df, df_seqs, oligos, sites, title, labels, colors, cluster_site = 'SUBP', metric='horn', method = 'ward'){
    # wrapper for all the functions defined above.
    labels_ordered <- OLIGO_DIST(df_seqs, oligos, labels)
    ordering <- match(labels_ordered, labels)
    labels <- labels[ordering]
    oligos <- oligos[ordering]
    colors <- colors[ordering]

    PLOT_OLIGO_PER_SITE_STACKBAR(df, oligos, sites = sites, title = title, labels= labels, colors = colors)
    samples_order <- PLOT_OLIGO_EVERYWHERE(df, oligos, cluster_site = cluster_site, sites = sites, method = method, title = title, labels = labels, colors = colors, metric = metric)
    for(s in oral_sites){
        PLOT_OLIGO_PERCENT(df, oligos, samples_order, site = s, title = paste(title, ' in ', s, sep=''))
    }
}


########################################################################################################
#  DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA DATA 
########################################################################################################

setwd("./v3v5/")
df_v3v5 <- as.data.frame(read.csv('environment.txt', header=TRUE, sep="\t"))
df_v3v5_seqs <- as.data.frame(read.csv('representative_sequences.txt', header=TRUE, sep="\t"))

setwd("../v1v3/")
df_v1v3 <- as.data.frame(read.csv('environment.txt', header=TRUE, sep="\t"))
df_v1v3_seqs <- as.data.frame(read.csv('representative_sequences.txt', header=TRUE, sep="\t"))

all_sites = c('BM', 'HP', 'KG', 'PT', 'SUBP', 'SUPP', 'SV', 'TD', 'TH', 'ST')
oral_sites = c('BM', 'HP', 'KG', 'PT', 'SUBP', 'SUPP', 'SV', 'TD', 'TH')
base_colors <- c('#FF0000', '#FFCC33', '#FFFF00', '#00FF00', '#009900', '#00FFFF', '#3399FF', '#0000CC', '#FF00FF', '#9933FF', '#009999', '#990099', '#999900', '#990000')

# where figures are going to be generated:
setwd("..")


########################################################################################################
# SAMPLE FIGURES
########################################################################################################

# MOST ABUNDANT NEISSERIA FIGURE V1V3
man_labels_v1v3 = c('Neisseria elongata 99.6%',
                    'Neisseria flavescens',
                    'Neisseria flavescens 98.8%',
                    'Neisseria flavescens 99.6%',
                    'Neisseria oralis',
                    'Neisseria polysaccharea / meningitidis)',
                    'Neisseria sicca / mucosa / flava / oralis',
                    'Neisseria subflava')

man_oligos_v1v3 = c('Betaproteobacteria_ATGCCGAGGCCGTTCAA',
                    'Betaproteobacteria_ATGCCCTGTTGACGCGA',
                    'Betaproteobacteria_ATGCCCTGCCGACGTGA',
                    'Betaproteobacteria_ATGCCCTGTCGACGCGA',
                    'Betaproteobacteria_ATGCCGCGTCATTCGGA',
                    'Betaproteobacteria_ATGCCCTGTTAGCGCGA',
                    'Betaproteobacteria_ATGCCGCGGCCCTTCGA',
                    'Betaproteobacteria_ATGCCATAGCCCTTCGA')

man_sites_v1v3 = c('KG', 'SUBP', 'PT')
man_colors_v1v3 <- c("#FF0000", "#FFCC33", "#FFFF00", "#00FF00", "#c2c2c2", "#00FFFF", "#009900", "#0000CC")


pdf('neisseria_v1v3.pdf', width=12, height=6)
STANDARD_FIG(df=df_v1v3,                               # use the data frame for v1v3 data
             df_seq=df_v1v3_seqs,                      # this is where the representative sequences are
             oligos= man_oligos_v1v3,                  # these are the oligotypes we are interested in
             sites = oral_sites,                       # include these oral sites
             title = 'most abundant Neisseria v1v3',   # put this as a title
             labels = man_labels_v1v3,                 # use these labels instead of raw oligotype names
             colors = man_colors_v1v3,                 # colors.
             cluster_site = 'TD',                      # cluster all patients based on this site
             metric = 'horn',                          # distance metric for clustering
             method = 'ward')                          # clustering algorithm
dev.off()



# STREP OLIGOS FIGURE V3V5
strep_oligos_v3v5 = c('Firmicutes_GGCTCATTACTTGGTTTCCCTTCTACGT',
                      'Firmicutes_GGCTCATTACTTGGTTTCTCTTCTACGT',
                      'Firmicutes_GGTTCGTTATTTGATTCCCCTTCTACGT',
                      'Firmicutes_GGTTCATTGTTTGACTCCCCTTCTACGT',
                      'Firmicutes_GGTTCATTGCTTGGCTTCCCTTCTACGT',
                      'Firmicutes_GGTTCGTTGCTTGGCTTCCCTTCTACGT',
                      'Firmicutes_GGTTCGTTGTTTGACTTCCCTTCTACGT',
                      'Firmicutes_GGTTCGTTGCTTGGCTTCGCTACTACGT')

strep_labels_v3v5 = c('S.peroris (728) / S.oralis (707) / S. mitis (677, 398) / S. sp. (431, 423, 071)',
                      'S. infantis (638) / S. sp. (486, 074, 061, 058)',
                      'S. salivarius (755) / S. vestibularis (021)',
                      'S. sinensis (767) / S. parasanguinis I (721) / S. parasanguinis II (411) / S. cristatus (578) / S. australis(073) / S. sp. (069, 067)',
                      'S. sanguinis (758) / S. agalactiae (537)',
                      'S. gordonii (622) / S. sp. (056)', 'S. intermedius (644) / S. constellatus (576)',
                      'S. mutans (686)')

strep_colors_v3v5 <- base_colors[1:length(strep_oligos_v3v5)]


pdf('strep_oligos_v3v5.pdf', width=12, height=6)
STANDARD_FIG(df=df_v3v5,
        df_seq=df_v3v5_seqs,
        oligos= strep_oligos_v3v5,
        sites = oral_sites,
        title = 'STREP SITE v3v5',
        labels = strep_labels_v3v5,
        colors = strep_colors_v3v5,
        cluster_site = 'SUBP',
        metric = 'horn',
        method = 'ward')
dev.off()


# STREP OLIGOS FIGURE V1V3
strep_oligos_v1v3 = c('Firmicutes_-CTCAGTGTTGCGAGTGGTAGTTCACACTTCG-',
                      'Firmicutes_-CTCAATGTTGCGAGTGGTAGTTCACACTTCG-',
                      'Firmicutes_-CTCTGTCCTCCGAGTGGTAGTTCACACGTCG-',
                      'Firmicutes_-CTCAGTGTTGCGAGTGGTAGTTCACACATCG-',
                      'Firmicutes_-CTCAGTCCTGCGAGTGGTAGTTCACACATCG-',
                      'Firmicutes_-CTCAGTGCTGCGGGTGGTAGTTCACACTTCG-',
                      'Firmicutes_-CTCTGTATTGCGGGTGGTAGTTCACACTTCG-',
                      'Firmicutes_-CTCTGTGCTGCGAGTGGTAGTTCATACCACG-',
                      'Firmicutes_-GTTAGTATTCCGTGTGGTAGTTCACACGTCG-',
                      'Firmicutes_-CTCAGTCCTGCGAGTGGTAGTTCACGCTTCG-',
                      'Firmicutes_-CTCAATCCTGCGAGTGGTAGTTCACACATCG-',
                      'Firmicutes_-CTCTGTCCTGCGAGTGGTAGTTCACACTTCG-')


strep_labels_v1v3 = c('S. mitis I / mitis II / australis / pneumoniae',
                      'S. oralis / infantis / mitis bv II',
                      'S. salivarius / S. vestibularis',
                      'S. parasanguinis I',
                      'S. parasanguinis II',
                      'S. sanguinis',
                      'S. gordonii',
                      'S. intermedius',
                      'S. mutans',
                      'S. peroris',
                      'S. sp. (HOT-065)',
                      'S. sp. (HOT-056)')

strep_colors_v1v3 <- base_colors[1:length(strep_oligos_v1v3)]

pdf('strep_oligos_v1v3.pdf', width=12, height=6)
STANDARD_FIG(df=df_v1v3,
        df_seq=df_v1v3_seqs,
        oligos= strep_oligos_v1v3,
        sites = oral_sites,
        title = 'STREP SITE v1v3',
        labels = strep_labels_v1v3,
        colors = strep_colors_v1v3,
        cluster_site = 'SUBP',
        metric = 'horn',
        method = 'ward')
dev.off()
