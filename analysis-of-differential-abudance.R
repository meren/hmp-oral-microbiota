#!/usr/bin/env Rscript
#

suppressPackageStartupMessages(library(compute.es))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

setwd("./v3v5")

data <- as.data.frame(read.csv("oligotyping-results.txt", header=TRUE, sep="\t", fileEncoding="UTF-8"))
percents <- as.data.frame(read.csv("percent-matrix-no-patient-id.txt", header=TRUE, sep="\t", fileEncoding="UTF-8"))

sites = c('BM', 'HP', 'KG', 'PT', 'SUBP', 'SUPP', 'ST', 'SV', 'TD', 'TH')
maps <- as.data.frame(expand.grid(rep(list(0:1), 10))) 

# blank theme for ggplot
no_margins <- theme(
        axis.line =         element_blank(),
        axis.text.x =       element_blank(),
        axis.text.y =       element_blank(),
        axis.ticks =        element_blank(),
        axis.title.x =      element_blank(),
        axis.title.y =      element_blank(),
        panel.background =  element_blank(),
        panel.border =      element_blank(),
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        plot.background =   element_blank(),
        plot.title =        element_blank(),
        plot.margin =       unit(c(0, 0, 0, 0), "lines")
)   


final_results <- data.frame(OLIGO=NA,
                            ABUNDANT_IN=NA,
                            W_P=NA,
                            W=NA,
                            K_P=NA,
                            K=NA,
                            MES_d=NA,
                            MES_var_d=NA,
                            MES_g=NA,
                            MES_var_g=NA)
final_results <- final_results[-1, ]


for(j in seq(1, nrow(data))){
    oligo <- as.character(data$OLIGO[j])
    oligo_df <- as.data.frame(data[data$OLIGO == oligo, ])
    
    cat('oligo: ', oligo, ' (', j, ')\n')
    
    l = c()
    for(site in sites){
        l <- append(l, as.numeric(as.character(oligo_df[[site]])))
    }
    

    oligo_dist_among_sites <- percents[colnames(percents) %in% c('samples', gsub('-', '.', oligo))]
    l = c()
    for(site in sites){
        l <- append(l, c(oligo_dist_among_sites[oligo_dist_among_sites$samples == site, ][2]))
    }    


    results <- data.frame(MAPPING=NA,
                          T_P=NA,
                          T=NA,
                          W_P=NA,
                          W=NA,
                          K_P=NA,
                          K=NA,
                          MES_d=NA,
                          MES_var_d=NA,
                          MES_g=NA,
                          MES_var_g=NA)
    
    names(results) <- c('MAPPING', 'T_P', 'T', 'W_P', 'W', 'K_P', 'K', 'MES_d', 'MES_var_d', 'MES_g', 'MES_var_g')
    results <- results[-1, ]
    
    for(i in seq(2, nrow(maps) - 1)){
        mapping <- as.numeric(as.vector(as.data.frame(maps[i, ])))
    
        d.test <- as.numeric(unlist(l[which(mapping %in% 0)]))
        d.ctrl <- as.numeric(unlist(l[which(mapping %in% 1)]))
        
        
        d.test[1] <- d.test[1] + 1e-12
        d.ctrl[1] <- d.ctrl[1] + 1e-12
        
        sd_t <- 0; if(length(d.test) > 1) sd_t <- sd(d.test)
        sd_c <- 0; if(length(d.ctrl) > 1) sd_c <- sd(d.ctrl)
    
        d.mes_results <- mes(mean(d.test), mean(d.ctrl), sd_t, sd_c, length(d.test), length(d.ctrl))
        mes_d <- as.numeric(d.mes_results$MeanDifference[1])
        mes_d_var <- as.numeric(d.mes_results$MeanDifference[2])
        mes_g <- as.numeric(d.mes_results$MeanDifference[3])
        mes_g_var <- as.numeric(d.mes_results$MeanDifference[4])

        w <- wilcox.test(d.test, d.ctrl)
        wilcox_p <- w$p.value
        wilcox_s <- w$statistic
        
        t <- t.test(d.test, d.ctrl, var.equal=TRUE)
        t_p <- t$p.value
        t_s <- t$statistic
        
        k <- kruskal.test(as.vector(l), g = as.vector(mapping))
        kruskal_p <- k$p.value
        kruskal_s <- as.numeric(k$statistic)
        
        test_result <- data.frame(i, t_p, t_s, wilcox_p, as.numeric(wilcox_s), kruskal_p, kruskal_s, mes_d, mes_d_var, mes_g, mes_g_var)
        names(test_result) <- c('MAPPING', 'T_P', 'T', 'W_P', 'W', 'K_P', 'K', 'MES_d', 'MES_var_d', 'MES_g', 'MES_var_g')
        results <- rbind(results, test_result)
    }
    
    ordered_results <- results[with(results, order(-T, -T_P)), ]
    if(nrow(ordered_results) == 0)
        next
    
    best_match <- ordered_results[1, ]
    best_map <- as.numeric(as.vector(as.data.frame(maps[best_match$MAPPING, ])))
    
    best_map_df <- data.frame(SITE=NA,
            GROUP=NA,
            ABUNDANCE=NA)
    best_map_df <- best_map_df[-1, ]
    
    for(site in sites[which(best_map %in% 0)])
        best_map_df <- rbind(best_map_df, data.frame(SITE=site, GROUP='G1', ABUNDANCE=as.numeric(as.vector(oligo_df[[site]]))))
    for(site in sites[which(best_map %in% 1)])
        best_map_df <- rbind(best_map_df, data.frame(SITE=site, GROUP='G2', ABUNDANCE=as.numeric(as.vector(oligo_df[[site]]))))


    final_results <- rbind(final_results,
                           data.frame(OLIGO=oligo,
                                      ABUNDANT_IN=paste(as.vector(sites[which(best_map %in% 0)]), collapse=","),
                                      T_P=best_match$T_P,
                                      T=best_match$T,
                                      W_P=best_match$W_P,
                                      W=best_match$W,
                                      K_P=best_match$K_P,
                                      K=best_match$K,
                                      MES_d=best_match$MES_d,
                                      MES_var_d=best_match$MES_var_d,
                                      MES_g=best_match$MES_g,
                                      MES_var_g=best_match$MES_var_g))
                                                
    # VISUALIZATION
    g <- ggplot(data=best_map_df, aes(SITE, ABUNDANCE, fill=GROUP)) 
    g <- g + geom_bar(stat='identity', position='stack', alpha=0.9, linetype=0, colour='black')
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) 
    g <- g + labs(x='', y='Percent Abundances', title = oligo) 


    label <- ''
    label <- paste(label, "Student's t; p: ", signif(best_match$T_P, digits = 2), ', t: ', signif(best_match$T, digits = 3), sep = "")    
    label <- paste(label, "\nWilcoxon test; p: ", signif(best_match$W_P, digits = 2), ', W: ', signif(best_match$W, digits = 3), sep = "")    
    label <- paste(label, "\nKruskal test; p: ", signif(best_match$K_P, digits = 2), ', Chi-S: ', signif(best_match$K, digits = 3), sep = "")    
    label <- paste(label, "\nMeans To Effect Size:", sep="")
    label <- paste(label, "\n    d: ", signif(best_match$MES_d, digits = 2), ", var(d): ", signif(best_match$MES_var_d, digits = 2), sep="")
    label <- paste(label, "\n    g: ", signif(best_match$MES_g, digits = 2), ", var(g): ", signif(best_match$MES_var_g, digits = 2), sep="")
    label <- paste(label, "\n\nsites: ", paste(sites, collapse = ' '), sep="")
    label <- paste(label, "\nmap: ", paste(best_map, collapse = ' '), sep="")


    p <- ggplot(data=data.frame())
    p <- p + annotate("text", x=2, y=40, label=label, size=6, hjust = 0)
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
    p <- p + no_margins
    p <- p + xlim(0, 100)
    p <- p + ylim(0, 50)
    
    pdf(paste("V_", oligo, '.pdf', sep=''), width=16, height=16)
    sidebysideplot <- grid.arrange(g, p, nrow=2)
    dev.off()
}


write.table(final_results, file = "analysis-of-differential-abundance.txt", sep = "\t", row.names = FALSE, quote = FALSE)
