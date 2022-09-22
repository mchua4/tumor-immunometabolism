pacman::p_load(edgeR, limma, RColorBrewwer, pheatmap, ggplot2, ggpubr, qvalue, 
               plyr, GSEABase, wesanderson, scales, grid, WebGestaltR, stringr)

################
# MDS function #
################
mds <- function(normDGE, metacol, title) {
  mcol <- as.factor(metacol)
  col <- rainbow(length(levels(mcol)), 1, 0.8, alpha=0.5)[mcol]
  plotMDS(normDGE, col=col, pch=16, cex=2)
  legend("center", fill=rainbow(length(levels(mcol)), 1, 0.8), legend=levels(mcol), horiz=F, btw="o", box.col="grey", xpd=TRUE)
  tile(main=title)
}

###################
# Camera function #
###################
Cameraplusplots <- function(contrast, genesetlist, vobject, design, catcolors, title) {
  tmp.list <- list()
  cam <- data.frame(matrix(ncol=5, nrow=0))
  for (i in 1:length(genesetlist)){
    cam.s <- camera(vobject, genesetlist[[i]], design, contrast=contrast, inter.gene.cor=0.01)
    tmp.list[[i]] <- cam.s
    names(tmp.list)[i] <- names(genesetlist)[i]
    tmp.list[[i]]$category <- names(tmp.list[i])
    colnames(cam) <- names(tmp.list[[i]])
    cam <- rbind.data.frame(cam, tmp.list[[i]])
    print(paste0("Gene set categories:", i))
  }
}
cam$neglogFDR <- log10(cam$FDR)
cam$dirNeglogFDR <- cam$neglogFDR
cam[(cam$Direction=="Down"), "dirNeglogFDR"] <- -cam[(cam$Direction=="Down"), "neglogFDR"]
grob <- grobTree(testGrob(c("UP", "DOWN"), x=c(0.94, 0.89). t=c(0.95, 0.05), gp=gpar(fontsize=13), hjust=0))
q <- ggplot(aes(x=cam$category, y=dirNeglogFDR, color=category), data=cam) + 
  geom_jitter(aes(size=NGenes, alpha=neglogFDR), pch=19, show.legend=F) + 
  scale_color_manual(values=catcolors) + scale_alpha_continuous(range=c(0.4, 1)) + 
  scale_size_continuous(range=c(4,16)) + geom_hline(yintercept=c(-1.3,1.3), color="red", alpha=0.5) + 
  labs(x="Gene set categories", y="-log10(FDR)", title=title) + geom_hline(yintercept=0) + 
  sclae_y_continous(limits=c(-10,10), oob=squish, labels=abs) theme_bw(14)+ 
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), axis.ticks.x=element_blac(),
        panel.grid.minor=element_blank(), panel.grid.major=element_blank()) + annotation_custom(grob)
print(q)
cam$geneSet <- row.names(cam)
cam10 <- as.data.frame(cam)
cam10 <- cam10[order(cam10$FDR),]
cam10 <- cam10[1:10,]
grob <- grobTree(textGrob(c("DOWN", "UP"), x=c(0.03, 0.09), y=0.025, hjust=0, col="grey60", gp=gpar(fontsize=9)))
g <- ggplot(aes(x=geneSet, y=dirNeglogFDR, fill=category), data=cam10) +
  geom_col() + geom_hline(yintercept=c(-1.3, 1.3), color="red", alpha=0.3) +
  aes(reorder(stringr::str_wrap(geneSet, 60), FDR), dirNeglogFDR) + xlab(NULL) + 
  geom_hline(y_intercept=0) + scale_y_continuous(limits=c(-10,10), oob=squish, labels=abs) +
  labs(y="-log10(FDR)", title=title) + scale_fill_manual(values=catcolors) + coord_flip() + 
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y=element_blank()) +
  annotation_custom(grab)
print(g)
return(cam)

#########################################
# Over-representation bar plot function #
#########################################
oraplot <- function(orares, catcolors, name){
  orares.n <- orares(order(orares$FDR),)
  orares.n <- orares.n[1:10,]
  orares.n$neglogFDR <- log10(orares.n$FDR)
  orares.n$geneSet <- gsub("_", "", title=paste0(name)) + 
    scale_fill_manual(values=catcolors) + coord_flip() + theme_bw(11)
  return(g)}

###########################
# compare log FC function #
###########################
power.compare.logFC <- function(sig1, sig2, N, N_other=c(2,4,6,8,10), alpha=0.05, n_test_20000){
  d <- seq(0,3,length.out=1000)
  alpha_multiple <- alpha/n_test
  df<- lapply(N_other/N, function(n_scale){
    sigSq <- (sig1^2 + sig2^2)/n_scale
    cutoff <- qnorm(alpha_multiple/2,0,sd = sqrt(sigSq), lower.tail=FALSE)
    p1 <- pnorm(-1*cutoff, d, sqrt(sigSq))
    p2 <- 1-pnorm(cutoff, d, sqrt(sigSq))
    data.frame(n_scale, power=p1+p2,d)
  })
  df <- do.call("rbind", df)
  ggplot(df, aes(d, power, color=as.factor(n_scale*N))) + 
    theme_bw(14) + geom_line()+
    scale_color_discrete("Samples") + 
    theme(aspect.ratio=1, plot.title=element_text(hjust=0.5)) + ylim(0,1) + 
    xlab(bquote(abs(logFC[observed] - logFC[expected]))) + 
    ggtitle("Power versus difference in logFC")
}
