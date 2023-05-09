basic_filter <- function(dartD, name, maxmisi = 50, mincalL = 0.80, mincalI = 0.50, minrep = 0.99, minmaf = 0.03, secd = TRUE, HWEF = TRUE, depthr = c(5,75), npopsHE = 6) {
    # dartD    - DArT genotype file
    # name  - prefix for output files (string)
    # maxmisi  - maximum proportion of missing data per individual
    # mincalL  - minimum call rate per locus, inverse of maximum missing data por locus
    # minrep  - minimum reproducibility
    # minmaf  - minimum minor allele frequency

    # secd  - remove secundaries
    # depthr  - min and maximum depth of coverage
    # HW  - filter for HWE
    # npopsHE  - Maximum number of population where a locus can be out of HWE
        my_results <- list()
        ind.missing <- function(gi, maxmis) {
                l <- list()
                x <- as.data.frame(rowSums(is.na(gi@tab))/dim(gi@tab)[2]*100)
                l$missing.data <- x
                l$hmi <- rownames(x[x > maxmis, drop=FALSE, 1])
                cat(paste0("There are ", length(l$hmi), " samples with > ", maxmis, "% missing data\n"))
                return(l)
        }
        fst.heatmap <- function(fstmat,order_list) {
                Fst_mat<- fstmat
                Fst_mat[upper.tri(Fst_mat)]<-Fst_mat[lower.tri(Fst_mat)]
                diag(Fst_mat)[is.na(diag(Fst_mat))] <- 0
                melted_cormat <- melt(Fst_mat, na.rm = TRUE)
                melted_cormat$Var1<-factor(melted_cormat$Var1, levels = order_list )
                melted_cormat$Var2<-factor(melted_cormat$Var2, levels = order_list )
                Fst_mat<- acast(melted_cormat, Var1 ~ Var2, value.var = "value")
                Fst_mat[upper.tri(Fst_mat)]<-NA
                melted_cormat <- melt(Fst_mat, na.rm = TRUE)
                return(melted_cormat)
        }
        dart.gi <- gl2gi(dartD)
        my_results$mdin <- ind.missing(dart.gi, maxmisi)
        dartD01 <- gl.drop.ind(dartD,ind.list = my_results$mdin$hmi,recalc = TRUE)
        my_results$dartD01 <- dartD01
        my_results$filstpstab <- data.frame(SNPs = nLoc(dartD01))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Dart raw SNPs")
        dartD02 <- gl.filter.monomorphs(dartD01)
        my_results$dartD02 <- dartD02
        my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD02))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Polymorphic SNPs")
        cat(paste0(nLoc(dartD02), " SNPs were polymorphic\n"))
        dartD03 <- gl.filter.reproducibility(dartD02, minrep)
        my_results$dartD03 <- dartD03
        my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD03))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Reprodicibility >= ", minrep)
        cat(paste0(nLoc(dartD02)," SNPs had higher reproducibility than ", minrep, "\n"))
        dartD04 <- gl.filter.callrate(dartD03, method = "loc", threshold = mincalL)
        my_results$dartD04 <- dartD04
        my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD04))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Call rate per locus >= ", mincalL)
        cat(paste0(nLoc(dartD04), " SNPs had higher call rate than ", mincalL, "\n"))
        dartD05 <- gl.filter.maf(dartD04,minmaf)
        my_results$dartD05 <- dartD05
        my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD05))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("MAF >= ", minmaf)
        cat(paste0(nLoc(dartD05), " SNPs had higher maf than ", minmaf, "\n"))
        dartD06<-gl.filter.rdepth(dartD05,depthr[1],depthr[2])
        my_results$dartD06 <- dartD06
        my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD06))
        rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Depth of coverage > ", depthr[1], " < ", depthr[2])
        cat(paste0(nLoc(dartD06), " SNPs had depth coverage between ", depthr[1], "and ", depthr[2],"\n"))
        if (HWEF){
                dartD07<-gl.filter.hwe(dartD06, n.pop.threshold = npopsHE)
                my_results$dartD07 <- dartD07
                my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD07))
                rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Out of HWE < ", npopsHE, " pops")
                cat(paste0(nLoc(dartD07), " SNPs are on HWE\n"))
        } else {
                dartD07<-dart06
        }
        if (secd){
                dartD08 <-gl.filter.secondaries(dartD07)
                my_results$dartD08 <- dartD08
                my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD08))
                rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("One SNP per read")
                cat(paste0(nLoc(dartD08), " SNPs are uniquic (one per read)\n"))
        } else {
                dartD08 <- dartD07
        }
        my_results$filstpstab$SNPs <-format(my_results$filstpstab, big.mark=",",scientific=FALSE)
        my_results$dartD08.gi <- gl2gi(dartD08)
        my_results$pophet <- gl.report.heterozygosity(dartD08)
        my_results$popdiv <- gl.report.diversity(dartD08)
        my_results$popfst <- gl.fst.pop(dartD08)
        my_results$popGD <-gl.dist.pop(dartD08)
        my_results$PCoA <- gl.pcoa(dartD08,4)
        my_results$dartD08.gi_Nmd<-tab(my_results$dartD08.gi, NA.method="mean")
        my_results$PCA<-dudi.pca(df = my_results$dartD08.gi_Nmd, scannf = FALSE, nf = 4)
        my_results$basdiv <- my_results$pophet[c("nInd","nLoc","polyLoc","Ho","He","FIS")]
        my_results$basdiv["nInd"]<-table(my_results$dartD08$pop)
        my_results$basdiv[c("Ho","He","FIS")]<-format(round(my_results$basdiv[c("Ho","He","FIS")],4),nsmall=4)
        my_results$basdiv[c("nLoc","polyLoc")]<-format(my_results$basdiv[c("nLoc","polyLoc")],big.mark=",",scientific=FALSE)
        my_results$fst4hm <- fst.heatmap(my_results$popfst$Fst, levels(dartD$pop))
        saveRDS(my_results, file=paste0(name, "_basicfilters.Rdata"))
        return(my_results)
}
