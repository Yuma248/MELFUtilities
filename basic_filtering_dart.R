##This scrip does basic filtering of dart data and create a table with filtering steps, basis statistic, Fst heat map, PCA plot and Admixture/Structure plot.
#' @name basic_filter
#' @rdname basic_filter
#' This function does basic filtering of dart SNPs.
#'
#' @param dartD DArT genotype file.
#' @param name prefix for output files (string).
#' @param maxmisi maximum proportion of missing data per individual.
#' @param mincalL minimum call rate per locus, inverse of maximum missing data por locus.
#' @param minrep minimum reproducibility.
#' @param minmaf minimum minor allele frequency.
#' @param maxsim max reletedness to remove duplicates.
#' @param secd remove secundaries.
#' @param depthr min and maximum depth of coverage.
#' @param HW filter for HWE.
#' @param npopsHE Maximum number of population where a locus can be out of HWE.
#' @param PDFplots If you want to create a PDF with basic results (filter steps table, diversity table, Fst heat map, PCA, Q plot).
#' @param nK maximum number of pop to test with cluster algorithm.
#' @return A description of the return value.
#' @examples
#' basic_filter(dartD, "Yuma", maxmisi = 50, mincalL = 0.80, mincalI = 0.50, minrep = 0.99, minmaf = 0.01, secd = TRUE, HWEF = TRUE, depthr = c(5,100), npopsHE = 2, maxsim = 0.8, PDFplots = TRUE, nK = 3)





basic_filter <- function(dartD, name, maxmisi = 50, mincalL = 0.80, mincalI = 0.50, minrep = 0.99, minmaf = 0.01, secd = TRUE, HWEF = TRUE, depthr = c(5,100), npopsHE = 2, maxsim = 0.8, PDFplots = TRUE, nK = 3) {
    # dartD    - DArT genotype file 
    # name  - prefix for output files (string)
    # maxmisi  - maximum proportion of missing data per individual 
    # mincalL  - minimum call rate per locus, inverse of maximum missing data por locus
    # minrep  - minimum reproducibility 
    # minmaf  - minimum minor allele frequency
    # maxsim  - max reletedness to remove duplicates
    # secd  - remove secundaries 
    # depthr  - min and maximum depth of coverage
    # HW  - filter for HWE 
    # npopsHE  - Maximum number of population where a locus can be out of HWE
    # PDFplots  - If you want to create a PDF with basic results (filter steps table, diversity table, Fst heat map, PCA)

	my_results <- list()
	ind.missing <- function(gi, maxmis) {
		l <- list()
		x <- as.data.frame(rowSums(is.na(gi@tab))/dim(gi@tab)[2]*100)
		l$missingdata <- x
		names(l$missingdata)[1]<-"%miss"		
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
	ind.dup <- function(rel, maxsim, missdat) {
		dup <- list()
		dup$dup <- data.frame(matrix(nrow = 0, ncol = 3)) 
		colnames(dup$dup) = c("Ind1.id","Ind2.id","Wang")
		dup$ind2rem <- c()
		for(i in 1:nrow(rel)) {
			if (rel$wang[i] < maxsim) {
				next
			}
			dup$dup[nrow(dup$dup) + 1,] = c(rel[i,"ind1.id"],rel[i,"ind2.id"], rel$wang[i] )
			if (which.max(c(missdat[match(rel[i,"ind1.id"],rownames(missdat)),"%miss"],missdat[match(rel[i,"ind2.id"], rownames(missdat)),"%miss"])) > 1){
				dup$ind2rem<- append(dup$ind2rem,rel[i,"ind2.id"])
			} else {
				dup$ind2rem<-append(dup$ind2rem,rel[i,"ind1.id"])
			}
		}
		cat(paste0("There are ", length(dup$ind2rem), " possible duplicates with > ", maxsim, " relatedness\n"))
		return (dup)
	}
	Qvsnmf<-function(project,dartDin ){
		myselK<-which.min(summary(project)$crossEntropy[2,])
		bestR <- which.min(cross.entropy(project, K = myselK))
		myorder <- order(dartDin$pop)
		mynamesort <- dartDin$ind.names[order(dartDin$pop)]
		mypopsort<- dartDin$pop[order(dartDin$pop)]
		myQv <- read.table(paste0(getwd(),"/dartD09.snmf/K",as.character(myselK),"/run",as.character(bestR),"/dartD09_r",as.character(bestR),".",as.character(myselK),".Q"))[myorder,]
		rownames(myQv)<-mynamesort
		myQv <- cbind(myQv,mypopsort)
		myQvmelte <- melt(rownames_to_column(myQv, "Ind"), value.name="Q")
		colnames(myQvmelte)[c(2,3)]<-c("Loc","Clusters")
		myQvmelte$Loc<-factor(myQvmelte$Loc, levels = levels(dartDin$pop))
		return(myQvmelte)
	}
	basic_filter_plots<-function(fresults,pref){
		FSY<-ggtexttable(fresults$filstpstab)
		DIVY<-ggtexttable(fresults$basdiv)
		hmY<-ggplot(data = fresults$fst4hm, aes(x = Var2, y = Var1, fill = value))+
			geom_tile(color = "white")+
			scale_fill_gradient(low = "#FFD100", high = "#7D1845",
			limit = c(min(fresults$fst4hm$value),max(fresults$fst4hm$value)), space = "Lab") +
			theme_minimal()+ xlab("") + ylab("")+
			theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1))+
			labs(fill= expression("F"["ST"]))+
			scale_y_discrete(limits = rev(levels(fresults$fst4hm$Var1))) +
			scale_x_discrete(limits = rev(levels(fresults$fst4hm$Var2))) +
			coord_fixed()
		PCAY<-ggplot(fresults$PCA$li, aes(x=Axis1, y=Axis2, shape=fresults$dartD09.gi$pop, color=fresults$dartD09.gi$pop))+
			geom_point(size = 4, alpha=0.5)+
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
			axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"))+
			xlab(paste0("PC1 (",round((fresults$PCA$eig[1]/sum(fresults$PCA$eig)*100),2),"%)"))+
			ylab(paste0("PC2 (",round((fresults$PCA$eig[2]/sum(fresults$PCA$eig)*100),2),"%)"))+
			scale_color_discrete("") +
			scale_shape_manual("",values = rep(c(15, 17,19),20)) +
			labs(color = "Pops")
		QPLOT<-ggplot(fresults$Q, aes(x = Ind, y = Q, fill = Clusters)) + geom_bar(stat = "identity", width = 1) + 
			scale_y_continuous(expand = c(0, 0)) +
			facet_grid(. ~ Loc, scales = "free_x", space = "free_x")+
			theme(legend.position = "none",axis.text.x=element_blank(), 
			axis.title.x=element_blank(),axis.ticks.x=element_blank(), panel.spacing.x = unit(0, "lines"), 
			panel.border = element_rect(color = "black", fill = NA), plot.margin = margin(0, 0.2, 0.2, 0, "cm"))
		pdf(paste0(pref,"_basicplot.pdf"), width = 12, height = 16)
		grid.arrange(FSY,DIVY,hmY, PCAY, QPLOT, ncol = 2, nrow = 3, layout_matrix = rbind(c(1,2), c(3,4), c(5,5)), 
			widths = c(0.01, 0.01), heights = c(0.01, 0.01, 0.01))
		dev.off()
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
	cat(paste0(nLoc(dartD06), " SNPs had depth coverage between ", depthr[1], " and ", depthr[2],"\n"))
	if (HWEF){
		dartD07<-gl.filter.hwe(dartD06, n.pop.threshold = npopsHE)
		my_results$dartD07 <- dartD07
		my_results$filstpstab <- rbind(my_results$filstpstab, nLoc(dartD07))
		rownames(my_results$filstpstab)[nrow(my_results$filstpstab)] <- paste0("Out of HWE < ", npopsHE, " pops")
		cat(paste0(nLoc(dartD07), " SNPs are on HWE\n"))
	} else {
		dartD07<-dartD06
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
	my_results$mdinD08 <- ind.missing(my_results$dartD08.gi, maxmisi)
	my_results$dartD08.rel <- gl2related(dartD08, save = FALSE)
	my_results$relest <- coancestry(my_results$dartD08.rel, wang =1)
	my_results$duplic <- ind.dup(my_results$relest$relatedness, maxsim, my_results$mdinD08$missingdata)
	dartD09 <- gl.filter.monomorphs(gl.drop.ind(dartD08,ind.list = my_results$duplic$ind2rem,recalc = TRUE))
	my_results$dartD09 <- dartD09
	my_results$dartD09.gi <- gl2gi(dartD09)
	my_results$pophet <- gl.report.heterozygosity(dartD09)
#	my_results$popdiv <- gl.report.diversity(dartD09)
	my_results$popfst <- gl.fst.pop(dartD09)
#	my_results$popGD <-gl.dist.pop(dartD09)
#	my_results$PCoA <- gl.pcoa(dartD09,4)
	my_results$dartD09.gi_Nmd<-tab(my_results$dartD09.gi, NA.method="mean")
	my_results$PCA<-dudi.pca(df = my_results$dartD09.gi_Nmd, scannf = FALSE, nf = 4)
	my_results$basdiv <- my_results$pophet[c("nInd","nLoc","polyLoc","Ho","He","FIS")]
	my_results$basdiv["nInd"]<-table(my_results$dartD09$pop)
	my_results$basdiv[c("Ho","He","FIS")]<-format(round(my_results$basdiv[c("Ho","He","FIS")],4),nsmall=4)
	my_results$basdiv[c("nLoc","polyLoc")]<-format(my_results$basdiv[c("nLoc","polyLoc")],big.mark=",",scientific=FALSE)
	my_results$fst4hm <- fst.heatmap(my_results$popfst$Fst, levels(dartD$pop))
	gl2geno(dartD09, outfile = "dartD09", outpath = getwd())
	my_results$snmfR <- snmf("dartD09.geno",K = 1:nK, CPU=10, entropy = TRUE, project = "new", repetitions = 3)
	my_results$Q <- Qvsnmf(my_results$snmfR, my_results$dartD09)
	saveRDS(my_results, file=paste0(name, "_basicfilters.Rdata"))
	if (PDFplots){
		basic_filter_plots(my_results, name)
	}
	return(my_results)
}
