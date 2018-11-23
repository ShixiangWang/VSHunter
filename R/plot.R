# Visualization Part ------------------------------------------------------
cnv_plotSignatures = function(Res = NULL, contributions = FALSE, color = NULL,
                              patient_order = NULL, font_size = 1.2, show_title = TRUE,
                              axis_lwd = 2, title_size = 1.3, show_barcodes = FALSE, barcode_size = 0.5,
                              yaxisLim = 0.3, ...) {
    # modify from https://github.com/PoisonAlien/maftools/blob/master/R/plotSignatures.R

    # input can be multipe objects from different functions
    # basically, we only use w and h matrix of NMF

    # assuming input is a basic NMF result
    if (inherits(Res, "NMFfitX1")) {
        #-- Signatures
        w = NMF::basis(Res)
        w = apply(w, 2, function(x) x / sum(x)) # Scale signatures (basis)
        colnames(w) = paste("Signature", 1:ncol(w), sep = "_")

        #-- Contribution (this mainly come from result of LCD, it is ok directly from NMF)
        h = NMF::coef(Res)
        h = apply(h, 2, function(x) x/sum(x)) # Scale contributions (coefs)

        rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
    } else if (is.list(Res)) {
        if (!all(c("signature", "exposure") %in% names(Res))) stop("signature and exposure elements must in input list.")
        w = Res$signature
        w = apply(w, 2, function(x) x / sum(x)) # Scale signatures (basis)

        if (is.null(colnames(w))) {
            colnames(w) = paste("Signature", 1:ncol(w), sep = "_")
        }

        h = Res$exposure$relative_exposure
        if (is.null(rownames(h))) {
            rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
        }
    }

    contrib = h

    if(contributions){
        contribt = t(contrib)
        #calculate sd
        if(!is.null(patient_order)){
            contribt = contribt[patient_order,] #order on user-specified ordering of the genomes
        }else{
            contribt = contribt[order(contribt[,ncol(contribt)]),] #order according to standard deviation
        }

        #contrib = t(contribt[,1:(ncol(contribt)-1)])
        contrib = t(contribt[,1:(ncol(contribt))])

        cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')

        if(show_barcodes){
            lo = layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
            par(mar = c(6, 4, 2, 1))
            b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
            axis(side = 1, at = b, labels = colnames(contrib), lwd = 2, cex.axis = barcode_size,
                 las = 2, line = 1, hadj = 0.8, font = 2, tick = FALSE)
            axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 2, las = 2, cex.axis = 0.9)
            mtext(text = "Signature exposures", side = 2, font = 2, cex = 1, line = 2.8)
            par(mar = rep(2, 4))
            plot.new()
            #par(mar = c(2, 3, 0, 0))
            legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
                   border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
                   cex = 1.2, pt.cex = 1.5, horiz = TRUE)
        }else{
            lo = layout(mat = matrix(data = c(1, 2), nrow = 2), heights = c(6, 2))
            par(mar = c(3, 4, 2, 1))
            b = barplot(contrib, axes = FALSE, horiz = FALSE, col = cols, border = NA, names.arg = rep("", ncol(contrib)))
            axis(side = 2, at = seq(0, 1, 0.25), lwd = 3, font = 2, las = 2, cex.axis = 0.9)
            mtext(text = "Signature exposure", side = 2, font = 2, cex = 1, line = 2.8)
            plot.new()
            par(mar = c(2, 3, 0, 0))
            legend(x = "left", legend = rownames(contrib), col = cols[1:nrow(contrib)],
                   border = NA, bty = "n", pch = 15, xpd = TRUE, ncol = 1,
                   cex = 1.2, pt.cex = 1.5, horiz = TRUE)
        }
    }else{

        plotData = as.data.frame(t(w))
        nsigs = nrow(plotData)

        if(is.null(color)){
            #color = c("blue","black","red","gray","green","maroon")
            color = c('coral4', 'lightcyan4', 'deeppink3', 'lightsalmon1', 'forestgreen', 'cornflowerblue')
        }

        cnames = colnames(plotData)
        len_c1 = sum(grepl("segsize", cnames))
        len_c2 = sum(grepl("bp10MB", cnames))
        len_c3 = sum(grepl("osCN", cnames))
        len_c4 = sum(grepl("changepoint", cnames))
        len_c5 = sum(grepl("copynumber", cnames))
        len_c6 = sum(grepl("bpchrarm", cnames))

        len_seqc = c(len_c1,
                     len_c2,
                     len_c3,
                     len_c4,
                     len_c5,
                     len_c6)

        colors = rep(color, times=len_seqc)

        len_seq = c(0L)
        for (i in 1:length(len_seqc)) {
            s = len_seqc[i] * 2 + len_seq[i]
            len_seq = c(len_seq, s)
        }

        par(mfrow = c(nsigs,1),oma = c(5,4,0,0) + 0.1, mar = c(0,0,2.5,0) + 0.1, las=1, tcl=-.25, font.main=4, xpd = NA)

        for(i in 1:nsigs){

            ae = rownames(plotData)[i]
            d = as.matrix(plotData[i,])
            if(is.na(yaxisLim)){
                bh = ceiling(max(d, na.rm = TRUE) * 10)/10 #Bar height
            }else{
                bh = 0.3
            }

            barplot(d, xaxt = "n", yaxt = "n", col = colors, beside = TRUE, ylim = c(-0.1, bh),
                    cex.main = 1, border = NA, font.axis = 2, font.lab = 2,
                    adj = 0.25, ...)
            if(show_title){
                title(main = ae, cex.main = title_size, line = 0, adj = 0.98)
            }

            #mtext(text = ae, side = 1, line = 2, font = 1, cex = 0.5, at = 0.3)
            axis(side = 2, at = seq(0, bh, 0.1),
                 pos = -2, las = 2, lwd = axis_lwd, hadj = 1.1,
                 font = 2, cex.axis = font_size)
            #abline(h = seq(0, 0.3, 0.1),lty=2,lwd=0.3, col = 'gray70')
            rect(xleft = len_seq, ybottom = -0.05,
                 xright = ncol(plotData)*2, ytop = -0.02, col = color, border = 'gray70')
            if(i == nsigs){
                text(labels = c("segsize","bp10MB","osCN","changepoint","copynumber","bpchrarm"),
                     y = rep(-0.1,6),x = len_seq[2:7]-len_seqc, cex = font_size,
                     font = 1.8, font.lab = 2, pos = 1.2, srt = 30)
            }
        }
    }


}


cnv_plotFeatureDistribution = function(features, ylab = "", ...) {
    features = lapply(features, function(x) {
        x$value = as.numeric(x$value)
        return(x)
        })

    p_1 = ggplot(data = features$segsize, aes(x = value)) +
        geom_line(stat="density") + labs(x = "Segment size", y = ylab) + theme_classic()
    p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(0, 7:9),
                                   labels = scales::trans_format("log10", scales::math_format(10^.x)))

    p_2 = ggplot(data = features$copynumber, aes(x = value)) +
        geom_line(stat="density") + labs(x = "Copy number", y = ylab) + theme_classic()
    p_3 = ggplot(data = features$changepoint, aes(x = value)) +
        geom_line(stat="density") + labs(x = "Copy number changepoint", y = ylab) + theme_classic()

    p_4 = ggplot(data = features$bp10MB, aes(x = value)) +
        geom_bar(stat = "count") + labs(x = "Breakpoint count per 10MB", y = ylab) + theme_classic()
    p_5 = ggplot(data = features$bpchrarm, aes(x = value)) +
        geom_bar(stat = "count") + labs(x = "Breakpoint count per chr arm", y = ylab) + theme_classic()
    p_6 = ggplot(data = features$osCN, aes(x = as.factor(value))) +
        geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) + theme_classic() + theme(text = element_text(color = "black"))

    p = cowplot::plot_grid(p_1, p_2, p_3, p_4, p_5, p_6,
                       nrow = 2, align = "hv", ...)
    p
}


cnv_plotMixComponents = function(features, components, ...) {


    cbPalette <- c(RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set1"),"black")
    plotNormDensity = function(value, matrix, xlab) {
        p = ggplot(data = data.frame(x = c(min(value), max(value))),
                   aes(x)) + ylab("")

        for (i in 1:ncol(matrix)) {
            p = p + stat_function(fun = stats::dnorm, n = 1000,
                                  args = list(
                                      mean = matrix[1, i],
                                      sd = matrix[2, i]),
                                  color = cbPalette[i])
        }

        p = p+xlab(xlab)+theme_classic()
        p
    }

    plotPoisDensity = function(value, lambda, xlab, max_value=10) {
        if (is.null(max_value)) {
            p = ggplot(data = data.frame(x = seq(min(value), max(value), length.out = 100)),
                       aes(x)) + ylab("")
        } else {
            p = ggplot(data = data.frame(x = seq(min(value), max_value, length.out = 100)),
                       aes(x)) + ylab("")
        }


        for (i in 1:length(lambda)) {
            p = p + stat_function(geom = "line",  n = 11, fun = stats::dpois,
                                  args = list(lambda = lambda[i]),
                                  color = cbPalette[i])
        }

        p = p+xlab(xlab)+theme_classic()
        p
    }



    features = lapply(features, function(x) {
        x$value = as.numeric(x$value)
        return(x)
    })

    # norm distribution
    comp_segsize = flexmix::parameters(components[["segsize"]])
    comp_copynumber = flexmix::parameters(components[["copynumber"]])
    comp_changepoint = flexmix::parameters(components[["changepoint"]])
    # pois distribution
    comp_bp10MB = flexmix::parameters(components[["bp10MB"]])
    comp_bpchrarm = flexmix::parameters(components[["bpchrarm"]])
    comp_osCN = flexmix::parameters(components[["osCN"]])

    # norm plots
    p_1 = plotNormDensity(features[["segsize"]]$value, comp_segsize, xlab = "Segment size")
    p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(0, 7:9),
            labels = scales::trans_format("log10", scales::math_format(10^.x)))
    p_2 = plotNormDensity(features[["copynumber"]]$value, comp_copynumber, xlab = "Copy number")
    p_3 = plotNormDensity(features[["changepoint"]]$value, comp_changepoint, xlab = "Copy-number changepoint")

    # pois plots
    p_4 = plotPoisDensity(features[["bp10MB"]]$value, comp_bp10MB, xlab = "Breakpoint count per 10MB", max_value = 10)
    p_5 = plotPoisDensity(features[["bpchrarm"]]$value, comp_bpchrarm, xlab = "Breakpoint count per arm", max_value = 50)
    p_6 = plotPoisDensity(features[["osCN"]]$value, comp_osCN, xlab = "Oscilating CN chain length")

    p = cowplot::plot_grid(p_1, p_2, p_3, p_4, p_5, p_6,
                           nrow = 2, align = "hv", ...)
    p
}


