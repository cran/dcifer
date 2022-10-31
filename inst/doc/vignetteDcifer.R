## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = ">"#,
#  fig.width = 7, fig.height = 4, fig.fullwidth = TRUE 
)

## ----setup--------------------------------------------------------------------
library(dcifer)
pardef <- par(no.readonly = TRUE)

## -----------------------------------------------------------------------------
sfile <- system.file("extdata", "MozParagon.csv", package = "dcifer")
dsmp <- readDat(sfile, svar = "sampleID", lvar = "locus", avar = "allele")
str(dsmp, list.len = 2)

## ---- eval = FALSE------------------------------------------------------------
#  dsmp <- formatDat(dlong, svar = "sampleID", lvar = "locus", avar = "allele")

## -----------------------------------------------------------------------------
# optionally, extract location information
meta <- unique(read.csv(sfile)[c("sampleID", "province")])
meta <- meta[match(names(dsmp), meta$sampleID), ]  # order samples as in dsmp

## -----------------------------------------------------------------------------
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)

## -----------------------------------------------------------------------------
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

## -----------------------------------------------------------------------------
afile  <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
afreq2 <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq") 

## ---- eval = FALSE------------------------------------------------------------
#  # alternatively, if allele frequencies are provided in an R data frame (aflong):
#  afreq2 <- formatAfreq(aflong, lvar = "locus", avar = "allele", fvar = "freq")

## -----------------------------------------------------------------------------
dsmp2  <- matchAfreq(dsmp, afreq2)

## -----------------------------------------------------------------------------
provinces <- c("Maputo", "Inhambane")
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

## ---- eval = FALSE------------------------------------------------------------
#  dres <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0,
#                 alpha = 0.05, nr = 1e3)

## -----------------------------------------------------------------------------
dres[9, 5, ]  

## ---- fig.width = 6, fig.height = 6-------------------------------------------
par(mar = c(3, 3, 1, 1))
nsmp  <- length(dsmp)
atsep <- cumsum(nsite)[-length(nsite)]

# create symmetric matrix
dmat <- dres[, , "estimate"]
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))] 

# determine significant, indices for upper triangle
alpha <- 0.05                          # significance level                    
isig <- which(dres[, , "p_value"] <= alpha, arr.ind = TRUE)
col_id <- rep(c("plum4", "lightblue4"), nsite)
plotRel(dmat, isig = isig[, 2:1], draw_diag = TRUE, alpha = alpha, idlab = TRUE,
        col_id = col_id)
abline(v = atsep, h = atsep, col = "gray45", lty = 5)

## ---- fig.height = 3.6--------------------------------------------------------
par(mfrow = c(1, 2), mar = c(1, 0, 1, 0.2))
plotRel(dres, draw_diag = TRUE, alpha = alpha)
mtext("p-values", 3, 0.2)
# p-values for upper triangle
pmat <- matrix(NA, length(dsmp), length(dsmp))
pmat[upper.tri(pmat)] <- t(log(dres[, , "p_value"]))[upper.tri(pmat)]
pmat[pmat == -Inf] <- min(pmat[is.finite(pmat)])#*1.2
plotRel(pmat, rlim = NULL, draw_diag = TRUE, col = hcl.colors(101, "Red-Purple"), 
        sig = FALSE, add = TRUE, col_diag = "white", border_diag = "gray45")
abline(v = atsep, h = atsep, col = "gray45", lty = 5)

# number of non-missing loci for upper triangle
dmiss <- lapply(dsmp, function(lst) sapply(lst, function(v) all(!v)))
nmat <- matrix(NA, nsmp, nsmp) 
for (jsmp in 2:nsmp) {
  for (ismp in 1:(jsmp - 1)) {
    nmat[ismp, jsmp] <- sum(!dmiss[[ismp]] & !dmiss[[jsmp]])
  }
}
nrng <- range(nmat, na.rm = TRUE)  #  
par(mar = c(1, 0.2, 1, 0))
plotRel(dres, draw_diag = TRUE, alpha = alpha)
mtext("number of loci", 3, 0.2)
coln <- hcl.colors(diff(nrng)*2.4, "Purple-Blue", rev = TRUE)[1:(diff(nrng) + 1)]
plotRel(nmat, rlim = NA, col = coln, add = TRUE, 
        draw_diag = TRUE, col_diag = "gray45", border_diag = "white")
par(pardef)

## ---- fit.width = 7, fig.height = 5.8-----------------------------------------
layout(matrix(1:2, 1), width = c(7, 1))
par(mar = c(1, 1, 2, 1))#, mgp = c(0, 0, 0.5))
plotRel(dmat, draw_diag = TRUE, isig = rbind(isig, isig[, 2:1]))
atclin <- cumsum(nsite) - nsite/2
abline(v = atsep, h = atsep, col = "gray45", lty = 5)
mtext(provinces, side = 3, at = atclin, line = 0.2)
mtext(provinces, side = 2, at = atclin, line = 0.2)
par(mar = c(1, 0, 2, 3))
plotColorbar()
par(pardef)

## ---- fig.width = 6.2, fig.height = 6-----------------------------------------
# horizontal colorbar (needed?)
par(mar = c(1, 1, 1, 3))
border_sig = "darkviolet"
plotRel(dres, draw_diag = TRUE, border_diag = border_sig, alpha = alpha, 
        border_sig = border_sig, lwd_sig = 2)
legend(32, 20, pch = 0, col = border_sig, pt.lwd = 2, pt.cex = 1.4, 
       box.col = "gray", legend = expression("significant at" ~~ alpha == 0.05))
text(1:nsmp + 0.3, 1:nsmp - 0.5, labels = names(dsmp), col = col_id, adj = 0,
     cex = 0.6, xpd = TRUE)
par(fig = c(0.25, 1, 0.81, 0.92), mar = c(1, 1, 1, 1), new = TRUE)
plotColorbar(at = c(0.2, 0.4, 0.6, 0.8), horiz = TRUE)
ncol <- 301
lines(c(0, ncol, ncol, 0, 0), c(0, 0, 1, 1, 0), col = "gray")
par(pardef)

## ---- eval = FALSE------------------------------------------------------------
#  # First, create a grid of r values to evaluate over
#  revals <- mapply(generateReval, 1:5, nr = c(1e3, 1e2, 32, 16, 12))

## -----------------------------------------------------------------------------
sig1 <- sig2 <- vector("list", nrow(isig))
for (i in 1:nrow(isig)) {
  sig1[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, Mmax = 5, 
                       equalr = FALSE, reval = revals)
}
for (i in 1:nrow(isig)) {
  sig2[[i]] <- ibdEstM(dsmp[isig[i, ]], coi[isig[i, ]], afreq, equalr = TRUE)
}
M1      <- sapply(sig1, function(r) sum(r > 0))  
M2      <- sapply(sig2, length)          
rtotal1 <- sapply(sig1, sum)    
rtotal2 <- sapply(sig2, sum)    

cor(M1, M2)                 
cor(rtotal1, rtotal2) 

## -----------------------------------------------------------------------------
samples <- names(dsmp)
sig <- as.data.frame(isig, row.names = FALSE)
sig[c("id1", "id2")]         <- list(samples[isig[, 1]], samples[isig[, 2]])
sig[c("M1", "M2")]           <- list(M1, M2)
sig[c("rtotal1", "rtotal2")] <- list(round(rtotal1, 3), round(rtotal2, 3))
head(sig)

## -----------------------------------------------------------------------------
i <- 18                              
pair <- dsmp[isig[i, ]]
coii <-  coi[isig[i, ]]

res1 <- ibdPair(pair, coii, afreq, M = M1[i], pval = TRUE, equalr = FALSE,  
                reval = revals[[M1[i]]])
res2 <- ibdPair(pair, coii, afreq, M = M2[i], pval = TRUE, equalr = TRUE,  
                confreg = TRUE, llik = TRUE, reval = revals[[1]])

res1$rhat                              # estimate with equalr = FALSE
rep(res2$rhat, M2[i])                  # estimate with equalr = TRUE

## ---- fig.width = 6, fig.height = 4.5-----------------------------------------
CI     <- range(res2$confreg)
llikCI <- max(res2$llik) - qchisq(1 - alpha, df = 1)/2
llrng  <- range(res2$llik[is.finite(res2$llik)])
yCI <- llrng + diff(llrng)*0.25
yLR <- (res2$llik[1] + max(res2$llik))/2
cols <- c("purple", "cadetblue3", "gray60")

par(mar = c(3, 2.5, 2, 0.1), mgp = c(1.5, 0.3, 0))
plot(revals[[1]], res2$llik, type = "l", xlab = "r", ylab = "log-likelihood", 
     yaxt = "n")#, tck = -0.01)
abline(v = res2$rhat, lty = 1, col = cols[1])
abline(h = c(max(res2$llik), llikCI, res2$llik[[1]]), lty = 2, col = cols[2])
abline(v = CI, col = cols[1], lty = 5)
arrows(CI[1], yCI, CI[2], yCI, angle = 20, length = 0.1, code = 3, 
       col = cols[3])
arrows(0.15, res2$llik[1], 0.15, max(res2$llik), angle = 20, length = 0.1, 
       code = 3, col = cols[3])
text(0.165, yLR, adj = 0, "0.5 LR statistic", col = cols[3])
text(mean(CI), yCI + 0.05*diff(llrng), "confidence interval", col = cols[3])
text(res2$rhat - 0.025, yLR, "MLE", col = cols[3], srt = 90)
par(pardef)

