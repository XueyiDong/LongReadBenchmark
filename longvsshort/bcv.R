library(edgeR)
library(ggplot2)
library(gridExtra)

DIR="/stornext/General/data/user_managed/grpu_mritchie_1/XueyiDong/long_read_benchmark"

s <- catchSalmon(file.path(DIR, "ONT/salmon_bs", list.files(file.path(DIR, "/ONT/salmon_bs"))))
dge <- DGEList(counts=s$counts, genes=s$annotation)

s.short <- catchSalmon(file.path(DIR, "illumina/salmon_bs", list.files(file.path(DIR, "illumina/salmon_bs"))))
dge.short <- DGEList(counts = s.short$counts, genes = s.short$annotation)

group <- rep(c("000", "100", "075", "050", "025"), rep(3, 5))
dge$samples$group <- dge.short$samples$group <- group

# BCV for pure and mixture, calculated together ----

dge.pure <- dge[, 1:6]
dge.mix <- dge[, 7:15]
dge.short.pure <- dge.short[, 1:6]
dge.short.mix <- dge.short[, 7:15]

filt1 <- filterByExpr(dge.pure)
filt2 <- filterByExpr(dge.mix)
dge.pure <- dge.pure[filt1 & filt2, ]
dge.mix <- dge.mix[filt1 & filt2, ]

filt1 <- filterByExpr(dge.short.pure)
filt2 <- filterByExpr(dge.short.mix)
dge.short.pure <- dge.short.pure[filt1 & filt2, ]
dge.short.mix <- dge.short.mix[filt1 & filt2, ]

dge.pure <- calcNormFactors(dge.pure)
dge.mix <- calcNormFactors(dge.mix)
dge.short.pure <- calcNormFactors(dge.short.pure)
dge.short.mix <- calcNormFactors(dge.short.mix)

design.pure <- model.matrix(~dge.pure$samples$group)
dge.pure <- estimateDisp(dge.pure, design.pure)
dge.short.pure <- estimateDisp(dge.short.pure, design.pure)

design.mix <- model.matrix(~dge.mix$samples$group)
dge.mix <- estimateDisp(dge.mix, design.mix)
dge.short.mix <- estimateDisp(dge.short.mix, design.mix)

dge.pure$common.dispersion
dge.mix$common.dispersion
dge.short.pure$common.dispersion
dge.short.mix$common.dispersion

## Plots ----

par(mfrow=c(1,2))
plotBCV(dge.pure)
plotBCV(dge.mix)

par(mfrow=c(1,2))
plotBCV(dge.short.pure)
plotBCV(dge.short.mix)


cor(dge.pure$tagwise.dispersion, dge.mix$tagwise.dispersion, method = "spearman")
cor(dge.short.pure$tagwise.dispersion, dge.short.mix$tagwise.dispersion, method = "spearman")

png("plots/BCV_tagwise_long.png")
plot(dge.pure$tagwise.dispersion, dge.mix$tagwise.dispersion, log = "xy", pch = 16, cex = .3,
     xlab = "pure samples", ylab = "in silico mixture samples", main = "ONT")
abline(0, 1, col = "red")
dev.off()

png("plots/BCV_tagwise_short.png")
plot(dge.short.pure$tagwise.dispersion, dge.short.mix$tagwise.dispersion, log = "xy", pch = 16, cex = .3,
     xlab = "pure samples", ylab = "in silico mixture samples", main = "Illumina")
abline(0, 1, col = "red")
dev.off()

# calculate BCV per group ----
## long ----
dge.group <- lapply(unique(group), function(x){
  tmp <- dge[, dge$samples$group == x]
  # tmp <- tmp[filterByExpr(tmp), ]
  # tmp <- calcNormFactors(tmp)
  # tmp <- estimateDisp(tmp)
  return(tmp)
})

filt <- lapply(dge.group, filterByExpr)
filt <- filt[[1]] & filt[[2]] & filt[[3]] & filt[[4]] & filt[[5]]
table(filt)

dge.group <- lapply(dge.group, function(x){
  tmp <- x[filt, ]
  tmp <- calcNormFactors(tmp)
  tmp <- estimateDisp(tmp)
})

disp <- data.frame(
  disp = sapply(dge.group, function(x){
    x$common.dispersion
  }, simplify = TRUE),
  group = unique(group)
)
disp$group <- factor(disp$group, levels = c("000", "025", "050", "075", "100"))
ggplot(disp, aes(x=group, y=disp))+
  geom_bar(stat="identity")

### plot trended BCV----
yrange <- sapply(dge.group, function(x){
  range(x$trended.dispersion)
}, simplify = TRUE)
yrange <- range(yrange)
xrange <- sapply(dge.group, function(x){
  range(x$AveLogCPM)
}, simplify = TRUE)
xrange <- range(xrange)

plotTagDisp <- function(x, col, lty){
  A <- x$AveLogCPM
  o <- order(A)
  lines(A[o], sqrt(x$trended.dispersion)[o], col = col, lwd = 2, lty = lty)
}

library(RColorBrewer)
col <- brewer.pal(5, "Set1")
lty <- rep(c(1, 2), c(2, 3))

lab <- sapply(dge.group, function(x){
  paste0(x$samples$group[1], ", common BCV=", round(sqrt(x$common.dispersion), 3))
}, simplify = TRUE)
o <- c(1, 5, 4, 3, 2)
pdf("plots/BCV_per_sample_long.pdf", height = 4, width = 8)
layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
par(mar = c(5, 4, 4, 0) + 0.1)
plot(1, 1, type = "n", ylim = sqrt(yrange), xlim = xrange, xlab = "Average log CPM", ylab = "Biological coefficient of variation")
for(i in 1:5){
  plotTagDisp(dge.group[[i]], col = col[i], lty = lty[i])
}
par(mar = c(5, 0, 4, 0) + 0.1)
plot(1,1, type="n", yaxt="n", xaxt="n", ylab="", xlab="", frame.plot=FALSE)
legend("topleft", legend = lab[o], text.col = col[o], col = col[o], lty = lty[o], lwd = 2)
dev.off()

### plot tagwise dispersion----
tagwise <- lapply(dge.group, function(x){
  tmp <- data.frame(
    tagwise.dispersion = x$tagwise.dispersion,
    group = x$samples$group[1]
  )
  return(tmp)
})
tagwise <- Reduce(rbind, tagwise)
tagwise$group <- factor(tagwise$group, levels = c("000", "025", "050", "075", "100"))
pdf("plots/tagwise_long.pdf", height = 4)
ggplot(tagwise, aes(x=group, y=sqrt(tagwise.dispersion), fill=group)) +
  geom_violin()+
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_fill_manual(values = col[c(1, 5, 4, 3, 2)]) +
  # scale_y_continuous(trans = "log") +
  labs(y = "Biological coefficient of variation") +
  theme_bw()
dev.off()

## short ----
dge.short.group <- lapply(unique(group), function(x){
  tmp <- dge.short[, dge.short$samples$group == x]
  # tmp <- tmp[filterByExpr(tmp), ]
  # tmp <- calcNormFactors(tmp)
  # tmp <- estimateDisp(tmp)
  return(tmp)
})

filt <- lapply(dge.short.group, filterByExpr)
filt <- filt[[1]] & filt[[2]] & filt[[3]] & filt[[4]] & filt[[5]]
table(filt)

dge.short.group <- lapply(dge.short.group, function(x){
  tmp <- x[filt, ]
  tmp <- calcNormFactors(tmp)
  tmp <- estimateDisp(tmp)
})
### plot trended BCV----
yrange <- range(sapply(dge.short.group, function(x){
  range(x$trended.dispersion)
}, simplify = TRUE))
xrange <- range(sapply(dge.short.group, function(x){
  range(x$AveLogCPM)
}, simplify = TRUE))
lab <- sapply(dge.short.group, function(x){
  paste0(x$samples$group[1], ", common BCV=", round(sqrt(x$common.dispersion), 3))
}, simplify = TRUE)
pdf("plots/BCV_per_sample_short.pdf", height = 4, width = 8)
layout(matrix(c(1, 2), ncol = 2), widths = c(2, 1))
par(mar = c(5, 4, 4, 0) + 0.1)
plot(1L, 1L, type = "n", ylim = sqrt(yrange), xlim = xrange, xlab = "Average log CPM", ylab = "Biological coefficient of variation")
for(i in 1:5){
  plotTagDisp(dge.short.group[[i]], col = col[i], lty = lty[i])
}
par(mar = c(5, 0, 4, 0) + 0.1)
plot(1,1, type="n", yaxt="n", xaxt="n", ylab="", xlab="", frame.plot=FALSE)
legend("topleft", legend = lab[o], text.col = col[o], col = col[o], lty = lty[o], lwd = 2)
dev.off()

### plot tagwise dispersion----
tagwise <- lapply(dge.short.group, function(x){
  tmp <- data.frame(
    tagwise.dispersion = x$tagwise.dispersion,
    group = x$samples$group[1]
  )
  return(tmp)
})
tagwise <- Reduce(rbind, tagwise)
tagwise$group <- factor(tagwise$group, levels = c("000", "025", "050", "075", "100"))
pdf("plots/tagwise_short.pdf", height = 4)
ggplot(tagwise, aes(x=group, y=sqrt(tagwise.dispersion), fill=group))  +
  geom_violin()+
  geom_boxplot(width = 0.2, outlier.colour = NA, alpha = 0) +
  scale_fill_manual(values = col[c(1, 5, 4, 3, 2)]) +
  labs(y = "Biological coefficient of variation") +
  theme_bw()
dev.off()
