### run dropletUtils to detect empty droplets
# conda activate project_lung_exercise_R
# based on https://github.com/CB-postech/2022_KOGO_workshop/blob/main/KOGO_QC1-DropletUtils.md

library(Seurat)
library(scater)
library(DropletUtils)
library(scran)

save_path = '/home/sjcho/projects/4-1BB/DA240060/outs/QC/1.dropletUtils'
rawsce.p1 <- read10xCounts('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool1/outs/multi/count/raw_feature_bc_matrix', type = "sparse", compressed = TRUE)
rawsce.p2 <- read10xCounts('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool2/outs/multi/count/raw_feature_bc_matrix', type = "sparse", compressed = TRUE)

##### Pool-1
counts.p1 <- counts(rawsce.p1)
counts.p1 <- counts.p1[1:(length(rownames(counts.p1)) - 6), ] # delete Antibody counts
br.out.p1 <- barcodeRanks(counts.p1)

png(paste0(save_path, '/DropletUtils_pool1.png'))
plot(br.out.p1$rank, br.out.p1$total, log = "xy", xlab = "Rank", ylab = "Total")

o <- order(br.out.p1$rank)
lines(br.out.p1$rank[o], br.out.p1$fitted[o], col = "red")

e.out <- emptyDrops(counts.p1, lower = 500)  ## Cells that have UMI counts lower than 100 are empty cells.
table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
is.cell <- e.out$FDR <= 0.05

print(sum(is.cell, na.rm=TRUE))
print(table(br.out.p1$rank == sum(is.cell, na.rm=TRUE)))

abline(h=250, col="purple", lty=2)
abline(h=500, col="purple", lty=2)
abline(h=750, col="purple", lty=2)

abline(h=min(br.out.p1$fitted[o], na.rm=TRUE), col="red", lty=2)
abline(h=metadata(br.out.p1)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out.p1)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))
dev.off()

colnames(rawsce.p1) = colData(rawsce.p1)$Barcode
rawsce.p1 <- rawsce.p1[,which(e.out$FDR <= 0.05)]
saveRDS(rawsce.p1, file = paste0(save_path, '/empty_filtered_pool1_lower500.rds'))

###### Pool-2
counts.p2 <- counts(rawsce.p2)
counts.p2 <- counts.p2[1:(length(rownames(counts.p2)) - 5), ] # delete Antibody counts
br.out.p2 <- barcodeRanks(counts.p2)

lower = 600
png(paste0(save_path, '/DropletUtils_pool2_lower', lower, '.png'))
plot(br.out.p2$rank, br.out.p2$total, log = "xy", xlab = "Rank", ylab = "Total")

o <- order(br.out.p2$rank)
lines(br.out.p2$rank[o], br.out.p2$fitted[o], col = "red")

e.out <- emptyDrops(counts.p2, lower = lower)  ## Cells that have UMI counts lower than 100 are empty cells.
table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
is.cell <- e.out$FDR <= 0.05

print(sum(is.cell, na.rm=TRUE))
print(table(br.out.p2$rank == sum(is.cell, na.rm=TRUE)))

abline(h=250, col="purple", lty=2)
abline(h=500, col="purple", lty=2)
abline(h=750, col="purple", lty=2)

abline(h=min(br.out.p2$fitted[o], na.rm=TRUE), col="red", lty=2)
abline(h=metadata(br.out.p2)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out.p2)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), legend=c("knee", "inflection", "FDR_0.05"))
dev.off()

colnames(rawsce.p2) = colData(rawsce.p2)$Barcode
rawsce.p2 <- rawsce.p2[,which(e.out$FDR <= 0.05)]
saveRDS(rawsce.p2, file = paste0(save_path, '/empty_filtered_pool2_lower', lower, '.rds'))
