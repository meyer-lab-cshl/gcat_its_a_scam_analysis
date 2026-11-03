library(Seurat)
library(Signac)
library(UpSetR)

# compute the number of fragmets per barcode
bwa.frag.count <- CountFragments("../peaks/pbmc_granulocyte_sorted_3k_noM_fragments.tsv.gz")
print(dim(bwa.frag.count))
print(head(bwa.frag.count))

# get a list of barcodes that meet the min frags
bwa.keep.bc <- bwa.frag.count[bwa.frag.count$frequency_count >= 1000, ]$CB
print(length(bwa.keep.bc))
print(head(bwa.keep.bc))


# compute the number of fragmets per barcode
cr.frag.count <- CountFragments("../../../og/fragments/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
print(dim(cr.frag.count))
print(head(cr.frag.count))

# get a list of barcodes that meet the min frags
cr.keep.bc <- cr.frag.count[cr.frag.count$frequency_count >= 1000, ]$CB
print(length(cr.keep.bc))
print(head(cr.keep.bc))
cr.keep.bc <- substr(cr.keep.bc, 1, 16)
print(head(cr.keep.bc))

print(length(intersect(bwa.keep.bc, cr.keep.bc)))

bwa.cr <- list(bwa = bwa.keep.bc, cr = cr.keep.bc)
pdf("bwa_cr_cells_upset.pdf", height = 5, width = 5)
upset(fromList(bwa.cr), order.by="freq", mb.ratio = c(0.55, 0.45))
dev.off()
