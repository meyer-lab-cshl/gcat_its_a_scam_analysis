library(Seurat)
library(Signac)


bwa.frag.file <- "pbmc_granulocyte_sorted_3k_noM_fragments.tsv.gz"
cr.frag.file <- "../../../og/fragments/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"


# BWA compute the number of fragmets per barcode
bwa.frag.count <- CountFragments("pbmc_granulocyte_sorted_3k_noM_fragments.tsv.gz")
print(dim(bwa.frag.count))
print(head(bwa.frag.count))

# BWA get a list of barcodes that meet the min frags
bwa.keep.bc <- bwa.frag.count[bwa.frag.count$frequency_count >= 1000, ]$CB
print(length(bwa.keep.bc))


# CR compute the number of fragmets per barcode
cr.frag.count <- CountFragments("../../../og/fragments/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz")
print(dim(cr.frag.count))
print(head(cr.frag.count))

# CR get a list of barcodes that meet the min frags
cr.keep.bc <- cr.frag.count[cr.frag.count$frequency_count >= 1000, ]$CB
print(length(cr.keep.bc))
cr.keep.bc <- substr(cr.keep.bc, 1, 16)

common.bc <- intersect(bwa.keep.bc, cr.keep.bc)
print(length(common.bc))



# BWA create a fragment object
bwa.frags <- CreateFragmentObject(path = bwa.frag.file, cells = common.bc)

# call peaks
peaks <- CallPeaks(bwa.frags, macs2.path = "macs3")
print(dim(peaks))
print(head(peaks))

# create cell-peak matrix
bwa.count.matrix <- FeatureMatrix(fragments = bwa.frags,
  features = peaks, cells = common.bc, sep = c(":", "-"))
print(dim(bwa.count.matrix))

# create chromatic assay
bwa.chrom.assay <- CreateChromatinAssay(counts = bwa.count.matrix,
  sep = c(":", "-"), fragments = bwa.frag.file,
  min.cells = 0, min.features = 0)

# create seruat assay
bwa <- CreateSeuratObject(counts = bwa.chrom.assay, assay = "peaks")
print(bwa)

SaveSeuratRds(bwa, file = "bwa_matrix.Rds")



# create a fragment object
print(head(common.bc))
common.bc <- paste0(common.bc, "-1")
print(head(common.bc))
cr.frags <- CreateFragmentObject(path = cr.frag.file, cells = common.bc)

# create cell-peak matrix
cr.count.matrix <- FeatureMatrix(fragments = cr.frags,
  features = peaks, cells = common.bc, sep = c(":", "-"))
print(dim(cr.count.matrix))

# create chromatic assay
cr.chrom.assay <- CreateChromatinAssay(counts = cr.count.matrix,
  sep = c(":", "-"), fragments = cr.frag.file,
  min.cells = 0, min.features = 0)

# create seruat assay
cr <- CreateSeuratObject(counts = cr.chrom.assay, assay = "peaks")
print(cr)

SaveSeuratRds(cr, file = "cr_matrix.Rds")
