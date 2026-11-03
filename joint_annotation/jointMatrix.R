library(Seurat)
library(Signac)


## read and process the gex data
gex.dir <- "/grid/meyer/home/kaliapp/data/10x_multiome/pbmc_3k_nogran/gex/mapped_reads/pbmc_granulocyte_sorted_3k_Solo.out/GeneFull_Ex50pAS/filtered/"
gex.matrix <- ReadMtx(mtx = paste0(gex.dir, "matrix.mtx"),
                       cells = paste0(gex.dir, "barcodes.tsv"),
                       features = paste0(gex.dir, "features.tsv"))

pbmc3k <- CreateSeuratObject(counts = gex.matrix, assay = "RNA")
pbmc3k

## read and process the atac data
bwa.frag.file <- "/grid/meyer/home/kaliapp/data/10x_multiome/pbmc_3k_nogran/atac/analysis/compare/paper_plots/pbmc_granulocyte_sorted_3k_noM_fragments.tsv.gz"

# BWA compute the number of fragmets per barcode
bwa.frag.count <- CountFragments(bwa.frag.file)
print(dim(bwa.frag.count))
print(head(bwa.frag.count))

# BWA get a list of barcodes that meet the min frags
bwa.keep.bc <- bwa.frag.count[bwa.frag.count$frequency_count >= 1000, ]$CB
print(length(bwa.keep.bc))

# get the common cells
common.bc <- intersect(bwa.keep.bc, colnames(pbmc3k))


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


# keep just the cells that have both
pbmc3k <- pbmc3k[, common.bc]
pbmc3k

# add in chromatic assay
pbmc3k[["ATAC"]] <- CreateChromatinAssay(
  counts = bwa.count.matrix,
  sep = c(":", "-"),
  fragments = bwa.frag.file
)
pbmc3k

## save the object
SaveSeuratRds(pbmc3k, file = "pbmc3k_joint_matrix.Rds")
