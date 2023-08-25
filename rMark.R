library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(matrixStats)
library(patchwork)
library(pheatmap)
library(Seurat)
library(RColorBrewer)
library(reshape2)

if (!require("tiledb", quietly = TRUE))
  remotes::install_github("TileDB-Inc/TileDB-R", force = TRUE)

if (!require("tiledbsc", quietly = TRUE))
  remotes::install_github("tiledb-inc/tiledbsc", force = TRUE)

library(tiledb)
library(tiledbsc)


# File path to TileDB array on local machine
tiledbURI <- "LiverDataRelease_TileDBArray/"


# AWS region of the URI
s3Region <- "us-west-2"

config <- tiledb::tiledb_config()
config["vfs.s3.region"] <- s3Region
ctx <- tiledb::tiledb_ctx(config)


directory
|-- __group
|   |-- __data_specific_alphanumeric_name_here
|   |-- __data_specific_alphanumeric_name_here
|-- __meta
|   |-- __data_specific_alphanumeric_name_here
|-- __tiledb_group.tdb
|-- soma_RNA
|   |-- X
|   |-- __group
|   |-- __meta
|   |-- __tiledb_group.tdb
|   |-- obs
|   |-- obsm
|   |-- obsp
|   |-- uns
|   |-- var
|   |-- varm
|   |-- varp
|-- uns
    |-- __group
    |-- __meta
    |-- __tiledb_group.tdb
    |-- commands



# read in SOMACollection
tiledb_scdataset <- tiledbsc::SOMACollection$new(uri = tiledbURI, 
                                                 verbose = FALSE)


names(tiledb_scdataset$somas)

## [1] "falsecode"      "negprobes"      "RNA"            "RNA_normalized"

names(tiledb_scdataset$somas$RNA$members)

## [1] "X"    "varp" "obs"  "varm" "obsp" "uns"  "obsm" "var"


# batch_mode parallelizes the readin, decreasing computation time
counts <- tiledb_scdataset$somas$RNA$X$members$counts$to_matrix(batch_mode = TRUE) 
dim(counts)

## [1]   1000 793318

counts[1:4,1:4]

## 4 x 4 sparse Matrix of class "dgTMatrix"
##       c_1_100_1000 c_1_100_1017 c_1_100_1028 c_1_100_1043
## AATK             1            1            2            1
## ABL1             .            .            .            .
## ABL2             .            .            .            .
## ACACB            .            .            .            .



norm <- tiledb_scdataset$somas$RNA_normalized$X$members$data$to_matrix(batch_mode = TRUE)
dim(norm)

## [1]   1000 793318

norm[1:4,1:4]

## 4 x 4 sparse Matrix of class "dgTMatrix"
##        c_1_100_1  c_1_100_10 c_1_100_100 c_1_100_1000
## AATK  -0.1429034 -0.06889529  -0.2645790    7.5780320
## ABL1  -0.1972264 -0.09509172  -0.3650734   -0.1790312
## ABL2  -0.1667603 -0.08039924  -0.3087212   -0.1513743
## ACACB -0.2951553 -0.14233400  -0.5460241   -0.2679370



metadata <- tiledb_scdataset$somas$RNA$obs$to_dataframe()
dim(metadata)

## [1] 793318     63

metadata[1:4,1:10]

##              RNA_pca_cluster_default RNA_pca_cluster_default.1 orig.ident nCount_RNA nFeature_RNA nCount_negprobes
## c_1_100_1                         16                        20          c        142           77                0
## c_1_100_10                        13                        15          c         33           26                0
## c_1_100_100                        3                         7          c        487          143                0
## c_1_100_1000                      13                        16          c        117           69                0
##              nFeature_negprobes nCount_falsecode nFeature_falsecode fov
## c_1_100_1                     0                0                  0 100
## c_1_100_10                    0                5                  5 100
## c_1_100_100                   0                1                  1 100
## c_1_100_1000                  0                2                  2 100



cellCoords <- tiledb_scdataset$somas$RNA$obs$to_dataframe(
  attrs = c("x_FOV_px", "y_FOV_px", "x_slide_mm", "y_slide_mm", 
            "slide_ID_numeric", "Run_Tissue_name", "fov"))
head(cellCoords)

##              x_FOV_px y_FOV_px x_slide_mm y_slide_mm slide_ID_numeric Run_Tissue_name fov
## c_1_100_1          42       36    8.70804    9.73368                1     NormalLiver 100
## c_1_100_10       2737       25    9.03144    9.73500                1     NormalLiver 100
## c_1_100_100      2888      409    9.04956    9.68892                1     NormalLiver 100
## c_1_100_1000     3457     3742    9.11784    9.28896                1     NormalLiver 100
## c_1_100_1001     1298     3763    8.85876    9.28644                1     NormalLiver 100
## c_1_100_1002     1001     3738    8.82312    9.28944                1     NormalLiver 100

ggplot(cellCoords, aes(x=x_slide_mm, y=y_slide_mm))+
  geom_point(alpha = 0.05, size = 0.01)+
  facet_wrap(~Run_Tissue_name)+
  coord_equal()+
  labs(title = "Cell coordinates in XY space")


transcriptCoords <- tiledb::tiledb_array(
  tiledb_scdataset$somas$RNA$obsm$members$transcriptCoords$uri,
  return_as="data.frame")[]
head(transcriptCoords)

##   slideID fov x_FOV_px y_FOV_px z_FOV_slice  target CellId     cell_id  CellComp
## 1       1  15     1375        9           1   HMGN2   1271 c_1_15_1271   Nuclear
## 2       1  15       25       10           6  MALAT1      1    c_1_15_1   Nuclear
## 3       1  15       25       10           5    XBP1      1    c_1_15_1   Nuclear
## 4       1  15       36       10           3 SLCO2B1      1    c_1_15_1   Nuclear
## 5       1  15      167       10           6   RPL32     23   c_1_15_23 Cytoplasm
## 6       1  15      201       10           7   APOA1     23   c_1_15_23 Cytoplasm


slide <- 1
fov <- 15

slideName <- unique(cellCoords$Run_Tissue_name[cellCoords$slide_ID_numeric == 
                                                 slide])

fovCoords <- cellCoords[cellCoords$slide_ID_numeric == slide & 
                          cellCoords$fov == fov,]
fovTranscriptCoords <- transcriptCoords[transcriptCoords$slideID == slide & 
                                          transcriptCoords$fov == fov,]

targetCounts <- table(fovTranscriptCoords$target)

targets <- names(targetCounts[which(targetCounts >= 2500 & 
                                      targetCounts <= 3000)])
fovTranscriptCoords <- fovTranscriptCoords[fovTranscriptCoords$target %in% 
                                             targets,]

ggplot(fovCoords, aes(x=x_FOV_px, y=y_FOV_px))+
  geom_point(alpha = 0.6, size = 0.1, color = "black")+
  geom_point(data = fovTranscriptCoords, 
             mapping = aes(x=x_FOV_px, 
                           y=y_FOV_px, 
                           color = target), 
             size = 0.3, alpha = 0.2)+
  theme_bw()+
  coord_equal()+
  guides(colour = guide_legend(override.aes = list(size=1,
                                                   alpha=1)))+
  labs(color = "RNA Target", title = paste0("RNA Transcripts in\n", 
                                            slideName, "\nFOV", fov))



pca <- tiledb_scdataset$somas$RNA$obsm$members$dimreduction_pca$to_matrix()
pca[1:4,1:4]

##                   PCA_1       PCA_2       PCA_3      PCA_4
## c_1_100_1    -1.6417182 -2.33161620 -2.99827387 -0.1049817
## c_1_100_10   -1.1434024 -0.08528228 -1.04553850  0.3209358
## c_1_100_100  -7.0135208  0.52200934  0.09413956  0.1917221
## c_1_100_1000 -0.2238483 -2.22115908 -0.78119275  1.5374492

neighbs <- tiledb_scdataset$somas$RNA$obsp$members$graph_RNA_pca_nn$to_seurat_graph()
neighbs[1:4,1:25]

## 4 x 25 sparse Matrix of class "dgCMatrix"
##                                                               
## c_1_100_1    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 . . . . .
## c_1_100_10   . . . . . . . . . . . . . . . . . . . . 1 1 1 1 1
## c_1_100_100  . . . . . . . . . . . . . . . . . . . . . . . . .
## c_1_100_1000 . . . . . . . . . . . . . . . . . . . . . . . . .



