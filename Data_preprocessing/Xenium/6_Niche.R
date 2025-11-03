library(Seurat);library(SeuratObject);library(BPCells)
library(dplyr);library(ggplot2);library(patchwork)
library(data.table);library(here);library(glmGamPoi)
library(tictoc);library(dbscan)
options(future.globals.maxSize = 50 * 1024^3)  
gc()
object=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/cell_state_integration/Seurat_cell_state_curated_localBPcell_Oct2025.rds")
object$cell_type[is.na(object$cell_state)]%>%table()
object$cell_state[is.na(object$cell_state)]=object$cell_type[is.na(object$cell_state)]
#threshold <- quantile(object$ADC_trafficking1[object$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC")], 0.75, na.rm = TRUE)
#object$cell_state[object$cell_state%in%c("Her2E_SC")&object$ADC_trafficking1>threshold]="ADC trafficking Her2E_SC"  
#object=subset(object,subset = orig.ident%in%c("1143_BL","1406_BL", "1408_BL"))
gc()
centroids=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/centroids_metadata.rds")
centroids=centroids[centroids$cell%in%Cells(object),]
all.equal(centroids$cell,Cells(object))
centroids$cell_state=object$cell_state

# adjust coordinate +20000
offset_tbl <- tibble(
  orig.ident = sort(unique(centroids$orig.ident)),
  y_offset   = (seq_along(sort(unique(centroids$orig.ident))) - 1L) * 20000
)

# merge
centroids <- centroids %>% left_join(offset_tbl, by = "orig.ident")
get_decimals <- function(x) {
  sapply(strsplit(format(x, scientific = FALSE), "\\."), function(parts) {
    if (length(parts) == 2) nchar(rtrim(parts[2], "0")) else 0
  })
}

rtrim <- function(x, pattern) sub(paste0(pattern, "+$"), "", x)

centroids$decimals <- get_decimals(centroids$y)
centroids$y <- round(centroids$y + centroids$y_offset, centroids$decimals)
centroids$y_offset <- NULL
centroids$decimals <- NULL
centroids$orig.ident=NULL
all.equal(centroids$cell,Cells(object))
segmentations.data <- list(
  "centroids" = CreateCentroids(centroids)
)
coords <- CreateFOV(
  coords = segmentations.data,
  type = "centroids",
  molecules = NULL,
  assay = "Xenium"
)
object[["fov"]] <- coords
ImageDimPlot(object, cols = "polychrome",fov = "fov", size = 0.75)

# neighorhoods analysis 
mat<- centroids[,c("x","y")]
mat<- as.matrix(mat)
rownames(mat)<- centroids[,"cell"]

# ===== radius analysis =====
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(dbscan)
library(rlang)
# Seurat 
suppressMessages(require(Seurat))

evaluate_radii <- function(
    mat,                      # numeric matrix, rownames = cell IDs; cols = coordinates (x,y) or higher-D
    object,                   # Seurat object，object$cell_state；Cells(object) rownames(mat)
    eps_values = seq(20, 100, by = 10),
    exclude_self = TRUE,     
    build_seurat = TRUE,    
    seed = 123
){
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required. Please install.packages('dbscan').")
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required. Please install.packages('tibble').")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Base 'stats' pkg not available?")
  }
  if (isTRUE(build_seurat) && !requireNamespace("Seurat", quietly = TRUE)) {
    warning("Package 'Seurat' not found. build_seurat will be set to FALSE.")
    build_seurat <- FALSE
  }
  
  stopifnot(is.matrix(mat))
  if (is.null(rownames(mat))) stop("mat must have rownames as cell IDs.")
  if (!all(Seurat::Cells(object) %in% rownames(mat))) {
    stop("Some Cells(object) are not present in rownames(mat). Please align.")
  }
  
  mat <- mat[Seurat::Cells(object), , drop = FALSE]
  
  # cell_state
  cs <- object$cell_state
  if (is.null(cs)) stop("object$cell_state not found. Please add a 'cell_state' column to meta.data or object@meta.data.")
  if (!is.factor(cs)) cs <- factor(cs)
  cs_levels <- levels(cs)
  if (length(cs_levels) == 0L) stop("cell_state has zero levels; please check your annotations.")
  names(cs) <- Seurat::Cells(object)  #
  
  # ---- Shannon （0-1）----
  .shannon_norm <- function(counts_vec) {
    s <- sum(counts_vec)
    if (s <= 0) return(0)
    p <- counts_vec / s
    p <- p[p > 0]
    H <- -sum(p * log(p))
    K <- length(counts_vec)
    if (K <= 1) return(0)
    H / log(K)
  }
  
  nn_obj_list <- list()
  summary_rows <- list()
  
  set.seed(seed)
  
  # ---- loop radius----
  for (eps in eps_values) {
    cat("===== Running eps =", eps, "=====\n")
    t0 <- proc.time()
    nn <- dbscan::frNN(x = mat, eps = eps, sort = FALSE)
    time_sec <- as.numeric((proc.time() - t0)["elapsed"])
    
    idx_list     <- nn$id
    neighbor_idx <- unlist(idx_list, use.names = FALSE)
    query_idx    <- rep.int(seq_len(nrow(mat)), lengths(idx_list))
    
    if (length(neighbor_idx) == 0L) {
      nn_mat <- matrix(
        0L,
        nrow = length(Seurat::Cells(object)),
        ncol = length(cs_levels),
        dimnames = list(Seurat::Cells(object), cs_levels)
      )
    } else {
      if (isTRUE(exclude_self)) {
        keep <- neighbor_idx != query_idx
        neighbor_idx <- neighbor_idx[keep]
        query_idx    <- query_idx[keep]
      }
      
      if (length(neighbor_idx) == 0L) {
        nn_mat <- matrix(
          0L,
          nrow = length(Seurat::Cells(object)),
          ncol = length(cs_levels),
          dimnames = list(Seurat::Cells(object), cs_levels)
        )
      } else {
        query_cell    <- rownames(mat)[query_idx]
        neighbor_cell <- rownames(mat)[neighbor_idx]
        
        neighbor_cs <- factor(cs[neighbor_cell], levels = cs_levels)
        
        # query_cell × cell_state matrices
        tab <- as.data.frame.matrix(
          stats::xtabs(~ query_cell + neighbor_cs, drop.unused.levels = FALSE)
        )
        
        # exports
        nn_mat <- matrix(
          0L,
          nrow = length(Seurat::Cells(object)),
          ncol = length(cs_levels),
          dimnames = list(Seurat::Cells(object), cs_levels)
        )
        rn <- intersect(rownames(nn_mat), rownames(tab))
        cn <- intersect(colnames(nn_mat), colnames(tab))
        if (length(rn) > 0 && length(cn) > 0) {
          nn_mat[rn, cn] <- as.matrix(tab[rn, cn, drop = FALSE])
        }
      }
    }
    
    # composition matrices
    row_tot  <- rowSums(nn_mat)
    prop_mat <- nn_mat
    nz <- row_tot > 0
    prop_mat[nz, ]  <- nn_mat[nz, , drop = FALSE] / row_tot[nz]
    prop_mat[!nz, ] <- 0
    
    # richness
    richness <- rowSums(nn_mat > 0)
    shannon_norm <- apply(nn_mat, 1, .shannon_norm)
    
    cell_metrics <- tibble::tibble(
      cell         = Seurat::Cells(object),
      neighbors    = as.integer(row_tot),
      richness     = as.integer(richness),
      shannon_norm = as.numeric(shannon_norm)
    )
    
    # summary
    sumrow <- tibble::tibble(
      eps                  = eps,
      time_sec             = time_sec,
      n_cells              = nrow(nn_mat),
      neighbors_mean       = mean(cell_metrics$neighbors),
      neighbors_median     = stats::median(cell_metrics$neighbors),
      neighbors_p5         = as.numeric(stats::quantile(cell_metrics$neighbors, 0.05, names = FALSE, type = 7)),
      neighbors_p95        = as.numeric(stats::quantile(cell_metrics$neighbors, 0.95, names = FALSE, type = 7)),
      frac_zero_neighbors  = mean(cell_metrics$neighbors == 0),
      richness_mean        = mean(cell_metrics$richness),
      shannon_mean         = mean(cell_metrics$shannon_norm)
    )
    summary_rows[[as.character(eps)]] <- sumrow
    
    # Seurat（features=cell_state, cells）
    seurat_obj <- NULL
    if (isTRUE(build_seurat)) {
      counts_for_seurat <- t(nn_mat)  # X features ，Y cells
      seurat_obj <- suppressMessages(
        Seurat::CreateSeuratObject(
          counts   = counts_for_seurat,
          meta.data = object@meta.data[Seurat::Cells(object), , drop = FALSE]
        )
      )
    }
    nn_obj_list[[paste0("eps_", eps)]] <- list(
      eps          = eps,
      nn           = nn,
      nn_mat       = nn_mat,
      prop_mat     = prop_mat,
      cell_metrics = cell_metrics,
      seurat       = seurat_obj
    )
  }
  
  summary_by_eps <- do.call(rbind, summary_rows)
  rownames(summary_by_eps) <- NULL
  
  list(
    summary = summary_by_eps,  
    details = nn_obj_list      
  )
}


eps_values=seq(10,100,10)
res <- evaluate_radii(
  mat = mat,                 
  object = object,           
  eps_values = eps_values,
  exclude_self = TRUE,
  build_seurat = TRUE
)
a=res$summary
saveRDS(a,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/evaluate_radii.rds")


a=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/evaluate_radii.rds")
library(ggplot2)
p=ggplot(a, aes(x = eps)) +
  geom_line(aes(y = richness_mean, color = "Richness")) +
  geom_point(aes(y = richness_mean, color = "Richness")) +
  geom_line(aes(y = shannon_mean * max(richness_mean), color = "Shannon (scaled)")) +
  geom_point(aes(y = shannon_mean * max(richness_mean), color = "Shannon (scaled)")) +
  scale_y_continuous(sec.axis = sec_axis(~ ./max(a$richness_mean),
                                         name = "Shannon mean (0–1)")) +
  scale_x_continuous(
    breaks = seq(10, 100, by = 10)
  )+
  theme_minimal(base_size = 13) +
  labs(x = "Radius (µm)",
       y = "Richness mean (# of cell types)",
       title = "Richness vs. Shannon across radius") +
  scale_color_manual(values = c("Richness"="#377eb8", "Shannon (scaled)"="#e41a1c"))
ggsave(p,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS10/Radius_Richness_Shannon.pdf")
rm(res)
gc()



#-----------------------------
# create VoltRon
#-----------------------------
library(VoltRon)
library(ComplexHeatmap)
coordinates <- centroids[, c("x", "y")]
row.names(coordinates) <- centroids$cell
metadata <- object@meta.data[, c("cell_type", "cell_state")]
stopifnot(all.equal(row.names(coordinates), row.names(metadata)))
vr_object <- formVoltRon(metadata = metadata, coords = coordinates)
gc()

r=40
Xen_data <- getSpatialNeighbors(vr_object, radius = r, method = "radius")
Xen_data <- getNicheAssay(Xen_data, label = "cell_state", graph.type = "radius")
vrMainFeatureType(Xen_data) <- "Niche"
Xen_data <- normalizeData(Xen_data, method = "CLR")
  

nclus_values  <- 5:12
for (k in nclus_values) {
    message("  Clustering with nclus = ", k)
    # K-means
    Xen_data <- getClusters(Xen_data, nclus = k, method = "kmeans", label = paste0("niche_clusters","_",k))
    p=vrHeatmapPlot(
      Xen_data,
      features = vrFeatures(Xen_data),
      show_row_names = TRUE,
      group.by = paste0("niche_clusters","_",k)
    )
    png(sprintf("E:Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche/VoltRon_heatmap_radius%d_k%d.png", r, k), width = 4000, height =2000)
    draw(p)                                   # render heatmap
    dev.off()                                 # close device
}

saveRDS(Xen_data,file="E:Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche/Xen_data_radius70.rds")
  



# radius 40, k=8 selected 
library(VoltRon)
library(ComplexHeatmap)
coordinates <- centroids[, c("x", "y")]
row.names(coordinates) <- centroids$cell
metadata <- object@meta.data[, c("cell_type", "cell_state")]
stopifnot(all.equal(row.names(coordinates), row.names(metadata)))
vr_object <- formVoltRon(metadata = metadata, coords = coordinates)
gc()

r=50
Xen_data <- getSpatialNeighbors(vr_object, radius = r, method = "radius")
Xen_data <- getNicheAssay(Xen_data, label = "cell_state", graph.type = "radius")
vrMainFeatureType(Xen_data) <- "Niche"
Xen_data <- normalizeData(Xen_data, method = "CLR")
Xen_data <- getClusters(Xen_data, nclus = 8, method = "kmeans", label = paste0("niche_clusters"))
library(circlize)
col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "#FFFFFF", "#B2182B")         # ColorBrewer 风格
)
vrHeatmapPlot(
    Xen_data,show_heatmap_legend=T,
    col = col_fun,
    features = vrFeatures(Xen_data),
    show_row_names = TRUE,
    group.by = paste0("niche_clusters")
  )
Xen_data@metadata@cell$niche_clusters=paste0("Niche",Xen_data@metadata@cell$niche_clusters)
res=Xen_data@metadata@cell
res$sampleID=str_extract(res$id.1, "^[^_]+_[^_]+")
res$cell.name=res$id.1
res$cell_id <- str_extract(res$id.1, "[^_]+-[0-9]+$")
res=res[,c("sampleID","cell_id","cell.name","niche_clusters","cell_state")]
res$niche_clusters <- recode(
  res$niche_clusters,
  "Niche1" = "Niche1_Plasma",
  "Niche2" = "Niche2_iCAF",
  "Niche3" = "Niche3_Stromal",
  "Niche4" = "Niche4_Immunoactive",
  "Niche5" = "Niche5_Luminal_tumor",
  "Niche6" = "Niche6_HER2E_tumor",
  "Niche7" = "Niche7_Basal_tumor",
  "Niche8" = "Niche8_Blood_vessel"
)

saveRDS(res,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/Figure6/Niche_meta.rds")


