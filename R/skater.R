# SKATER clustering of the unit graph, run in its own R process because rgeoda can crash the caller.
# Usage: Rscript R/skater.R <input.rds> <output.rds>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
suppressPackageStartupMessages({ library(spdep); library(rgeoda); library(data.table) })

inp <- readRDS(args[1])
attr(inp$nb_aug, "region.id") <- as.character(seq_len(inp$n_grids))
gal <- tempfile(fileext = ".gal")
spdep::write.nb.gal(inp$nb_aug, gal, oldstyle = FALSE)
w <- rgeoda::read_gal(gal, as.character(seq_len(inp$n_grids)))
unlink(gal)
saveRDS(data.table(grid_id = inp$grid_id, idx_st_cluster = rgeoda::skater(inp$nc, w, inp$co_df)$Clusters), args[2])
