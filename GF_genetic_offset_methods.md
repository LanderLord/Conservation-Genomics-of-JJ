#Notes:as we analyzed several brown bear (Ursus arctos) populations across the Tibetan Plateau and the Asian mountains, for simplicity, this script will just use all sampled populations as an example.
#This script used for simulating Gradient Forest and three genetic offset metrics (local, forward and backward offset)
#Install the dependencies package
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
rm(list=ls())
library(gradientForest)#version 0.1-24
library(raster); library(dplyr); library(fields); library(geosphere)
library(foreach); library(doParallel)
#Load pop allele frequency (individual x SNP, genotype coded 0/1/2)
geno <- read.delim("417snps_genotype_matrix.tsv", sep=",", header=TRUE, check.names=FALSE)
alfreq <- t(as.matrix(geno[,-1]))/2
colnames(alfreq) <- geno$SNP_ID
#Load sampling coordinates and align sample order
sample_info <- read.csv("sample_geo.csv", stringsAsFactors=FALSE)
sample_info <- sample_info[match(rownames(alfreq), sample_info$SampleID),]
coord <- as.matrix(sample_info[,c("lon","lat")])
#Load bio-climatic variables, in this study, we used elev, bio2, bio8, bio13, and bio18
vec <- c("elev","bio_8","bio_18","bio_2","bio_13")
load_stack <- function(path){
  files <- list.files(path, pattern="\\.tif$", full.names=TRUE)
  stk <- stack(files[match(paste0(vec,".tif"), basename(files))])
  names(stk) <- vec; stk}
bio_cur <- load_stack("current")
bio_fut <- resample(load_stack("BCC-CSM2-MR_ssp126_2041-2060"),
                    bio_cur, method="bilinear")
#Extract climate at sampling sites and drop missing
env_cur <- extract(bio_cur, coord)
env_fut <- extract(bio_fut, coord)
keep <- complete.cases(env_cur)
alfreq <- alfreq[keep,]; env_cur <- env_cur[keep,]; env_fut <- env_fut[keep,]
#Standardize predictors, parameters reused for grid prediction
env_cur_s <- scale(env_cur)
ctr <- attr(env_cur_s,"scaled:center"); scl <- attr(env_cur_s,"scaled:scale")
#GF analysis
gf <- gradientForest(cbind(as.data.frame(env_cur_s), as.data.frame(alfreq)),
                     predictor.vars=vec, response.vars=colnames(alfreq),
                     ntree=500, maxLevel=1, corr.threshold=0.5,
                     compact=TRUE, nbin=201)
saveRDS(gf,"gf_model.rds")
gf
#Predict the current and future genetic composition on ~50000 grids
set.seed(123)
grid_pts <- sampleRandom(bio_cur[[1]], size=50000, sp=TRUE)
grid_df <- as.data.frame(coordinates(grid_pts)); colnames(grid_df) <- c("lon","lat")
grid_cur <- scale(extract(bio_cur, grid_df), center=ctr, scale=scl)
grid_fut <- scale(extract(bio_fut, grid_df), center=ctr, scale=scl)
Trns_grid <- cbind(grid_df, predict(gf, as.data.frame(grid_cur)))
Trns_grid_fut <- cbind(grid_df, predict(gf, as.data.frame(grid_fut)))
gf_cols <- setdiff(colnames(Trns_grid), c("lon","lat"))
#Local genetic offset
Trns_grid$local_offset <- sapply(1:nrow(Trns_grid), function(i){
  rdist(Trns_grid[i,gf_cols], Trns_grid_fut[i,gf_cols])})
#Forward and reverse genetic offset
require(parallel); require(doParallel)
cl <- makeCluster(detectCores()-1); registerDoParallel(cl)
#Forward metric
#x:current transformed genetic composition
#y:future transformed genetic composition of the same grid
forward <- foreach(i=1:nrow(Trns_grid), .packages=c("fields","geosphere"),
                   .combine=rbind) %dopar% {
  onePop <- Trns_grid[i,]
  #nearest future genetic composition, then geographic distance
  d <- rdist(onePop[,gf_cols], Trns_grid_fut[,gf_cols])
  k <- which.min(d)
  data.frame(lon=onePop$lon, lat=onePop$lat, forwardOffset=d[k],
             geoDist=distGeo(onePop[,c("lon","lat")], Trns_grid_fut[k,c("lon","lat")]),
             lon2=Trns_grid_fut$lon[k], lat2=Trns_grid_fut$lat[k])}
#Backward metric
#x:future transformed genetic composition
#y:current transformed genetic composition of the same grid
backward <- foreach(i=1:nrow(Trns_grid_fut), .packages=c("fields","geosphere"),
                    .combine=rbind) %dopar% {
  onePop <- Trns_grid_fut[i,]
  d <- rdist(onePop[,gf_cols], Trns_grid[,gf_cols])
  k <- which.min(d)
  data.frame(lon=onePop$lon, lat=onePop$lat, reverseOffset=d[k],
             geoDist=distGeo(onePop[,c("lon","lat")], Trns_grid[k,c("lon","lat")]),
             lon2=Trns_grid$lon[k], lat2=Trns_grid$lat[k])}
stopCluster(cl)
#Save outliers
write.csv(Trns_grid[,c("lon","lat","local_offset")], "GF_Local_Genetic_Offset.csv", row.names=FALSE)
write.csv(forward, "GF_Forward_Offset.csv", row.names=FALSE)
write.csv(backward, "GF_Backward_Offset.csv", row.names=FALSE)
