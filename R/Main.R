library(terra)
library(sf)
library(leastcostpath)
library(cppRouting)
library(spatstat)
library(maptools)
library(gtools)
library(splines)
library(dplyr)
library(tmap)
library(ggplot2)

source("./R/Functions.R")
set.seed(NULL)
set.seed(1)

agg <- 2

ext <- terra::ext(150000, 500000, 100000, 400000)

dem2 <- terra::rast("./Data/DEM/RW_100m.tif")

rivers <- readRDS("./Data/Rivers/EU_RN_5.RDS")
rivers <- sf::st_transform(rivers, sf::st_crs(dem2))
rivers <- sf::st_crop(rivers, ext)

polygony <- st_make_grid(ext, square = T, cellsize = c(10000, 10000)) %>%
  st_sf()

polygony_vect <- terra::vect(polygony)
polygony_vect2 <- terra::as.points(polygony_vect)

polygony_vect <- terra::extract(dem2, polygony_vect2, cells = TRUE)
polygony_vect <- polygony_vect[!is.na(polygony_vect$HP40),]
cells <- polygony_vect$cell
cells <- unique(cells)

slope_cs <- leastcostpath::create_slope_cs(x = dem2, cost_function = "herzog", neighbours = 8)
river_cs <- create_slope_cs(x = dem2, 
                             cost_function = function(x) { ifelse(x >= 0, 0.6/3.6, 2.5/3.6) }, neighbours = 8)
slope_cs2 <- leastcostpath::replace_values(x = slope_cs, y = river_cs, sf = rivers)
slope_cs3 <- leastcostpath::update_values(x = slope_cs, sf = rivers, FUN = function(j) { replace(x = j, values = 0)})

create_nc_surface(cs = slope_cs, cells = cells, dem = dem2, filename = paste0("RW_NC_NR_", 50*agg, "m"), ncores = 25)
create_nc_surface(cs = slope_cs2, cells = cells, dem = dem2, filename = paste0("RW_NC_R_", 50*agg, "m"), ncores = 25)
create_nc_surface(cs = slope_cs3, cells = cells, dem = dem2, filename = paste0("RW_NC_BR_", 50*agg, "m"), ncores = 25)

nc_nr <- terra::rast(paste0("./Output/Natural_Corridors/RW_NC_NR_", 50*agg, "m.tif"))
nc_r <- terra::rast(paste0("./Output/Natural_Corridors/RW_NC_R_", 50*agg, "m.tif"))
nc_br <- terra::rast(paste0("./Output/Natural_Corridors/RW_NC_BR_", 50*agg, "m.tif"))

window_size <- seq(3, 400, 2)
window_size <- window_size[c(1:75, seq(80, 180, 20))]

for(i in window_size) {
  print(paste0(i, " of ", max(window_size)))

  win_matrix <- terra::focalMat(x = nc_nr, max(terra::res(nc_nr))*i, "circle")
  
  nc_nr2 <- terra::focal(x = nc_nr, w = win_matrix, fun = "sum", na.rm = TRUE, na.policy = "omit")
  nc_nr2 <- normalise_raster(nc_nr2)
  terra::writeRaster(nc_nr2, filename = paste0("./Output/Natural_Corridors/RW_NC_NR_", max(terra::res(nc_nr2))*i, "m.tif"), overwrite = TRUE)
   
  nc_r2 <- terra::focal(x = nc_r, w = win_matrix, fun = "sum", na.rm = TRUE, na.policy = "omit")
  nc_r2 <- normalise_raster(nc_r2)
  terra::writeRaster(nc_r2, filename = paste0("./Output/Natural_Corridors/RW_NC_R_", max(terra::res(nc_r2))*i, "m.tif"), overwrite = TRUE)
  
  nc_br2 <- terra::focal(x = nc_br, w = win_matrix, fun = "sum", na.rm = TRUE, na.policy = "omit")
  nc_br2 <- normalise_raster(nc_br2)
  terra::writeRaster(nc_br2, filename = paste0("./Output/Natural_Corridors/RW_NC_BR_", max(terra::res(nc_br2))*i, "m.tif"), overwrite = TRUE)

}

nc_nr2_filepaths <- list.files("./Output/Natural_Corridors/", pattern = ".tif", full.names = TRUE)
nc_nr2_filepaths <- nc_nr2_filepaths[grep("RW_NC_NR", nc_nr2_filepaths)]
nc_nr2_filepaths <- gtools::mixedsort(nc_nr2_filepaths)

nc_r2_filepaths <- list.files("./Output/Natural_Corridors/", pattern = ".tif", full.names = TRUE)
nc_r2_filepaths <- nc_r2_filepaths[grep("RW_NC_R", nc_r2_filepaths)]
nc_r2_filepaths <- gtools::mixedsort(nc_r2_filepaths)

nc_br2_filepaths <- list.files("./Output/Natural_Corridors/", pattern = ".tif", full.names = TRUE)
nc_br2_filepaths <- nc_br2_filepaths[grep("RW_NC_BR", nc_br2_filepaths)]
nc_br2_filepaths <- gtools::mixedsort(nc_br2_filepaths)

dem2_ext <- dem2
dem2_ext[!is.na(dem2_ext)] <- 1
dem2_ext_polygon <- terra::as.polygons(dem2_ext, aggregate = TRUE)
dem2_ext_polygon <- sf::st_as_sf(dem2_ext_polygon)
dem2_ext_polygon <- sf::st_cast(dem2_ext_polygon, "POLYGON")
dem2_ext_polygon <- dem2_ext_polygon[c(18,84,85),]
dem2_ext_polygon <- sf::st_crop(dem2_ext_polygon, terra::ext(150000, 400000, 160000, 400000))
dem2_ext_window <- as.owin(dem2_ext_polygon)

tribes <- sf::st_read("./Data/Sites/tribes.gpkg")
tribes <- tribes[tribes$Name %in% c("DECEANGLI", "ORDOVICES", "SILURES"),]
sites <- st_read("./Data/Sites/Sites_Wales.csv", options=c("X_POSSIBLE_NAMES=X","Y_POSSIBLE_NAMES=Y"), crs = sf::st_crs(27700))
sites <- sites[sites$Type == "Campaign base/Fort",]
sites_coords <- sf::st_coordinates(sites)

sites_ppp <- ppp(x = sites_coords[,1], y = sites_coords[,2], window= dem2_ext_window)
unitname(sites_ppp) <- "metres"
sites_ppp <- spatstat.geom::rescale(X = sites_ppp, s = 1000, unitname = "km")

res_val = rep(NA, length(nc_nr2_filepaths))
AIC_vals <- rep(NA, length(nc_nr2_filepaths))
r2_vals <- rep(NA, length(nc_nr2_filepaths))

nc_nr_df <- data.frame(res_val, AIC_vals, r2_vals)
nc_r_df <- data.frame(res_val, AIC_vals, r2_vals)
nc_br_df <- data.frame(res_val, AIC_vals, r2_vals)

# quadrature every 0.5km (500m)
sites_ppp2 <- quadscheme(data = sites_ppp, method = "grid", eps = 0.5)

for(rast_index in 1:length(nc_nr2_filepaths)) {
  
  print(rast_index)
  
  nc_nr2_rast <- terra::rast(nc_nr2_filepaths[rast_index])
  nc_r2_rast <- terra::rast(nc_r2_filepaths[rast_index])
  nc_br2_rast <- terra::rast(nc_br2_filepaths[rast_index])
  
  nc_nr2_rast2 <- maptools::as.im.RasterLayer(as(nc_nr2_rast, "Raster"))
  unitname(nc_nr2_rast2) <- "metres"
  nc_nr2_rast2 <- spatstat.geom::rescale(X = nc_nr2_rast2, s = 1000, unitname = "km")
  
  nc_r2_rast2 <- maptools::as.im.RasterLayer(as(nc_r2_rast, "Raster"))
  unitname(nc_r2_rast2) <- "metres"
  nc_r2_rast2 <- spatstat.geom::rescale(X = nc_r2_rast2, s = 1000, unitname = "km")
  
  nc_br2_rast2 <- maptools::as.im.RasterLayer(as(nc_br2_rast, "Raster"))
  unitname(nc_br2_rast2) <- "metres"
  nc_br2_rast2 <- spatstat.geom::rescale(X = nc_br2_rast2, s = 1000, unitname = "km")
  
  nc_nr2_ppm <- ppm(sites_ppp2 ~ bs(nc), data = list(nc = nc_nr2_rast2), use.gam=TRUE)
  nc_r2_ppm <- ppm(sites_ppp2 ~ bs(nc), data = list(nc = nc_r2_rast2), use.gam=TRUE)
  nc_br2_ppm <- ppm(sites_ppp2 ~ bs(nc), data = list(nc = nc_br2_rast2), use.gam=TRUE)
  
  nc_nr_df[rast_index, 3] <- pseudoR2(nc_nr2_ppm)*100
  nc_nr_df[rast_index, 2] <- AIC(nc_nr2_ppm)
  nc_nr_df[rast_index, 1] <- c(max(res(dem2)), window_size*max(res(dem2)))[rast_index]
    
  nc_r_df[rast_index, 3] <- pseudoR2(nc_r2_ppm)*100
  nc_r_df[rast_index, 2] <- AIC(nc_r2_ppm)
  nc_r_df[rast_index, 1] <- c(max(res(dem2)), window_size*max(res(dem2)))[rast_index]
  
  nc_br_df[rast_index, 3] <- pseudoR2(nc_br2_ppm)*100
  nc_br_df[rast_index, 2] <- AIC(nc_br2_ppm)
  nc_br_df[rast_index, 1] <- c(max(res(dem2)), window_size*max(res(dem2)))[rast_index]
  
}

nc_nr_df$model <- "Natural corridors no rivers"
nc_r_df$model <- "Natural corridors with rivers"
nc_br_df$model <- "Natural corridors with rivers as barriers"

nc_df <- rbind(nc_nr_df, nc_r_df, nc_br_df)
colnames(nc_df)[4] <- "Scenario"
write.csv(nc_df, "./Output/nc_AIC.csv", row.names = FALSE)

### Figure 1 ######

figure1A <- tm_shape(dem2, raster.downsample = FALSE) + 
  tm_raster(style = "pretty", palette = grey.colors(5), legend.is.portrait = TRUE, title = "Elevation (m)", n = 5, legend.show = TRUE, alpha = 1) +
  tm_shape(rivers) + 
  tm_lines(col = "#9ECAE1", lty = 1, lwd = 1) + 
  tm_shape(sites[sites$Type == "Campaign base/Fort",]) + 
  tm_dots(col = "white", shape = 22, size = 0.1, border.col = "black") + 
  tm_shape(sites[sites$Name %in% c("Clyro", "Rhyn Park"),]) + 
  tm_dots(col = "red", shape = 22, size = 0.2, border.col = "black") + 
  tm_shape(sites[sites$Name %in% c("Clyro", "Rhyn Park"),]) + 
  tm_text("Name", size = 0.6, just = "top", ymod = -0.25, col = "white") + 
  tm_shape(tribes) + 
  tm_text("Name", size = 0.6, just = "top", col = "white", fontface = "bold") + 
  tm_layout(frame = TRUE, panel.label.bg.color = NA) + 
  tm_layout(legend.position = c("left", "top"))

tmap::tmap_save(tm = figure1A, filename = "./Output/Figures/Figure1A.png", width = 5, height = 5, dpi = 300)

### Figure 2 ######

nc_nr <- terra::rast(paste0("./Output/Natural_Corridors/RW_NC_NR_", 50*agg, "m.tif"))
nc_r <- terra::rast(paste0("./Output/Natural_Corridors/RW_NC_R_", 50*agg, "m.tif"))

nc_nr <- normalise_raster(nc_nr)
nc_r <- normalise_raster(nc_r)

cs_nc_nr <- leastcostpath::rasterise(slope_cs)
cs_nc_r <- leastcostpath::rasterise(slope_cs2)

cs_nc_nr <- normalise_raster(cs_nc_nr)
cs_nc_r <- normalise_raster(cs_nc_r)

names(nc_nr) <- "A"
names(cs_nc_nr) <- "B"

figure2A <- tm_shape(cs_nc_nr, raster.downsample = FALSE) + 
  tm_raster(style = "quantile", palette = "cividis", legend.is.portrait = TRUE, title = "Conductance", n = 5, labels = c("Low", rep("", 3), "High")) + 
  tm_layout(frame = TRUE, panel.label.bg.color = NA) + 
  tm_layout(legend.position = c("left", "top"), 
            main.title = "A", 
            main.title.position = "left",
            title.position = c('left', 'top'))

polygony_vect3 <- polygony_vect2
polygony_vect3 <- sf::st_as_sf(polygony_vect3)
sf::st_crs(polygony_vect3) <- sf::st_crs(dem2)
polygony_vect3 <- polygony_vect3[!is.na(terra::extract(dem2, polygony_vect3)$HP40),]

figure2B <- tm_shape(nc_nr, raster.downsample = FALSE) + 
  tm_raster(style = "fisher", palette = "cividis", legend.is.portrait = TRUE, title = "Density of LCPs crossing each cell", n = 10, labels = c("Low", rep("", 8), "High")) + 
  tm_shape(polygony_vect3) + 
  tm_dots(col = "grey85") + 
  tm_layout(frame = TRUE, panel.label.bg.color = NA) + 
  tm_layout(legend.position = c("left", "top"), 
            main.title = "B", 
            main.title.position = "left",
            title.position = c('left', 'top'))

Figure2 <- tmap::tmap_arrange(figure2A, figure2B, ncol = 2)
tmap::tmap_save(tm = Figure2, filename = "./Output/Figures/Figure2.png", width = 10, height = 5, dpi = 300)

### Figure 3 ######

Figure3A_rasts <- (c(terra::rast(nc_nr2_filepaths[5]), terra::rast(nc_nr2_filepaths[25]), terra::rast(nc_nr2_filepaths[50]), terra::rast(nc_nr2_filepaths[75])))
names(Figure3A_rasts) <- c("900m", "4,900m", "9,900m", "14,900m")
Figure3B_rasts <- (c(terra::rast(nc_r2_filepaths[5]), terra::rast(nc_r2_filepaths[25]), terra::rast(nc_r2_filepaths[50]), terra::rast(nc_r2_filepaths[75])))
names(Figure3B_rasts) <- c("900m", "4,900m", "9,900m", "14,900m")
Figure3C_rasts <- (c(terra::rast(nc_br2_filepaths[5]), terra::rast(nc_br2_filepaths[25]), terra::rast(nc_br2_filepaths[50]), terra::rast(nc_br2_filepaths[75])))
names(Figure3C_rasts) <- c("900m", "4,900m", "9,900m", "14,900")

Figure3A <- tm_shape(Figure3A_rasts, raster.downsample = FALSE) + 
  tm_raster(style = "fisher", palette = "cividis", legend.is.portrait = FALSE, title = "Density of LCPs\ncrossing cell", n = 10, labels = c("Low", rep("", 8), "High")) + 
  tm_facets(ncol = 4) + 
  tm_layout(main.title = "A", 
            main.title.position = "left",
            title.position = c('left', 'top'),
            legend.outside = T,
            legend.outside.position = "bottom")

Figure3B <- tm_shape(Figure3B_rasts, raster.downsample = FALSE) + 
  tm_raster(style = "fisher", palette = "cividis", legend.is.portrait = FALSE, title = "Density of LCPs\ncrossing cell", n = 10, labels = c("Low", rep("", 8), "High")) + 
  tm_facets(ncol = 4) + 
  tm_layout(main.title = "B", 
            main.title.position = "left",
            title.position = c('left', 'top'),
            legend.outside = T,
            legend.outside.position = "bottom")

Figure3C <- tm_shape(Figure3C_rasts, raster.downsample = FALSE) + 
  tm_raster(style = "fisher", palette = "cividis", legend.is.portrait = FALSE, title = "Density of LCPs crossing cell", n = 10, labels = c("Low", rep("", 8), "High")) + 
  tm_facets(ncol = 4) + 
  tm_layout(main.title = "C", 
            main.title.position = "left",
            title.position = c('left', 'top'),
            legend.outside = T,
            legend.outside.position = "bottom")

Figure3 <- tmap::tmap_arrange(Figure3A, Figure3B, Figure3C, ncol = 1)
tmap::tmap_save(tm = Figure3, filename = "./Output/Figures/Figure3.png", width = 10, height = 10, dpi = 300)

##### FIGURE 4 ######

Figure4A <- tm_shape(terra::rast(nc_nr2_filepaths[6]), raster.downsample = FALSE) + 
  tm_raster(style = "fisher", palette = "cividis", legend.is.portrait = TRUE, title = "Density of LCPs crossing each cell", n = 10, labels = c("Low", rep("", 8), "High")) + 
  tm_shape(sites[sites$Type == "Campaign base/Fort",]) + 
  tm_dots(col = "red", shape = 22, size = 0.5, border.col = "black") + 
  tm_layout(main.title = "A B", 
            main.title.position = "left",
            title.position = c('left', 'top'))

tmap::tmap_save(tm = Figure4A, filename = "./Output/Figures/Figure4A.png", width = 10, height = 10, dpi = 300)

nc_df$AIC_vals_diffs <- nc_df$AIC_vals - min(nc_df$AIC_vals)

Figure4B <- ggplot(nc_df) + 
  geom_vline(xintercept = unique(nc_df$res_val)[which.min(nc_df$AIC_vals)]/1000,lty = 2, colour = "grey40") +
  geom_line(aes(x = res_val/1000, y = AIC_vals_diffs, colour = Scenario, group = Scenario), linewidth = 1) +
  geom_point(aes(x = res_val/1000, y = AIC_vals_diffs, group = Scenario), size = 0.75, shape=21, fill = "black") +
  scale_x_continuous(labels = scales::comma, breaks = seq(0, 37, 2), expand = c(0,1)) +
  scale_colour_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) + 
  scale_fill_manual(values = c("#4daf4a", "#377eb8", "#e41a1c")) + 
  labs(x = "Scale of Window (km)", y = "Δ Akaike information criterion (AIC) score") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  guides(colour= guide_legend(nrow=3,byrow=TRUE, title.position = "top")) + 
  coord_flip()

ggplot2::ggsave(plot = Figure4B, filename = "./Output/Figures/Figure4B.png", dpi = 300, width = 6, height = 10)
