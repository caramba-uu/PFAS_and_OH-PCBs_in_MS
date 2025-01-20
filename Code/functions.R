#####################
#
# cell_fun_bottom - cell function for complex heatmap bottom triangel
#
#####################
cell_fun_bottom = function(j, i, x, y, w, h, fill){
  if(i == j & as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
    grid.rect(x, y, w, h, gp = gpar(fill = "gray", col = "gray"))
  }else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
  }
  if (!is.na(fcP[i, j]) & fcP[i, j]  < 0.001 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
    grid.text("***", x, y, vjust = 0.7, gp = gpar(fontsize = 6))
  }else if (!is.na(fcP[i, j]) & fcP[i, j]  < 0.01 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
    grid.text("**", x, y, vjust = 0.7, gp = gpar(fontsize = 6))
  } else if (!is.na(fcP[i, j]) & fcP[i, j]  <= 0.05 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
    grid.text("*", x, y, vjust = 0.7, gp = gpar(fontsize = 6))
  }
}

#####################
#
# cell_fun_bottom_r - cell function for complex heatmap bottom triangel including r values in the cells
#
#####################
cell_fun_bottom_r <- function(j, i, x, y, w, h, fill){
  if(i == j & as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
    grid.rect(x, y, w, h, gp = gpar(fill = "gray", col = "gray"))
  }else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    grid.text(round(fc[i,j], digits = 2), x, y, vjust = 0.7, gp = gpar(fontsize = 6))
  }
}

#######################
#
# forestPlotAllComp <- Forestplot of all compounds (PFAS)
#
######################
forestPlotAllComp <- function(res, digits = 2){
  data <- tibble(Compound = rownames(res),
                 p.value = res$p.value,
                 mean = res$`hazard ratio`,
                 lower = res$HR.confint.lower,
                 upper = res$HR.confint.upper,
                 N = res$N)
  
  # Number of digits shown 
  data$mean <- round(data$mean, digits = digits)
  data$lower <- round(data$lower, digits = digits)
  data$upper <- round(data$upper, digits = digits)
  
  for(i in 1:nrow(data)){
    if(data$p.value[i] < 0.01){
      data$p.value[i] <- format(as.numeric(data$p.value[i]), scientific = T)
    }else{
      data$p.value[i] <- round(as.numeric(data$p.value[i]), digits = digits)
    }
  }
  
  data$conf.int <- NA
  for(i in 1:nrow(data)){
    data$conf.int[i] <- paste0(data$lower[i], "-", data$upper[i])
  }
  
  data <- rbind(rep(NA,7), data)
  
  data$colour <- append(rep(c("white", "gray95"), nrow(data)/2), "white")
  p <- ggplot(data, aes(x = mean, y = Compound, xmin = lower, xmax = upper)) +
    geom_hline(aes(yintercept = Compound, colour = colour), size = 7) + 
    geom_pointrange(shape = 22, fill = "black") +
    geom_vline(xintercept = 1, linetype = 3) +
    xlab("Hazard ratio (95% CI)") +
    ylab("") +
    theme_classic() +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(data$Compound)) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  
  fplottable  <- data[-1,c(1,6,2)]
  fplottable <- rbind(c("Compound", "N", "P"),
                      fplottable)
  fplottable$`Hazard ratio (95% CI)` <- c("Hazard ratio",
                                          unlist(lapply(data$Compound[-1], function(x) paste0(data$mean[which(data$Compound == x)], " (", data$lower[which(data$Compound == x)], "-", data$upper[which(data$Compound == x)], ")"))))
  fplottable$color <- c("white", unlist(lapply(data$Compound[-1], function(x) data$colour[which(data$Compound == x)])))
  fplottable$Compound <- factor(fplottable$Compound, levels = fplottable$Compound)
  
  
  data_table <- ggplot(data = fplottable, aes(y = Compound)) +
    geom_hline(aes(yintercept = Compound, colour = color), size = 5) +
    geom_text(aes(x = 0, label = Compound), hjust = 0) +
    geom_text(aes(x = 4, label = p.value)) +
    geom_text(aes(x = 10, label = `Hazard ratio (95% CI)`), hjust = 1) +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(fplottable$Compound)) +
    theme_void() + 
    theme(plot.margin = margin(5, 0, 35, 0))                                       
  
  return(grid.arrange(as.grob(data_table), as.grob(p), ncol = 2))
}   

########################
#
# heatmapPlot - plot complex heatmap
#
########################
heatmapPlot <- function(fc, fcP, title, legendName = expression("log"[2]~"FC"), showLegend = FALSE, rowClust, minVal = -0.5, maxVal = 0.5, width_size = 8, annotate){
  result <- Heatmap(fc,
                    col = colorRamp2(c(minVal, 0, maxVal), c("blue", "white", "red")),
                    cell_fun = function(j, i, x, y, w, h, fill){
                      if(!is.na(fcP[i, j]) & fcP[i, j] <= 0.001){
                        grid.text("***", x, y, vjust = 0.7, gp = gpar(fontsize = 14))
                      }else if(!is.na(fcP[i, j]) & fcP[i, j] <= 0.01){
                        grid.text("**", x, y, vjust = 0.7,  gp = gpar(fontsize = 14))
                      }else if(!is.na(fcP[i, j]) & fcP[i, j] <= 0.05){
                        grid.text("*", x, y, vjust = 0.7,  gp = gpar(fontsize = 14))
                      }},
                    column_names_gp = grid::gpar(fontsize = 12),
                    row_names_gp = grid::gpar(fontsize = 10),
                    cluster_columns = FALSE,
                    cluster_rows = rowClust,
                    column_title = title,
                    heatmap_legend_param = list(title = legendName, title_gp = grid::gpar(fontface = "bold")),
                    width = width_size,
                    left_annotation = annotate,
                    show_heatmap_legend = showLegend)
  
  return(result)
}

#######################
#
# heatmapBySex <- Heatmap plot divided in males and females
#
########################
heatmapBySex <- function(fc.Male, fcP.Male, fc.Female, fcP.Female, row.names, filename, ha, width = 16){
  fc <- fc.Male
  fc <- apply(fc, 2, function(x) as.numeric(x))
  rownames(fc) <- row.names
  colnames(fc) <-  c("MS vs MC", "RRMS vs MC", "PMS vs MC")
  fcP <- fcP.Male
  fcP <- apply(fcP, 2, function(x) as.numeric(x))
  rownames(fcP) <- row.names
  colnames(fcP) <- c("MS vs MC", "RRMS vs MC", "PMS vs MC")
  
  hxp_MS_male <- heatmapPlot(fc = fc, fcP = fcP, title = "Males", rowClust = FALSE, maxVal = 1, minVal = -1,  showLegend = TRUE, annotate = ha)
  
  fc <- fc.Female
  fc <- apply(fc, 2, function(x) as.numeric(x))             
  rownames(fc) <- row.names
  colnames(fc) <-  c("MS vs MC", "RRMS vs MC", "PMS vs MC")
  fcP <- fcP.Female
  fcP <- apply(fcP, 2, function(x) as.numeric(x))
  rownames(fcP) <- row.names
  colnames(fcP) <- c("MS vs MC", "RRMS vs MC", "PMS vs MC")
  
  hxp_MS_female <- heatmapPlot(fc = fc, fcP = fcP, title = "Females", rowClust = FALSE, maxVal = 1, minVal = -1,  showLegend = FALSE, annotate = ha)
  
  png(file=paste0(filename, ".png"),width = 8, height = 8, units = "in", res = 300)
  draw(hxp_MS_male + hxp_MS_female)
  dev.off()
  
  pdf(file=paste0(filename, ".pdf"), width = width)
  draw(hxp_MS_male + hxp_MS_female)
  dev.off()
}