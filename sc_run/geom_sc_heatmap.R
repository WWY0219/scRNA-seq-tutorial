library(ggplot2)
library(grid)
library(dplyr)
library(Seurat)
library(ggdendro)
library(ggh4x)
library(reshape2)
library(tidyr)
library(tibble)
library(rlang)
ggname <- getFromNamespace("ggname","ggplot2")
 
fetch_data <- function(object = NULL,
                       reduction = "umap",
                       features = NULL,
                       featuresAnno = 0,
                       pct.exp.var = "seurat_clusters",
                       group.vars = NULL,
                       slot = c("data","counts")){
    slot <- match.arg(slot,c("data","counts"))
    options(warn=-1)
    options(dplyr.summarise.inform = FALSE)
    reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))
    colnames(reduc) <- paste0("Dim",1:2)
    meta <- object@meta.data
    ident <- data.frame(ident = Seurat::Idents(object))
    meta <- cbind(ident,meta)
    pc12 <- cbind(reduc, meta)
    pc12$cell <- rownames(pc12)
    if(!is.null(group.vars)){
        pc12 <- reshape2::melt(pc12,
                               id.vars = setdiff(x = colnames(pc12),y = group.vars),
                               variable.name = "group.vars",
                               value.name = "group.value")
    }
    if(!is.null(features)){
        geneExp <- Seurat::FetchData(object = object, vars = features,slot = slot)
        mer <- cbind(pc12, geneExp)
        megredf <- reshape2::melt(
            mer,
            id.vars = colnames(pc12),
            variable.name = "gene_name",
            value.name = "value"
        )
        ave_exp <- megredf %>% dplyr::group_by(.data[["gene_name"]],.data[[pct.exp.var]]) %>%
            dplyr::summarise(mean_exp = mean(value),median_exp = median(value)) %>%
            dplyr::ungroup()
        gene_pct <- megredf %>% dplyr::group_by(.data[["gene_name"]],.data[[pct.exp.var]]) %>%
            dplyr::summarise(pct = Seurat::PercentAbove(value,threshold = 0)*100) %>%
            dplyr::ungroup()
        fanno <- data.frame(gene_name = features,featureAnno = featuresAnno)
        megredf <- megredf %>%
            dplyr::left_join(y = gene_pct,by = c("gene_name",pct.exp.var)) %>%
            dplyr::left_join(y = ave_exp,by = c("gene_name",pct.exp.var)) %>%
            dplyr::left_join(y = fanno,by = c("gene_name"),
                             relationship = "many-to-many")
        gs <- table(features) > 1
        dup_genes <- names(gs)[gs]
        if(length(dup_genes) > 0){
            uni_f <- dup_genes
            lapply(seq_along(uni_f), function(x){
                tmp <- subset(megredf,gene_name %in% uni_f[x])
                anno_f <- unique(tmp$featureAnno)
                # check anno types length
                if(length(anno_f) > 1){
                    lapply(seq_along(anno_f), function(x){
                        tmp2 <- subset(tmp,featureAnno %in% anno_f[x])
                        tmp2$gene_name <- paste0(tmp2$gene_name,"_",x)
                        return(tmp2)
                    }) %>% Reduce("rbind",.) -> add_name
                    return(add_name)
                }else{
                    return(tmp)
                }
            }) %>% Reduce("rbind",.) %>% dplyr::arrange(.data[["featureAnno"]]) -> dup_df
            # rbind
            megredf <- rbind(dup_df,subset(megredf,gene_name %in% setdiff(features,dup_genes)))
        }
        # order
        od <- unique(megredf[,c("featureAnno","gene_name")])
        megredf$gene_name <- factor(megredf$gene_name,levels = od$gene_name)
        megredf <- megredf %>% dplyr::group_by(.data[["gene_name"]]) %>%
            dplyr::arrange(.data[["value"]])
    }else{
        megredf <- pc12
    }
    return(megredf)
}
ggscplot <- function(data = NULL,
                     mapping = aes(),
                     object = NULL,
                     reduction = "umap",
                     features = NULL,
                     featuresAnno = 0,
                     pct.exp.var = "seurat_clusters",
                     slot = "data",
                     gene.order = NULL,
                     group.vars = NULL,
                     t = 0.1,r = 0.1,b = 0.1,l = 0.1,
                     environment = parent.frame(),
                     ...) {
    if(length(mapping) == 0){
        mapping <- aes(x = Dim1,y = Dim2)
    }
    if(!is.null(object)){
        data <- fetch_data(object = object,reduction = reduction,
                           features = features,featuresAnno = featuresAnno,
                           pct.exp.var = pct.exp.var,
                           group.vars = group.vars,
                           slot = slot)
        # add levels
        if(!is.null(gene.order) & !is.null(features)){
            data$gene_name <- factor(data$gene_name,levels = gene.order)
        }
    }
    p <- ggplot(data = data,
                mapping = mapping,
                ...) +
        theme(plot.margin = margin(t = t,r = r,b = b,l = l,unit = "npc")) +
        coord_cartesian(clip = "off")
    class(p) <- c("ggscplot", class(p))
    return(p)
}
facet_new <- function(facet_col = NULL,
                      facet_row = NULL,
                      scales = "free_x",
                      strip.col = NULL,
                      space = "free_x",
                      no.xlabel = FALSE,
                      no.ylabel = FALSE,
                      title.hjust = 0,
                      t = 0.025,r = 0.1,b = 0.025,l = 0.1,
                      low.col = "grey90",high.col = "red",
                      x.angle = 0,x.label.hjust = 0.5,
                      legend.position = NULL,
                      ...){
    if(!is.null(facet_col)){
        facet_col <- vars(!!!rlang::ensyms(facet_col))
    }
    if(!is.null(facet_row)){
        facet_row <- vars(!!!rlang::ensyms(facet_row))
    }
    if(no.xlabel == TRUE){
        ele.tick.x <- element_blank()
        ele.text.x <- element_blank()
    }else{
        ele.tick.x <- element_line()
        ele.text.x <- element_text(angle = x.angle,hjust = ifelse(x.angle == 45,1,x.label.hjust))
    }
    if(no.ylabel == TRUE){
        ele.tick.y <- element_blank()
        ele.text.y <- element_blank()
    }else{
        ele.tick.y <- element_line()
        ele.text.y <- element_text()
    }
    list(ggh4x::facet_grid2(rows = facet_row,cols = facet_col,
                            scales = scales,
                            strip = ggh4x::strip_themed(background_x =
                                                            ggh4x::elem_list_rect(fill = strip.col),
                                                        background_y =
                                                            ggh4x::elem_list_rect(fill = strip.col),
                            ),
                            space = space,
                            ...),
         theme_bw() +
             theme(axis.ticks.x = ele.tick.x,
                   axis.text.x = ele.text.x,
                   axis.ticks.y = ele.tick.y,
                   axis.text.y = ele.text.y,
                   axis.text = element_text(colour = "black"),
                   strip.text = element_text(face = "bold.italic"),
                   strip.placement = "outside",
                   legend.position = legend.position,
                   legend.background = element_blank(),
                   plot.margin = margin(t = t,r = r,b = b,l = l,unit = "npc"),
                   panel.grid = element_blank()),
         scale_fill_gradient(low = low.col,high = high.col),
         guides(size = guide_legend(title = "Fraction of cells\nin group (%)",
                                    title.hjust = title.hjust),
                fill = guide_colorbar(title = "Mean expression\n in group",
                                      title.hjust = title.hjust,
                                      frame.colour = "black",
                                      ticks.colour = "black"))
    )
}
create_annosegment <- function(data = NULL,panel_params = NULL,
                               branch.side = c("top","right","bottom","left"),
                               branch.shift = 0,
                               branch.len = 0.9,
                               branch.height = 0.05,
                               branch.lwd = 0.5,
                               branch.label.rot = NULL,
                               branch.label.hjust = NULL,
                               branch.label.size = 10,
                               branch.label.shift = 0.2){
  branch.side <- match.arg(branch.side,c("top","right","bottom","left"))
  # generate celltype anno grobs
  cell_type <- unique(data$celltype)
  lapply(seq_along(cell_type), function(ii){
    if(branch.side %in% c("top","bottom")){
      tmp <- subset(data,celltype %in% cell_type[ii],select = c(x,celltype)) %>% unique()
      seg <- data.frame(x = min(tmp$x) - branch.len*0.5,
                        xend = max(tmp$x) + branch.len*0.5,
                        celltype = cell_type[ii])
    }else{
      tmp <- subset(data,celltype %in% cell_type[ii],select = c(y,celltype)) %>% unique()
      seg <- data.frame(x = min(tmp$y) - branch.len*0.5,
                        xend = max(tmp$y) + branch.len*0.5,
                        celltype = cell_type[ii])
    }
    return(seg)
  }) %>% Reduce("rbind",.) -> seg_df
  # ============================================================================
  # generate segments grob
  # check position
  if(branch.side == "top"){
    vp.width = 1
    vp.height = branch.height
    xscale = panel_params$x.range
    yscale = c(0,1)
    vp.x = 0.5
    vp.y = 1 + branch.shift
    vp.just = "bottom"
    x0 = seg_df$x
    x1 = seg_df$xend
    y0 = 1; y1 = 1
    vx0 = c(seg_df$x,seg_df$xend)
    vx1 = c(seg_df$x,seg_df$xend)
    vy0 = 0; vy1 = 1
    labelx = (x0 + x1)*0.5
    labely = vp.y + branch.label.shift
    label.rot = 90
    label.hjust = 0
  }else if(branch.side == "right"){
    vp.width = branch.height
    vp.height = 1
    xscale = c(0,1)
    yscale = panel_params$y.range
    vp.x = 1 + branch.shift
    vp.y = 0.5
    vp.just = "left"
    x0 = 1; x1 = 1
    y0 = seg_df$x
    y1 = seg_df$xend
    vx0 = 0; vx1 = 1
    vy0 = c(seg_df$x,seg_df$xend)
    vy1 = c(seg_df$x,seg_df$xend)
    labelx = vp.x + branch.label.shift
    labely = (y0 + y1)*0.5
    label.rot = 0
    label.hjust = 0
  }else if(branch.side == "bottom"){
    vp.width = 1
    vp.height = branch.height
    xscale = panel_params$x.range
    yscale = c(0,1)
    vp.x = 0.5
    vp.y = 0 - branch.shift
    vp.just = "top"
    x0 = seg_df$x
    x1 = seg_df$xend
    y0 = 0; y1 = 0
    vx0 = c(seg_df$x,seg_df$xend)
    vx1 = c(seg_df$x,seg_df$xend)
    vy0 = 0; vy1 = 1
    labelx = (x0 + x1)*0.5
    labely = vp.y - branch.label.shift
    label.rot = 90
    label.hjust = 1
  }else{
    vp.width = branch.height
    vp.height = 1
    xscale = c(0,1)
    yscale = panel_params$y.range
    vp.x = 0 - branch.shift
    vp.y = 0.5
    vp.just = "right"
    x0 = 0
    x1 = 0
    y0 = seg_df$x
    y1 = seg_df$xend
    vx0 = 0; vx1 = 1
    vy0 = c(seg_df$x,seg_df$xend)
    vy1 = c(seg_df$x,seg_df$xend)
    labelx = vp.x - branch.label.shift
    labely = (y0 + y1)*0.5
    label.rot = 0
    label.hjust = 1
  }
  vp <- viewport(width = vp.width,height = vp.height,
                 xscale = xscale,
                 yscale = yscale,
                 x = vp.x,y = vp.y,just = vp.just)
  h_seg <- segmentsGrob(x0 = x0,x1 = x1,
                        y0 = y0,y1 = y1,
                        gp = gpar(lwd = branch.lwd),
                        default.units = "native",
                        vp = vp)
  v_seg <- segmentsGrob(x0 = vx0,
                        x1 = vx1,
                        y0 = vy0,y1 = vy1,
                        gp = gpar(lwd = branch.lwd),
                        default.units = "native",
                        vp = vp)
  if(!is.null(branch.label.rot)) label.rot <- branch.label.rot
  if(!is.null(label.hjust)) label.hjust <- branch.label.hjust
  label_grob <- textGrob(label = seg_df$celltype,
                         x = labelx,
                         y = labely,
                         default.units = "native",
                         gp = gpar(fontsize = branch.label.size,
                                   fontface = "bold"),
                         rot = label.rot,hjust = label.hjust,
                         vp = vp)
  # ============================================================================
  grobs <- list(h_seg = h_seg,
                v_seg = v_seg,
                label_grob = label_grob)
}
create_dendrogrob <- function(data = NULL,
                              panel_params = NULL,
                              exp_mat = NULL,
                              tree.type = c("rectangle","triangle"),
                              add.tree.y = FALSE,
                              tree.y.side = "left",
                              tree.y.width = 0.06,
                              add.tree.x = FALSE,
                              tree.x.side = "bottom",
                              tree.x.height = 0.06,
                              tree.y.label.hjust = 1,
                              tree.x.label.hjust = 1,
                              tree.x.shift = 0,
                              tree.y.shift = 0,
                              new.ylabel.x = -0.05,
                              new.ylabel.width = 0.025,
                              new.ylabel.size = 8,
                              new.ylabel.face = "plain",
                              new.ylabel.rot = 0,
                              new.xlabel.y = -0.05,
                              new.xlabel.height = 0.025,
                              new.xlabel.size = 8,
                              new.xlabel.face = "plain",
                              new.xlabel.rot = 90){
  tree.type <- match.arg(tree.type,c("rectangle","triangle"))
  # long to matrix
  if(is.null(exp_mat)){
    exp_mat <- data %>%
      dplyr::select(x,y,exp) %>%
      tidyr::spread(key = x,value = exp) %>%
      tibble::column_to_rownames(var = "y")
  }else{
    exp_mat <- exp_mat
  }
  # =====================================================
  # add y axis tree
  if(add.tree.y == TRUE){
    # hclust for row
    dendy = as.dendrogram(hclust(dist(exp_mat)))
    # get tree data
    dend_dfy <- ggdendro::dendro_data(model = dendy,type = tree.type)$segments
    dend_dfy$y1 <- max(dend_dfy$y,dend_dfy$yend) - dend_dfy$yend
    dend_dfy$yend1 <- max(dend_dfy$y,dend_dfy$yend) - dend_dfy$y
    dend_dfy$x1 <- max(dend_dfy$x,dend_dfy$xend) - dend_dfy$xend + 1
    dend_dfy$xend1 <- max(dend_dfy$x,dend_dfy$xend) - dend_dfy$x + 1
    # check side
    if(tree.y.side == "right"){
      tree.y.x0 = dend_dfy$y
      tree.y.x1 = dend_dfy$yend
      tree.y.y0 = dend_dfy$x1
      tree.y.y1 = dend_dfy$xend1
      tree.y.just = "left"
      tree.y.x = 1 + tree.x.shift
      tree.y.label.hjust = 1
    }else{
      tree.y.x0 = dend_dfy$y1
      tree.y.x1 = dend_dfy$yend1
      tree.y.y0 = dend_dfy$x1
      tree.y.y1 = dend_dfy$xend1
      tree.y.just = "right"
      tree.y.x = 0 - tree.x.shift
      tree.y.label.hjust = 0
    }
    tree.groby <- segmentsGrob(x0 = tree.y.x0,x1 = tree.y.x1,
                               y0 = tree.y.y0,y1 = tree.y.y1,
                               default.units = "native",
                               vp = viewport(width = tree.y.width,height = 1,
                                             xscale = range(dend_dfy$y,dend_dfy$yend),
                                             yscale = panel_params$y.range,
                                             x = tree.y.x,y = 0.5,just = tree.y.just))
    # reorder data
    order_y <- order.dendrogram(dendy)
    yl <- unique(data$y)
    lapply(seq_along(yl), function(ii){
      tmp <- subset(data,y %in% yl[ii])
      idx <- match(yl[ii],order_y)
      tmp$y <- idx
      return(tmp)
    }) %>% Reduce("rbind",.) -> data
    # re-order coord y limit according to hclust order
    raw_limits <- panel_params$y$limits
    new_limits <- raw_limits[order.dendrogram(dendy)]
    yaxis_grob <- textGrob(label = new_limits,
                           x = 0,y = 1:length(new_limits),
                           gp = gpar(fontsize = new.ylabel.size,
                                     fontface = new.ylabel.face),
                           rot = new.ylabel.rot,
                           hjust = tree.y.label.hjust,
                           default.units = "native",
                           vp = viewport(width = 0.025,height = 1,
                                         xscale = c(0,1),
                                         yscale = panel_params$y.range,
                                         x = new.ylabel.x,y = 0.5))
  }else{
    tree.groby <- zeroGrob()
    yaxis_grob <- zeroGrob()
  }
  # =============================================================================
  # add x axis tree
  if(add.tree.x == TRUE){
    # hclust for row
    dendx = as.dendrogram(hclust(dist(t(exp_mat))))
    # get tree data
    dend_dfx <- ggdendro::dendro_data(model = dendx,type = tree.type)$segments
    dend_dfx$y1 <- max(dend_dfx$y,dend_dfx$yend) - dend_dfx$yend
    dend_dfx$yend1 <- max(dend_dfx$y,dend_dfx$yend) - dend_dfx$y
    dend_dfx$x1 <- max(dend_dfx$x,dend_dfx$xend) - dend_dfx$xend + 1
    dend_dfx$xend1 <- max(dend_dfx$x,dend_dfx$xend) - dend_dfx$x + 1
    # check side
    if(tree.x.side == "top"){
      tree.x.x0 = dend_dfx$x1
      tree.x.x1 = dend_dfx$xend1
      tree.x.y0 = dend_dfx$y
      tree.x.y1 = dend_dfx$yend
      tree.x.just = "bottom"
      tree.x.y = 1 + tree.y.shift
      tree.x.label.hjust = 1
    }else{
      tree.x.x0 = dend_dfx$x
      tree.x.x1 = dend_dfx$xend
      tree.x.y0 = dend_dfx$y1
      tree.x.y1 = dend_dfx$yend1
      tree.x.just = "top"
      tree.x.y = 0 - tree.y.shift
      tree.x.label.hjust = 0
    }
    tree.grobx <- segmentsGrob(x0 = tree.x.x0,x1 = tree.x.x1,
                               y0 = tree.x.y0,y1 = tree.x.y1,
                               default.units = "native",
                               vp = viewport(width = 1,height = tree.x.height,
                                             xscale = panel_params$x.range,
                                             yscale = range(dend_dfx$y,dend_dfx$yend),
                                             x = 0.5,y = tree.x.y,just = tree.x.just))
    # reorder data
    order_x <- order.dendrogram(dendx)
    xl <- unique(data$x)
    lapply(seq_along(xl), function(ii){
      tmp <- subset(data,x %in% xl[ii])
      idx <- match(xl[ii],order_x)
      tmp$x <- idx
      return(tmp)
    }) %>% Reduce("rbind",.) -> data
    # re-order coord y limit according to hclust order
    raw_limits <- panel_params$x$limits
    new_limits <- raw_limits[order.dendrogram(dendx)]
    xaxis_grob <- textGrob(label = new_limits,
                           x = 1:length(new_limits),y = 0,
                           gp = gpar(fontsize = new.xlabel.size,
                                     fontface = new.xlabel.face),
                           rot = new.xlabel.rot,
                           hjust = tree.x.label.hjust,
                           default.units = "native",
                           vp = viewport(width = 1,height = new.xlabel.height,
                                         xscale = panel_params$x.range,
                                         yscale = c(0,1),
                                         x = 0.5,y = new.xlabel.y))
  }else{
    tree.grobx <- zeroGrob()
    xaxis_grob <- zeroGrob()
  }
  # ============================================================================
  grobs <- list(data = data,
                tree.groby = tree.groby,
                yaxis_grob = yaxis_grob,
                tree.grobx = tree.grobx,
                xaxis_grob = xaxis_grob)
}
calc_density <- function(data = NULL,
                         trim = TRUE,
                         bw = "nrd0",
                         adjust = 1,
                         kernel = "gaussian",
                         n = 512){
  # calculate density
  density_data <- stats::density(data$exp,
                                 bw = bw, adjust = adjust,
                                 kernel = kernel,
                                 n = n)
 
  # data range
  range_data <- range(data$exp)
 
  # whether tirm head and tail
  if(trim == TRUE){
    # data frame
    new_daframe <- data.frame(vio_y = density_data$x,vio_x = density_data$y)
 
    # trim head and tail
    new_daframe$vio_y <- dplyr::case_when(
      new_daframe$vio_y < range_data[1] ~ range_data[1],
      new_daframe$vio_y > range_data[2] ~ range_data[2],
      TRUE ~ new_daframe$vio_y)
 
  }else{
    # data frame
    new_daframe <- data.frame(vio_y = density_data$x,vio_x = density_data$y)
  }
 
  # add x y
  new_daframe$x <- data$x[1]
  new_daframe$y <- data$y[1]
 
  # add median exp
  new_daframe$meadian_exp <- median(data$exp)
  new_daframe
}
 
geom_sc_heatmap <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        linejoin = "mitre",
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        add.tree = FALSE,
                        tree.type = "rectangle",
                        add.tree.y = FALSE,
                        tree.y.side = "left",
                        tree.y.width = 0.06,
                        add.tree.x = FALSE,
                        tree.x.side = "bottom",
                        tree.x.height = 0.06,
                        tree.y.label.hjust = 1,
                        tree.x.label.hjust = 1,
                        new.ylabel.x = -0.08,
                        new.ylabel.width = 0.025,
                        new.ylabel.size = 8,
                        new.ylabel.face = "plain",
                        new.ylabel.rot = 0,
                        new.xlabel.y = -0.08,
                        new.xlabel.height = 0.025,
                        new.xlabel.size = 8,
                        new.xlabel.face = "plain",
                        new.xlabel.rot = 90,
                        tree.x.shift = 0,
                        tree.y.shift = 0,
                        # =======================
                        branch.side = "top",
                        branch.shift = 0,
                        branch.height = 0.05,
                        branch.len = 0.8,
                        branch.lwd = 0.5,
                        branch.label.rot = NULL,
                        branch.label.hjust = NULL,
                        branch.label.size = 10,
                        branch.label.shift = 0.2) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSCtile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(
      linejoin = linejoin,
      na.rm = na.rm,
      add.tree = add.tree,
      tree.type = tree.type,
      add.tree.y = add.tree.y,
      tree.y.side = tree.y.side,
      tree.y.width = tree.y.width,
      add.tree.x = add.tree.x,
      tree.x.side = tree.x.side,
      tree.x.height = tree.x.height,
      tree.y.label.hjust = tree.y.label.hjust,
      tree.x.label.hjust = tree.x.label.hjust,
      new.ylabel.x = new.ylabel.x,
      new.ylabel.width = new.ylabel.width,
      new.ylabel.size = new.ylabel.size,
      new.ylabel.face = new.ylabel.face,
      new.ylabel.rot = new.ylabel.rot,
      new.xlabel.y = new.xlabel.y,
      new.xlabel.height = new.xlabel.height,
      new.xlabel.size = new.xlabel.size,
      new.xlabel.face = new.xlabel.face,
      new.xlabel.rot = new.xlabel.rot,
      tree.x.shift = tree.x.shift,
      tree.y.shift = tree.y.shift,
      # =======================
      branch.side = branch.side,
      branch.shift = branch.shift,
      branch.height = branch.height,
      branch.len = branch.len,
      branch.lwd = branch.lwd,
      branch.label.rot = branch.label.rot,
      branch.label.hjust = branch.label.hjust,
      branch.label.size = branch.label.size,
      branch.label.shift = branch.label.shift,
      ...
    )
  )
}
 
GeomSCtile <- ggproto("GeomSCtile", ggplot2::GeomRect,
                      required_aes = c("x", "y"),
 
                      extra_params = c("na.rm"),
 
                      default_aes = aes(colour = NA, fill = "grey35", linewidth = 0.5, linetype = 1,
                                        alpha = NA,
                                        cluster = NULL,exp = NULL,celltype = NULL),
 
                      # These aes columns are created by setup_data(). They need to be listed here so
                      # that GeomRect$handle_na() properly removes any bars that fall outside the defined
                      # limits, not just those for which x and y are outside the limits
                      non_missing_aes = c("xmin", "xmax", "ymin", "ymax"),
 
                      setup_data = function(data, params) {
                        data$width <- data$width %||% params$width %||% resolution(data$x, FALSE)
                        data$height <- data$height %||% params$height %||% resolution(data$y, FALSE)
 
                        transform(data,
                                  xmin = x - width / 2,  xmax = x + width / 2,  width = NULL,
                                  ymin = y - height / 2, ymax = y + height / 2, height = NULL
                        )
                      },
 
                      draw_panel = function(self, data, panel_params, coord,
                                            lineend = "butt", linejoin = "mitre",
                                            add.tree = FALSE,
                                            tree.type = "rectangle",
                                            add.tree.y = FALSE,
                                            tree.y.side = "left",
                                            tree.y.width = 0.06,
                                            add.tree.x = FALSE,
                                            tree.x.side = "bottom",
                                            tree.x.height = 0.06,
                                            tree.y.label.hjust = 1,
                                            tree.x.label.hjust = 1,
                                            new.ylabel.x = -0,
                                            new.ylabel.width = 0.025,
                                            new.ylabel.size = 8,
                                            new.ylabel.face = "plain",
                                            new.ylabel.rot = 0,
                                            new.xlabel.y = -0,
                                            new.xlabel.height = 0.025,
                                            new.xlabel.size = 8,
                                            new.xlabel.face = "plain",
                                            new.xlabel.rot = 90,
                                            tree.x.shift = 0,
                                            tree.y.shift = 0,
                                            # =======================
                                            branch.side = "top",
                                            branch.shift = 0,
                                            branch.height = 0.05,
                                            branch.len = 0.8,
                                            branch.lwd = 0.5,
                                            branch.label.rot = NULL,
                                            branch.label.hjust = NULL,
                                            branch.label.size = 10,
                                            branch.label.shift = 0.2){
                        # =====================================================
                        # remove dulicate expressions
                        data <- unique(data)
 
                        # =====================================================
                        # add dendrogram to row or col
                        if(add.tree == TRUE){
                          grobs <- create_dendrogrob(data = data,
                                                     panel_params = panel_params,
                                                     tree.type = tree.type,
                                                     add.tree.y = add.tree.y,
                                                     tree.y.side = tree.y.side,
                                                     tree.y.width = tree.y.width,
                                                     add.tree.x = add.tree.x,
                                                     tree.x.side = tree.x.side,
                                                     tree.x.height = tree.x.height,
                                                     tree.y.label.hjust = tree.y.label.hjust,
                                                     tree.x.label.hjust = tree.x.label.hjust,
                                                     tree.x.shift = tree.x.shift,
                                                     tree.y.shift = tree.y.shift,
                                                     new.ylabel.x = new.ylabel.x,
                                                     new.ylabel.width = new.ylabel.width,
                                                     new.ylabel.size = new.ylabel.size,
                                                     new.ylabel.face = new.ylabel.face,
                                                     new.ylabel.rot = new.ylabel.rot,
                                                     new.xlabel.y = new.xlabel.y,
                                                     new.xlabel.height = new.xlabel.height,
                                                     new.xlabel.size = new.xlabel.size,
                                                     new.xlabel.face = new.xlabel.face,
                                                     new.xlabel.rot = new.xlabel.rot)
 
                          tree.groby <- grobs$tree.groby
                          yaxis_grob <- grobs$yaxis_grob
                          tree.grobx <- grobs$tree.grobx
                          xaxis_grob <- grobs$xaxis_grob
 
                          data <- grobs$data
                        }else{
                          tree.groby <- zeroGrob()
                          yaxis_grob <- zeroGrob()
                          tree.grobx <- zeroGrob()
                          xaxis_grob <- zeroGrob()
                        }
 
                        # =======================================================
                        # annotation for genes
                        if("celltype" %in% colnames(data)){
                          # generate celltype anno grobs
                          seg_grobs <- create_annosegment(data = data,
                                                          panel_params = panel_params,
                                                          branch.side = branch.side,
                                                          branch.shift = branch.shift,
                                                          branch.height = branch.height,
                                                          branch.len = branch.len,
                                                          branch.lwd = branch.lwd,
                                                          branch.label.rot = branch.label.rot,
                                                          branch.label.hjust = branch.label.hjust,
                                                          branch.label.size = branch.label.size,
                                                          branch.label.shift = branch.label.shift)
 
                          h_seg <- seg_grobs$h_seg
                          v_seg <- seg_grobs$v_seg
                          label_grob <- seg_grobs$label_grob
                        }else{
                          h_seg <- zeroGrob()
                          v_seg <- zeroGrob()
                          label_grob <- zeroGrob()
                        }
 
                        # =========================================================================
                        coords <- coord$transform(data, panel_params)
 
                        rect_grob <- rectGrob(coords$xmin, coords$ymax,
                                              width = coords$xmax - coords$xmin,
                                              height = coords$ymax - coords$ymin,
                                              default.units = "native",
                                              just = c("left", "top"),
                                              gp = gpar(
                                                col = coords$colour,
                                                fill = alpha(coords$fill, coords$alpha),
                                                lwd = coords$linewidth * .pt,
                                                lty = coords$linetype,
                                                linejoin = linejoin,
                                                lineend = lineend))
 
                        ggname("geom_scTile",grid::gTree(children = gList(rect_grob,
                                                                          h_seg,v_seg,label_grob,
                                                                          tree.groby,yaxis_grob,
                                                                          tree.grobx,xaxis_grob)))
 
                      },
 
                      draw_key = draw_key_polygon
)
