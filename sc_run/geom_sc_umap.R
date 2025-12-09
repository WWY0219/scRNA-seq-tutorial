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
library(lifecycle)
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
geom_sc_umap <- function(mapping = NULL, data = NULL,
                          stat = "identity", position = "identity",
                          ...,
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          label.gp = gpar(),add_label = TRUE,
                          add_legend = FALSE,
                          ncol = 1,vgap = 0.5,hgap = 1,
                          point.size = 2,point.label.size = 10,
                          point.label.col = "black",
                          lgd_x = 1.05,lgd_y = 0.5,
                          lgd_width = 0.1,lgd_height = 0.9,
                          fontsize = 10,fontface = "plain") {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSCpoint2,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(
      na.rm = na.rm,
      label.gp = label.gp,add_label = add_label,
      add_legend = add_legend,
      ncol = ncol,vgap = vgap,hgap = hgap,
      point.size = point.size,point.label.size = point.label.size,
      point.label.col = point.label.col,
      lgd_x = lgd_x,lgd_y = lgd_y,
      lgd_width = lgd_width,lgd_height = lgd_height,
      fontsize = fontsize,fontface = fontface,
      ...
    )
  )
}
 
#' ggproto for GeomSCpoint2
#' @format NULL
#' @usage NULL
#' @export
GeomSCpoint2 <- ggproto("GeomSCpoint2", Geom,
                        required_aes = c("x", "y"),
                        non_missing_aes = c("size", "shape", "colour"),
                        default_aes = aes(
                          shape = 19, colour = "black", size = 1, fill = "black",
                          alpha = NA, stroke = 0.5,cluster = NULL,cluster_anno = NULL
                        ),
 
                        draw_panel = function(self, data, panel_params, coord, na.rm = FALSE,
                                              label.gp = gpar(),add_label = TRUE,
                                              add_legend = FALSE,
                                              ncol = 1,vgap = 0.5,hgap = 1,
                                              point.size = 2,point.label.size = 10,
                                              point.label.col = "black",
                                              lgd_x = 1.05,lgd_y = 0.5,
                                              lgd_width = 0.1,lgd_height = 0.9,
                                              fontsize = 10,fontface = "plain") {
 
                          # ==================================================================
                          if(add_legend == TRUE){
                            if(!("cluster_anno" %in% colnames(data) & "cluster_anno" %in% colnames(data))){
                              message("Please supply cluster_anno mapping variable.")
                              stop()
                            }
 
                            ld_df <- subset(data,select = c(cluster,cluster_anno,fill,colour,alpha)) %>%
                              unique() %>% dplyr::arrange(cluster_anno)
 
                            vp = viewport(x = lgd_x,y = lgd_y,just = "left",
                                          height = lgd_height,width = lgd_width)
 
                            legend_grob <- legendGrob2(labels = ld_df$cluster,
                                                       labels.point = ld_df$cluster_anno,
                                                       labels.point.gp = gpar(fontface = "bold",
                                                                              col = point.label.col,
                                                                              fontsize = point.label.size),
                                                       ncol = ncol,
                                                       pch = unique(data$shape),
                                                       vgap = unit(vgap, "lines"),
                                                       hgap = unit(hgap, "lines"),
                                                       gp = gpar(fontsize = fontsize,
                                                                 fontface = fontface,
                                                                 col = alpha(ld_df$colour, ld_df$alpha),
                                                                 fill = alpha(ld_df$fill, ld_df$alpha)),
                                                       point.size = point.size,
                                                       vp = vp)
 
                            lgd <- legend_grob
                          }else{
                            lgd <- zeroGrob()
                          }
 
                          # ==================================================================
 
                          if (is.character(data$shape)) {
                            data$shape <- translate_shape_string(data$shape)
                          }
 
                          coords <- coord$transform(data, panel_params)
 
                          # calculate cluster centers
                          if("cluster" %in% colnames(coords)){
                            centers <- coords %>%
                              dplyr::group_by(cluster) %>%
                              dplyr::summarise(x = median(x = x), y = median(x = y))
 
                            label_grob <- textGrob(x = centers$x,y = centers$y,
                                                   label = centers$cluster,
                                                   gp = label.gp,
                                                   default.units = "npc")
 
                            if(add_label == FALSE){
                              label_grob <- zeroGrob()
                            }
                          }else{
                            label_grob <- zeroGrob()
                          }
 
                          # =======================================================
                          # point grob
                          stroke_size <- coords$stroke
                          stroke_size[is.na(stroke_size)] <- 0
 
                          point_grob <- pointsGrob(
                            coords$x, coords$y,
                            pch = coords$shape,
                            default.units = "npc",
                            gp = gpar(
                              col = alpha(coords$colour, coords$alpha),
                              fill = alpha(coords$fill, coords$alpha),
                              # Stroke is added around the outside of the point
                              fontsize = coords$size * .pt + stroke_size * .stroke / 2,
                              lwd = coords$stroke * .stroke / 2
                            )
                          )
 
                          ggname("geom_scPoint2",
                                 grid::gTree(children = gList(point_grob,label_grob,
                                                              lgd)))
                        },
 
                        draw_key = draw_key_point
)
theme_sc <- function(x.line.len = 0.25,x.label = "UMAP 1",
                     y.line.len = 0.25,y.label = "UMAP 2",
                     label.shift = 0.5,fontsize = 10,
                     border.col = NA,add.arrow = TRUE,
                     t = 0.1,r = 0.1,b = 0.1,l = 0.1,
                     ...) {
  if(add.arrow == TRUE){
    axis.line.x.bottom <- element_line2(x.line.ed = x.line.len,y.line.st = 1,
                                        label = x.label,label.shift = label.shift,
                                        fontsize = fontsize,
                                        arrow = arrow(length = unit(0.3,"cm"),
                                                      type = "closed",
                                                      ends = "last"))
 
    axis.line.y.left <- element_line2(y.line.ed = y.line.len,x.line.st = 1,
                                      label = y.label,label.shift = label.shift,
                                      fontsize = fontsize,
                                      arrow = arrow(length = unit(0.3,"cm"),
                                                    type = "closed",
                                                    ends = "last"))
  }else{
    axis.line.x.bottom <- element_blank()
    axis.line.y.left <- element_blank()
  }
 
  list(theme_classic(...) %+replace%
         theme(
           # axis.line = element_blank(),
           axis.line.x.bottom = axis.line.x.bottom,
           axis.line.y.left = axis.line.y.left,
 
           # white background and dark border
           panel.background = element_rect(fill = "white", colour = NA),
           panel.border     = element_rect(fill = NA, colour = border.col),
           # contour strips to match panel contour
           strip.background = element_rect(fill = "grey90", colour = "black"),
           strip.text = element_text(face = "bold.italic"),
 
           axis.ticks = element_blank(),
           axis.text = element_blank(),
           axis.title = element_blank(),
           panel.grid = element_blank(),
           plot.margin = margin(t = t,r = r,b = b,l = l,unit = "npc"),
           panel.spacing = unit(0.5,"cm"),
 
           complete = F
         )
       )
 
}
element_line2 <- function(colour = NULL, linewidth = NULL, linetype = NULL,
                          lineend = NULL, color = NULL, arrow = NULL,
                          inherit.blank = FALSE, size = lifecycle::deprecated(),
                          x.line.st = 0,x.line.ed = 1,
                          y.line.st = 0,y.line.ed = 1,
                          label.shift = -0.1,label = "dim",
                          fontface = "bold.italic",fontsize = 8) {
 
  if (lifecycle::is_present(size)) {
    deprecate_soft0("3.4.0", "element_line(size)", "element_line(linewidth)")
    linewidth <- size
  }
 
  if (!is.null(color))  colour <- color
  if (is.null(arrow)) arrow <- FALSE
  structure(
    list(colour = colour, linewidth = linewidth, linetype = linetype, lineend = lineend,
         arrow = arrow, inherit.blank = inherit.blank,
         x.line.st = x.line.st,x.line.ed = x.line.ed,
         y.line.st = y.line.st,y.line.ed = y.line.ed,
         label.shift = label.shift,label = label,
         fontface = fontface,fontsize = fontsize
    ),
    class = c("element_line2","element_line", "element")
  )
}
 
 
 
 
#' @importFrom ggplot2 element_grob
#' @method element_grob element_line2
#' @export
element_grob.element_line2 <- function(element,...,
                                       x = 0:1, y = 0:1,
                                       colour = NULL, linewidth = NULL, linetype = NULL, lineend = NULL,
                                       default.units = "npc", id.lengths = NULL, size = deprecated()) {
 
  if (lifecycle::is_present(size)) {
    deprecate_soft0("3.4.0", "element_grob.element_line(size)", "element_grob.element_line(linewidth)")
    linewidth <- size
  }
 
  # The gp settings can override element_gp
  gp <- gpar(
    col = colour, fill = colour,
    lwd = linewidth, lty = linetype, lineend = lineend
  )
  element_gp <- gpar(
    col = element$colour, fill = element$colour,
    lwd = linewidth, lty = element$linetype,
    lineend = element$lineend
  )
  arrow <- if (is.logical(element$arrow) && !element$arrow) {
    NULL
  } else {
    element$arrow
  }
 
  line.grob <- polylineGrob(
    x = c(element$x.line.st,element$x.line.ed),
    y = c(element$y.line.st,element$y.line.ed),
    # x = x,y = y,
    default.units = default.units,
    gp = ggplot2:::modify_list(element_gp, gp),
    id.lengths = id.lengths, arrow = arrow
  )
 
  # ============================================================================
  # labels
 
 
  if(TRUE){
  if(element$y.line.st == element$y.line.ed){
    label.x = (element$x.line.st + element$x.line.ed)*0.5
    label.y = element$label.shift
    rot = 0
 
    hjust=0.5
    vjust=element$label.shift+1
  }else{
    label.y = (element$y.line.st + element$y.line.ed)*0.5
    label.x = element$label.shift
    rot = 90
 
    hjust=0.5
    vjust=-element$label.shift
  }
}
 
  #print(c(element$label,label.x,label.y,hjust,vjust))
  label.grob <- textGrob(label = element$label,
                         x = label.x,y = label.y,
                         hjust = hjust,vjust=vjust,
                         rot = rot,
                         gp = gpar(fontface = element$fontface,fontsize = element$fontsize),
                         default.units = default.units)
 
  gTree(children = gList(line.grob,label.grob))
}
legendGrob2 <- function(labels = NULL,
                        labels.point = NULL,
                        labels.point.gp = gpar(),
                        nrow = NULL, ncol = NULL, byrow = FALSE,
                        do.lines = has.lty || has.lwd, lines.first = TRUE,
                        hgap = unit(1, "lines"), vgap = unit(1, "lines"),
                        default.units = "lines",
                        point.size = 2,
                        pch = NULL, gp = gpar(), vp = NULL){
  # ============================================================================
  ## Type checking on arguments; labels: character, symbol or expression:
  labels <- as.graphicsAnnot(labels)
  labels <- if(is.character(labels)) as.list(labels) else as.expression(labels)
  nkeys <- if(is.call(labels)) 1 else length(labels)
  if(nkeys == 0) return(nullGrob(vp=vp))
  if (!is.unit(hgap))
    hgap <- unit(hgap, default.units)
  if (length(hgap) != 1) stop("'hgap' must be single unit")
  if (!is.unit(vgap))
    vgap <- unit(vgap, default.units)
  if (length(vgap) != 1) stop("'vgap' must be single unit")
 
  # ============================================================================
  ## nrow, ncol
  miss.nrow <- missing(nrow)
  miss.ncol <- missing(ncol)
  if(miss.nrow && miss.ncol) {ncol <- 1; nrow <- nkeys} # defaults to 1-column legend
  else if( miss.nrow && !miss.ncol) nrow <- ceiling(nkeys / ncol)
  else if(!miss.nrow &&  miss.ncol) ncol <- ceiling(nkeys / nrow)
  if(nrow < 1) stop("'nrow' must be >= 1")
  if(ncol < 1) stop("'ncol' must be >= 1")
  if(nrow * ncol < nkeys)
    stop("nrow * ncol < #{legend labels}")
 
  # ============================================================================
  ## pch, gp
  if(has.pch <- !missing(pch) && length(pch) > 0) pch <- rep_len(pch, nkeys)
  if(doGP <- length(nmgp <- names(gp)) > 0) {
    if(has.lty  <-  "lty" %in% nmgp) gp$lty  <- rep_len(gp$lty, nkeys)
    if(has.lwd  <-  "lwd" %in% nmgp) gp$lwd  <- rep_len(gp$lwd, nkeys)
    if(has.col  <-  "col" %in% nmgp) gp$col  <- rep_len(gp$col,  nkeys)
    if(has.fill <- "fill" %in% nmgp) gp$fill <- rep_len(gp$fill, nkeys)
  } else {
    gpi <- gp
    if(missing(do.lines)) do.lines <- FALSE
  }
 
  # ============================================================================
  ## main
  u0 <- unit(0, "npc")
  u1 <- unit(1, "char")
  ord <- if(lines.first) 1:2 else 2:1
  fg <- frameGrob(vp = vp)    # set up basic frame grob (for packing)
  for (i in seq_len(nkeys)) {
    if(doGP) {
      gpi <- gp
      if(has.lty)    gpi$lty <- gp$lty[i]
      if(has.lwd)    gpi$lwd <- gp$lwd[i]
      if(has.col)    gpi$col <- gp$col[i]
      if(has.fill) gpi$fill<- gp$fill[i]
    }
    if(byrow) {
      ci <- 1+ (i-1) %%  ncol
      ri <- 1+ (i-1) %/% ncol
    } else {
      ci <- 1+ (i-1) %/% nrow
      ri <- 1+ (i-1) %%  nrow
    }
 
    # ============================================================================
    ## borders; unit.c creates a 4-vector of borders (bottom, left, top, right)
    vg <- if(ri != nrow) vgap else u0
    symbol.border <- unit.c(vg, u0, u0, 0.5 * hgap)
    text.border   <- unit.c(vg, u0, u0, if(ci != ncol) hgap else u0)
 
    # ============================================================================
    ## points/lines grob:
    plGrob <- if(has.pch && do.lines){
      gTree(children = gList(linesGrob(0:1, 0.5, gp=gpi),
                             pointsGrob(0.5, 0.5, default.units="npc",
                                        size = unit(point.size, "char"),
                                        pch=pch[i], gp=gpi),
                             textGrob(label = labels.point[[i]],
                                      x = 0.5,y = 0.5,
                                      default.units="npc",
                                      gp = labels.point.gp))[ord])
    }else if(has.pch){
      gTree(children = gList(pointsGrob(0.5, 0.5,default.units="npc",
                                        size = unit(point.size, "char"),
                                        pch=pch[i], gp=gpi),
                             textGrob(label = labels.point[[i]],
                                      x = 0.5,y = 0.5,
                                      default.units="npc",
                                      gp = labels.point.gp)))
    }else if(do.lines){
      linesGrob(0:1, 0.5, gp=gpi)
    }else{
      nullGrob() # should not happen...
    }
 
    fg <- packGrob(fg, plGrob,
                   col = 2*ci-1, row = ri, border = symbol.border,
                   width = u1, height = u1, force.width = TRUE)
 
    # ============================================================================
    ## text grob: add the labels
    gpi. <- gpi
    gpi.$col <- "black" # maybe needs its own 'gp' in the long run (?)
    fg <- packGrob(fg, textGrob(labels[[i]], x = 0, y = 0.5,
                                just = c("left", "centre"), gp=gpi.),
                   col = 2*ci, row = ri, border = text.border)
  }
  fg
}
