###########################################################################################
### SAMSARA - Figure 6: Networks                                                        ###
###########################################################################################

###########################################################################################
### Gets a network with edges and nodes correctly classified
get_nw <- function(edgesx, nodes) {
  nodes <- nodes %>%
    mutate(id = 1:n())
  # Update edges to use these IDs
  edges <- edgesx
  edges <- edges %>%
    left_join(nodes, by = c("from" = "nodenames")) %>%
    rename(from_id = id) %>%
    left_join(nodes, by = c("to" = "nodenames")) %>%
    rename(to_id = id) %>%
    select(from = from_id, to = to_id)
  edges$type <- edgesx$both 
  graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = "id")
  return(graph)
}


#adj <- basic_adj 
#sprsc <- basic_sprsc

get_dataprep <- function(adj,sprsc,eucm,intm,dset,repl) {
  ###########################################################################################
  sp_nms <- seq(1,25,1)
  ### Specify data
  adjd <- adj[adj$d == dset,]
  sprscd <- sprsc[sprsc$d == dset,]
  adjr <- adjd[adjd$repl == repl,1:25]
  sprscr <- sprscd[sprscd$repl == repl,1:3]
  eucmr <- eucm[eucm$repl == repl, 1:25]
  intmr <- intm[intm$repl == repl, 1:25]
  ###########################################################################################
  ### Wrangle
  colnames(adjr) <- sp_nms
  rownames(adjr) <- sp_nms
  rownames(sprscr) <- sp_nms
  intmr <- intmr[,1:25]
  eucmr <- eucmr[,1:25]
  rownames(intmr) <- sp_nms
  rownames(eucmr) <- sp_nms
  eucmr[eucmr > pos_thr] <- 0 # & eucmr < neg_thr
  eucmr[eucmr == 0] <- NA
  eucmr[eucmr < pos_thr] <- 'pos_env'
  #eucmr[eucmr > neg_thr & eucmr != 'pos_env'] <- 'neg_env'
  intmr[intmr >= 0] <- NA
  intmr[intmr < 0] <- 'pos_int'
  #intmr[intmr > 0 & intmr != 'pos_int'] <- 'neg_int'
  ###########################################################################################
  ### Wrangle data
  adjr[lower.tri(adjr)] <- 0
  adjr <- melt(as.matrix(adjr))
  eucmr <- melt(as.matrix(eucmr))
  intmr <- melt(as.matrix(intmr))
  colnames(adjr) <- c('from','to','value')
  adjr$euc <- eucmr$value
  adjr$int <- intmr$value
  adjr$both <- paste0(adjr$euc,'x',adjr$int)
  ###########################################################################################
  ### Node data
  sprscx <- as.data.frame(melt(as.matrix(sprscr)))
  sprscx[,1:2] <- sapply(sprscx[,1:2],as.character)
  sprscx <- sprscx[sprscx$value > 0,]
  duplsp <- as.character(sprscx$Var1[duplicated(sprscx$Var1)])
  nodesx <- as.data.frame(cbind(name = sp_nms))
  nodes <- sprscx[match(nodesx$name,sprscx$Var1),]
  nodes$Var1 <- sp_nms
  nodes <- nodes[,1:2]
  colnames(nodes) <- c('nodenames','description')
  if (length(duplsp) > 0) {
    nodes[nodes$nodenames %in% duplsp,]$description <- 'mult_rsc'
  }
  nodes[is.na(nodes$description),]$description <- 'none'
  return(list(adj = adjr, node = nodes))
}


###########################################################################################
plot_nw_pos_allthree <- function(basic_adj,basic_sprsc,nullint_adj,nullint_sprsc,nullrsc_adj,nullrsc_sprsc,eucm,intm,dset,repl) {
  
  basic <- get_dataprep(basic_adj,basic_sprsc,eucm,intm,dset,repl)
  nullint <- get_dataprep(nullint_adj,nullint_sprsc,eucm,intm,dset,repl)
  nullrsc <- get_dataprep(nullrsc_adj,nullrsc_sprsc,eucm,intm,dset,repl)
  ###########################################################################################
  basic_adjp <- basic$adj[basic$adj$value > 0,]
  nullint_adjp <- nullint$adj[nullint$adj$value > 0,]
  nullrsc_adjp <- nullrsc$adj[nullrsc$adj$value > 0,]
  ###
  basic_nods <- basic$node
  nullint_nods <- nullint$node
  nullrsc_nods <- nullrsc$node
  ###########################################################################################
  ### Get network
  basic_gp <- get_nw(basic_adjp, basic_nods)
  nullint_gp <- get_nw(nullint_adjp, nullint_nods)
  nullrsc_gp <- get_nw(nullrsc_adjp, nullrsc_nods)
  
  ###########################################################################################
  ### Stop. Plottin time
  set_graph_style(plot_margin = margin(1,1,1,1))
  basic_plt <- ggraph(basic_gp, layout = 'kk') + 
    geom_edge_link(aes(colour = type), width = 1) + 
    geom_node_point(aes(colour = description, size = 3)) +
    geom_node_text(aes(label = nodenames), size = 2.9)
  nullint_plt <- ggraph(nullint_gp, layout = 'kk') + 
    geom_edge_link(aes(colour = type), width = 1) + 
    geom_node_point(aes(colour = description, size = 3)) +
    geom_node_text(aes(label = nodenames), size = 2.9)
  nullrsc_plt <- ggraph(nullrsc_gp, layout = 'kk') + 
    geom_edge_link(aes(colour = type), width = 1) + 
    geom_node_point(aes(colour = description, size = 3)) +
    geom_node_text(aes(label = nodenames), size = 2.9)
  nullrsc_plt + basic_plt + nullint_plt  
}




###########################################################################################
plot_nw_posneg <- function(adj,sprsc,eucm,intm,dset,repld) {
  ###########################################################################################
  ### Specify data
  adjd <- adj[adj$d == dset,]
  sprscd <- sprsc[sprsc$d == dset,]
  adjr <- adjd[adjd$repl == repl,1:25]
  sprscr <- sprscd[sprscd$repl == repl,1:3]
  eucmr <- eucm[eucm$repl == repl, 1:25]
  intmr <- intm[intm$repl == repl, 1:25]
  ###########################################################################################
  ### Wrangle
  sp_nms <- colnames(adjr)
  rownames(adjr) <- sp_nms
  rownames(sprscr) <- sp_nms
  intmr <- intmr[,1:25]
  eucmr <- eucmr[,1:25]
  rownames(intmr) <- sp_nms
  rownames(eucmr) <- sp_nms
  eucmr[eucmr > pos_thr & eucmr < neg_thr] <- 0
  eucmr[eucmr == 0] <- NA
  eucmr[eucmr < pos_thr] <- 'pos_env'
  eucmr[eucmr > neg_thr & eucmr != 'pos_env'] <- 'neg_env'
  intmr[intmr == 0] <- NA
  intmr[intmr < 0] <- 'pos_int'
  intmr[intmr > 0 & intmr != 'pos_int'] <- 'neg_int'
  ###########################################################################################
  ### Wrangle data
  adjr[lower.tri(adjr)] <- 0
  adjr <- melt(as.matrix(adjr))
  eucmr <- melt(as.matrix(eucmr))
  intmr <- melt(as.matrix(intmr))
  colnames(adjr) <- c('from','to','value')
  adjr$euc <- eucmr$value
  adjr$int <- intmr$value
  adjr$both <- paste0(adjr$euc,'x',adjr$int)
  ###########################################################################################
  adjp <- adjr[adjr$value > 0,]
  adjn <- adjr[adjr$value < 0,]
  adjr <- adjr[adjr$value != 0,]
  ###########################################################################################
  ### Node data
  sprscx <- as.data.frame(melt(as.matrix(sprscr)))
  sprscx[,1:2] <- sapply(sprscx[,1:2],as.character)
  sprscx <- sprscx[sprscx$value > 0,]
  duplsp <- as.character(sprscx$Var1[duplicated(sprscx$Var1)])
  nodesx <- as.data.frame(cbind(name = sp_nms))
  nodes <- sprscx[match(nodesx$name,sprscx$Var1),]
  nodes$Var1 <- sp_nms
  nodes <- nodes[,1:2]
  colnames(nodes) <- c('nodenames','description')
  if (length(duplsp) > 0) {
    nodes[nodes$nodenames %in% duplsp,]$description <- 'mult_rsc'
  }
  nodes[is.na(nodes$description),]$description <- 'none'
  #if (sum(unique(as.character(nodes$description)) %in% c('R1','R2','R3','mult_rsc')) == 0) {
  #  nodes$description <- 'none'
  #}
  ###########################################################################################
  ### Get network
  gp <- get_nw(adjp, nodes)
  gn <- get_nw(adjn, nodes)
  ###########################################################################################
  ### Stop. Plottin time
  set_graph_style(plot_margin = margin(1,1,1,1))
  p_plt <- ggraph(gp, layout = 'kk') + 
    geom_edge_link(aes(colour = type)) + 
    geom_node_point(aes(colour = description, size = 3)) +
    geom_node_text(aes(label = nodenames))
  n_plt <- ggraph(gn, layout = 'kk') + 
    geom_edge_link(aes(colour = type)) + 
    geom_node_point(aes(colour = description, size = 3)) +
    geom_node_text(aes(label = nodenames))
  p_plt + n_plt
}




