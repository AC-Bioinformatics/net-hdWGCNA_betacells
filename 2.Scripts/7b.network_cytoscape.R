# Cytoscape view of one module network from 7.network_visualisation.R, built with RCy3.
# RCy3 talks to a running Cytoscape (3.10+) through its REST API, so Cytoscape must be open
# on the machine running this script.
# Inputs (3.Results/12.Network/<module>/):
#   <module>_cytoscape_edges.txt  - WGCNA exportNetworkToCytoscape edges (TOM weights)
#   <module>_node_attributes.csv  - kME, connectivity, hub flag and Hallmark label per gene
# Style: node size = scaled kWithin, node colour = kME, hub genes with a thick border,
# edge width / opacity = TOM weight. Saves a .cys session and a PDF image.

# Path configuration: folder holding 3.Results
Root <- "~/network"
mod <- "brown"
In_dir <- file.path(Root, "3.Results", "12.Network", mod)

library(RCy3)
cytoscapePing()

# Edges: WGCNA columns fromNode / toNode / weight -> RCy3 source / target / weight
edges <- read.delim(file.path(In_dir, paste0(mod, "_cytoscape_edges.txt")))
edges <- data.frame(source = edges$fromNode, target = edges$toNode,
                    interaction = "TOM", weight = edges$weight)

# Nodes: every module gene (genes without edges stay as unconnected nodes)
nodes <- read.csv(file.path(In_dir, paste0(mod, "_node_attributes.csv")))
nodes <- data.frame(id = nodes$gene_name, nodes[, setdiff(colnames(nodes), "gene_name")])

net_title <- paste(mod, "module network")
createNetworkFromDataFrames(nodes, edges, title = net_title, collection = "hdWGCNA modules")

# Visual style
mod_hex <- rgb(t(col2rgb(mod)), maxColorValue = 255)
style <- paste0(mod, "_module_style")
createVisualStyle(
  style,
  defaults = list(NODE_SHAPE = "ELLIPSE", NODE_BORDER_PAINT = "#333333",
                  NODE_LABEL_FONT_SIZE = 10, EDGE_STROKE_UNSELECTED_PAINT = "#999999"),
  mappings = list(
    mapVisualProperty("node label", "id", "p"),
    mapVisualProperty("node size", "kWithin_scaled", "c", c(0, 1), c(15, 70)),
    mapVisualProperty("node fill color", "kME", "c",
                      range(nodes$kME), c("#FFFFFF", mod_hex)),
    mapVisualProperty("node border width", "hub", "d", c("hub", "other"), c(6, 1)),
    mapVisualProperty("edge width", "weight", "c", range(edges$weight), c(0.5, 4)),
    mapVisualProperty("edge transparency", "weight", "c", range(edges$weight), c(60, 220))
  )
)
setVisualStyle(style)

# Force-directed layout weighted by TOM
layoutNetwork("force-directed edgeAttribute=weight")

saveSession(file.path(In_dir, paste0(mod, "_network.cys")))
exportImage(file.path(In_dir, paste0(mod, "_network")), type = "PDF")
message("Saved: ", In_dir)
