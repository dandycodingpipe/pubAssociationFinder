#Visualization
#' Bio-concept co-occurence network visualization
#'
#' This function allows you visualize the literature landscape as an interactive network and explore how bioconcepts directly or indirectly relate to one another. 
#' @name systox
#' @param data stricly-bioCAF association rules data-frame
#' @param style The network layout style. We recommend "layout_nicely" or "layout_in_circle" (For additional styles see https://igraph.org/r/doc/layout_.html)
#' @param class The classification display option ("all","bme","anatomy", "organism","diseases","chemicals","techniques","other","unlab")

#' @keywords Network, co-occurence network, exposome, exposomics
#' @export
#' @examples
#' concepts <- PubTator("tobacco smoke toxicity")
#' results <- bioCAF(concepts, "all", 0.003, "C:/Users/Chris/OneDrive/2023/Systox/venvJune19", "en_core_web_lg")
#' viz <- systox(results, "layout_nicely")
#' 
#' #export option as html
#' library(htmlwidgets)
#' saveWidget(viz, file = "bioTM_hyperthyroidism.html")

systox <- function(data, style){
data <- data[,1:2]
#data <- data[c(94:103),]
data$LHS <- gsub("[{}]", "", data$LHS)
data$RHS <- gsub("[{}]", "", data$RHS)

cells <- paste(data$LHS, data$RHS, sep = ",")

matches <- str_extract_all(cells, "\\b(g_|d_|c_)\\w+")  # Extract gene names



data$genes <- sapply(matches, paste, collapse = " ")

unique_genes <- unique(unlist(matches))

co_occurrence_matrix <- matrix(0, nrow = length(unique_genes), ncol = length(unique_genes))
rownames(co_occurrence_matrix) <- colnames(co_occurrence_matrix) <- unique_genes

for (i in 1:nrow(data)) {
  genes <- unlist(matches[i])
  for (gene1 in genes) {
    for (gene2 in genes) {
      if (gene1 != gene2) {
        co_occurrence_matrix[gene1, gene2] <- co_occurrence_matrix[gene1, gene2] + 1
      }
    }
  }
}

co_occurrence_network <- co_occurrence_matrix

library(igraph)

require(visNetwork, quietly = TRUE)

# Stack overflow on combinig igraph and visNetwork

 test.gr <- graph_from_adjacency_matrix(co_occurrence_network, mode = "undirected", weighted = T, diag = FALSE)
  
 test.visn <- toVisNetworkData(test.gr)

 test.visn$edges$value <- 1+log10(test.visn$edges$weight)

 # Count the number of unique edges for each "from" observation
 edge_counts <- test.visn$edges %>%
   group_by(from) %>%
   summarize(unique_edges = n_distinct(to))
 
 edge_counts2 <- test.visn$edges %>%
   group_by(to) %>%
   summarize(unique_edges = n_distinct(from))
 
 frequency <- full_join(edge_counts, edge_counts2, by = c("from" = "to"))
 frequency$unique_edges.x <-  ifelse(is.na(frequency$unique_edges.x), 0, frequency$unique_edges.x)
 frequency$unique_edges.y <- ifelse(is.na(frequency$unique_edges.y), 0, frequency$unique_edges.y)
  
 frequency$count <- frequency$unique_edges.x + frequency$unique_edges.y
 frequency <- frequency %>% arrange(from)
 
 test.visn$nodes <- test.visn$nodes %>%
   arrange(id)
 
 not_in_frequency <- test.visn$nodes$id[!test.visn$nodes$id %in% frequency$from]
 test.visn$nodes <- test.visn$nodes[!test.visn$nodes$id %in% not_in_frequency, ]
 
 test.visn$nodes$frequency <- frequency$count
 

 
 
 node_size <- 1
 edge_size <- 0.001
 
 #test.visn$nodes$title <- paste0("<p><a href='https://www.uniprot.org/uniprotkb?query=", substr(test.visn$nodes$id, 3, nchar(test.visn$nodes$id)), "'>",substr(test.visn$nodes$id, 3, nchar(test.visn$nodes$id)), "</a></p>")
 
 test.visn$nodes$title <- ifelse(substr(test.visn$nodes$id, 1, 2) == "d_",
                                 paste0("<p>", test.visn$nodes$id, " (Connections: ", test.visn$nodes$frequency, ")<br>disease information not available</p>"),
                                 ifelse(substr(test.visn$nodes$id, 1, 2) == "g_",
                                        paste0("<p><a href='https://www.uniprot.org/uniprotkb?query=", substr(test.visn$nodes$id, 3, nchar(test.visn$nodes$id)), "'>", substr(test.visn$nodes$id, 3, nchar(test.visn$nodes$id)), "</a>", " (Connections: ", test.visn$nodes$frequency, ")</p>"),
                                        ifelse(substr(test.visn$nodes$id, 1, 2) == "c_",
                                               paste0("<p>", test.visn$nodes$id, " (Connections: ", test.visn$nodes$frequency, ")<br>chemical information not available</p>"), "tomato")
                                 )
 )
 
 
 test.visn$nodes$color <- ifelse(substr(test.visn$nodes$id,1,2) == "d_",
                          "#FF1D15",
                          ifelse(substr(test.visn$nodes$id,1,2) == "g_",
                                 "#3EC300",
                                 ifelse(substr(test.visn$nodes$id,1,2) == "c_",
                                        "#337CA0", "grey"
                                 )
                          )
 )
 
 test.visn$nodes$group <- ifelse(substr(test.visn$nodes$id,1,2) == "d_",
                                 "Disease",
                                 ifelse(substr(test.visn$nodes$id,1,2) == "g_",
                                        "Gene",
                                        ifelse(substr(test.visn$nodes$id,1,2) == "c_",
                                               "Chemical", "tomato"
                                        )
                                 )
 )
 
 
 
output <- visNetwork(test.visn$nodes, test.visn$edges, randomSeed = 1124, height = "900px", width = "100%", main = "bioTM Co-occurence Network", submain = paste(length(test.visn$nodes$id), "unique concepts") ) %>% 
                                    visIgraphLayout(layout = style, type = "full") %>% 
                                    visGroups(groupname = "Chemical", color = "#337CA0") %>%
                                    visGroups(groupname = "Disease", color = "#FF1D15") %>%  
                                    visGroups(groupname = "Gene", color = "#3EC300") %>% 
                                    visLegend() %>%
                                    visOptions(selectedBy = "group", highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T) 


return(output)
}
 