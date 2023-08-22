## ----setup, include = FALSE, echo = FALSE, message = F, warning = F-----------
library(igraph)
library(knitr)
library(bookdown)
library(htmlTable)
library(iFraudSimulator)
options(table_counter = T)
# knit_hooks$set(plot = function(x, options) {
#   paste('<figure><figcaption>', options$fig.cap, '</figcaption><img src="',
#         opts_knit$get('base.url'), paste(x, collapse = '.'),
#         '"></figure>',
#         sep = '')
# })
FixTable <- function(x, ...) {
  htmlTable::htmlTable(x, rnames = F,
                           css.cell = rbind(rep("background: lightgrey; padding-left: .5em; padding-right: .2em;", 
                               times = ncol(x)),
                           matrix("", 
                                  ncol = ncol(x), 
                                  nrow = nrow(x))),
                           ...)
}

## ----SampleNetwork, echo = FALSE, fig.cap = "Sample network", fig.width = 7.5, fig.align='center'----
WebGraph = data.frame(from = c('C3','C3', 'C4', 'P2','P3','P3', 'C1','C5', 'P1','P4'),
                       to = c('P2', 'P3', 'P3','C1','C1', 'C5','P1','P4', 'C2', 'C2'  ))
g        = graph_from_data_frame(WebGraph, directed = FALSE)
pos      = cbind(c(1,1,2,2,2.5,3,4,4,5),c(3,1,4,2,3,0.1,2,0.1,1 ))
igraph::V(g)$color = sapply(igraph::V(g)$name, function(x) {
  if (grepl("C4", x)) {
    "red"
  } else if (grepl("C", x)) {
    "green"
  } else {
    "grey"
  }
  })
igraph::V(g)$shape = sapply(igraph::V(g)$name, function(x) {
  if (grepl("C", x)) {
    "circle"
  } else {
    "square"
  }
})
set.seed(63)
pos = layout_with_dh(g, coords = pos, weight.edge.lengths = edge_density(g) / 1e5)
plot.igraph(g, edge.label = NA, edge.color = 'black', layout = pos, 
            vertex.label = igraph::V(g)$name, vertex.color = igraph::V(g)$color, 
            vertex.label.color = 'black', vertex.size = 50, main = "Social network",
            vertex.shape = igraph::V(g)$shape)
legend("topleft", c("Unknown/ non-fraudulent claim", "Fraudulent claim", "Party (PH, broker, ...)"),
       pch = c(rep(21, 2), 22), col = "black", bty = "n", pt.bg = c("green", "red", "grey"),
       pt.cex = 2, y.intersp = 1.5,)

## ----echo = F, results = 'asis'-----------------------------------------------
NetwLabel = data.frame(
  'start node' = c('P2', 'P3', 'P3', 'C1', 'C1', 'C5', 'P1', 'P4', 'C2', 'C2'),
  'end node' = c('C3', 'C3', 'C4', 'P2', 'P3', 'P3', 'C1', 'C5', 'P1', 'P4'),
  stringsAsFactors = F,
  check.names = F
)
NetwLabel[grepl("C", NetwLabel$`start node`), 1:2] = NetwLabel[grepl("C", NetwLabel$`start node`), 2:1] 
NetwLabel = NetwLabel[order(NetwLabel$`start node`), ]
rownames(NetwLabel) = NULL
print(FixTable(NetwLabel, caption = 'Edges network', label = 'tab:Edges'))

## ----echo = F, results = 'asis'-----------------------------------------------
c0  = c(rep(0, 3), 1, 0) 
FrM = data.frame("Fraud indicator" = c0, check.names = F, row.names = NULL)
print(FixTable(FrM, caption = 'Fraud indicator', label = 'tab:FraudInd'))


## ---- echo = T, message = F, warning = F, results = 'hide'--------------------
NetwLabel = data.frame(
  startNode = c('P2', 'P3', 'P3', 'C1', 'C1', 'C5', 'P1', 'P4', 'C2', 'C2'),
  endNode = c('C3', 'C3', 'C4', 'P2', 'P3', 'P3', 'C1', 'C5', 'P1', 'P4'),
  stringsAsFactors = F
)
NetwLabel[grepl("C", NetwLabel$startNode), 1:2] = NetwLabel[grepl("C", NetwLabel$startNode), 2:1] 
NetwLabel = NetwLabel[order(NetwLabel$startNode), ]
NetwLabel$FraudInd = sapply(NetwLabel$endNode, function(x)
  if (x == "C4")
    1
  else
    0)
NetwLabel$startNode = as.numeric(gsub("P", "", NetwLabel$startNode))
NetwLabel$endNode   = as.numeric(gsub("C", "", NetwLabel$endNode))
Results  = BiRankFr(NetwLabel, data.frame(FraudInd = c0))

## ---- echo = F, results = 'asis', message = F, warning = F--------------------
FrSc     = Results$ResultsClaims
FrSc[, 2:4] = lapply(FrSc[, 2:4], round, digits = 3) 
colnames(FrSc) = c("Claim", "Fraud score", "Std score", "Scaled score")
print(FixTable(FrSc, caption = 'Fraud scores', label = 'tab:FraudScores'))


