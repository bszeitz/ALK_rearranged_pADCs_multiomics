##############################
#
# Color list for heatmaps
# Author: Bea Szeitz
#
##############################


colorlist <- list(
  "Pathology" = 
    c("normal"="white",
      "solid"="#377EB8",
      #"unannotated"="grey",
      "papillary" ="#4DAF4A",
      "tubular" ="#984EA3"),
  "Mucin score" = 
    c("0"="white",
      "1" ="#D9F0A3",
      "2"="#78C679",
      "3" ="#238443"),
  "Stroma score" = 
    c("0"="white",
      "1" ="#FED976",
      "2"="#FD8D3C",
      "2.5" ="#E31A1C",
      "3" ="black"),
  "pROI cluster" = 
    c("1"="#b2df8a",
      "2" ="#33a02c",
      "3"="#1f78b4",
      "4" ="#a6cee3",
      "5" ="#fb9a99",
      "6" = "#e31a1c",
      "unannotated" = "grey"),
  "Immune score" = 
    c("0"="white",
      "1" ="yellow",
      "2"="orange",
      "3" ="red"),
  "TIL (%)" = colorRamp2(c(0, 10, 20, 30 ,40), c("white","yellow", "orange","red","black")),
  #"Immune score" = colorRamp2(c(0,1,2,3), c("white","yellow", "orange","red")),
  "OS (years)" = colorRamp2(c(2, 6, 14), c("#FFD124","#00AFC1", "#006778")),
  "Case" = 
    c( 
      "Case1" = "gold",
      "Case2"="#A65628",
      "Case3" = "#FF7F00",
      "Case4" = "#F781BF",
      "Case5" = "darkblue",
      "Case6" = "#a1045a",
      "Case7" = "lightgreen")
)
