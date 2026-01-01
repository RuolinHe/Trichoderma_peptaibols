library(readxl)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")
path <- "D:/课题组/zhiyuan_Lab/4-Peptaibols_NRPS/data/Peptaibol肽段馏分组成.xlsx"
haircolors <- c("Aib","Ala","Gln","Glu","Glu(Me)","Gly","Iva","Leu","Lxx","Phe","Pro","Tyr","Val","Vxx") #↓按照matlab里unique_substrate的顺序
groupColors <- c( "#FF3744", "#FFFF00","#7EE1A8","#CFE9B0","#B8E001","#3C7747","#6eb4f6","#3844C6","#296EC1","#6C5BFF","#D3C1B3","#7A5D8B","#727171","#000000")
## origin data
m2 <- read_excel(path,sheet = "tabulate", range = "B2:O15",col_names = FALSE)
m3 <-matrix(unlist(m2), ncol = 14, nrow = 14)
dimnames(m3) <- list(have = haircolors,
                     prefer = haircolors)
# 导出到网页html再用https://www.sejda.com/html-to-pdf转成pdf
# Build the chord diagram:
chorddiag(m3, groupColors = groupColors, groupnamePadding = 27,
          groupnameFontsize = 24,  showTicks = TRUE,
          tickInterval = 5, ticklabelFontsize = 20,)


# 下面是旧的版本
path <- "D:/课题组/zhiyuan_Lab/4-Peptaibols_NRPS/data/trichohypolin structure.xlsx"
haircolors <- c("Aib","Ala","Gln","Glu","Glu(Me)","Gly","Lxx","Phe","Tyr","Pro","Vxx") #↓按照matlab里unique_substrate的顺序
groupColors <- c( "#FF3744", "#F8DB3F", "#7EE1A8","#CFE9B0","#9EB483","#3C7747","#296EC1","#6C5BFF","#7A5D8B","#D3C1B3","#000000")

## origin data
m2 <- read_excel(path,sheet = "tabulate", range = "B2:L12",col_names = FALSE)
m3 <-matrix(unlist(m2), ncol = 11, nrow = 11)
dimnames(m3) <- list(have = haircolors,
                     prefer = haircolors)
# Build the chord diagram:
chorddiag(m3, groupColors = groupColors, groupnamePadding = 27,
          groupnameFontsize = 24,  showTicks = TRUE,
          tickInterval = 5, ticklabelFontsize = 20,)



