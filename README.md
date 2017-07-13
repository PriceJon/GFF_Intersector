# GFF_Intersector
An R program capable of intersecting a .GFF file with a chosen flanking region and large files containg genomic co-ordinates from multiple experiments and visualising the intersections genome wide as well as producing individual region plots and quantification of the intersections. 

# Prerequsites

R programming language

The program also relies on a numer of R libraries, although install proceedures are included you can check before and install manually. 

```{r}
  library(shiny)

  library(ggrepel)

  library(gplots)

  library(RColorBrewer)

  library(GenomicRanges)

  library(DT)

  library(ggplot2)

  library(ggbio)
```

# Running

```{r}
library(shiny)

shiny::runGitHub('GFF_Intersector','PriceJon')
```


Currently only GFF files are accepted along with another file containing at least chromsome, start, stop and a unique ID columns contained within the first 10 columns. 

Instructions are available whilst using the program. 


# Updates

This is an active development and new features will be added monthly.

- Quantification of bases intersected in the upstream/downstream/genic/IG etc regions. 

- The ability to separate different experiment based on their ID and analyse spearately

- The ability to download high quality PDFs





