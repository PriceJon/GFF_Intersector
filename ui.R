###########################################
###########################################
# SETUP
###########################################
###########################################


#---------------------------------------
#load libraries
options(shiny.maxRequestSize=200*1024^2) 

if (!require("shiny")) {
  install.packages("shiny", dependencies = TRUE)
  library(shiny)
}


#---------------------------------------
#load funcions for GFF

getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {  
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}



gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}



###########################################
###########################################
# UI
###########################################
###########################################

ui <- fluidPage(
  ###########################################
  # 1 get the input files from a list, render the go button.
  ###########################################
  headerPanel(
    HTML('<img src="./gene_model.png", height="140px"    
         style="float:right"/>','<img src="region.png", height="140px"    
         style="float:left"/>','<center><em><p style="color:rgb(8, 48, 100)">GFF Intersector </p></em><br></center>' )),
  
  br(),
  br(),
  hr(),
  br(),
  br(),
  br(),
  br(),
  # fluidRow(
  #  column(12, align="center",
  #  headerPanel("Coordinate Intersector", windowTitle = "COINT1.0")
  #  )
  #    ),
  br(),
  p("GFF Intersector has been designed to intersect a set of genomic cooridinates (abtained from any experiment) with a set of genes in GFF format. The input for the program is two files; 1) A GFF format file contaning exon information for the genes of interest, It is best here to use a subsetted GFF although you can use a whole genome GFF (it just may take some time), 2) A file contaning the genomic locations you wish to intersect, this file should contain chromsome, start, end and a unique ID column. When these files are ready, find them using the pull downs below and press read in. The program then runs in 3 stages; 1) Loading and setup, Here you can load in your data, examine it and tell the program which columns contains the useful information. 2) Global analysis, here you can examine the intersections on a global level. 3) Region level analysis, the last section allows you to explore particular regions of interest."),
  br(),
  
  fluidRow(
    column(6,align="center",
           fileInput(inputId = "file1",multiple = F,label = "Select GFF")),
    column(6,align="center",
           fileInput(inputId = "file2",multiple = F,label = "Select co-ordinate file")
    )),
  actionButton(inputId = "go",
               label = "Read In"),
  ###########################################
  # 2 after pressing go print the paths specified and
  #render the tables nicely so one can examine them
  ###########################################
  verbatimTextOutput("path"),
  br(),
  hr(),
  hr(),
  h3("1. Loading and setup"),
  br(),
  dataTableOutput('mytable'),
  dataTableOutput('mytable2'),
  
  
  ###########################################
  # 3 now render the column choices, (not reactive)
  ###########################################
  DT::dataTableOutput('foo'),
  br(),
  p("Examine the tables above and choose the correct columns for Chromsome, start, end and ID. When selections are complete you can examine the two boxes below to check that the selection is correct. Then choose the required flanking region for the genes in the GFF below, choose a filename then intersect."),
  ###########################################
  # 4 this bit renders the choices in real time
  ###########################################
  
  fluidRow(
    column(3,   sliderInput("flank", "Flanking bp for GFF:", 
                            min=0, max=10000, value=2000)
    ),
    column(3,
           textInput("file", label = ("Intersection Analysis"), 
                     value = "Filename here . . "),
           actionButton(inputId = "go4",
                        label = "Find intersections and write to table")
    ),
    column(3,
           h4("Example from file1"),
           verbatimTextOutput('start1')
    ),
    column(3,
           h4("Example from file2"),
           verbatimTextOutput('start2')
    )
  ),
  br(),
  hr(),
  hr(),
  br(),
  h3("2. Global Intersections"),
  p("Global intersection are displayed here, the top plot shows a histogram with genes in blue and your regions in green. Intersections are shown as red bar. By highlighting on the top plot you can render a zoomed in plot below (if labels are confused choose a smaller region)"),
  verbatimTextOutput('intstats'),
  plotOutput("lrp1",   brush = "plot_brush"),
  plotOutput("lrp2"),
  #test for brush : uncomment if needed:
  #verbatimTextOutput("info"),
  
  ###########################################
  # 5 now once is complete choose analyse to intersect the coordinates
  ###########################################
  # actionButton(inputId = "go2",
  #     label = "If above looks correct, find out if there are ANY intersections"),
  ###########################################
  # 6 and print the stats
  ###########################################
  # verbatimTextOutput("numbers"),
  
  ###########################################
  # 8 simple test, are there intersects?
  ###########################################
  # verbatimTextOutput("ranges1"),
  #  verbatimTextOutput("ranges2"),
  
  ###########################################
  # 8 simple test, are there intersects?
  ###########################################
  # actionButton(inputId = "go3",
  #               label = "Give me the Regions of Intersections"),
  
  ###########################################
  # 9 properly
  ###########################################
  #verbatimTextOutput("funions"),
  
  #allow me to print these nicely
  br(),
  hr(),
  hr(),
  br(),
  h3("3. Individual Regions"),
  p("The below table has been written to file and contains all the intersect regions identified in this dataset. The next plot then shows how many of your regions were intersected within each region. This can be useful to identify the most promenant regions. Then you can choose a particular region to plot below that."),
  #verbatimTextOutput("check"),
  dataTableOutput('mytableInts'),
  #verbatimTextOutput("check"),
  plotOutput("allplot"),
  ###########################################
  # 10 analys individual intersections
  ###########################################
  hr(),
  fluidRow(
    column(4,
           uiOutput("regionchoice")
    ),
    column(5,
           textInput("pdf", label = (""), 
                     value = "save.pdf"),
           actionButton(inputId = "gopdf",
                        label = "Save PDF")
    )),
  # verbatimTextOutput("check"),
  plotOutput("rangeplot")
  
  
  ###########################################
  #UI end
  ###########################################
    ) #<---- end


