###########################################
###########################################
# SETUP
###########################################
###########################################


#---------------------------------------
#load libraries
options(shiny.maxRequestSize=200*1024^4) 

if (!require("shiny")) {
  install.packages("shiny", dependencies = TRUE)
  library(shiny)
}

if (!require("ggrepel")) {
  install.packages("ggrepel", dependencies = TRUE)
  library("ggrepel")
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("GenomicRanges")) {
  install.packages("GenomicRanges", dependencies = TRUE)
  library("GenomicRanges")
}
if (!require("DT")) {
  install.packages("DT", dependencies = TRUE)
  library("DT")
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library("ggplot2")
}
if (!require("ggbio")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggbio")
  library("ggbio")
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
# server
###########################################
###########################################

server <- function(input, output) {
  ###########################################
  #this is where the outputs are made
  ###########################################
  
  ###########################################
  # 1 get the file paths and make two variables, 
  # then render the variables for printing
  ###########################################
  file1 = eventReactive(input$go,{
    input$file1$name
  })
  file2 = eventReactive(input$go,{
    input$file2$name
  })
  output$path = renderPrint({
    print(file1())
    file2()
  })
  
  ###########################################
  #2 render the tables for examination
  ###########################################
  #get the files read in to the app
  table1 = eventReactive(input$go,{
    gff <- gffRead(input$file1$datapath)
    #make GFF usable
    gff = gff[gff[,3] == "exon",]
    gff$id = sapply(gff$attributes,FUN = getAttributeField,"Parent")
    gff
    
  })
  table2 = eventReactive(input$go,{
    read.table(input$file2$datapath)
  })
  
  
  
  output$mytable2 = renderDataTable({
    read.table(input$file2$datapath)
  })
  output$mytable = renderDataTable({
    table1()
  })
  
  
  
  ###########################################
  # 3 now render the column choices, (not reactive)
  ###########################################
  m = matrix(as.character(1:10),nrow = 8, ncol = 10, byrow = TRUE, dimnames = list(c("start1","end1","chr1","ID1","start2","end2","chr2","ID2"),seq(1,10)))
  
  colname = c("start1","end1","chr1","ID1","start2","end2","chr2","ID2")
  row.names(m) = colname
  for (i in seq_len(nrow(m))) {
    m[i, ] = sprintf(
      '<input type="radio" name="%s" value="%s"/>',
      colname[i], m[i, ]
    )
  }
  output$foo =
    DT::renderDataTable(m, escape = FALSE, selection = 'none', server = FALSE,
                        options = list(dom = 't', paging = FALSE, ordering = FALSE),
                        callback = JS("table.rows().every(function(i, tab, row) {
                                      var $this = $(this.node());
                                      $this.attr('id', this.data()[0]);
                                      $this.addClass('shiny-input-radiogroup');
});
                                      Shiny.unbindAll(table.table().node());
                                      Shiny.bindAll(table.table().node());"
    )
    )
  ###########################################
  #4 this bit renders the choices in real time
  ###########################################
  output$sel = renderPrint({
    str(sapply(colname, function(i) input[[i]]))
  })
  output$start1 = renderPrint({
    print("chr       start         end ")
    paste(table1()[,as.numeric(input[["chr1"]])][1], table1()[,as.numeric(input[["start1"]])][1],table1()[,as.numeric(input[["end1"]])][1],table1()[,as.numeric(input[["ID1"]])][1],sep = "     ")
  })
  output$start2 = renderPrint({
    print("chr       start         end ")
    paste(table2()[,as.numeric(input[["chr2"]])][1], table2()[,as.numeric(input[["start2"]])][1],table2()[,as.numeric(input[["end2"]])][1],table2()[,as.numeric(input[["ID2"]])][1],sep = "     ")
  })
  
  
  ###########################################
  # 5 now once is complete choose analyse to intersect the coordinates
  ###########################################
  
  
  
  ###########################################
  # 6 and print the stats
  ###########################################
  #get the numbers of corrdinates in each file
  number1 = eventReactive(input$go4,{
    nrow(table1())
  })
  number2 = eventReactive(input$go4,{
    nrow(table2())
  })
  output$numbers = renderPrint({
    paste("Corrdinates in file 1: ",number1()," and in file 2: ",number2())
  })
  
  ###########################################
  # 7 save ranges in ranges()$gr1 and ranges()$gr2
  ###########################################
  ranges = eventReactive(input$go4,{
    gr1 <-
      GRanges(seqnames =
                Rle(as.numeric(regmatches(table1()[,as.numeric(input[["chr1"]])],regexpr("\\d{1,2}",table1()[,as.numeric(input[["chr1"]])])))),
              ranges =
                IRanges(table1()[,as.numeric(input[["start1"]])], end = table1()[,as.numeric(input[["end1"]])],names = table1()[,as.numeric(input[["ID1"]])]),
              strand =
                Rle(strand(rep("*",as.numeric(number1())))),
              sample2 = table1()[,as.numeric(input[["ID1"]])])
    
    
    gr2 <-
      GRanges(seqnames =
                Rle(regmatches(table2()[,as.numeric(input[["chr2"]])],regexpr("\\d{1,2}",table2()[,as.numeric(input[["chr2"]])]))),
              ranges =
                IRanges(table2()[,as.numeric(input[["start2"]])], end = table2()[,as.numeric(input[["end2"]])],names = table2()[,as.numeric(input[["ID2"]])]),
              strand =
                Rle(strand(rep("*",as.numeric(number2())))),
              sample2 = table2()[,as.numeric(input[["ID2"]])])
    grl <- GRangesList("gr1" = gr1, "gr2" = gr2)
    
  })
  
  ###########################################
  # 8 simple test, are there intersects?
  ###########################################
  output$ranges1 = renderPrint({
    print("Intersect of two files")
    GenomicRanges::intersect(ranges()$gr1,ranges()$gr2)
  })
  output$ranges2 = renderPrint({
    print("Union of two files")
    GenomicRanges::union(ranges()$gr1,ranges()$gr2)
  })
  
  
  ###########################################
  # 9 FIND THE ACTUAL OVERLAPS
  ###########################################
  unions = eventReactive(input$go4, {
    gr1 = ranges()$gr1
    start(gr1) <- start(gr1) - input$flank
    end(gr1) <- end(gr1) + input$flank
    withProgress(message = 'Intersecting .. ..', value = 0, {
      incProgress(1/4, detail = paste("Doing part 1"))
      uni =    GenomicRanges::union(gr1,ranges()$gr2)
      incProgress(2/4, detail = paste("Doing part 2"))
      int = GenomicRanges::intersect(gr1,ranges()$gr2)
      incProgress(3/4, detail = paste("Doing part 3"))
      i1 = findOverlaps(uni,int)
      incProgress(4/4, detail = paste("Doing part 4"))
      final_unions = uni[queryHits(i1),]
      final_unions = reduce(final_unions)
      final_unions
    })
  })
  
  # output$funions = renderPrint({
  #   print(paste("Final Regions of Intersections = ",length(unions())))
  #   unions()
  # })
  
  
  
  ###########################################
  # 9 write th intersects to file
  ###########################################
  #gte the intersections in a dataframe file
  file = eventReactive(input$go4,{
    if(nchar(input$file) < 6 ){
      x = "Please enter a file name grater that 4 letters"
      x
    }else{
      finaldf = c()
      for (i in 1:length(unions())){
        grx = unions()[i,]
        grx1 = as.data.frame(grx)[1:3]
        grx1$file = "region"
        grx1$region = i
        grx1$names = paste("intersectregion",i,sep="")
        
        if1 = findOverlaps(ranges()$gr1,grx)
        file1ints = ranges()$gr1[queryHits(if1),]
        f11 = as.data.frame(file1ints,row.names = NULL)[1:3]
        f11$file = "file1"
        f11$region = i
        f11$names = names(file1ints)
        
        if2 = findOverlaps(ranges()$gr2,grx)
        file2ints = ranges()$gr2[queryHits(if2),]
        f21 = as.data.frame(file2ints)[1:3]
        f21$file = "file2"
        f21$region = i
        f21$names = names(file2ints)
        
        finaldf = rbind(finaldf,grx1,f11,f21)
        
      }
      write.table(finaldf,file = input$file,quote = F,row.names = F, sep = "\t")
      #x = "File is written"
      finaldf = as.data.frame(finaldf,row.names = NULL)
      finaldf
    }
  })
  output$check = renderPrint({
    file()
  })
  output$mytableInts = renderDataTable({
    file()
  })
  
  
  ###########################################
  # 10 Global analysis
  ###########################################
  
  
  #get stats of intersections
  output$intstats = renderText({
    table = file()
    regions = nrow(table[table$file == "region",])
    genes = length(unique(table[table$file == "file1",]$names))
    region =  nrow(table[table$file == "file2",])
    paste("Global Intersection Stats:
          Number of Intersect Regions: ",regions, "Number of genes intersected: ",genes,"Number of regions intersecting these genes: ",region)
    
    
  })
  
  
  
  
  output$lrp1 = renderPlot({
    
    
    df <- data.frame(seqnames=seqnames(ranges()$gr1),
                     starts=start(ranges()$gr1)-1,
                     file = rep("file1",length(ranges()$gr1)))
    df2 <- data.frame(seqnames=seqnames(ranges()$gr2),
                      starts=start(ranges()$gr2)-1,
                      file = rep("file1",length(ranges()$gr2)))
    df3 = file()
    
    #df = data.frame(x = rnorm(100), x2 = rnorm(100, mean=2))
    
    # g = ggplot(df, aes(x)) + geom_histogram( aes(x = x, y = ..density..),
    #                                         binwidth = diff(range(df$x))/30, fill="blue") + 
    #   geom_histogram( aes(x = x2, y = -..density..), binwidth = diff(range(df$x))/30, fill= "green")
    # print(g)
    # pdf("region_plot.pdf")
    # ggplot() + geom_histogram(data = df,aes(x=starts, y = ..density..),alpha = 0.4,fill ="slateblue",binwidth = .1e6 )+
    #   geom_histogram(data = df2,aes(x=starts, y = -..density..),fill = "darkgreen",alpha = 0.5 ,binwidth = .1e6)+
    #   geom_vline(data = df3,aes(xintercept=start), colour = "darkred")+
    #   #geom_histogram(data = df3,aes(x=start, y = 1),fill = "darkred",alpha = 0.4 ,binwidth = .1e6)+
    #   #geom_histogram(data = df3,aes(x=start, y = -1),fill = "darkred",alpha = 0.4 ,binwidth = .1e6)+
    #   facet_grid(seqnames ~ .) +
    #   theme_bw()+
    #   theme(panel.grid = element_blank())
    # dev.off()
    # 
    ggplot() + geom_histogram(data = df,aes(x=starts, y = ..density..),alpha = 0.4,fill ="slateblue",binwidth = .1e6 )+
      geom_histogram(data = df2,aes(x=starts, y = -..density..),fill = "darkgreen",alpha = 0.5 ,binwidth = .1e6)+
      geom_vline(data = df3,aes(xintercept=start), colour = "darkred")+
      #geom_histogram(data = df3,aes(x=start, y = 1),fill = "darkred",alpha = 0.4 ,binwidth = .1e6)+
      #geom_histogram(data = df3,aes(x=start, y = -1),fill = "darkred",alpha = 0.4 ,binwidth = .1e6)+
      facet_grid(seqnames ~ .) +
      theme_bw()+
      theme(panel.grid = element_blank())
    
    
    
    
  })
  
  #this just checks to see if he brush is working, output canbe toggled by comments on UI
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    xy_range_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("xmin=", round(e$xmin, 1), " xmax=", round(e$xmax, 1), 
             " ymin=", round(e$ymin, 1), " ymax=", round(e$ymax, 1))
    }
    
    paste0(
      
      "brush: ", xy_range_str(input$plot_brush)
    )
  })
  
  #This renders the second plot of the zoomed region
  output$lrp2 = renderPlot({
    
    
    df3 = file()
    
    df3 = df3[df3$start > input$plot_brush$xmin & df3$end < input$plot_brush$xmax,]
    df3 = df3[df3$file == "region",] 
    
    gr1 <-
      GRanges(seqnames =
                Rle( df3$seqnames),
              ranges =
                IRanges(df3$start,end = df3$end),
              strand =
                Rle("*",length(df3$start)),
              sample2 = df3$names,
              namess = "Intersect region",
              level = 0.5)
    
    
    
    grg = ranges()$gr1[start(ranges()$gr1) > input$plot_brush$xmin & start(ranges()$gr1) < input$plot_brush$xmax,]
    grd = ranges()$gr2[start(ranges()$gr2) > input$plot_brush$xmin & start(ranges()$gr2) < input$plot_brush$xmax,]
    
    mcols(grg)$namess = "Gene"
    mcols(grd)$namess = "User region"
    
    mcols(grg)$level = 1
    mcols(grd)$level = 1.5
    
    gL <- GRangesList(gl = gr1, g1 = grg, g2 = grd)
    
    
    
    
    txt = gL
    txt$g1 = txt$g1[!(duplicated(mcols(txt$g1)$sample2))]
    
    ggplot(txt) + geom_rect(aes(colour = namess,fill = namess) ) +
      #  geom_text(aes(label = sample2, x = start, y = level))+
      geom_text_repel(aes(label = sample2, x = start,y = level),
                      color = "gray20",
                      force = 10)+
      geom_rect(data = gL, aes(colour = namess,fill = namess) ) +
      xlim(input$plot_brush$xmin,input$plot_brush$xmax)+
      facet_grid(seqnames ~ .)+
      theme_bw()+
      theme(legend.position="bottom")
    
  })
  
  
  
  plotall = eventReactive(input$go4,{
    table = file()
    regions = unique(table$region)
    df = c()
    for(i in regions){
      line = c(as.numeric(i),as.numeric(nrow(table[table$region == i & table$file == "file2",])) )
      df = rbind(df,line)
    }
    colnames(df) = c("region","number")
    row.names(df) = seq(1:nrow(df))
    as.data.frame(df)
  })
  
  
  output$check = renderPrint({
    plotall()
  })
  
  
  output$allplot = renderPlot({
    ggplot(plotall(),aes(x = number,y = region)) + geom_point() +xlab("Number of regions intersecting Gene")+ ylab("Region")+ theme_bw()+
      geom_text(aes(label=ifelse(number>1,as.character(region),'')),hjust=0,vjust=0)
  })
  
  ###########################################
  # 10 analyse regions indivudally
  ###########################################
  
  
  
  #get choice
  output$regionchoice <- renderUI({
    cities <- seq(1:length(unions()))
    selectInput("regionChoice", "Choose Region to Analyse", cities)
  })
  
  #make plot of choice region
  output$rangeplot = renderPlot({
    table= file()[file()$region == input$regionChoice, ]
    #table
    x = GRanges(seqnames =
                  Rle(table$seqnames),
                ranges =
                  IRanges(table$start,table$end,names = table$names),
                strand =
                  Rle(rep("*",nrow(table))),
                sample = table$names,
                sample2 = sub("_exon.*","",table$names,perl = T))
    diff = (( max(table$end ) - min(table$start )) /100) * 10
    labels = table[table$file == "file2,"]
    #  start(x) <- start(x) + 1000
    #   end(x) <- end(x) - 1000
    print(ggplot(x) + geom_alignment(stat = "stepping",aes(group = sample2, fill = sample2, colour = sample2),group.selfish = F)  + xlim( ( min(table$start ) -diff),( max(table$end ) +diff))  + theme_bw())
    #ggplot(x) + geom_chevron( aes(color = rep(names(x),2)))+theme_bw()
  })
  
  #save pdf doesnt work
  savepdf = eventReactive(input$gopdf,{
    table= file()[file()$region == input$regionChoice, ]
    #table
    x = GRanges(seqnames =
                  Rle(table$seqnames),
                ranges =
                  IRanges(table$start,table$end,names = table$names),
                strand =
                  Rle(rep("*",nrow(table))),
                sample = table$names)
    pdf(file = input$pdf)
    print(ggplot(x) + geom_alignment(stat = "stepping",aes(group = sample),group.selfish = F) + theme_bw())
    dev.off()
  })
  ###########################################
  #server end
  ###########################################
  }# <-- end


