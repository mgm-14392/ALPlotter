library(shiny)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_111')
library(rJava)
library(rcdk)
library(fingerprint)

shinyServer(function(input, output){
  

output$sasplot <- renderPlot({
  inFile <- input$file1
  
  if (is.null(inFile))
    return(NULL)
  
  compounds <-  read.csv(inFile$datapath, header = TRUE)
  mols <- parse.smiles(as.vector(compounds[,1]))
  
  rad <- as.numeric(input$dep)
  ABX <- as.numeric(input$ablineX)
  ABY <- as.numeric(input$ablineY)
  ACTIV <- as.numeric(input$act)
  
  
  if(input$FP =="ECFP"){fps <- lapply(mols, get.fingerprint,
  type='circular', fp.mode='bit', 
  depth=rad)} else if(input$FP =="Pubchem"){fps <- lapply(mols,
  get.fingerprint, type='pubchem', fp.mode='bit',
  depth=rad)} else if(input$FP == "MACCS"){fps <- lapply(mols,
  get.fingerprint, type='maccs', fp.mode='bit',
  depth=rad)}
  
  sim.matrix <- fp.sim.matrix(fps, method='tanimoto')
  rownames(sim.matrix) <- rownames(compounds)
  colnames(sim.matrix) <- rownames(compounds)
  
  
  activity.column.1 <- ACTIV
  act.matrix <- matrix(nrow=nrow(compounds), ncol=nrow(compounds), 
                       dimnames=list(rownames(compounds),rownames(compounds)))
  max.values <- matrix(nrow=nrow(compounds), ncol=nrow(compounds), 
                       dimnames=list(rownames(compounds),rownames(compounds)))
  for(i in 1:nrow(compounds)){
    act.matrix[,i] <- abs(compounds[,activity.column.1] - compounds[i,activity.column.1])
    t <- cbind(compounds[,activity.column.1], rep(compounds[i,activity.column.1], nrow(compounds)))
    max.values[,i] <- apply(t, 1, max)
  }
  
  #CONVERT MATRIX TO LIST
  matrix2list <- function(matrix) {
    matrix[!upper.tri(matrix)] = NA
    matrix = as.data.frame(as.table(matrix))
    matrix = as.matrix(matrix[!is.na(matrix[,3]),])
    row.names(matrix) = paste(matrix[,1], matrix[,2], sep="_")
    matrix = as.matrix(matrix[,-1])
    matrix = as.matrix(matrix[,-1])
    return(matrix)
  }
  
  sim.list <- matrix2list(sim.matrix)
  colnames(sim.list) <- c("Similarity")
  act.list <- matrix2list(act.matrix)
  colnames(act.list) <- c("Activity.Difference")
  max.value.list <- matrix2list(max.values)
  colnames(max.value.list) <- c("Max.Values")
  
  sim.list <- as.matrix(merge(sim.list, act.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  sim.list <- as.matrix(merge(sim.list, max.value.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  
  #CALCULATE SALI
  abs.difference <- abs(as.numeric(sim.list[,2]))
  distance <- 1-as.numeric(sim.list[,1])
  distance[distance == 0] <- 0.0001
  SALI <- abs.difference/distance
  sim.list <- cbind(sim.list,SALI)
  
  #PLOT SAS MAP
  no.data <- which(is.na(sim.list[,2]))
  if(length(no.data) > 0){sas.data <- sim.list[-no.data,]} else {sas.data <- sim.list}
  sas.data3 <- as.data.frame(sas.data)
  sas.data3$Similarity <- as.numeric(as.character(sas.data3$Similarity))
  sas.data3$Activity.Difference <- as.numeric(as.character(sas.data3$Activity.Difference))
  sas.data3$SALI <- as.numeric(as.character(sas.data3$SALI))
  sas.data3$Max.Values <- as.numeric(as.character(sas.data3$Max.Values))
  
  
  if(input$col =="SALI"){COL1 <- sas.data3[,4]} else if(input$col =="Max.Activity"){COL1 <- sas.data3[,3] }
  
  
   color.gradient <- function(x, colors=c("green","yellow","red"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  
  plot(sas.data3$Similarity, sas.data3$Activity.Difference, main="SAS map", xlab ="Structural Similarity", ylab="Activity Difference",
       col=color.gradient(COL1), pch = 20, cex=1.5)
  abline(v=ABX, h=ABY, lty=2)

  
  
  SASdata <- as.data.frame(sim.list)
  SASdata$Similarity <- as.numeric(as.character(SASdata[,1]))
  SASdata$Activity.Difference <- as.numeric(as.character(SASdata[,2]))
  
  b <- subset(SASdata,  Similarity > ABX)
  d  <- subset(SASdata, Similarity < ABX)
  UNDX <- subset(SASdata, Similarity = ABX)
  UNDY <- subset(SASdata, Similarity = ABY)
  C1 <- subset(d, Activity.Difference > ABY )
  R.HOPPING <- subset(d, Activity.Difference < ABY )
  ACTIV.CLIFF <- subset(b, Activity.Difference > ABY )
  SM.SAR <- subset(b, Activity.Difference < ABY )
  
  output$brush_info <- renderPrint({
    brushedPoints(sas.data3, input$plot_brush, xvar = "Similarity", yvar = "Activity.Difference")
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste("SASdata", '.csv', sep='') },
    content = function(file) {
      write.csv(SASdata, file)
    })
  
  output$ACdata <- downloadHandler(
    filename = function() { paste("ActivityCliffdata", '.csv', sep='') },
    content = function(file) {
      write.csv(ACTIV.CLIFF, file)
    })
  
  output$SSARdata <- downloadHandler(
    filename = function() { paste("SmoothSARdata", '.csv', sep='') },
    content = function(file) {
      write.csv(SM.SAR, file)
    })
  
  output$RHata <- downloadHandler(
    filename = function() { paste("Scaffoldhoppingdata", '.csv', sep='') },
    content = function(file) {
      write.csv(R.HOPPING, file)
    })
  
  output$C1data <- downloadHandler(
    filename = function() { paste("C1data", '.csv', sep='') },
    content = function(file) {
      write.csv(C1, file)
    })
  
  output$UX <- downloadHandler(
    filename = function() { paste("UNDXdata", '.csv', sep='') },
    content = function(file) {
      write.csv(UNDX, file)
    })
  
  output$UY <- downloadHandler(
    filename = function() { paste("UNDYdata", '.csv', sep='') },
    content = function(file) {
      write.csv(UNDY, file)
    })
  
  output$down <- downloadHandler(
    filename = function(){
    paste("sasmap", '.tiff', sep='')  
    },
    content = function(file){
      tiff(file, width= 900, height= 800, units = "px", compression = "lzw")
      plot(sas.data3$Similarity, sas.data3$Activity.Difference, main="SAS map", xlab ="Structural Similarity", ylab="Activity Difference",
           col=color.gradient(COL1), pch = 20, cex=1.5, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5)
      abline(v=ABX, h=ABY, lty=2)
      usr <- par("usr")
      text(usr[2], usr[3],"Generated with Activity Landscape Plotter, DIFACQUIM", cex = 0.9, col = "azure4", adj = c(1,0))
      dev.off()
    }
  )
  
})

################################################
output$sasplot2 <- renderPlot({
  inFile3 <- input$file3
  
  if (is.null(inFile3)) return(NULL)
  
  compounds3 <-  read.csv(inFile3$datapath, header = TRUE)
  mols3 <- parse.smiles(as.vector(compounds3[,1]))
  
  rad3 <- as.numeric(input$dep3)
  ABX2 <- as.numeric(input$ablineX2)
  ABY2 <- as.numeric(input$ablineY2)
  ACTIV4 <- as.numeric(input$act4)
  
  
  if(input$FP3 =="ECFP"){fps <- lapply(mols3, get.fingerprint,
                                      type='circular', fp.mode='bit', 
                                      depth=rad3)} else if(input$FP3 == "Pubchem"){fps <- lapply(mols3,
                                                                                                get.fingerprint, type='pubchem', fp.mode='bit',
                                                                                                depth=rad3)} else if(input$FP3 == "MACCS"){fps <- lapply(mols3,
                                                                                                                                                       get.fingerprint, type='maccs', fp.mode='bit',
                                                                                                                                                    depth=rad3)}
  
  sim.matrix <- fp.sim.matrix(fps, method='tanimoto')
  rownames(sim.matrix) <- rownames(compounds3)
  colnames(sim.matrix) <- rownames(compounds3)
  
  
  activity.column.1 <- ACTIV4
  act.matrix <- matrix(nrow=nrow(compounds3), ncol=nrow(compounds3), 
                       dimnames=list(rownames(compounds3),rownames(compounds3)))
  for(i in 1:nrow(compounds3)){
    act.matrix[,i] <- abs(compounds3[,activity.column.1] - compounds3[i,activity.column.1])
  }
  
  #CONVERT MATRIX TO LIST
  matrix2list <- function(matrix) {
    matrix[!upper.tri(matrix)] = NA
    matrix = as.data.frame(as.table(matrix))
    matrix = as.matrix(matrix[!is.na(matrix[,3]),])
    row.names(matrix) = paste(matrix[,1], matrix[,2], sep="_")
    matrix = as.matrix(matrix[,-1])
    matrix = as.matrix(matrix[,-1])
    return(matrix)
  }
  
  sim.list <- matrix2list(sim.matrix)
  colnames(sim.list) <- c("Similarity")
  act.list <- matrix2list(act.matrix)
  colnames(act.list) <- c("Activity.Difference")
  
  sim.list <- as.matrix(merge(sim.list, act.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  
  #PLOT SAS MAP
  no.data <- which(is.na(sim.list[,2]))
  if(length(no.data) > 0){sas.data <- sim.list[-no.data,]} else {sas.data <- sim.list}
  sas.data4 <- as.data.frame(sas.data)
  sas.data4$Similarity <- as.numeric(as.character(sas.data4$Similarity))
  sas.data4$Activity.Difference <- as.numeric(as.character(sas.data4$Activity.Difference))
 
  b <- as.matrix.data.frame(sas.data4)
  dcols <- densCols(b[,1],b[,2], 
                    colramp=colorRampPalette(c("snow3", "red")))
  plot(b[,1],b[,2], col = dcols, pch = 20, 
       main = "Density SAS map", xlab = "Structural similarity",
       ylab = "Activity Difference", cex= 1.5)
  abline(v=ABX2, h= ABY2, lty=2)
 
  
  b <- subset(sas.data4,  Similarity > ABX2)
  d  <- subset(sas.data4, Similarity < ABX2)
  UNDX <- subset(sas.data4, Similarity = ABX2)
  UNDY <- subset(sas.data4, Similarity = ABY2)
  C1 <- subset(d, Activity.Difference > ABY2)
  R.HOPPING <- subset(d, Activity.Difference < ABY2)
  ACTIV.CLIFF <- subset(b, Activity.Difference > ABY2)
  SM.SAR <- subset(b, Activity.Difference < ABY2)
  
  output$brush_info2 <- renderPrint({
    brushedPoints(sas.data4, input$plot_brush2, xvar = "Similarity", yvar = "Activity.Difference")
  })
  
  output$downloadData2 <- downloadHandler(
    filename = function() { paste("SASdata", '.csv', sep='') },
    content = function(file) {
      write.csv(sas.data4, file)
    })
  
  output$ACdata2 <- downloadHandler(
    filename = function() { paste("ActivityCliffdata", '.csv', sep='') },
    content = function(file) {
      write.csv(ACTIV.CLIFF, file)
    })
  
  output$SSARdata2 <- downloadHandler(
    filename = function() { paste("SmoothSARdata", '.csv', sep='') },
    content = function(file) {
      write.csv(SM.SAR, file)
    })
  
  output$RHata2 <- downloadHandler(
    filename = function() { paste("Scaffoldhoppingdata", '.csv', sep='') },
    content = function(file) {
      write.csv(R.HOPPING, file)
    })
  
  output$C1data2 <- downloadHandler(
    filename = function() { paste("C1data", '.csv', sep='') },
    content = function(file) {
      write.csv(C1, file)
    })
  
  output$UX2 <- downloadHandler(
    filename = function() { paste("UNDXdata", '.csv', sep='') },
    content = function(file) {
      write.csv(UNDX, file)
    })
  
  output$UY2 <- downloadHandler(
    filename = function() { paste("UNDYdata", '.csv', sep='') },
    content = function(file) {
      write.csv(UNDY, file)
    })
  
  output$down2 <- downloadHandler(
    filename = function(){
      paste("Densitysasmap", '.tiff', sep='')  
    },
    content = function(file){
      tiff(file, width= 900, height= 800, units = "px", compression = "lzw")
      b <- as.matrix.data.frame(sas.data4)
      dcols <- densCols(b[,1],b[,2], 
                        colramp=colorRampPalette(c("snow3", "red")))
      plot(b[,1],b[,2], col = dcols, pch = 20, 
           main = "Density SAS map", xlab = "Structural similarity",
           ylab = "Activity Difference", cex= 1.5,
           cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5)
      abline(v=ABX2, h= ABY2, lty=2)
      usr <- par("usr")
      text(usr[2], usr[3], "Generated with Activity Landscape Plotter, DIFACQUIM", cex = 0.9, col = "azure4", adj = c(1,0))
      dev.off()
    }
  )
  
})

#################################################################

output$DADmap <-renderPlot({
  
  inFile2 <- input$file2
  
  if (is.null(inFile2))
    return(NULL)
  
  compounds2 <-  read.csv(inFile2$datapath, header = TRUE)
  
  rad2 <- as.numeric(input$dep2)
  
  mols2 <- parse.smiles(as.vector(compounds2[,1]))
  if(input$FP2 =="ECFP"){fps2 <- lapply(mols2, get.fingerprint,
  type='circular', fp.mode='bit', 
  depth=rad2)} else if(input$FP2 == "Pubchem"){fps2 <- lapply(mols2,
  get.fingerprint, type='pubchem', fp.mode='bit',
  depth=rad2)} else if(input$FP2 == "MACCS"){fps2 <- lapply(mols2,
  get.fingerprint, type='maccs', fp.mode='bit',
  depth=rad2)}
  
  sim.matrix2 <- fp.sim.matrix(fps2, method='tanimoto')
  rownames(sim.matrix2) <- rownames(compounds2)
  colnames(sim.matrix2) <- rownames(compounds2)
  
  
  ACTIV1 <- as.numeric(input$act2)
  ACTIV2 <- as.numeric(input$act3)
  MDADX <- as.numeric(input$minlinex)
  MAXDADx <- as.numeric(input$maxlinx)
  MDADY <- as.numeric(input$minliney)
  MAXDADY <- as.numeric(input$maxliney)
  
  activity.column.1 <- ACTIV1
  activity.column.2 <- ACTIV2
  
  act.matrix.1 <- matrix(nrow=nrow(compounds2), ncol=nrow(compounds2), 
                         dimnames=list(rownames(compounds2),rownames(compounds2)))
  diferen2 <- matrix(nrow=nrow(compounds2), ncol=nrow(compounds2), 
                     dimnames=list(rownames(compounds2),rownames(compounds2)))
   for(i in 1:nrow(compounds2)){
    act.matrix.1[,i] <- compounds2[,activity.column.1] - compounds2[i,activity.column.1]
  }
  
  act.matrix.2 <- matrix(nrow=nrow(compounds2), ncol=nrow(compounds2), 
                         dimnames=list(rownames(compounds2),rownames(compounds2)))
  for(i in 1:nrow(compounds2)){
    act.matrix.2[,i] <- compounds2[,activity.column.2] - compounds2[i,activity.column.2]
  }
  
  
  diferen <- abs(compounds2[,activity.column.1] - compounds2[,activity.column.2])
  diferen <-as.data.frame(diferen)
  for(i in 1:nrow(diferen)){
    diferen2[,i] <- abs(diferen[,1] - diferen[i,1])
  }
  
  #CONVERT MATRIX TO LIST
  matrix2list <- function(matrix) {
    matrix[!upper.tri(matrix)] = NA
    matrix = as.data.frame(as.table(matrix))
    matrix = as.matrix(matrix[!is.na(matrix[,3]),])
    row.names(matrix) = paste(matrix[,1], matrix[,2], sep="_")
    matrix = as.matrix(matrix[,-1])
    matrix = as.matrix(matrix[,-1])
    return(matrix)
  }
  
  sim.list2 <- matrix2list(sim.matrix2)
  colnames(sim.list2) <- c("Similarity")
  act.list.1 <- matrix2list(act.matrix.1)
  colnames(act.list.1) <- c("Activity.Difference.1")
  act.list.2 <- matrix2list(act.matrix.2)
  colnames(act.list.2) <- c("Activity.Difference.2")
  selectivity.list <- matrix2list(diferen2)
  colnames(selectivity.list) <- c("Selectivity")
  
  
  sim.list2 <- as.matrix(merge(sim.list2, act.list.1, by="row.names", all = TRUE))
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
  sim.list2 <- merge(sim.list2, act.list.2, by="row.names", all = TRUE)
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
  sim.list2 <- as.matrix(merge(sim.list2, selectivity.list, by="row.names", all = TRUE))
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
  
  
  #PLOT DAD MAP
  no.data.i <- which(is.na(sim.list2[,2]))
  no.data.j <- which(is.na(sim.list2[,3]))
  no.data <- union(no.data.i,no.data.j)
  if(length(no.data) > 0){dad.data <- sim.list2[-no.data,]} else {dad.data <- sim.list2}
  dad.data <- as.data.frame(dad.data)
  dad.data$Activity.Difference.1 <- as.numeric(as.character(dad.data$Activity.Difference.1))
  dad.data$Activity.Difference.2 <- as.numeric(as.character(dad.data$Activity.Difference.2))
  dad.data$Similarity <- as.numeric(as.character(dad.data$Similarity))
  dad.data$Selectivity <- as.numeric(as.character(dad.data$Selectivity))
  
  if(input$colp =="Similarity"){COL3 <- dad.data[,1]} else if(input$colp =="Selectivity"){COL3 <- dad.data[,4] }
  
  color.gradient <- function(x, colors=c("green","yellow","red"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  plot(dad.data$Activity.Difference.1, dad.data$Activity.Difference.2, main="DAD map", xlab =(colnames(dad.data)[2]), ylab=(colnames(dad.data)[3]),
       col=color.gradient(COL3), pch = 20, cex= 1.5)
  abline(v=c(MDADX,MAXDADx), h=c(MDADY,MAXDADY), lty=2)

 
  
  Z1s <- subset(dad.data, Activity.Difference.2 > MAXDADY & Activity.Difference.1 > MAXDADx)
  Z2i <- subset(dad.data, Activity.Difference.2 < MDADY & Activity.Difference.1 > MAXDADx)
  Z2s <- subset(dad.data, Activity.Difference.2 > MAXDADY & Activity.Difference.1 < MDADX)
  Z1i <- subset(dad.data, Activity.Difference.2 < MDADY & Activity.Difference.1 < MDADX)
  Z3s <- subset(dad.data, Activity.Difference.2 > MAXDADY & Activity.Difference.1 < MAXDADx & Activity.Difference.1 > MDADX)
  Z3i <- subset(dad.data, Activity.Difference.1 < MAXDADx & Activity.Difference.1 > MDADX & Activity.Difference.2 < MDADX)
  Z4l <- subset(dad.data, Activity.Difference.2 < MAXDADY & Activity.Difference.2 > MDADY & Activity.Difference.1 < MDADX)
  Z4r <- subset(dad.data, Activity.Difference.2 < MAXDADY & Activity.Difference.2 > MDADY & Activity.Difference.1 > MAXDADx)
  Z5 <- subset(dad.data, Activity.Difference.2 < MAXDADY & Activity.Difference.2 > MDADY & Activity.Difference.1 < MAXDADx & Activity.Difference.1 > MDADX)
 
  output$brush_info1 <- renderPrint({
    brushedPoints(dad.data, input$plot_brush1, xvar = "Activity.Difference.1", yvar = "Activity.Difference.2")
  })
   
  output$downloadDADData <- downloadHandler(
    filename = function() { paste("DADdata", '.csv', sep='') },
    content = function(file) {
      write.csv(dad.data, file)
    })
  
  output$ZON1S <- downloadHandler(
    filename = function() { paste("Z1sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z1s, file)
    })
  
  output$ZON2I <- downloadHandler(
    filename = function() { paste("Z2idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z2i, file)
    })
  
  output$ZON2S <- downloadHandler(
    filename = function() { paste("Z2sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z2s, file)
    })
  
  output$ZON1I <- downloadHandler(
    filename = function() { paste("Z1idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z1i, file)
    })
  
  output$ZON3S <- downloadHandler(
    filename = function() { paste("Z3sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z3s, file)
    })
  
  output$ZON3I <- downloadHandler(
    filename = function() { paste("Z3idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z3i, file)
    })
  
  output$ZON4L <- downloadHandler(
    filename = function() { paste("Z4ldata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z4l, file)
    })
  
  output$ZON4R <- downloadHandler(
    filename = function() { paste("Z4rdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z4r, file)
    })
  
  output$ZON5 <- downloadHandler(
    filename = function() { paste("Z5data", '.csv', sep='') },
    content = function(file) {
      write.csv(Z5, file)
    })
  
  output$down3 <- downloadHandler(
    filename = function(){
      paste("DADmap", '.tiff', sep='')  
    },
    content = function(file){
      tiff(file, width= 900, height= 800, units = "px", compression = "lzw")
      
      color.gradient <- function(x, colors=c("green","yellow","red"), colsteps=100) {
        return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
      }
      plot(dad.data$Activity.Difference.1, dad.data$Activity.Difference.2, main="DAD map", xlab =(colnames(dad.data)[2]), ylab=(colnames(dad.data)[3]),
           col=color.gradient(COL3), pch = 20, cex= 1.5, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5)
      abline(v=c(MDADX,MAXDADx), h=c(MDADY,MAXDADY), lty=2)
      usr <- par("usr")
      text(usr[2], usr[3],"Generated with Activity Landscape Plotter, DIFACQUIM", cex = 0.9, col = "azure4", adj = c(1,0))
      dev.off()
    }
  )
  
})

#################################################################################################
 

output$DADmap2 <-renderPlot({
  
  inFile4 <- input$file4
  
  if (is.null(inFile4))
    return(NULL)
  
  compounds4 <-  read.csv(inFile4$datapath, header = TRUE)
  
  mols4 <- parse.smiles(as.vector(compounds4[,1]))
  fps2 <- lapply(mols4, get.fingerprint, type='maccs', fp.mode='bit')
  sim.matrix2 <- fp.sim.matrix(fps2, method='tanimoto')
  rownames(sim.matrix2) <- rownames(compounds4)
  colnames(sim.matrix2) <- rownames(compounds4)
  
  
  ACTIV5 <- as.numeric(input$act5)
  ACTIV6 <- as.numeric(input$act6)
  MDADX2 <- as.numeric(input$mlx)
  MAXDADx2 <- as.numeric(input$maxlx)
  MDADY2 <- as.numeric(input$mley)
  MAXDADY2 <- as.numeric(input$maxly)
  
  activity.column.1 <- ACTIV5
  activity.column.2 <- ACTIV6
  
  act.matrix.1 <- matrix(nrow=nrow(compounds4), ncol=nrow(compounds4), 
                         dimnames=list(rownames(compounds4),rownames(compounds4)))
  for(i in 1:nrow(compounds4)){
    act.matrix.1[,i] <- compounds4[,activity.column.1] - compounds4[i,activity.column.1]
  }
  
  act.matrix.2 <- matrix(nrow=nrow(compounds4), ncol=nrow(compounds4), 
                         dimnames=list(rownames(compounds4),rownames(compounds4)))
  for(i in 1:nrow(compounds4)){
    act.matrix.2[,i] <- compounds4[,activity.column.2] - compounds4[i,activity.column.2]
  }
  
  
  #CONVERT MATRIX TO LIST
  matrix2list <- function(matrix) {
    matrix[!upper.tri(matrix)] = NA
    matrix = as.data.frame(as.table(matrix))
    matrix = as.matrix(matrix[!is.na(matrix[,3]),])
    row.names(matrix) = paste(matrix[,1], matrix[,2], sep="_")
    matrix = as.matrix(matrix[,-1])
    matrix = as.matrix(matrix[,-1])
    return(matrix)
  }
  
  sim.list2 <- matrix2list(sim.matrix2)
  colnames(sim.list2) <- c("Similarity")
  act.list.1 <- matrix2list(act.matrix.1)
  colnames(act.list.1) <- c("Activity.Difference.1")
  act.list.2 <- matrix2list(act.matrix.2)
  colnames(act.list.2) <- c("Activity.Difference.2")
  
  
  sim.list2 <- as.matrix(merge(sim.list2, act.list.1, by="row.names", all = TRUE))
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
  sim.list2 <- merge(sim.list2, act.list.2, by="row.names", all = TRUE)
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
  
  
  #PLOT DAD MAP
  no.data.i <- which(is.na(sim.list2[,2]))
  no.data.j <- which(is.na(sim.list2[,3]))
  no.data <- union(no.data.i,no.data.j)
  if(length(no.data) > 0){dad.data2 <- sim.list2[-no.data,]} else {dad.data2 <- sim.list2}
  dad.data2$Activity.Difference.1 <- as.numeric(as.character(dad.data2$Activity.Difference.1))
  dad.data2$Activity.Difference.2 <- as.numeric(as.character(dad.data2$Activity.Difference.2))
  dad.data2 <- dad.data2[,-1]
  
  b <- as.matrix.data.frame(dad.data2)
  dcols <- densCols(b[,1],b[,2], 
                    colramp=colorRampPalette(c("snow3", "red")))
  plot(b[,1],b[,2], col = dcols, pch = 20, 
       main = "Density DAD map", xlab = "Activity Difference 1",
       ylab = "Activity Difference 2", cex= 1.5)
  abline(v=c(MDADX2,MAXDADx2), h=c(MDADY2,MAXDADY2), lty=2)

   
  Z1s <- subset(dad.data2, Activity.Difference.2 > MAXDADY2 & Activity.Difference.1 > MAXDADx2)
  Z2i <- subset(dad.data2, Activity.Difference.2 < MDADY2 & Activity.Difference.1 > MAXDADx2)
  Z2s <- subset(dad.data2, Activity.Difference.2 > MAXDADY2 & Activity.Difference.1 < MDADX2)
  Z1i <- subset(dad.data2, Activity.Difference.2 < MDADY2 & Activity.Difference.1 < MDADX2)
  Z3s <- subset(dad.data2, Activity.Difference.2 > MAXDADY2 & Activity.Difference.1 < MAXDADx2 & Activity.Difference.1 > MDADX2)
  Z3i <- subset(dad.data2, Activity.Difference.1 < MAXDADx2 & Activity.Difference.1 > MDADX2 & Activity.Difference.2 < MDADX2)
  Z4l <- subset(dad.data2, Activity.Difference.2 < MAXDADY2 & Activity.Difference.2 > MDADY2 & Activity.Difference.1 < MDADX2)
  Z4r <- subset(dad.data2, Activity.Difference.2 < MAXDADY2 & Activity.Difference.2 > MDADY2 & Activity.Difference.1 > MAXDADx2)
  Z5 <- subset(dad.data2, Activity.Difference.2 < MAXDADY2 & Activity.Difference.2 > MDADY2 & Activity.Difference.1 < MAXDADx2 & Activity.Difference.1 > MDADX2)
  
  output$brush_info3 <- renderPrint({
    brushedPoints(dad.data2, input$plot_brush3, xvar = "Activity.Difference.1", yvar = "Activity.Difference.2")
  })
  
  output$downloadDADData3 <- downloadHandler(
    filename = function() { paste("DADdata", '.csv', sep='') },
    content = function(file) {
      write.csv(dad.data2, file)
    })
  
  output$ZON1S2 <- downloadHandler(
    filename = function() { paste("Z1sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z1s, file)
    })
  
  output$ZON2I2 <- downloadHandler(
    filename = function() { paste("Z2idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z2i, file)
    })
  
  output$ZON2S2 <- downloadHandler(
    filename = function() { paste("Z2sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z2s, file)
    })
  
  output$ZON1I2 <- downloadHandler(
    filename = function() { paste("Z1idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z1i, file)
    })
  
  output$ZON3S2 <- downloadHandler(
    filename = function() { paste("Z3sdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z3s, file)
    })
  
  output$ZON3I2 <- downloadHandler(
    filename = function() { paste("Z3idata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z3i, file)
    })
  
  output$ZON4L2 <- downloadHandler(
    filename = function() { paste("Z4ldata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z4l, file)
    })
  
  output$ZON4R2 <- downloadHandler(
    filename = function() { paste("Z4rdata", '.csv', sep='') },
    content = function(file) {
      write.csv(Z4r, file)
    })
  
  output$ZON52 <- downloadHandler(
    filename = function() { paste("Z5data", '.csv', sep='') },
    content = function(file) {
      write.csv(Z5, file)
    })
  
  output$down4 <- downloadHandler(
    filename = function(){
      paste("DensityDADmap", '.tiff', sep='')  
    },
    content = function(file){
      tiff(file, width= 900, height= 800, units = "px", compression = "lzw")
      b <- as.matrix.data.frame(dad.data2)
      dcols <- densCols(b[,1],b[,2], 
                        colramp=colorRampPalette(c("snow3", "red")))
      plot(b[,1],b[,2], col = dcols, pch = 20, 
           main = "Density DAD map", xlab = "Activity Difference 1",
           ylab = "Activity Difference 2", cex= 1.5, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5)
      abline(v=c(MDADX2,MAXDADx2), h=c(MDADY2,MAXDADY2), lty=2)
      usr <- par("usr")
      text(usr[2], usr[3], "Generated with Activity Landscape Plotter, DIFACQUIM", cex = 0.9, col = "azure4", adj = c(1,0))
      
      dev.off()
    }
  )
  
})

})