library(shiny)
library(rJava)
library(rcdk)
library(fingerprint)
library(KernSmooth)
library(ggplot2)
library(pryr)
library(aqfig)
library(geoR)
library(Cairo)

shinyServer(function(input, output){
  

output$sasplot <- renderPlot({
  
  validate(need(input$file1 !="", "Introduce a comma (,) or tab delimited file. The last column in your file MUST have the IDs of your compounds."))
  
  inFile <- input$file1
  
  if (is.null(inFile))
    return(NULL)
  
  compounds <-  read.csv(inFile$datapath, header = TRUE)
  if(ncol(compounds) == 1){compounds <-  read.csv(inFile$datapath, sep="\t", header = TRUE)}
 
  if(ncol(compounds) == 6){
    colnames(compounds)[6] <- "ID"}
  else if(ncol(compounds) != 6){
    names(compounds)[length(names(compounds))]<-"ID"}
  
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
  rownames(sim.matrix) <- compounds$ID
  colnames(sim.matrix) <- compounds$ID
  
  smiles <- sim.matrix 
  rownames(smiles) <- compounds[,1]
  colnames(smiles) <- compounds[,1]
  
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
  
  rownames(act.matrix) <- compounds$ID
  colnames(act.matrix) <- compounds$ID
  rownames(max.values) <- compounds$ID
  colnames(max.values) <- compounds$ID
  
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
  smiles <- matrix2list(smiles)
  smiles[,1] <- row.names(smiles)
  row.names(smiles) <- row.names(sim.list)
  smiles <- t(apply(smiles, 1, function(x){unlist(strsplit(x, "_"))}))
  colnames(smiles) <- c("SMILES1", "SMILES2")
  
  sim.list <- as.matrix(merge(smiles, sim.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  sim.list <- as.matrix(merge(sim.list, act.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  sim.list <- as.matrix(merge(sim.list, max.value.list, by="row.names", all = TRUE))
  row.names(sim.list) <- sim.list[,1]
  sim.list <- sim.list[,-1]
  
  #CALCULATE SALI
  abs.difference <- abs(as.numeric(sim.list[,4]))
  distance <- 1-as.numeric(sim.list[,3])
  distance[distance == 0] <- 0.01
  SALI <- abs.difference/distance
  sim.list <- cbind(sim.list,SALI)
  
  #PLOT SAS MAP
  no.data <- which(is.na(sim.list[,4]))
  if(length(no.data) > 0){sas.data <- sim.list[-no.data,]} else {sas.data <- sim.list}
  sas.data3 <- as.data.frame(sas.data)
  sas.data3$Similarity <- as.numeric(as.character(sas.data3$Similarity))
  sas.data3$Activity.Difference <- as.numeric(as.character(sas.data3$Activity.Difference))
  sas.data3$SALI <- as.numeric(as.character(sas.data3$SALI))
  sas.data3$Max.Values <- as.numeric(as.character(sas.data3$Max.Values))
  sas.data4 <-sas.data3[,-1]
  sas.data4 <-sas.data4[,-1]
  
  
  color.gradient <- function(x, colors=c("#00FF00","#FFFF00","#FFE200",
                                         "#FFC600","#FFAA00","#FF8D00","#FF7100","#FF5500","#FF3800",
                                         "#FF1C00","#FF0000"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  
  color.gradient2 <- function(x, colors=c("green","yellow","red"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }
  
  dm <- as.matrix.data.frame(sas.data4)
  dcols <-densCols(x=dm[,1],y=dm[,2],
                   colramp = colorRampPalette(c("snow3","rosybrown1","red")))
  
  if(input$col =="SALI"){
    sl%<a-%(plot(sas.data4$Similarity, sas.data4$Activity.Difference, main="SAS map",
                 xlab ="Structural Similarity", ylab="Activity Difference",
         col=color.gradient(sas.data4[,4]), pch = 20,cex=1.5,
         cex.axis=1.2,cex.lab=1.3,cex.main=1.5) +
    abline(v=ABX, h=ABY, lty=2))
    par(mar=c(5,6,5,5)+0.1)
    sl
    par(ps=12,cex=1.7,cex.main=1)
    vertical.image.legend(zlim = range(sas.data4[,4]), col = color.gradient(sort(unique(round(sas.data4[,4],1)))))
  }
  else if(input$col =="Max.Activity"){
    m%<a-%(plot(sas.data4$Similarity, sas.data4$Activity.Difference, main="SAS map", 
                  xlab ="Structural Similarity", ylab="Activity Difference",
                  col=color.gradient2(sas.data4[,3]), pch = 20,cex=1.5,
                cex.axis=1.2,cex.lab=1.3,cex.main=1.5) +
             abline(v=ABX, h=ABY, lty=2))
    par(mar=c(5,6,5,5)+0.1)
    m
    par(ps=12,cex=1.6,cex.main=1)
    vertical.image.legend(zlim = range(sas.data4[,3]), col = color.gradient2(sort(unique(round(sas.data4[,3],1)))))
    }
  else if(input$col =="Density"){
    a%<a-%(plot(dm[,1],dm[,2],col=dcols,pch=20, 
                  main="Density SAS map",xlab="Structural Similarity",
                  ylab="Activity Difference",cex=1.5,
                  cex.axis=1.2,cex.lab=1.3,cex.main=1.5) +
               abline(v=ABX, h=ABY, lty=2, col="black"))
    a
    }
  
  output$down<-downloadHandler(
    filename = function(){
      paste("sasmap",'.tiff',sep='')
    },
    content = function(file){
      Cairo(file,type="tiff",width = 900,height = 800,units="px",dpi="auto")
      if (input$col=="SALI"){
        par(mar=c(5,6,5,5)+0.1)
        sl
        par(ps=12,cex=1.7,cex.main=1)
        vertical.image.legend(zlim = range(sas.data4[,4]), col = color.gradient(sort(unique(round(sas.data4[,4],1)))))
        
      }else if(input$col == "Max.Activity"){
        par(mar=c(5,6,5,5)+0.1)
          m
        par(ps=12,cex=1.6,cex.main=1)
        vertical.image.legend(zlim = range(sas.data4[,3]), col = color.gradient(sort(unique(round(sas.data4[,3],1)))))
        
      }else if(input$col == "Density"){
        a
        }
      dev.off()
    })

  SASdata <- as.data.frame(sas.data)
  SASdata$Similarity <- as.numeric(as.character(SASdata[,3]))
  SASdata$Activity.Difference <- as.numeric(as.character(SASdata[,4]))
  
  b <- subset(SASdata,  Similarity > ABX)
  d  <- subset(SASdata, Similarity < ABX)
  UNDX <- subset(SASdata, Similarity = ABX)
  UNDY <- subset(SASdata, Similarity = ABY)
  C1 <- subset(d, Activity.Difference > ABY )
  R.HOPPING <- subset(d, Activity.Difference < ABY )
  ACTIV.CLIFF <- subset(b, Activity.Difference > ABY )
  SM.SAR <- subset(b, Activity.Difference < ABY )
  
  output$brush_info <- renderPrint({
    brushedPoints(sas.data4, input$plot_brush, xvar = "Similarity", yvar = "Activity.Difference")
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
  
})

##################################################################################
#################################    DAD MAP    ##################################
##################################################################################

output$DADmap <-renderPlot({
  
  validate(need(input$file2 !="", "Introduce a comma (,) or tab delimited file. The last column in your file MUST have the IDs of your compounds."))

  
  inFile2 <- input$file2
  
  if (is.null(inFile2))
    return(NULL)
  
  compounds2 <-  read.csv(inFile2$datapath, header = TRUE)
  if(ncol(compounds2) == 1){compounds2 <-  read.csv(inFile2$datapath, sep="\t", header = TRUE)}
   
  if(ncol(compounds2) == 6){
    colnames(compounds2)[6] <- "ID"}
  else if(ncol(compounds2) != 6){
    names(compounds)[length(names(compounds2))]<-"ID"}
  
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
  rownames(sim.matrix2) <- compounds2$ID
  colnames(sim.matrix2) <- compounds2$ID
  
  smiles2 <- sim.matrix2 
  rownames(smiles2) <- compounds2[,1]
  colnames(smiles2) <- compounds2[,1]
  

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
  
  rownames(act.matrix.1) <- compounds2$ID
  colnames(act.matrix.1) <- compounds2$ID
  rownames(act.matrix.2) <- compounds2$ID
  colnames(act.matrix.2) <- compounds2$ID
  rownames(diferen2) <- compounds2$ID
  colnames(diferen2) <- compounds2$ID
  
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
  smiles2 <- matrix2list(smiles2)
  smiles2[,1] <- row.names(smiles2)
  row.names(smiles2) <- row.names(sim.list2)
  smiles2 <- t(apply(smiles2, 1, function(x){unlist(strsplit(x, "_"))}))
  colnames(smiles2) <- c("SMILES1", "SMILES2")
  
  sim.list2 <- as.matrix(merge(smiles2, sim.list2, by="row.names", all = TRUE))
  row.names(sim.list2) <- sim.list2[,1]
  sim.list2 <- sim.list2[,-1]
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
  no.data.i <- which(is.na(sim.list2[,4]))
  no.data.j <- which(is.na(sim.list2[,5]))
  no.data <- union(no.data.i,no.data.j)
  if(length(no.data) > 0){dad.data <- sim.list2[-no.data,]} else {dad.data <- sim.list2}
  dad.data <- as.data.frame(dad.data)
  dad.data$Activity.Difference.1 <- as.numeric(as.character(dad.data$Activity.Difference.1))
  dad.data$Activity.Difference.2 <- as.numeric(as.character(dad.data$Activity.Difference.2))
  dad.data$Similarity <- as.numeric(as.character(dad.data$Similarity))
  dad.data$Selectivity <- as.numeric(as.character(dad.data$Selectivity))
  dad.data2 <-dad.data[,-1]
  dad.data2 <-dad.data2[,-1]
  
  color.gradient <- function(x, colors=c("green","yellow","red"), colsteps=100) {
    return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  }

  dmd <- as.matrix.data.frame(dad.data2)
  dcolsd <-densCols(x=dmd[,2],y=dmd[,3],
                   colramp = colorRampPalette(c("snow3","rosybrown1","red")))
  
  if(input$colp =="Similarity"){
    sm%<a-%(plot(dad.data2$Activity.Difference.1, dad.data2$Activity.Difference.2, main="DAD map", 
                 xlab =(colnames(dad.data2)[2]), ylab=(colnames(dad.data2)[3]),
                   col=color.gradient(dad.data2[,1]),pch = 20,cex= 1.5,cex.axis=1.2,cex.lab=1.2,cex.main=1.5) +
              abline(v=c(MDADX,MAXDADx), h=c(MDADY,MAXDADY), lty=2)) 
    par(mar=c(5,6,5,5)+0.1)
    sm
    par(ps=12,cex=1.6,cex.main=1)
    vertical.image.legend(zlim = range(dad.data2[,1]), col = color.gradient(sort(unique(round(dad.data2[,1],1)))))
    
  }
  else if(input$colp =="Selectivity"){
    se%<a-%(plot(dad.data2$Activity.Difference.1, dad.data2$Activity.Difference.2, main="DAD map", 
                   xlab =(colnames(dad.data2)[2]), ylab=(colnames(dad.data2)[3]),
                   col=color.gradient(dad.data2[,4]),pch = 20,cex= 1.5,cex.axis=1.2,cex.lab=1.2,cex.main=1.5) +
                abline(v=c(MDADX,MAXDADx), h=c(MDADY,MAXDADY), lty=2))
    par(mar=c(5,6,5,5)+0.1)
    se
    par(ps=12,cex=1.6,cex.main=1)
    vertical.image.legend(zlim = range(dad.data2[,4]), col = color.gradient(sort(unique(round(dad.data2[,4],1)))))
  }
  else if(input$colp =="Density"){
    ds%<a-%(plot(dmd[,2],dmd[,3], col = dcolsd, pch = 20, 
                  main = "Density DAD map", xlab = "Activity Difference 1",
                  ylab = "Activity Difference 2",cex= 1.5,cex.axis=1.2,cex.lab=1.2,cex.main=1.5) +
             abline(v=c(MDADX,MAXDADx), h=c(MDADY,MAXDADY), lty=2))
    ds
  }
  
  output$down3<-downloadHandler(
    filename = function(){
      paste("DADmap",'.tiff',sep='')
    },
    content = function(file){
      Cairo(file,type="tiff",width = 900,height = 800,units="px",dpi="auto")
      if (input$colp=="Similarity"){
        par(mar=c(5,6,5,5)+0.1)
        sm
        par(ps=12,cex=1.6,cex.main=1)
        vertical.image.legend(zlim = range(dad.data2[,1]), col = color.gradient(sort(unique(round(dad.data2[,1],1)))))
        
      } else if(input$colp == "Selectivity"){
        par(mar=c(5,6,5,5)+0.1)
          se
        par(ps=12,cex=1.6,cex.main=1)
        vertical.image.legend(zlim = range(dad.data2[,4]), col = color.gradient(sort(unique(round(dad.data2[,4],1)))))
        
      } else if(input$colp == "Density"){
          ds
        }
      dev.off()
    })
  
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
    brushedPoints(dad.data2, input$plot_brush1, xvar = "Activity.Difference.1", yvar = "Activity.Difference.2")
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
  
})


})

