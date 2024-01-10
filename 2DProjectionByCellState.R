# Codes for two-dimensional projection of patient malignant cells by cell state
# Can be run in R

library(ggplot2)

# Read rds file contained results of cell state signature scores (can be obtained using example codes in "MalignantCellStateAssignment.R"), cell state assignment (can be obtained using example codes in "MalignantCellStateAssignment.R") and connectivity scores (can be obtained using example codes in "ConnectivitySignatureScore.R").
df <- readRDS("Patient.MalignantCell.Meta.rds")

# Calculate X and Y axis based on cell state signature scores
xyl <- list()
for(i in 1:nrow(df)){
  am <- df[i,]
  ac <- max(am$MES1,am$MES2,am$AC)
  oc <- max(am$NPC1,am$NPC2,am$OPC)
  y <- ac - oc
  if(y > 0){
    x <- am$AC - max(am$MES1,am$MES2)
  }else{
    x <- am$OPC - max(am$NPC1,am$NPC2)
  }
  xy <- c(x,y)
  xyl[[i]] <- xy
}
xydf <- matrix(unlist(xyl),ncol=2,byrow=T)
xydf <- as.data.frame(xydf)
colnames(xydf) <- c('MEStoAC','OPCtoAC')
df$MEStoAC <- xydf$MEStoAC
df$OPCtoAC <- xydf$OPCtoAC

# Visualize cells in 2D plot
Color.ConnectivityScore <- c('#e66101','#5e3c99')
Color.CellState <- c('#ff7f00','#e31a1c','#fb9a99','#33a02c','#6a3d9a','#cab2d6','#1f78b4','#a6cee3')
names(Color.CellState) <- c('AC','MES1','MES2','OPC','NPC1','NPC2','G1_S','G2_M')

# Color by connectivity signature scores
ggplot(df, aes(x=MEStoAC, y=OPCtoAC, color=ConnectivityScore)) + 
  geom_point(size=0.01) +
  scale_color_gradient2(high = Color.ConnectivityScore[1],mid = '#e0e0e0',low=Color.ConnectivityScore[2],midpoint = 0) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white",colour = 'black'),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(breaks=c(0)) +
  scale_y_continuous(breaks=c(0)) +
  coord_fixed()


# Color by cell states
ggplot(df, aes(x=MEStoAC, y=OPCtoAC, color=CellState)) + 
  geom_point(size=0.01) +
  scale_color_manual(values = Color.CellState) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  theme(panel.background = element_rect(fill = "white",colour = 'black'),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(breaks=c(0)) +
  scale_y_continuous(breaks=c(0)) +
  coord_fixed()
