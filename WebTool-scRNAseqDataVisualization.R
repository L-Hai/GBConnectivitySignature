library(shiny)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(Matrix)

choice.gene <- readRDS('choice.gene.rds')

meta.pdx <- readRDS('PDX.Meta.rds')
Matrix.PDX <- readRDS('PDX.Counts.rds')

meta.patient <- readRDS('Patient.Meta.rds')
Matrix.Patient <- readRDS('Patient.Counts.rds')

ind <- meta.patient$CellType == 'Malignant'
meta.malignant <- meta.patient[ind,]
# Matrix.Malignant <- Matrix.Patient[,ind]

ui <- navbarPage(
    "GBMConnectivity",
    tabPanel("Introduction",
             mainPanel(
               br(),
               p("Welcome to the interactive webtool for the connectivity in glioblastoma!"),
               p(
                 "It enables users to visualize and further explore the scRNA-seq datasets from paper:"
               ),
               p(strong("A clinically applicable connectivity signature for glioblastoma includes the tumor network driver CHI3L1")),
               br(),
               p("There are two datasets available:"),
               p(
                 strong("PDGCL: "),
                 "35K cells from three patient derived glioma cell line (PDGCL) mouse models. The cells from each PDGCL mouse was FACS sorted by SR101 to seperate TM-connected and TM-unconnected cells."
               ),
               p(
                 strong("Patient: "),
                 "213K cells from 21 glioblastoma patient samples, including 172K maligant cells and 41K nonmalignant cells."
               ),
               p("Please go to next tab for further exploring data!"),
               br(),
               p("The data is available in EGA database under the accession number EGAS00001007611.")
             )),
    tabPanel("Metadata",
             mainPanel(
               tabsetPanel(
                 tabPanel(
                   "PDGCL",
                   selectInput(
                     "Meta.PDGCL",
                     "Select value",
                     c(
                       "ConnectivityScore_SingleCell",
                       "ConnectivityScore_Bulk",
                       "PDGCL",
                       "Group",
                       "CellState"
                     ),
                     selected = "ConnectivityScore_SingleCell"
                   ),
                   plotOutput("UMAP.PDGCL")
                 ),
                 tabPanel(
                   "PatientMalignant",
                   selectInput(
                     "Meta.Malignant",
                     "Select value",
                     c(
                       "ConnectivityScore_SingleCell",
                       "ConnectivityScore_Bulk",
                       "Sample",
                       "CellState"
                     ),
                     selected = "ConnectivityScore_SingleCell"
                   ),
                   plotOutput("UMAP.Malignant")
                 ),
                 tabPanel(
                   "Patient",
                   selectInput(
                     "Meta.Patient",
                     "Select value",
                     c("Sample", "CellType", "Cluster"),
                     selected = "CellType"
                   ),
                   plotOutput("UMAP.Patient")
                 )
               )
             )),
    tabPanel(
      "GeneInPDGCL",
      sidebarPanel(
        selectInput("Gene.PDGCL", "Select gene", choice.gene, selected = "CHI3L1"),
        radioButtons(
          "Zero",
          "Zero expression cells:",
          c("Keep" = FALSE,
            "Remove" = TRUE),
          selected = FALSE
        ),
        actionButton('B.G.PDGCL', ' Plot '),
      ),
      mainPanel(tabsetPanel(
        tabPanel("UMAP",
                 plotOutput("UMAPplotSep")),
        tabPanel("InSortedGroup",
                 fluidRow(
                   column(6,plotOutput("Boxplot")),
                   column(6,plotOutput("Barplot"))
                 )),
        tabPanel("InCellState",
                 plotOutput("Boxplot.CellState")),
        tabPanel(
          "Correlation",
          radioButtons(
            "scatBulk",
            "Correlate with connectivity score derived in:",
            c("Single Cell" = FALSE,
              "Bulk" = TRUE),
            selected = FALSE
          ),
          plotOutput("Scatterplot")
        )
      ))
    ),
    tabPanel(
      "GeneInPatientMalignant",
      sidebarPanel(
        selectInput("Gene.Malignant", "Select gene", choice.gene, selected = "CHI3L1"),
        radioButtons(
          "Zero.Malignant",
          "Zero expression cells:",
          c("Keep" = FALSE,
            "Remove" = TRUE),
          selected = FALSE
        ),
        actionButton('B.G.Malignant', ' Plot ')
      ),
      mainPanel(tabsetPanel(
        tabPanel("UMAP",
                 plotOutput("UMAP.Gene.Maligant")),
        tabPanel(
          "InConnectivityScoreGroup",
          fluidRow(
            column(6,plotOutput("Boxplot.Group.Malignant")),
            column(6,plotOutput("Barplot.Group.Malignant"))
          )
        ),
        tabPanel("InCellState",
                 plotOutput("Boxplot.Malignant")),
        tabPanel(
          "Correlation",
          radioButtons(
            "scatBulk.Malignant",
            "Correlate with connectivity score derived in:",
            c("Single Cell" = FALSE,
              "Bulk" = TRUE),
            selected = FALSE
          ),
          plotOutput("Scatterplot.Malignant")
        )
      ))
    ),
    tabPanel(
      "GeneInPatient",
      sidebarPanel(
        selectInput("Gene.Patient", "Select gene", choice.gene, selected = "EGFR"),
        radioButtons(
          "Zero.Patient",
          "Zero expression cells:",
          c("Keep" = FALSE,
            "Remove" = TRUE),
          selected = FALSE
        ),
        actionButton('B.G.Patient', ' Plot ')
      ),
      mainPanel(tabsetPanel(
        tabPanel("UMAP",
                 plotOutput("UMAP.Gene.Patient")),
        tabPanel("InCellType",
                 plotOutput("Boxplot.Patient"))
      ))
    )
 )

server <- function(input, output, session) {
  
  output$UMAP.PDGCL <- renderPlot({
    df <- as.data.frame(meta.pdx[, c('UMAP_1', 'UMAP_2')])
    df$Legend <- meta.pdx[, input$Meta.PDGCL]
    p <-
      ggplot(df,
             aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Legend'))
    if (input$Meta.PDGCL == 'ConnectivityScore_Bulk' |
        input$Meta.PDGCL == 'ConnectivityScore_SingleCell') {
      p <-
        p + scale_color_gradientn(colours = rev(brewer.pal(11, 'RdBu')))
    } else{
      colorn <- length(unique(df$Legend))
      if (colorn <= 9) {
        colorcode <- brewer.pal(colorn, "Set1")
        p <- p + scale_color_manual(values = colorcode)
      } else{
        colorcode <- colorRampPalette(brewer.pal(9, "Set1"))(colorn)
        p <- p + scale_color_manual(values = colorcode)
      }
    }
    p <-
      p + geom_point(size = 0.05) + theme_classic() + guides(colour = guide_legend(override.aes = list(size =
                                                                                                         2))) + coord_fixed()
    return(p)
  })
  
  output$UMAP.Malignant <- renderPlot({
    df <- as.data.frame(meta.malignant[, c('UMAP_1_Mal', 'UMAP_2_Mal')])
    df$Legend <- meta.malignant[, input$Meta.Malignant]
    p <-
      ggplot(df,
             aes_string(x = 'UMAP_1_Mal', y = 'UMAP_2_Mal', color = 'Legend'))
    if (input$Meta.Malignant == 'ConnectivityScore_Bulk' |
        input$Meta.Malignant == 'ConnectivityScore_SingleCell') {
      p <-
        p + scale_color_gradientn(colours = rev(brewer.pal(11, 'RdBu')))
    } else{
      colorn <- length(unique(df$Legend))
      if (colorn <= 9) {
        colorcode <- brewer.pal(colorn, "Set1")
        p <- p + scale_color_manual(values = colorcode)
      } else{
        colorcode <- colorRampPalette(brewer.pal(9, "Set1"))(colorn)
        p <- p + scale_color_manual(values = colorcode)
      }
    }
    p <-
      p + geom_point(size = 0.05) + theme_classic() + guides(colour = guide_legend(override.aes = list(size =
                                                                                                         2))) + coord_fixed()
    return(p)
  })
  
  output$UMAP.Patient <- renderPlot({
    df <- as.data.frame(meta.patient[, c('UMAP_1', 'UMAP_2')])
    df$Legend <- meta.patient[, input$Meta.Patient]
    p <-
      ggplot(df,
             aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Legend'))
    colorn <- length(unique(df$Legend))
    if (colorn <= 9) {
      colorcode <- brewer.pal(colorn, "Set1")
      p <- p + scale_color_manual(values = colorcode)
    } else{
      colorcode <- colorRampPalette(brewer.pal(9, "Set1"))(colorn)
      p <- p + scale_color_manual(values = colorcode)
    }
    p <-
      p + geom_point(size = 0.05) + theme_classic() + guides(colour = guide_legend(override.aes = list(size =
                                                                                                         2))) + coord_fixed()
    p
  })
  
  observeEvent(input$B.G.PDGCL, {
    df <- data.frame(
      Gene = Matrix.PDX[input$Gene.PDGCL,],
      Group = meta.pdx$Group,
      ConnectivityScore_Bulk = meta.pdx$ConnectivityScore_Bulk,
      ConnectivityScore_SingleCell = meta.pdx$ConnectivityScore_SingleCell,
      CellState = meta.pdx$CellState
    )
    dfU <- meta.pdx
    dfU$Gene <- Matrix.PDX[input$Gene.PDGCL,]
    if (input$Zero) {
      df <- df[df$Gene > 0,]
    }
    box <-
      ggboxplot(df,
                x = "Group",
                y = "Gene",
                color = "Group") + stat_compare_means(label.y = 6) + coord_fixed(1 / 2)
    dfa <- data.frame(Group=c('Con','Uncon'),
                      Mean=c(mean(df$Gene[df$Group == 'Connected']),
                             mean(df$Gene[df$Group == 'Unconnected'])))
    bar <-
      ggbarplot(dfa,
                x = "Group",
                y = "Mean",
                color = "Group") + coord_fixed(2)
    gg <- ggplot(dfU, aes(x = UMAP_1, y = UMAP_2, color = Gene)) +
      geom_point(size = 0.1) +
      scale_color_gradientn(colours = brewer.pal(9, 'Oranges')) +
      theme_classic() +
      coord_fixed() +
      facet_wrap( ~ Group)
    
    output$UMAPplotSep <- renderPlot({
      return(gg)
    })
    output$Boxplot <- renderPlot({
      box
    })
    output$Barplot <- renderPlot({
      bar
    })
    
    output$Boxplot.CellState <- renderPlot({
      box <-
        ggboxplot(df,
                  x = "CellState",
                  y = "Gene",
                  color = "CellState") + coord_fixed(1 / 2) + stat_compare_means(label.y = 6)
      return(box)
    })
    
    output$Scatterplot <- renderPlot({
      if (input$scatBulk) {
        scat <- ggscatter(
          df,
          x = "Gene",
          y = "ConnectivityScore_Bulk",
          add = "reg.line",
          conf.int = FALSE,
          cor.coef = TRUE,
          cor.method = "pearson",
          color = '#fdae6b',
          size = 0.01,
          add.params = list(color = "#e6550d")
        ) + coord_fixed()
      } else{
        scat <-
          ggscatter(
            df,
            x = "Gene",
            y = "ConnectivityScore_SingleCell",
            add = "reg.line",
            conf.int = FALSE,
            cor.coef = TRUE,
            cor.method = "pearson",
            color = '#fdae6b',
            size = 0.01,
            add.params = list(color = "#e6550d")
          ) + coord_fixed()
      }
      return(scat)
    })
  })
  
  observeEvent(input$B.G.Malignant, {
    df <- meta.malignant
    df$Gene <- Matrix.Patient[input$Gene.Malignant, ind]
    
    if (input$Zero.Malignant) {
      df <- df[df$Gene > 0,]
    }
    
    output$UMAP.Gene.Maligant <- renderPlot({
      p <-
        ggplot(df,
               aes_string(x = 'UMAP_1_Mal', y = 'UMAP_2_Mal', color = 'Gene'))
      p <-
        p + scale_color_gradientn(colours = brewer.pal(9, 'Oranges'))
      p <-
        p + geom_point(size = 0.05) + theme_classic() + guides(colour = guide_legend(override.aes = list(size =
                                                                                                           2))) + coord_fixed()
      return(p)
    })
    
    output$Boxplot.Malignant <- renderPlot({
      box <-
        ggboxplot(df,
                  x = "CellState",
                  y = "Gene",
                  color = "CellState") + coord_fixed(1 / 2) + stat_compare_means(label.y = 6)
      return(box)
    })
    
    output$Boxplot.Group.Malignant <- renderPlot({
      box <-
        ggboxplot(df,
                  x = "ConnectivityScore_Group",
                  y = "Gene",
                  color = "ConnectivityScore_Group") + coord_fixed(1 / 2) + stat_compare_means(label.y = 6)
      return(box)
    })
    output$Barplot.Group.Malignant <- renderPlot({
      dfa <- data.frame(Group=c('Q1','Q2','Q3','Q4'),
                        Mean=c(mean(df$Gene[df$ConnectivityScore_Group == 'Q1']),
                               mean(df$Gene[df$ConnectivityScore_Group == 'Q2']),
                               mean(df$Gene[df$ConnectivityScore_Group == 'Q3']),
                               mean(df$Gene[df$ConnectivityScore_Group == 'Q4'])))
      bar <-
        ggbarplot(dfa,
                  x = "Group",
                  y = "Mean",
                  color = "Group") + coord_fixed(2)
      return(bar)
    })
    
    output$Scatterplot.Malignant <- renderPlot({
      if (input$scatBulk.Malignant) {
        scat <- ggscatter(
          df,
          x = "Gene",
          y = "ConnectivityScore_Bulk",
          add = "reg.line",
          conf.int = FALSE,
          cor.coef = TRUE,
          cor.method = "pearson",
          color = '#fdae6b',
          size = 0.01,
          add.params = list(color = "#e6550d")
        ) + coord_fixed()
      } else{
        scat <-
          ggscatter(
            df,
            x = "Gene",
            y = "ConnectivityScore_SingleCell",
            add = "reg.line",
            conf.int = FALSE,
            cor.coef = TRUE,
            cor.method = "pearson",
            color = '#fdae6b',
            size = 0.01,
            add.params = list(color = "#e6550d")
          ) + coord_fixed()
      }
      return(scat)
    })
  })
  
  
  observeEvent(input$B.G.Patient, {
    df <- meta.patient
    df$Gene <- Matrix.Patient[input$Gene.Patient, ]
    
    if (input$Zero.Patient) {
      df <- df[df$Gene > 0,]
    }
    
    output$UMAP.Gene.Patient <- renderPlot({
      p <-
        ggplot(df,
               aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Gene'))
      p <-
        p + scale_color_gradientn(colours = brewer.pal(9, 'Oranges'))
      p <-
        p + geom_point(size = 0.05) + theme_classic() + guides(colour = guide_legend(override.aes = list(size =
                                                                                                           2))) + coord_fixed()
      return(p)
    })
    
    output$Boxplot.Patient <- renderPlot({
      box <-
        ggboxplot(df,
                  x = "CellType",
                  y = "Gene",
                  color = "CellType") + coord_fixed(1 / 2.5) + stat_compare_means(label.y = 6)
      return(box)
    })
    
  })
  #=========END
}
shinyApp(ui, server)
