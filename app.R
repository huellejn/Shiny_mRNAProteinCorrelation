library(SummarizedExperiment)
library(shiny)
library(shinydashboard)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr) # scatter plots
#library(plotly)

# Load data
source("config.R")


rnaFPKMlog <- log2(assay_rna+0.1)

rowdata <- rowdata %>% arrange(HGNCApprovedSymbol)
coldata <- coldata %>% arrange(ProjectID, PatientID, Proteomics.SampleID, Batch, TMT.Channel)

# Function to calculate correlation coefficient
lm_eqn <- function(df, colx, coly){
  x <- df[[colx]]
  y <- df[[coly]]
  m <- lm(y ~ x, df);
  eq <- substitute(italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


# UI
ui <- dashboardPage(
  
  dashboardHeader(title = "RNA and protein expression in the MASTER cohort", titleWidth = '100%'),
  
  dashboardSidebar(collapsed = T),
  
  dashboardBody(
    
    # Hide the sidebar
    tags$script("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';"),
    
    fluidRow(
      
      box(
        title = "Correlation plot", status = "primary", height = "626px", width = 5,
        plotOutput("plotCorrelation", height = "570px")
      ),
      
      box(
        title = "Genes", status = "success", height = "626px", width = 7,
        DT::dataTableOutput("tableRowdata")
      )
      
    ),
    
    fluidRow(
      
      box(title = "Samples", status = "warning", width = 8,
          DT::dataTableOutput("tableColdata")
      ),
      
      box(title = "Expression values", width = 4,
          DT::dataTableOutput("tableExpression"))
    )
    
  )
)



# Server
server <- function(input, output, session) {
  
  # Table rowdata
  output$tableRowdata <- DT::renderDataTable({
    
    rowdata %>%
      select(`HGNC Symbol` = HGNCApprovedSymbol, `RNASeq Symbol` = Gene, `Proteomics Symbol` = Protein, `Ensembl Gene ID` = EnsemblGeneID) 
    
  }, 
  selection = "single", 
  rownames = FALSE, 
  filter = 'top',
  options = 
    list(
      scrollX = T,
      autoWidth = T
      ))
  
  
  # Get selected gene (HGNC)
  rct_rowdata_selected <- reactive({
    
    rowdata[input$tableRowdata_rows_selected,]
    
  })

  
  # Get RNA and protein expression of selected gene
  rct_gene_selected <- reactive({
    rct_rowdata_selected()$RNA.EnsemblGeneID
  })
  
  rct_protein_selected <- reactive({
    rct_rowdata_selected()$Protein
  })
  
  rct_exp <- reactive({

    gene_selected <- rct_gene_selected()
    protein_selected <- rct_protein_selected()
    
    # Subset RNAseq data for selected gene
    sub_rna <- rnaFPKMlog[gene_selected,]

    # Subset proteomics assay for selected gene
    sub_protein <- assay_proteomics[protein_selected, ]

    # Combine data
    df <- coldata
    df$RNA.Expression <- sub_rna[match(df$RNA.SampleID,names(sub_rna))]
    df$Protein.Expression <- sub_protein[match(df$Proteomics.SampleID,names(sub_protein))]
    
    # Add sample selection status for coloring
    df$selected_samples <- "not selected"
    
    # Sort table
    df <- arrange(df, Sample)

    return(df)
    
  })
  
  
  # Plot correlation
  output$plotCorrelation <- renderPlot({
    
    # Draw plot if gene is selected
    req(input$tableRowdata_rows_selected)
    
    row_selected <- rct_rowdata_selected()
    df <- rct_exp()
    
    # Calculate position of the lm coefficient label
    df_noNA <- df[!is.na(df$Protein.Expression) & ! is.na(df$RNA.Expression),]
    maxx <- max(df_noNA$RNA.Expression, na.rm=T)
    minx <- min(df_noNA$RNA.Expression, na.rm=T)
    rangex <- abs(maxx - minx)
    rangex10 <- rangex*0.1
    
    maxy <- max(df_noNA$Protein.Expression, na.rm=T)
    miny <- min(df_noNA$Protein.Expression, na.rm=T)
    rangey <- abs(maxy - miny)
    rangey10 <- rangey*0.1
    
    # Generate the plot
    p <- ggplot(df_noNA, aes(RNA.Expression, Protein.Expression)) +
      geom_point(color = "black") +
      theme_bw() +
      geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
      labs(title = row_selected$HGNCApprovedSymbol, x = "RNA Expression [log2(FPKM+0.1)]", y = "Protein Expression [abundance]") +
      geom_text(x = maxx-rangex10, y = maxy-rangey10, label = lm_eqn(df_noNA, colx="RNA.Expression", coly = "Protein.Expression"), parse = TRUE) 
    
    # Highlight data for selected samples
    if(length(input$tableColdata_rows_selected) > 0) {
      
      idx <- input$tableColdata_rows_selected
      samples_selected <- coldata[idx,"SampleID"]
      
      df_noNA$selected_sample[df_noNA$SampleID %in% samples_selected] <- "selected"
      
      p <- p + 
        geom_point(data = df_noNA[df_noNA$selected_sample=="selected",], color = "orange")
      
    }

    #return(ggplotly(p=p))
    return(p)
      
    })
  
  
  # Sample table
  output$tableColdata <- DT::renderDataTable({
    
    coldata %>%
      mutate(`Proteomics batch` = ifelse(is.na(Batch), NA, paste(Batch, TMT.Channel, sep = "/ch")) ) %>%
      select(Project = ProjectID, Patient = PatientID, `ICD-O-3 Topology` = ICD.O3.Topology, `ICD-O-3 Morphology` = ICD.O3.Morphology, Sample, `RNA sample name` = RNA.SampleID, `Protein sample name` = Proteomics.SampleID, `Proteomics batch` )
    
  }, 
  selection = "multiple", 
  rownames = FALSE, 
  filter = 'top',  
  options = list(
    scrollX = T,
    lengthMenu = c(10, 25, 50, 100, 150),
    pageLength = 150
    ))
  
  
  # Expression table
  output$tableExpression <- DT::renderDataTable({
    
    # Draw plot if gene is selected
    req(input$tableRowdata_rows_selected)
    
    tbl  <- rct_exp()
    
    tbl <- tbl %>%
      mutate(
        # Round expression values
        RNA.Expression = round(RNA.Expression, 4),
        Protein.Expression = round(Protein.Expression, 4),
      ) %>%
      select(`Sample ID` = SampleID, `RNA expression [log2(FPKM)]` = "RNA.Expression", `Protein expression [Intensity]` = "Protein.Expression")
    
    if(length(input$tableColdata_rows_selected) > 0) {
      idx <- input$tableColdata_rows_selected
      samples_selected <- coldata[idx,"SampleID"]
      tbl <- filter(tbl, `Sample ID` %in% samples_selected)
    }
    
    return(tbl)
    
  }, 
  rownames = FALSE, 
  filter = 'top', options = 
    list(
      scrollX = TRUE,
      lengthMenu = c(10, 25, 50, 100, 150),
      pageLength = 150
      ))
  
}


shinyApp(ui, server)
