library(shiny)
library(colourpicker)
library(DT)
library(ggplot2)
library(igraph)
library(bslib)
library(gplots)
library(heatmaply)
library(shinycssloaders)
library(plotly)
library(scales)
library(tidyverse)
library(beeswarm)
library(vioplot)






ui <- fluidPage(
  titlePanel("BF591 Final Project - bsd112"),
  p(" Post-mortem Huntington’s Disease prefrontal cortex compared with neurologically healthy controls"),
  tabsetPanel(
    # Tab 1: Sample Information Exploration
    tabPanel("Sample Information",
             sidebarLayout(
               sidebarPanel(
                 fileInput("sample_file", "Upload Sample Information (CSV)"),
                 # Add other input controls as needed
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            # Summary of the table with distinct values
                            withSpinner(dataTableOutput("summary_table"))
                   ),
                   tabPanel("Table",
                            # Data table displaying sample information
                            withSpinner(DTOutput("sample_table"))
                   ),
                   tabPanel("Plots",
                            # No static plotOutput("hist_plot") here
                            uiOutput("hist_tabs")
                   )
                 )
               )
             )
    ),
    # Tab 2: Counts Matrix Exploration
    tabPanel("Counts Matrix",
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_file", "Upload Counts Matrix (CSV)"),
                 sliderInput("variance_threshold", "Percentile of Variance", min = 0, max = 100, value = 80),
                 sliderInput("nonzero_samples_threshold", "Samples with Non-Zero Counts", min = 0, max = 100, value = 60)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            # Summary of the filtering effects
                            withSpinner(DTOutput("filter_summary_table"))
                   ),
                   tabPanel("Diagnostic Plots",
                            # Diagnostic scatter plots
                            withSpinner(plotOutput("scatter_median_variance")),
                            withSpinner(plotOutput("scatter_median_zeros"))
                   ),
                   tabPanel("Heatmap",
                            # Heatmap of counts after filtering
                            withSpinner(plotOutput("clustered_heatmap"))
                   ),
                   tabPanel("PCA Plot",
                            selectInput("PCAX", "Select PCA for x axis",
                                        choices = paste0('PC', 1:15)),
                            selectInput("PCAY", "Select PCA for y axis",
                                        choices = paste0('PC', 1:15), selected = "PC2"),
                            withSpinner(plotOutput("CPCA"))
                   ),
                 )
               )
             )
    ),
    # Tab 3: Differential Expression
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 fileInput("de_file", "Upload Differential Expression Results (CSV)"),
                 
                 radioButtons("x_axis", "Choose the column for the x-axis",
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                              selected = "log2FoldChange"),
                 radioButtons("y_axis", "Choose the column for the y-axis",
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"),
                              selected = "padj"),
                 colourInput("base", "Base point color", value = "#22577A"),
                 colourInput("highlight", "Highlight point color", value = "#FFCF56"),
                 sliderInput("slider", "Select the magnitude of the p adjusted coloring:",
                             min = -100, max = 100, value = -10),
                 actionButton("plot_button", "Plot", icon = icon("thumbs-up"), style = "width:100%;"),
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Sortable table",
                            DTOutput("de_table")
                   ),
                   tabPanel("Plot",
                    
                            withSpinner(plotOutput("volcano"))  # Add this line for the volcano plot
                          
                   )
                 )
               )
             )
    ),
    # Tab 4: Individual Gene Expression(s)
    tabPanel("Individual Gene Expression(s)",
             fluidPage(
               # Input controls
               sidebarLayout(
                 sidebarPanel(
                   # Input for normalized counts matrix
                   fileInput("countsFile", "Choose Normalized Counts CSV File",
                             accept = c(".csv")),
                   
                   # Input for sample information matrix
                   fileInput("sampleInfoFile", "Choose Sample Information CSV File",
                             accept = c(".csv")),
                   
                   # Input for selecting categorical field
                   selectInput("categoryField", "Select Categorical Field",
                               choices = c("PMI",	"Age.of.Death",	"RIN",	"mRNA-Seq.reads",	
                               "Age.of.Onset", "Duration",	"CAG")),
                   # Input for selecting plot tpe
                   selectInput("plotType", "Select plot type",
                               choices = c("Barplot", "Boxplot", "Beeswarmplot", "Violinplot")),
                   
                   # Input for searching and selecting genes
                   textInput("geneSearch", "Search and Select Gene (For example: ENSG00000000003.10", 
                             value = "ENSG00000000003.10"),
                   
                   
                   # Button to generate the plot
                   actionButton("generatePlotBtn", "Generate Plot")
                 ),
                 
                 # Display plot
                 mainPanel(
                   withSpinner(plotOutput("geneExpressionPlot"))
                 )
               )
             )
    )
    
             
             
    
  )
)



server <- function(input, output, session) {

# Set the file size limit to 10 megabytes
  options(shiny.maxRequestSize = 100*1024^2)
  
  
  
  # Server logic for Sample Information Exploration
  # Reactive function to read the uploaded CSV file
  sample_data <- reactive({
    req(input$sample_file)
    read.csv(input$sample_file$datapath)
  })
  
  
  # Summary table output
  output$summary_table <- renderDataTable({
    summary_df <- data.frame(
      #Column_Name = sample_data(),
      Type = sapply(sample_data(), class),
      Mean_or_Distinct_Values = sapply(sample_data(), function(col) {
        if (is.numeric(col)) {
          paste0(mean(col), " (+/- ", sd(col), ")")
        } else {
          paste0(length(unique(col)), " distinct values")
        }
      })
    )
    datatable(summary_df)
  })
  
  # Data table output
  output$sample_table <- renderDT({
    datatable(sample_data(), options = list(order = list(1, 'asc')))
  })
  
  # Histogram tabs output
  output$hist_tabs <- renderUI({
    numeric_columns <- names(sample_data())[sapply(sample_data(), is.numeric)]
    
    tabs <- lapply(numeric_columns, function(col) {
      tabPanel(
        col,
        plotOutput(paste0("hist_plot_", col))
      )
    })
    
    do.call(tabsetPanel, tabs)
  })
  
  # Histogram plots output
  observe({
    numeric_columns <- names(sample_data())[sapply(sample_data(), is.numeric)]
    
    lapply(numeric_columns, function(col) {
      output[[paste0("hist_plot_", col)]] <- renderPlot({
        hist(sample_data()[[col]], main = col, xlab = col)
      })
    })
  })
  
  

  
# Server logic for Counts
  
  # Function to read counts matrix
  counts_data <- reactive({
    req(input$counts_file)
    read.csv(input$counts_file$datapath, row.names = 1, check.names = FALSE)
  })
  
  # Function to filter genes based on variance and non-zero samples thresholds
  filtered_genes <- reactive({
    counts_matrix <- counts_data()
    
    # Filter genes based on variance threshold
    variance_threshold <- quantile(apply(counts_matrix, 1, var), input$variance_threshold / 100)
    genes_passing_variance <- rownames(counts_matrix)[apply(counts_matrix, 1, var) >= variance_threshold]
    
    # Filter genes based on non-zero samples threshold
    nonzero_samples_threshold <- input$nonzero_samples_threshold / 100 * ncol(counts_matrix)
    genes_passing_nonzero_samples <- rownames(counts_matrix)[apply(counts_matrix > 0, 1, sum) >= nonzero_samples_threshold]
    
    # Combined list of genes passing both thresholds
    selected_genes <- intersect(genes_passing_variance, genes_passing_nonzero_samples)
    
    # Create a filtered counts matrix
    filtered_counts_matrix <- counts_matrix[selected_genes, ]
    
    # Ensure that filtered_counts_matrix is a numeric matrix
    filtered_counts_matrix <- as.matrix(filtered_counts_matrix)
    
    # Handle missing values (replace NA with 0 or another appropriate value)
    filtered_counts_matrix[is.na(filtered_counts_matrix)] <- 0
    
    # Create a data frame with filtering summary
    filter_summary <- data.frame(
      Metric = c(
        "Number of Samples",
        "Total Number of Genes",
        "Number of Genes Passing Variance Filter",
        "Percentage of Genes Passing Variance Filter",
        "Number of Genes Passing NonZero Samples Filter",
        "Percentage of Genes Passing NonZero Samples Filter"
    
      ),
      Value = c(
        ncol(counts_matrix),
        nrow(counts_matrix),
        length(genes_passing_variance),
        length(genes_passing_variance) / nrow(counts_matrix) * 100,
        length(genes_passing_nonzero_samples),
        length(genes_passing_nonzero_samples) / nrow(counts_matrix) * 100
        
      )
    )
    
    return(list(filtered_counts_matrix = filtered_counts_matrix, filter_summary = filter_summary))
  })
  
  # Output for filter summary table
  output$filter_summary_table <- renderDT({
    datatable(filtered_genes()$filter_summary, options = list(paging = FALSE, searching = FALSE, dom = ''), 
              colnames = c("Metric", "Value"))
  })
  
  # Function to create diagnostic scatter plots
  plot_median_vs_variance <- function(data, selected_genes) {
    library(ggplot2)
    
    medians <- apply(data, 1, median)
    variances <- apply(data, 1, var)
    
    # Create a data frame with medians, variances, and gene inclusion status
    plot_data <- data.frame(medians = medians, variances = variances, Included = rownames(data) %in% selected_genes)
    
    # Plot median vs variance with color differentiation
    median_vs_variance_plot <- ggplot(plot_data, aes(x = medians, y = variances, color = Included)) +
      geom_point() +
      scale_y_continuous(trans = 'log10') +
      scale_x_continuous(trans = 'log10') +
      labs(title = "Median Count vs Variance", x = "Median Count", y = "Variance") +
      scale_color_manual(values = c("lightblue", "darkblue"))
    
    return(median_vs_variance_plot)
  }
  
  # Function to generate diagnostic scatter plot for median count vs number of zeros
  plot_median_vs_zeros <- function(data, selected_genes) {
    library(ggplot2)
    
    medians <- apply(data, 1, median)
    zeros <- apply(data == 0, 1, sum)
    
    # Create a data frame with medians, number of zeros, and gene inclusion status
    plot_data <- data.frame(medians = medians, zeros = zeros, Included = rownames(data) %in% selected_genes)
    
    # Add a small constant to medians before taking the logarithm
    plot_data$medians <- log10(plot_data$medians + 1e-10)
    
    # Plot median vs number of zeros with color differentiation
    median_vs_zeros_plot <- ggplot(plot_data, aes(x = medians, y = zeros, color = Included)) +
      geom_point() +
      scale_y_continuous(trans = 'log10') +
      scale_x_continuous(trans = 'log10') +
      labs(title = "Median Count vs Number of Zeros", x = "Median Count", y = "Number of Zeros") +
      scale_color_manual(values = c("lightblue", "darkblue"))
    
    return(median_vs_zeros_plot)
  }
  
  # Output for scatter plots
  output$scatter_median_variance <- renderPlot({
    req(filtered_genes()$filtered_counts_matrix) # Ensure that the filtered counts matrix is available
    counts_matrix <- counts_data()
    
    plot_median_vs_variance(counts_matrix, rownames(filtered_genes()$filtered_counts_matrix))
  })
  
  output$scatter_median_zeros <- renderPlot({
    req(filtered_genes()$filtered_counts_matrix) # Ensure that the filtered counts matrix is available
    counts_matrix <- counts_data()
    plot_median_vs_zeros(counts_matrix, rownames(filtered_genes()$filtered_counts_matrix))
  })
  
  
  
  # Function to create heatmap
  output$clustered_heatmap <- renderPlot({
    filtered_data <- filtered_genes()
    
    # Optionally log-transform the counts matrix for better visualization
    log_transformed_counts <- log2(filtered_data$filtered_counts_matrix + 1)
    
    # Create a clustered heatmap using heatmap.2
    heatmap_object <- heatmap.2(
      log_transformed_counts,
      trace = "none",
    )
    
    return(heatmap_object)
  })
  
  
  
  # Function to create PCA plot
  plot_pca <- function(counts_data, PCAX, PCAY, title = "") {
    counts_no_gene <- counts_data[, -1]
    
    pca_result <- prcomp(t(counts_no_gene), center = TRUE)
    
    pca_data <- as.data.frame(pca_result$x[, 1:15])
    
    pca_plot <- ggplot(data = pca_data, aes(x = pca_data[, PCAX], y = pca_data[, PCAY])) +
      geom_point() +
      labs(title = title, x = paste("Scores of", PCAX), y = paste("Scores of", PCAY))
    
    return(pca_plot)
  }
  
  # Reactive expression to update PCA plot
  output$CPCA <- renderPlot({
    req(input$counts_file, input$PCAX, input$PCAY)
    counts_data <- read.csv(input$counts_file$datapath, row.names = 1, check.names = FALSE)
    plot_pca(counts_data, input$PCAX, input$PCAY, title = "PCA Plot")
  })
  
  
  
  
# Server logic for Differential Expression
  de_data <- reactive({
    req(input$de_file)
    read.csv(input$de_file$datapath)
  })
  
  # Differential expression table output
  output$de_table <- renderDT({
    datatable(de_data())
  })
  
  # Volcano plot output
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    
    dataf <- na.omit(dataf)
    
    vol_plot <- ggplot(dataf, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = (padj < 1 * 10^slider))) +
      geom_point() +
      scale_color_manual(values = c(color1, color2)) +
      labs(
        title = "Volcano Plot",
        x = x_name,
        y = paste0("-log10(", y_name, ")"),
        color = "padj" 
      ) +
      theme_minimal()+
      theme(legend.position = "bottom")
    
    return(vol_plot)
  }

  
  # Button functionality for volcano plot
  output$volcano <- renderPlot({
    input$plot_button
    isolate({
      volcano_plot(dataf = de_data(), input$x_axis, input$y_axis, 
                   input$slider, input$base, input$highlight)
    })
  })
  
  
  
# Server logic for Individual Gene Expression(s)

  # Function to read counts matrix
  counts_data_ige <- reactive({
    req(input$countsFile)
    read.csv(input$countsFile$datapath)
  })
  
  # Reactive function to read the uploaded CSV file
  sample_data_ige <- reactive({
    req(input$sampleInfoFile)
    read.csv(input$sampleInfoFile$datapath)
  })
  

  
  # Function to process gene input and generate merged data
  combined_data <- reactive({
    counts <- counts_data_ige()
    transposed_counts <- data.frame(t(counts[-1]))  # Exclude the first column if it's an identifier
    names(transposed_counts) <- counts$Column1
    
    sample <- sample_data_ige()
    
    # Merge data
    merged_data <- merge(sample, transposed_counts, by.x = "Sample.ID", by.y = "row.names", all.x = TRUE)
    
    # Set the row names as Sample.ID and remove the unnecessary column
    rownames(merged_data) <- merged_data$Sample.ID
    merged_data <- merged_data[, -1]
    
    # Return the merged data
    return(merged_data)
  })
  
 
  output$geneExpressionPlot <- renderPlot({
    req(input$generatePlotBtn)
    
    # Retrieve merged data
    merged_data <- combined_data()
    
    # Filter data based on gene search term
    selected_gene <- input$geneSearch
    selected_gene_data <- merged_data[, grepl(selected_gene, names(merged_data))]
    
    # Check if the selected categorical field is valid
    req(input$categoryField %in% colnames(merged_data))
    
    # Generate plot based on selected plot type
    plot_title <- paste("Plot of", selected_gene, "vs.", input$categoryField)
    
    if (input$plotType == "Barplot") {
      # Barplot
      barplot(height = selected_gene_data, names.arg = merged_data[[input$categoryField]],
              xlab = input$categoryField, ylab = selected_gene, main = plot_title)
    } else if (input$plotType == "Boxplot") {
      # Boxplot
      boxplot(selected_gene_data ~ merged_data[[input$categoryField]],
              xlab = input$categoryField, ylab = selected_gene, main = plot_title)
    } else if (input$plotType == "Beeswarmplot") {
      # Beeswarm plot
      beeswarm(selected_gene_data ~ merged_data[[input$categoryField]],
               pch = 16, col = "blue", xlab = input$categoryField,
               ylab = selected_gene, main = plot_title)
    } else if (input$plotType == "Violinplot") {
      # Violin plot using ggplot2
      ggplot(merged_data, aes(x = factor(merged_data[[input$categoryField]]), y = selected_gene)) +
        geom_violin() +
        xlab(input$categoryField) +
        ylab(selected_gene) +
        ggtitle(plot_title)
    }
  })
  

  

  
  
  
}


# Run the application
shinyApp(ui = ui, server = server)