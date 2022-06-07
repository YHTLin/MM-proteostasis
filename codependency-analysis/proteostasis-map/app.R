### ShinyApp accompanying proteostasis.R

#############################################
# Step 1 - Set working directory and load libraries
#############################################
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-proteostais/proteostasis-map/')
library(shinythemes)
library(shiny)
library(plotly)
library(Hmisc)
library(dendsort)
library(dplyr)
library(stringr)
#library(DT)  # For interactive data tables


#############################################
# Step 2 - Load data files from proteostasis.R
#############################################
load("./Data/proteostasis-processed-lite.RData")


#############################################
# Step 2a - Include heme and solid tumor filters
#############################################
## Separate cell lines by heme and solid tumors (05/24/2019)
meta$tissue = str_extract(meta$CCLE_Name, "_.*")
meta$tissue = sub("^_", "", meta$tissue)
solid_meta = filter(meta, tissue != "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")
heme_meta = filter(meta, tissue == "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE")


#############################################
# Step 3 - Data visualization functions
#############################################
corr_heatmap = function(df, corr_type, dist_method, clus_method, 
                        gene_row, gene_col, window_size) {
  # df = data frame of processed CRISPR or expression data from depmap (proCrispr or proExp)
  # corr_type = character of length 1 indicating "pearson" or "spearman"
  # dist_method = character indicating type of distance calculations 
  #    ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
  # clus_method = character indicating type of clustering method
  #    ("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  # gene_row = character indicating the name of the gene for row search
  # gene_col = character indicating the name of the gene for column search
  # window_size = numeric dictating number of genes above and below the target for display
  
  # Construct correlation matrix
  cor_mat = rcorr(t(df), type = corr_type)
  
  # Cluter matrix
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  rowSort = sort_hclust(hclust(dist(cor_mat$r, method = dist_method),
                               method = clus_method))   # Find row distance
  colSort = sort_hclust(hclust(dist(t(cor_mat$r), method = dist_method),
                               method = clus_method))   # Find column distance
  
  cor_mat$r = cor_mat$r[rowSort$labels[rowSort$order],
                        colSort$labels[rev(colSort$order)]]
  
  # Apply filter
  if (gene_row != "") {
    gene_index = which(rownames(cor_mat$r) == gene_row)
    win_range = seq(gene_index - window_size, gene_index + window_size)
    win_range = win_range[win_range > 0 & win_range <= nrow(cor_mat$r)]
    cor_mat$r = cor_mat$r[win_range, ]
  }
  
  if (gene_col != "") {
    gene_index = which(colnames(cor_mat$r) == gene_col)
    win_range = seq(gene_index - window_size, gene_index + window_size)
    win_range = win_range[win_range > 0 & win_range <= ncol(cor_mat$r)]
    cor_mat$r = cor_mat$r[, win_range]
  }
  
  # Plot interactive heatmap
  plot_ly(x = colnames(cor_mat$r),
          y = rownames(cor_mat$r),
          z = cor_mat$r,
          type = "heatmap",
          colorscale = colorRampPalette(c("blue", "white", "red"))(300)
  ) %>%
    colorbar(limits = c(-1, 1))
}


corr_heatmap2 = function(df, corr_type, dist_method, clus_method,
                         gene_row, gene_col, window_size){
  # df = data frame of processed CRISPR or expression data from depmap (proCrispr or proExp)
  # corr_type = character of length 1 indicating "pearson" or "spearman"
  # dist_method = character indicating type of distance calculations 
  #    ("euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski")
  # clus_method = character indicating type of clustering method
  #    ("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
  # gene_row = character indicating the name of the gene for row search
  # gene_col = character indicating the name of the gene for column search
  # window_size = numeric dictating number of genes above and below the target for display
  
  # Construct correlation matrix
  cor_mat = rcorr(t(df), type = corr_type)
  cor_mat$r = -cor_mat$r[grepl("_CRISPR", rownames(cor_mat$r)),   # A negative has been inserted to reflect correlation between expression and SENSITIVITY
                         grepl("_EXP", colnames(cor_mat$r))]
  #cor_mat$P = cor_mat$P[grepl("_CRISPR", rownames(cor_mat$P)),
  #                      grepl("_EXP", colnames(cor_mat$P))]
  
  # Filter out rows or columns with NA
  if (sum(apply(cor_mat$r, 1, function(x) all(is.nan(x)))) > 0) {
    filt = apply(apply(cor_mat$r, 1, function(x) all(is.nan(x))))
    cor_mat$r = cor_mat$r[!filt, ]
    #cor_mat$P = cor_mat$P[!filt, ]
  }
  if (sum(apply(cor_mat$r, 2, function(x) all(is.nan(x)))) > 0) {
    filt = apply(cor_mat$r, 2, function(x) all(is.nan(x)))
    cor_mat$r = cor_mat$r[, !filt]
    #cor_mat$P = cor_mat$P[, !filt]
  }
  
  # Cluster matrix
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  rowSort = sort_hclust(hclust(dist(cor_mat$r, method = dist_method),
                               method = clus_method))   # Find row distance
  colSort = sort_hclust(hclust(dist(t(cor_mat$r), method = dist_method),
                               method = clus_method))   # Find column distance
  
  cor_mat$r = cor_mat$r[rowSort$labels[rowSort$order],
                        colSort$labels[rev(colSort$order)]]
  
  # Apply filter
  if (gene_row != "") {
    gene_index = which(rownames(cor_mat$r) == gene_row)
    win_range = seq(gene_index - window_size, gene_index + window_size)
    win_range = win_range[win_range > 0 & win_range <= nrow(cor_mat$r)]
    cor_mat$r = cor_mat$r[win_range, ]
  }
  
  if (gene_col != "") {
    gene_index = which(colnames(cor_mat$r) == gene_col)
    win_range = seq(gene_index - window_size, gene_index + window_size)
    win_range = win_range[win_range > 0 & win_range <= ncol(cor_mat$r)]
    cor_mat$r = cor_mat$r[, win_range]
  }
  
  # Plot interactive matrix
  plot_ly(x = colnames(cor_mat$r),
          y = rownames(cor_mat$r),
          z = cor_mat$r,
          type = "heatmap",
          colorscale = colorRampPalette(c("blue", "white", "red"))(300)
  ) %>%
    colorbar(limits = c(-1, 1))
}


scatter_plot = function(df, geneX, geneY) {
  # df = data frame of processed CRISPR or expression data from depmap (proCrispr or proExp)
  # geneX = character indicating gene to plot on the x-axis
  # geneY = character indicating gene to plot on the y-axis
  df = as.data.frame(t(df))
  df = df %>%
    mutate(DepMap_ID = rownames(df)) %>%
    mutate(GENEx = df[[geneX]]) %>%
    mutate(GENEy = df[[geneY]]) %>%
    left_join(meta, by = "DepMap_ID") %>%
    select(GENEx, GENEy, CCLE_Name, Primary.Disease)
  df$Primary.Disease[!grepl("Myeloma|Lymphoma|Leukemia", df$Primary.Disease)] = "Solid Tumor"
  
  # Linear fit
  fit <- lm(GENEy ~ GENEx, data = df)
  
  # Plot scatter
  plot_ly(type = 'scatter', mode = 'markers') %>%
    add_trace(
      x = df$GENEx[df$Primary.Disease != "Solid Tumor"],
      y = df$GENEy[df$Primary.Disease != "Solid Tumor"],
      opacity = 0.7,
      color = df$Primary.Disease[df$Primary.Disease != "Solid Tumor"],
      text = df$CCLE_Name[df$Primary.Disease != "Solid Tumor"],
      colors = "Set1",
      marker = list(
        size = 12)
    ) %>%
    add_trace(
      x = df$GENEx[df$Primary.Disease == "Solid Tumor"],
      y = df$GENEy[df$Primary.Disease == "Solid Tumor"],
      opacity = 0.2,
      color = df$Primary.Disease[df$Primary.Disease == "Solid Tumor"],
      text = df$CCLE_Name[df$Primary.Disease == "Solid Tumor"],
      marker = list(
        color = "gray",
        size = 7)
    ) %>%
    add_lines(x = ~GENEx, y = fitted(fit), data = df, 
              line = list(color = "orange"),
              showlegend = FALSE) %>%
    layout(
      title = paste0("Number of Cell Lines: ", nrow(df), "  \\\  ",
                     "Equation: y = ", round(summary(fit)$coefficients[2, 1], 3), 
                     "x + (", round(summary(fit)$coefficients[1, 1], 3), ")  \\\  ",
                     "R^2: ", round(summary(fit)$r.squared, 4)),
      xaxis = list(title = geneX), 
      yaxis = list(title = geneY))
}


#############################################
# Step 4 - ShinyApp code
#############################################
# USER INTERFACEs ------------------------------------------
ui = navbarPage(
  "Proteostasis Network Maps",
  theme = shinytheme("sandstone"),
  #theme = shinytheme("simplex"),
  
  
  # CRISPR PAGE ----------------------------------------------
  tabPanel("CRISPR",
           h3("CRISPR Correlation Analysis"),
           
           hr(),
           
           wellPanel(
             h4("Specify Clustering and Filtering Parameters"),
             
             br(),
             
             fluidRow(
               # Input correlation measure
               column(4, selectInput("CRISPR_corr_type",
                                     label = "Correlation Type",
                                     choices = list("Pearson" = "pearson",
                                                    "Spearman" = "spearman"),
                                     selected = "pearson")
               ),
               # Input distance measure
               column(4, selectInput("CRISPR_dist_method", 
                                     label = "Distance Method", 
                                     choices = list("Euclidean" = "euclidean",
                                                    "Maximum" = "maximum", 
                                                    "Manhattan" = "manhattan", 
                                                    "Canberra" = "canberra"), 
                                     selected = "euclidean")
               ),
               # Input clustering method
               column(4, selectInput("CRISPR_clus_method", 
                                     label = "Cluster Method", 
                                     choices = list("Complete" = "complete",
                                                    "Ward" = "ward.D",
                                                    "Average" = "average",
                                                    "Centroid" = "centroid"), 
                                     selected = "complete")
                      )
               ),
             
             fluidRow(
               # Input Gene for Row Search
               column(4, selectizeInput("CRISPR_select_gene_row",
                                        "Gene Search Filter (row)",
                                        choices = sort(rownames(proCrispr)),
                                        options = list(
                                          placeholder = 'Please select an option below',
                                          onInitialize = I('function() { this.setValue(""); }')
                                        )),
                      
                      helpText("Note: clear filter to reset heatmap.")
                      ),
               
               # Input Gene for Column Search
               column(4, selectizeInput("CRISPR_select_gene_col",
                                        "Gene Search Filter (column)",
                                        choices = sort(rownames(proCrispr)),
                                        options = list(
                                          placeholder = 'Please select an option below',
                                          onInitialize = I('function() { this.setValue(""); }')
                                        ))
                      ),
               
               # Input Gene for Window Size
               column(4, sliderInput("CRISPR_select_window", 
                                     "Window Size",
                                     min = 0, max = 100,
                                     value = 10))
             )
           ),
           
           tabsetPanel(
             tabPanel("All Lines",
                      plotlyOutput(outputId = "CRISPR_all_heatmap", width = "1100px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR Gene Dependency Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("CRISPR_all_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proCrispr))
                                                   )),
                          # Input Gene 2
                          column(4, selectizeInput("CRISPR_all_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proCrispr))
                                                   ))
                        ),
                        
                        plotlyOutput(outputId = "CRISPR_all_scatter", width = "1100px", height = "600px")
                      )),
             tabPanel("Myeloma Only",
                      plotlyOutput(outputId = "CRISPR_MM_heatmap", width = "900px", height = "500px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR Gene Dependency Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("CRISPR_MM_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proCrispr))
                                                   )),
                          # Input Gene 2
                          column(4, selectizeInput("CRISPR_MM_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proCrispr))
                                                   ))
                          ),
                        
                        plotlyOutput(outputId = "CRISPR_MM_scatter", width = "900px", height = "500px")
                        )),
             tabPanel("Heme Only",
                      plotlyOutput(outputId = "CRISPR_heme_heatmap", width = "900px", height = "500px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR Gene Dependency Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("CRISPR_heme_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proCrispr))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("CRISPR_heme_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proCrispr))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "CRISPR_heme_scatter", width = "900px", height = "500px")
                      )),
             tabPanel("Solid Tumor Only",
                      plotlyOutput(outputId = "CRISPR_solid_heatmap", width = "900px", height = "500px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR Gene Dependency Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("CRISPR_solid_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proCrispr))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("CRISPR_solid_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proCrispr))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "CRISPR_solid_scatter", width = "900px", height = "500px")
                      ))
           ),
           
           hr(),
           
           helpText("Please direct any inquiries to Arun Wiita (arun.wiita@ucsf.edu).")
  ),
  
  
  # EXPRESSION PAGE ----------------------------------------------
  tabPanel("Gene Expression",
           h3("Gene Expression Correlation Analysis"),
           
           hr(),
           
           wellPanel(
             h4("Specify Clustering and Filtering Parameters"),
             
             br(),
             
             fluidRow(
               # Input correlation measure
               column(4, selectInput("EXP_corr_type",
                                     label = "Correlation Type",
                                     choices = list("Pearson" = "pearson",
                                                    "Spearman" = "spearman"),
                                     selected = "pearson")
               ),
               # Input distance measure
               column(4, selectInput("EXP_dist_method", 
                                     label = "Distance Method", 
                                     choices = list("Euclidean" = "euclidean",
                                                    "Maximum" = "maximum", 
                                                    "Manhattan" = "manhattan", 
                                                    "Canberra" = "canberra"), 
                                     selected = "euclidean")
               ),
               # Input clustering method
               column(4, selectInput("EXP_clus_method", 
                                     label = "Cluster Method", 
                                     choices = list("Complete" = "complete",
                                                    "Ward" = "ward.D",
                                                    "Average" = "average",
                                                    "Centroid" = "centroid"), 
                                     selected = "complete")
               )
             ),
             
             fluidRow(
               # Input Gene for Row Search
               column(4, selectizeInput("EXP_select_gene_row",
                                        "Gene Search Filter (row)",
                                        choices = sort(rownames(proExp)),
                                        options = list(
                                          placeholder = 'Please select an option below',
                                          onInitialize = I('function() { this.setValue(""); }')
                                        )),
                      
                      helpText("Note: clear filter to reset heatmap.")
               ),
               
               # Input Gene for Column Search
               column(4, selectizeInput("EXP_select_gene_col",
                                        "Gene Search Filter (column)",
                                        choices = sort(rownames(proExp)),
                                        options = list(
                                          placeholder = "Please select an option below",
                                          onInitialize = I('function() { this.setValue(""); }')
                                        ))
               ),
               
               # Input Gene for Window Size
               column(4, sliderInput("EXP_select_window", 
                                     "Window Size",
                                     min = 0, max = 100,
                                     value = 10))
             )
           ),
           
           tabsetPanel(
             tabPanel("All Lines", 
                      plotlyOutput(outputId = "EXP_all_heatmap", width = "1100px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("Gene Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("EXP_all_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proExp))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("EXP_all_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proExp))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "EXP_all_scatter", width = "1100px", height = "600px")
                      )
             ), 
             tabPanel("Myeloma Only",
                      plotlyOutput(outputId = "EXP_MM_heatmap", width = "1100px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("Gene Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("EXP_MM_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proExp))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("EXP_MM_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proExp))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "EXP_MM_scatter", width = "900px", height = "500px")
                      )),
             tabPanel("Heme Only",
                      plotlyOutput(outputId = "EXP_heme_heatmap", width = "1100px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("Gene Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("EXP_heme_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proExp))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("EXP_heme_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proExp))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "EXP_heme_scatter", width = "900px", height = "500px")
                      )),
             tabPanel("Solid Tumor Only",
                      plotlyOutput(outputId = "EXP_solid_heatmap", width = "1100px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("Gene Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("EXP_solid_scatter_geneX",
                                                   "Scatter X:",
                                                   choices = sort(rownames(proExp))
                          )),
                          # Input Gene 2
                          column(4, selectizeInput("EXP_solid_scatter_geneY",
                                                   "Scatter Y:",
                                                   choices = sort(rownames(proExp))
                          ))
                        ),
                        
                        plotlyOutput(outputId = "EXP_solid_scatter", width = "900px", height = "500px")
                      ))
           ),
           
           hr(),
           
           helpText("Please direct any inquiries to Arun Wiita (arun.wiita@ucsf.edu).")
  ),
  
  
  # CRISPR vs Expression PAGE ----------------------------------------------
  tabPanel("CRISPR vs Gene Expression",
           h3("CRISPR-vs-Expression Correlation Analysis"),
           
           hr(),
           
           wellPanel(
             h4("Specify Clustering and Filtering Parameters"),
             
             br(),
             
             fluidRow(
               # Input correlation measure
               column(4, selectInput("COMBO_corr_type",
                                     label = "Correlation Type",
                                     choices = list("Pearson" = "pearson",
                                                    "Spearman" = "spearman"),
                                     selected = "pearson")
               ),
               # Input distance measure
               column(4, selectInput("COMBO_dist_method", 
                                     label = "Distance Method", 
                                     choices = list("Euclidean" = "euclidean",
                                                    "Maximum" = "maximum", 
                                                    "Manhattan" = "manhattan", 
                                                    "Canberra" = "canberra"), 
                                     selected = "euclidean")
               ),
               # Input clustering method
               column(4, selectInput("COMBO_clus_method", 
                                     label = "Cluster Method", 
                                     choices = list("Complete" = "complete",
                                                    "Ward" = "ward.D",
                                                    "Average" = "average",
                                                    "Centroid" = "centroid"), 
                                     selected = "complete")
               )),
               
               fluidRow(
                 # Input Gene for Row Search
                 column(4, selectizeInput("COMBO_select_gene_row",
                                          "Gene Search Filter (row)",
                                          choices = sort(grep("_CRISPR", rownames(combo), value = T)),
                                          options = list(
                                            placeholder = 'Please select an option below',
                                            onInitialize = I('function() { this.setValue(""); }')
                                          )),
                        
                        helpText("Note 1: Clear filter to reset heatmap."),
                        helpText("Note 2: Heatmap scores reflect correlation between expression and
                                 sensitivity instead effect scores.")
                 ),
                 
                 # Input Gene for Column Search
                 column(4, selectizeInput("COMBO_select_gene_col",
                                          "Gene Search Filter (column)",
                                          choices = sort(grep("_EXP", rownames(combo), value = T)),
                                          options = list(
                                            placeholder = 'Please select an option below',
                                            onInitialize = I('function() { this.setValue(""); }')
                                          ))
                 ),
                 
                 # Input Gene for Window Size
                 column(4, sliderInput("COMBO_select_window", 
                                       "Window Size",
                                       min = 0, max = 100,
                                       value = 10))
               )
           ),
           
           tabsetPanel(
             tabPanel("All Lines", 
                      plotlyOutput(outputId = "COMBO_all_heatmap", width = "1200px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR-Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("COMBO_all_scatter_geneX",
                                                   "Scatter X (Expression):",
                                                   choices = sort(grep("_EXP", rownames(combo), value = T)))
                          ),
                          # Input Gene 2
                          column(4, selectizeInput("COMBO_all_scatter_geneY",
                                                   "Scatter Y (CRISPR):",
                                                   choices = sort(grep("_CRISPR", rownames(combo), value = T)))
                          )
                        ),
                        
                        plotlyOutput(outputId = "COMBO_all_scatter", width = "1100px", height = "600px")
                      )
             ), 
             
             tabPanel("Myeloma Only",
                      plotlyOutput(outputId = "COMBO_MM_heatmap", width = "1200px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR-Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("COMBO_MM_scatter_geneX",
                                                   "Scatter X (Expression):",
                                                   choices = sort(grep("_EXP", rownames(combo), value = T)))
                          ),
                          # Input Gene 2
                          column(4, selectizeInput("COMBO_MM_scatter_geneY",
                                                   "Scatter Y (CRISPR):",
                                                   choices = sort(grep("_CRISPR", rownames(combo), value = T)))
                          )
                        ),
                        
                        plotlyOutput(outputId = "COMBO_MM_scatter", width = "1100px", height = "600px")
                      )
             ),
             
             tabPanel("Heme Only",
                      plotlyOutput(outputId = "COMBO_heme_heatmap", width = "1200px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR-Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("COMBO_heme_scatter_geneX",
                                                   "Scatter X (Expression):",
                                                   choices = sort(grep("_EXP", rownames(combo), value = T)))
                          ),
                          # Input Gene 2
                          column(4, selectizeInput("COMBO_heme_scatter_geneY",
                                                   "Scatter Y (CRISPR):",
                                                   choices = sort(grep("_CRISPR", rownames(combo), value = T)))
                          )
                        ),
                        
                        plotlyOutput(outputId = "COMBO_heme_scatter", width = "1100px", height = "600px")
                      )
             ),
             tabPanel("Solid Tumors Only",
                      plotlyOutput(outputId = "COMBO_solid_heatmap", width = "1200px", height = "700px"),
                      
                      hr(),
                      
                      wellPanel(
                        h4("CRISPR-Expression Scatter Plot"),
                        
                        fluidRow(
                          # Input Gene 1
                          column(4, selectizeInput("COMBO_solid_scatter_geneX",
                                                   "Scatter X (Expression):",
                                                   choices = sort(grep("_EXP", rownames(combo), value = T)))
                          ),
                          # Input Gene 2
                          column(4, selectizeInput("COMBO_solid_scatter_geneY",
                                                   "Scatter Y (CRISPR):",
                                                   choices = sort(grep("_CRISPR", rownames(combo), value = T)))
                          )
                        ),
                        
                        plotlyOutput(outputId = "COMBO_solid_scatter", width = "1100px", height = "600px")
                      )
             )
           ),
           
           hr(),
           
           helpText("Please direct any inquiries to Arun Wiita (arun.wiita@ucsf.edu).")
  ),
  
  # DOWNLOAD ---------------------------------------------------------------
  tabPanel("Download",
           h3("Source Data"),
           
           hr(),
           
           fluidRow(
             column(8,
               p(HTML("CRISPR deletion screen and gene expression data sets were downloaded
                      from the <a href='https://depmap.org/portal/download/'>Cancer Dependency Map</a> 
                      (19Q1 release). The raw data were filtered to include genes from a curated
                      list of 441 members known to regulate protein 
                      homeostasis. The processed data tables used for visualization are available 
                      for download.")))),
           br(),
           
           fluidRow(column(4, downloadButton("DL_CRISPR", "Download Filtered CRISPR Data"))),
           
           br(),
           
           fluidRow(column(4, downloadButton("DL_EXP", "Download Filtered Expression Data"))),
           
           br(),
           
           fluidRow(column(4, downloadButton("DL_PROTEOSTASIS", "Download 441 Proteostasis Genes"))),
           
           br(),
           
           fluidRow(column(4, downloadButton("DL_META", "Download Cell Line Metadata"))),
           
           hr(),
           
           helpText("Please direct any inquiries to Arun Wiita (arun.wiita@ucsf.edu).")
  )
)



# SERVER ------------------------------------------------------------------
server = function(input, output, session) {
  
  # CRISPR PAGE --------------------------------------------------
  output$CRISPR_all_heatmap = renderPlotly({ 
    corr_heatmap(proCrispr, 
                 input$CRISPR_corr_type, 
                 input$CRISPR_dist_method, 
                 input$CRISPR_clus_method,
                 input$CRISPR_select_gene_row,
                 input$CRISPR_select_gene_col,
                 input$CRISPR_select_window)})
  
  output$CRISPR_MM_heatmap = renderPlotly({ 
    corr_heatmap(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
                 input$CRISPR_corr_type, 
                 input$CRISPR_dist_method, 
                 input$CRISPR_clus_method,
                 input$CRISPR_select_gene_row,
                 input$CRISPR_select_gene_col,
                 input$CRISPR_select_window)})
  
  output$CRISPR_heme_heatmap = renderPlotly({ 
    corr_heatmap(proCrispr[names(proCrispr) %in% heme_meta$DepMap_ID], 
                 input$CRISPR_corr_type, 
                 input$CRISPR_dist_method, 
                 input$CRISPR_clus_method,
                 input$CRISPR_select_gene_row,
                 input$CRISPR_select_gene_col,
                 input$CRISPR_select_window)})
  
  output$CRISPR_solid_heatmap = renderPlotly({ 
    corr_heatmap(proCrispr[names(proCrispr) %in% solid_meta$DepMap_ID], 
                 input$CRISPR_corr_type, 
                 input$CRISPR_dist_method, 
                 input$CRISPR_clus_method,
                 input$CRISPR_select_gene_row,
                 input$CRISPR_select_gene_col,
                 input$CRISPR_select_window)})
  
  output$CRISPR_all_scatter = renderPlotly({
    scatter_plot(proCrispr, 
                 input$CRISPR_all_scatter_geneX, 
                 input$CRISPR_all_scatter_geneY)
  })
  
  output$CRISPR_MM_scatter = renderPlotly({
    scatter_plot(proCrispr[names(proCrispr) %in% MM_meta$DepMap_ID], 
                 input$CRISPR_MM_scatter_geneX, 
                 input$CRISPR_MM_scatter_geneY)
  })
  
  output$CRISPR_heme_scatter = renderPlotly({
    scatter_plot(proCrispr[names(proCrispr) %in% heme_meta$DepMap_ID], 
                 input$CRISPR_heme_scatter_geneX, 
                 input$CRISPR_heme_scatter_geneY)
  })
  
  output$CRISPR_solid_scatter = renderPlotly({
    scatter_plot(proCrispr[names(proCrispr) %in% solid_meta$DepMap_ID], 
                 input$CRISPR_solid_scatter_geneX, 
                 input$CRISPR_solid_scatter_geneY)
  })
  
  
  # EXPRESSION PAGE ----------------------------------------------
  output$EXP_all_heatmap = renderPlotly({ 
    corr_heatmap(proExp, 
                 input$EXP_corr_type, 
                 input$EXP_dist_method, 
                 input$EXP_clus_method,
                 input$EXP_select_gene_row,
                 input$EXP_select_gene_col,
                 input$EXP_select_window)})
  
  output$EXP_MM_heatmap = renderPlotly({ 
    corr_heatmap(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
                 input$EXP_corr_type, 
                 input$EXP_dist_method, 
                 input$EXP_clus_method,
                 input$EXP_select_gene_row,
                 input$EXP_select_gene_col,
                 input$EXP_select_window)})
  
  output$EXP_heme_heatmap = renderPlotly({ 
    corr_heatmap(proExp[names(proExp) %in% heme_meta$DepMap_ID], 
                 input$EXP_corr_type, 
                 input$EXP_dist_method, 
                 input$EXP_clus_method,
                 input$EXP_select_gene_row,
                 input$EXP_select_gene_col,
                 input$EXP_select_window)})
  
  output$EXP_solid_heatmap = renderPlotly({ 
    corr_heatmap(proExp[names(proExp) %in% solid_meta$DepMap_ID], 
                 input$EXP_corr_type, 
                 input$EXP_dist_method, 
                 input$EXP_clus_method,
                 input$EXP_select_gene_row,
                 input$EXP_select_gene_col,
                 input$EXP_select_window)})
  
  output$EXP_all_scatter = renderPlotly({
    scatter_plot(proExp, 
                 input$EXP_all_scatter_geneX, 
                 input$EXP_all_scatter_geneY)
  })
  
  output$EXP_MM_scatter = renderPlotly({
    scatter_plot(proExp[names(proExp) %in% MM_meta$DepMap_ID], 
                 input$EXP_MM_scatter_geneX, 
                 input$EXP_MM_scatter_geneY)
  })
  
  output$EXP_heme_scatter = renderPlotly({
    scatter_plot(proExp[names(proExp) %in% heme_meta$DepMap_ID], 
                 input$EXP_heme_scatter_geneX, 
                 input$EXP_heme_scatter_geneY)
  })
  
  output$EXP_solid_scatter = renderPlotly({
    scatter_plot(proExp[names(proExp) %in% solid_meta$DepMap_ID], 
                 input$EXP_solid_scatter_geneX, 
                 input$EXP_solid_scatter_geneY)
  })
  
  # CRISPR vs Expression PAGE ------------------------------------
  select_gene_col = reactive( {   # The following genes show no change in expression in myeloma lines, so correlation cannot be calculated
    validate(
      need(!(input$COMBO_select_gene_col %in% c("PDILT_EXP", "DNAJB3_EXP")), 
           paste0("'", input$COMBO_select_gene_col, "'", " is not in heatmap. Please select another gene.")))
    
    input$COMBO_select_gene_col
    })
  
  output$COMBO_all_heatmap = renderPlotly({ 
    corr_heatmap2(combo, 
                  input$COMBO_corr_type, 
                  input$COMBO_dist_method, 
                  input$COMBO_clus_method,
                  input$COMBO_select_gene_row,
                  input$COMBO_select_gene_col,
                  input$COMBO_select_window)})
  
  output$COMBO_MM_heatmap = renderPlotly({ 
    corr_heatmap2(combo[names(combo) %in% MM_meta$DepMap_ID], 
                  input$COMBO_corr_type, 
                  input$COMBO_dist_method, 
                  input$COMBO_clus_method,
                  input$COMBO_select_gene_row,
                  select_gene_col(),
                  input$COMBO_select_window)})
  
  output$COMBO_heme_heatmap = renderPlotly({ 
    corr_heatmap2(combo[names(combo) %in% heme_meta$DepMap_ID], 
                  input$COMBO_corr_type, 
                  input$COMBO_dist_method, 
                  input$COMBO_clus_method,
                  input$COMBO_select_gene_row,
                  input$COMBO_select_gene_col,
                  input$COMBO_select_window)})
  
  output$COMBO_solid_heatmap = renderPlotly({ 
    corr_heatmap2(combo[names(combo) %in% solid_meta$DepMap_ID], 
                  input$COMBO_corr_type, 
                  input$COMBO_dist_method, 
                  input$COMBO_clus_method,
                  input$COMBO_select_gene_row,
                  input$COMBO_select_gene_col,
                  input$COMBO_select_window)})
  
  output$COMBO_all_scatter = renderPlotly({
    scatter_plot(combo, 
                 input$COMBO_all_scatter_geneX, 
                 input$COMBO_all_scatter_geneY)
  })
  
  output$COMBO_MM_scatter = renderPlotly({
    scatter_plot(combo[names(combo) %in% MM_meta$DepMap_ID], 
                 input$COMBO_MM_scatter_geneX, 
                 input$COMBO_MM_scatter_geneY)
  })
  
  output$COMBO_heme_scatter = renderPlotly({
    scatter_plot(combo[names(combo) %in% heme_meta$DepMap_ID], 
                 input$COMBO_heme_scatter_geneX, 
                 input$COMBO_heme_scatter_geneY)
  })
  
  output$COMBO_solid_scatter = renderPlotly({
    scatter_plot(combo[names(combo) %in% solid_meta$DepMap_ID], 
                 input$COMBO_solid_scatter_geneX, 
                 input$COMBO_solid_scatter_geneY)
  })
  
  # DOWNLOAD PAGE -------------------------------------------------------
  output$DL_CRISPR <- downloadHandler(
    filename = function() {"filtered_DepMap_19Q1_CRISPR.csv"},
    content = function(file) {
      write.csv(proCrispr, file, row.names = FALSE)
    })
  output$DL_EXP <- downloadHandler(
    filename = function() {"filtered_DepMap_19Q1_EXPRESSION_TPM.csv"},
    content = function(file) {
      write.csv(proExp, file, row.names = FALSE)
    })
  output$DL_PROTEOSTASIS <- downloadHandler(
    filename = function() {"441_proteostasis_genes.csv"},
    content = function(file) {
      write.csv(proteasome, file, row.names = FALSE, col.names = FALSE)
    })
  output$DL_META <- downloadHandler(
    filename = function() {"filtered_DepMap_19Q1_CELL_LINE_METADATA.csv"},
    content = function(file) {
      write.csv(meta, file, row.names = FALSE)
    })
}

shinyApp(ui = ui, server = server)


#############################################
# Step 5 - Deploy App!
#############################################
# Copy and run code in console
#library(rsconnect)
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Computational_tools/random_r_scripts/arun-proteostais/proteostasis-map/')
#deployApp()
