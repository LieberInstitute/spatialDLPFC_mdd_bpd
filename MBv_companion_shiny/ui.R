suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(markdown)
  library(SpatialExperiment)
})
set.seed(1234)

spe <- readRDS("data/MBv-shiny_pseudobulk-spe_both-annotations.rds")
source("custom-plot_utils_Jun-17.R")

navbarPage("MBv Pseudobulk Companion App",
           tabPanel("Project Overview",
                    layout_columns(
                      card(includeMarkdown("appdescription.md"),
                           layout_columns(
                             card(includeMarkdown("landingpage.md")),
                             card(
                                  img(src="graphical_abstract-1.png", height=600, width=532)
                                  ),
                             gap = 50
                             )
                           ),
                    )
           ),
           tabPanel("Pseudobulk Expression Plots",
                    # Define the sidebar with one input
                    sidebarPanel(
                      selectizeInput("gene_name", "Gene name:", 
                                     choices=sort(rownames(spe)),
                                     #selected="SST"),
                                     options = list(
                                       placeholder = 'Start typing to select from list of available genes',
                                       onInitialize = I('function() { this.setValue(""); }')
                                     )),
                      p("Use the bar to select a gene to plot. Gene selection is possible for 21,080 genes that passed edgeR filtering for both domain", tags$sub("SP", .noWS = "outside"), " and domain", tags$sub("CT", .noWS = "outside"), " datasets."),
                      p("Plots present pseudobulk data, where RNAseq counts are aggregated for each annotated domain within each donor. Pseudobulk counts were then normalized and transformed (logcount expression, ", tags$em("y", .noWS="outside"), "-axis). Pseudobulk processing occurred separately for the domain", tags$sub("SP", .noWS = "outside"), " and domain", tags$sub("CT", .noWS = "outside"), " annotation strategies."),
                      p("Each point represents a single pseudobulk sample corresponding to each domain identified within each donor. Points are colored based on diagnosis group (DX). Crossbar indicates the mean +/- SD.")
                    ),
                    # Create a spot for the barplot
                    mainPanel(
                      card(h3("Whole-tissue DE model"),
                           p("Facets separate results from the two annotation strategies and contain ", tags$em("F", .noWS="outside"), "-test statistics. If ", tags$em("F", .noWS="outside"), "-test adj. ",tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks."),
                           plotOutput("whole.tissue")),
                      card(h3("Domain-restricted DE model: domain", tags$sub("SP", .noWS = "outside")),
                           p("Facets separate donors by biological sex. Subtitle contains ", tags$em("F", .noWS="outside"), "-test statistics for domain", tags$sub("SP", .noWS = "outside"), " and will be bolded if ", tags$em("F", .noWS="outside"), "-test results indicated a significant difference across diagnosis-by-sex domains (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks."),
                           plotOutput("domain.sp")),
                      card(h3("Domain-restricted DE model: domain", tags$sub("CT", .noWS = "outside")),
                           p("Facets separate donors by biological sex. Subtitle contains ", tags$em("F", .noWS="outside"), "-test statistics for domain", tags$sub("CT", .noWS = "outside"), " and will be bolded if ", tags$em("F", .noWS="outside"), "-test results indicated a significant difference across diagnosis-by-sex domains (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks."),
                           plotOutput("domain.ct"))
                    )
           ),
           tabPanel("Pseudobulk Methods Details",
                    layout_columns(
                      card(
                        includeMarkdown("methods.md"),
                        fill=FALSE, max_height = 600
                        ),
                     card(
                       img(src="pseudobulk-pipeline.png", height=600, width=536),
                     ),
                     gap = 50
                    )
           )
)
