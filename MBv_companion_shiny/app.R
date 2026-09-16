suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(markdown)
  library(SpatialExperiment)
})
set.seed(1234)

spe <- readRDS("data/MBv-shiny_pseudobulk-spe_both-annotations.rds")
source("custom-plot_utils_Jun-17.R")

myitems = list(
  accordion_panel("Step 1: Spot-level annotation", includeMarkdown("markdown/methods_1.md")),
  accordion_panel("Step 2: Pseudobulk processing", includeMarkdown("markdown/methods_2.md")),
  accordion_panel("Step 3: Differential expression (DE) pre-processing", includeMarkdown("markdown/methods_3.md")),
  accordion_panel("Step 4: DE modeling", includeMarkdown("markdown/methods_4.md")),
  accordion_panel("FAQ: Why do some gene expression differences look significant, but DE results indicate that they are not?", includeMarkdown("markdown/methods_5.md")),
  accordion_panel("FAQ: Why are the lowest pseudobulk expression values not equal to zero?", includeMarkdown("markdown/methods_6.md"))
)

ui <- page_navbar(title = "MBv Pseudobulk Companion App",
                 tabPanel("Study Design",
                          p("This is a companion app for ", tags$a(href="http://www.google.com", target="_blank", "Thompson et al. 2026"), " that allows users to easily check the expression of their favorite gene across our diagnosis-by-sex groups and domain annotations."),
                          tags$ul(
                            tags$li(tags$strong("Study Design", .noWS="outside"), ": Below you will find a quick refresher on our study design that provides helpful context for interpreting gene expression plots."),
                            tags$li("In the ", tags$strong("Pseudobulk Expression Plots", .noWS="outside"), " tab you can view violin plots of our pseudobulked spatial gene expression data. Plot subtitles contain the results of the gene-level omnibus ", tags$em("F", .noWS="outside"), "-test. Diagnosis-by-sex group contrasts are annotated to indicate statistically significant ", tags$em("t", .noWS="outside"), "-test results, where appropriate"),
                            tags$li("In the ", tags$strong("Pseudobulk Methods Details", .noWS="outside"), " tab we provide additional information on our pseudobulk processing pipeline that may answer common questions.")
                          ),
                          card(
                            layout_columns(
                              includeMarkdown(landingpage),
                              img(src="https://github.com/LieberInstitute/Posit_Connect_shiny_apps/blob/devel/dlpfc/DLPFC_MBv_companion/code/www/graphical_abstract-1.png?raw=TRUE", height=633, width=550),
                              gap = 50
                            )
                          )
                 ),
                 tabPanel("Pseudobulk Expression Plots",
                          sidebarLayout(
                            sidebarPanel(
                              selectizeInput("gene_name", "Gene name:",
                                             choices=sort(rownames(spe)),
                                             #selected="SST"),
                                             options = list(
                                               placeholder = 'Start typing to select from list of available genes',
                                               onInitialize = I('function() { this.setValue(""); }')
                                             )),
                              p("Use the bar above to select a gene to plot (e.g. SST). Gene selection is possible for 21,080 genes that include about 7,000 genes lowly expressed in SRT but highly expressed in snRNAseq from the dlPFC. Only the subset of genes included in DE modeling have statistics available (see 'Pseudobulk Methods Details' tab for more information)."),
                              p("Plots present pseudobulk data, where spot-level RNAseq counts are aggregated for each annotated domain within each donor. Pseudobulk counts were then normalized and transformed (logcount expression, ", tags$em("y", .noWS="outside"), "-axis). Pseudobulk processing occurred separately for the domain", tags$sub("SP", .noWS = "outside"), " and domain", tags$sub("CT", .noWS = "outside"), " annotation strategies."),
                              p("Each point represents a single pseudobulk sample corresponding to each domain identified within each donor. Points are colored based on diagnosis group (DX). Crossbar indicates the mean +/- SD."),
                              p("For a complete list of statistics for all tested genes please see our ", tags$a(href = "https://www.google.com", target = "_blank", "Supplementary Tables"), ".")
                            ),
                            mainPanel(
                              card(h3("Whole-tissue DE model"),
                                   plotOutput("whole.tissue"),
                                   helpText("Facets separate results from the two annotation strategies and contain ", tags$em("F", .noWS="outside"), "-test statistics. If ", tags$em("F", .noWS="outside"), "-test adj. ",tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks.")),
                              card(h3("Domain-restricted DE model: domain", tags$sub("SP", .noWS = "outside")),
                                   plotOutput("domain.sp"),
                                   helpText("Facets separate donors by biological sex. Subtitle contains ", tags$em("F", .noWS="outside"), "-test statistics for domain", tags$sub("SP", .noWS = "outside"), " and will be bolded if ", tags$em("F", .noWS="outside"), "-test results indicated a significant difference across diagnosis-by-sex domains (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks."),),
                              card(h3("Domain-restricted DE model: domain", tags$sub("CT", .noWS = "outside")),
                                   plotOutput("domain.ct"),
                                   helpText("Facets separate donors by biological sex. Subtitle contains ", tags$em("F", .noWS="outside"), "-test statistics for domain", tags$sub("CT", .noWS = "outside"), " and will be bolded if ", tags$em("F", .noWS="outside"), "-test results indicated a significant difference across diagnosis-by-sex domains (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " < 0.05, significant diagnosis-by-sex ", tags$em("t", .noWS="outside"), "-test results are shown on plot with black asterisks (adj. ", tags$em("p", .noWS="outside"), " < 0.05). If the ", tags$em("F", .noWS="outside"), "-test adj. ", tags$em("p", .noWS="outside"), " > 0.05, any ", tags$em("t", .noWS="outside"), "-test results with an adj. ", tags$em("p", .noWS="outside"), " < 0.05 are indicated with grey asterisks."),)
                            )
                          )
                 ),
                 tabPanel("Pseudobulk Methods Details",
                          p("For a full description of project methods, please see our ", tags$a(href="http://www.google.com", target="_blank", "publication in Cell", .noWS = "outside"), ". Visit the project ", tags$a(href="https://github.com/LieberInstitute/spatialDLPFC_mdd_bpd/tree/devel", target="_blank", "Github"), " for complete documentation of the code used in this analysis."),
                          layout_columns(
                            accordion(!!!myitems, id = "acc"),
                            img(src="https://github.com/LieberInstitute/Posit_Connect_shiny_apps/blob/devel/dlpfc/DLPFC_MBv_companion/code/www/pseudobulk-pipeline.png?raw=TRUE", height=616, width=550),
                            gap = 100
                          )
                 )
)

server <- function(input, output, session) {

  output$whole.tissue  <- renderPlot({
    if(input$gene_name=="" | is.null(input$gene_name)) {
      EMPTY("whole-tissue")
    } else {
      WHOLE_TISSUE(spe, input$gene_name)
    }
  })

  output$domain.sp <- renderPlot({
    if(input$gene_name=="" | is.null(input$gene_name)) {
      EMPTY("domain-restricted")
    } else {
      DOMAIN_RESTRICTED(spe, input$gene_name, "domain-SP")
    }
  })

  output$domain.ct <- renderPlot({
    if(input$gene_name=="" | is.null(input$gene_name)) {
      EMPTY("domain-restricted")
    } else {
      DOMAIN_RESTRICTED(spe, input$gene_name, "domain-CT")
    }
  })

}

shinyApp(ui = ui, server = server)
