function(input, output, session) {
  spe <- readRDS("data/MBv-shiny_pseudobulk-spe_both-annotations.rds")
  
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
