# initial files

WT_GENERATOR <- createCustomTable(WT_SUMMARY, className = "WholeTissueDE", fullName="Whole-tissue DE results")
DR_GENERATOR <- createCustomTable(DR_SUMMARY, className = "DomainRestrictedDE", fullName="Domain-restricted DE results")

initial <- list()

initial[["RowDataTable1"]] <- RowDataTable(Selected="SNAP25", Search="SNAP25", PanelWidth=12L)

initial[["FeatureAssayPlot1"]] <- FeatureAssayPlot(DataBoxOpen=FALSE, PanelWidth=6L, 
                                                   XAxis = "Column data", XAxisColumnData = "dx.sex",
                                                   YAxisFeatureDynamicSource = T, YAxisFeatureSource = "RowDataTable1",
                                                   FacetColumnBy = "Column data", FacetColumnByColData = "annotation",
                                                   ColorBy = "Column data", ColorByColumnData = "diagnosis",
                                                   FontSize = 1.2, HoverInfo = FALSE)

initial[["WholeTissueDE1"]] <- WT_GENERATOR(PanelWidth=6L, RowSelectionSource = "RowDataTable1", RowSelectionDynamicSource = T)

initial[["FeatureAssayPlot2"]] <- FeatureAssayPlot(DataBoxOpen=FALSE, PanelWidth=12L, 
                                                   XAxis = "Column data", XAxisColumnData = "dx.sex",
                                                   YAxisFeatureDynamicSource = T, YAxisFeatureSource = "RowDataTable1",
                                                   FacetColumnBy = "Column data", FacetColumnByColData = "consensus_domains",
                                                   FacetRowBy = "Column data", FacetRowByColData = "annotation",
                                                   ColorBy = "Column data", ColorByColumnData = "diagnosis",
                                                   FontSize = 1.2, HoverInfo = FALSE)

initial[["DomainRestrictedDE1"]] <- DR_GENERATOR(annot_name="domain-SP", PanelWidth=6L, 
                                                 RowSelectionSource = "RowDataTable1", RowSelectionDynamicSource = T)

initial[["DomainRestrictedDE2"]] <- DR_GENERATOR(annot_name = "domain-CT", PanelWidth=6L, 
                                                 RowSelectionSource = "RowDataTable1", RowSelectionDynamicSource = T)



ecm <- ExperimentColorMap(global_continuous = viridis::viridis,
                          colData = list(diagnosis = function(n) c("NTC"= "#7f7f7f", "MDD"= "#FB8861", "BD"= "#9260b2"),
                                         domain_annotation = function(n) c("domain-SP: L1"="#cfa45c", "domain-SP: L2"="#5D9940", "domain-SP: L3/4"="#5095CD", 
                                                                           "domain-SP: L5"="#ddc94e", "domain-SP: L6"="#E45C5F",  "domain-SP: WM"="#D1C4B0",
                                                                           "domain-CT: M/V"="#911223", "domain-CT: Ast"="#cfa45c", "domain-CT: L2/3"="#088F8F", "domain-CT: L4"="#c2cfcf", 
                                                                           "domain-CT: Inb"="#9377AC", "domain-CT: L5"="#ddc94e", "domain-CT: L6"="#E45C5F", "domain-CT: Olg"="#D1C4B0")))


mytour =  data.frame(element=c("#Welcome",".navbar-static-top","#RowDataTable1","#RowDataTable1","#FeatureAssayPlot1","#WholeTissueDE1","#WholeTissueDE1","#FeatureAssayPlot2","#DomainRestrictedDE1","#DomainRestrictedDE2"),
                     intro=c("This page contains pseudobulk data summaries for the Thompson et al. 2026 SRT data.",
                             "To return to this tour at any point, select the question mark icon and 'Click me for a quick tour'.",
                             "Use this table to select features to generate pseudobulk violin plots (<i>FeatureAssayPlots</i>) and populate differential expression tables (<i>WholeTissueDE</i> and <i>DomainRestrictedDE</i>).",
                             "Use the Search bar in the top right (currently containing 'SNAP25') to filter for a regex string within any of the table columns.<br><br>Alternatively or additionally, the table columns can be used to filter.",
                             "The selected feature in the row data table will appear in these pseudobulk <i>FeatureAssayPlots</i>.<br><br>This specific plot displays all pseudobulk samples for domain-SP and domain-CT annotations (facet columns), reflecting the dx-sex comparisons made by whole-tissue DE models.",
                             "These differential expression (DE) results tables will show all of the features from the row data table search criteria (see both <i>SNAP25</i> and <i>SNAP25-AS1</i>). To restrict DE to a single feature use the column filters.",
                             "This specific table summarizes the whole-tissue DE model results for both domain-SP and domain-CT, including the omnibus F-test statistics and post-hoc t-test significance for dx-sex group contrasts.",
                             "This plot stratifies pseudobulk samples to the different domains-SP and domains-CT (facet columns), reflecting comparisons made by domain-restricted DE models.<br><br>Consensus groups are used to facet out domains, highlighting the similarities/differences provided by the domain-SP and domain-CT annotations.",
                             "This table summarizes the domain-restricted DE model results for domain-SP annotations, including the omnibus F-test statistics and post-hoc t-test significance for dx-sex-domain group contrasts.",
                             "This table summarizes the domain-restricted DE model results for domain-CT annotations, including the omnibus F-test statistics and post-hoc t-test significance for dx-sex-domain group contrasts."))
