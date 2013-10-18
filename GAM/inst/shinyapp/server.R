library(shiny)
library(data.table)
library(igraph)
library(GAM)

data(kegg.mouse.network)
data(kegg.human.network)
networks <- list(
    "Mouse musculus"=kegg.mouse.network,
    "Homo sapiens"=kegg.human.network)

data(gene.id.map)
data(met.id.map)

heinz.py <- "/usr/local/lib/heinz/heinz.py"
solver <- heinz.solver(heinz.py=heinz.py, timeLimit=60)

renderGraph <- function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)
    
    function() {
        val <- func()
        if (is.null(val)) {
            return(list(nodes=list(), links=list()));
        }
        module2list(val)
    }
}

necessary.de.fields <- c("ID", "pval", "logFC")

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
    
    geneDEInput <- reactive({
        if (is.null(input$geneDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        res <- data.table(read.table(input$geneDE$datapath, sep="\t", header=T, stringsAsFactors=F))
        if (!all(necessary.de.fields %in% names(res))) {
            stop(paste0("Genomic differential expression data should contain at least these fields: ", 
                        paste(necessary.de.fields, collapse=", ")))
        }
        res
    })
    
    geneIdsType <- reactive({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        res <- getIdsType(data$ID, gene.id.map)
        if (length(res) != 1) {
            stop("Can't determine type of IDs for genes")
        }
        res
    })
    
    output$geneIdsType <- reactive({
        geneIdsType()
        
    })
    
    output$geneDETable <- renderTable({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(head(data[order(pval)]))
    })
    
    
    metDEInput <- reactive({
        if (is.null(input$metDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        res <- data.table(read.table(input$metDE$datapath, sep="\t", header=T, stringsAsFactors=F))
        if (!all(necessary.de.fields %in% names(res))) {
            stop(paste0("Metabolic differential expression data should contain at least these fields: ", 
                        paste(necessary.de.fields, collapse=", ")))
        }
        res
    })
    
    metIdsType <- reactive({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        res <- getIdsType(data$ID, met.id.map)
        if (length(res) != 1) {
            stop("Can't determine type of IDs for metabolites")
        }
        res
    })
    
    output$metIdsType <- reactive({
        metIdsType()
    })
    
    output$metDETable <- renderTable({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(head(data[order(pval)]))
    })
    
    esInput <- reactive({
        input$preprocess
        network <- networks[[isolate(input$network)]]
        gene.de <- isolate(geneDEInput())
        gene.ids <- isolate(geneIdsType())
        met.de <- isolate(metDEInput())
        met.ids <- isolate(metIdsType())
        if (is.null(gene.de) && is.null(met.de)) {
            return(NULL)
        }
        
        reactions.as.edges = isolate(input$reactionsAs) == "edges"
        collapse.reactions = isolate(input$collapseReactions)
        use.rpairs = isolate(input$useRpairs)
        
        makeExperimentSet(
            network=network,
            met.de=met.de, gene.de=gene.de,
            met.ids=met.ids, gene.ids=gene.ids,
            reactions.as.edges=reactions.as.edges,
            collapse.reactions=collapse.reactions,
            use.rpairs=use.rpairs,
            plot=F)
    })
    
    output$networkSummary <- renderPrint({
        es <- esInput()
        es$subnet
    })
    
    output$showModulePanel <- reactive({
        if (!is.null(esInput())) {
            return(paste0("mp = $('#module-panel'); mp[0].scrollIntoView();",
                          paste(sample(1:20, 10, replace=T), collapse="")))
        }
        # return("mp = $('#module-panel'); mp.hide();")
        return("")
    })
    
    rawModuleInput <- reactive({
        print(input$find)
        met.fdr <- isolate(input$metFDR)
        gene.fdr <- isolate(input$geneFDR)
        absent.met.score=isolate(input$absentMetScore)
        absent.rxn.score=isolate(input$absentRxnScore)
        
        es <- isolate(esInput())
        
        if (is.null(es)) {
            return(NULL)
        }

        findModule(es,
                    met.fdr=met.fdr,
                    gene.fdr=gene.fdr,
                    absent.met.score=absent.met.score,
                    absent.rxn.score=absent.rxn.score,
                    solver=solver)
    })
    
    moduleInput <- reactive({
        module <- rawModuleInput()
        if (is.null(module)) {
            return(NULL)
        }
        es <- isolate(esInput())
        
        if (es$reactions.as.edges) {
            if (isolate(input$useRpairs)) {
                if (input$addTransPairs) {
                    module <- addTransEdges(module, es)
                }
            }
        } else {
            if (input$addMetabolitesForReactions) {
                module <- addMetabolitesForReactions(module, es)
            }
            if (input$addInterconnections) {
                module <- addInterconnections(module, es)
            }
            module <- addNormLogFC(module)
            
            if (input$removeHangingNodes) {
                module <- removeHangingNodes(module)
            }
            
            if (input$removeSimpleReactions) {
                module <- removeSimpleReactions(module, es)
            }
            module <- expandReactionNodeAttributesToEdges(module)
        }
            
        module
    })
    
    output$moduleSummary <- renderPrint({
        moduleInput()
    })
    
    output$module <- renderGraph({
        moduleInput()
    })
    
    output$downloadXGMML <- downloadHandler(
        filename = "module.xgmml",
        content = function(file) {
            saveModuleToXgmml(moduleInput(), "module", file)
        })
})
