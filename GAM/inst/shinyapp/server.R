library(shiny)
library(data.table)
library(igraph)
library(GAM)
library(markdown)

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

vector2html <- function(v) {
    res.md <- paste("*", names(v), ": ", v, "\n", collapse="")
    return(markdownToHTML(text=res.md, fragment.only=T))
}

renderJs <- function(expr, env=parent.frame(), quoted=FALSE) {
    # Convert the expression + environment into a function
    func <- exprToFunction(expr, env, quoted)
    
    function() {
        val <- func()
        paste0(val, ";", 
               paste(sample(1:20, 10, replace=T), collapse=""))
    }
}

toJsLiteral <- function(x) {
    if (is(x, "numeric")) {
        return(as.character(x))
    } else if (is(x, "character")) {
        return(shQuote(x));
    } else if (is(x, "logical")) {
        return(if (x) "true" else "false")
    } else {
        stop(paste0("can't convert ", x, " to JS literal"))
    }
}

makeJsAssignments  <- function(...) {
    args <- list(...)
    values <- sapply(args, toJsLiteral)
    paste0(names(values), " = ", values, ";\n", collapse="")
}

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
    
    output$geneDESummary <- reactive({
        gene.de <- geneDEInput()
        ids.type <- geneIdsType()
        if (is.null(gene.de)) {
            return("There is no genomic data")
        }
        
        vector2html(c(
            "length" = nrow(gene.de),
            "ID type" = ids.type
            ))
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
    
    output$metDESummary <- reactive({
        met.de <- metDEInput()
        ids.type <- metIdsType()
        if (is.null(met.de)) {
            return("There is no metabolic data")
        }
        
        vector2html(c(
            "length" = nrow(met.de),
            "ID type" = ids.type
            ))
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
    
    output$networkSummary <- reactive({
        es <- esInput()
        net <- es$subnet
        if (is.null(net)) {
            return("There is no built network")
        }
        
        vector2html(c(
            "number of nodes" = length(V(net)),
            "number of edges" = length(E(net))
            ))
    })
    
    output$postProcessingOptions <- reactive({
        es <- esInput()
        makeJsAssignments(
            networkHasReactionsAsNodes = !is.null(es) && !es$reactions.as.edges,
            networkHasReactionsAsEdges = !is.null(es) && es$reactions.as.edges
            )        
    })
    
    output$showModulePanel <- renderJs({
        if (!is.null(esInput())) {
            return("mp = $('#module-panel'); mp[0].scrollIntoView();")
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
    
    output$moduleSummary <- reactive({
        module <- moduleInput()
        if (is.null(module)) {
            return("There is no module yet")
        }
        
        vector2html(c(
            "number of nodes" = length(V(module)),
            "number of edges" = length(E(module))
            ))
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
