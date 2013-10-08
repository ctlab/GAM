library(shiny)
library(data.table)

heinz.py <- "~/lib/heinz/heinz.py"

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

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
    
    geneDEInput <- reactive({
        if (is.null(input$geneDE)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        data.table(read.table(input$geneDE$datapath, sep="\t", header=T))
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
        
        data.table(read.table(input$metDE$datapath, sep="\t", header=T))
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
        gene.ids <- gene.ids[isolate(input$geneIds)]
        met.de <- isolate(metDEInput())
        met.ids <- met.ids[isolate(input$metIds)]
        if (is.null(gene.de) || is.null(met.de)) {
            return(NULL)
        }
        
        makeExperimentSet(network=network, met.de=met.de, gene.de=gene.de, met.ids=met.ids, gene.ids=gene.ids, plot=F)
    })
    
    output$networkSummary <- renderPrint({
        es <- esInput()
        es$subnet
    })
    
    moduleInput <- reactive({
        input$find
        met.fdr <- isolate(input$metFDR)
        gene.fdr <- isolate(input$geneFDR)
        absent.met.score=isolate(input$absentMetScore)
        
        es <- isolate(esInput())
        
        if (is.null(es)) {
            return(NULL)
        }

        findModule(es,
                    met.fdr=met.fdr,
                    gene.fdr=gene.fdr,
                    absent.met.score=absent.met.score,
                    heinz.py=heinz.py)
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
