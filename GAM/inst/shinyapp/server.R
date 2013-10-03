library(shiny)
library(data.table)

heinz.py <- "~/lib/heinz/heinz.py"

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {
    
    geneDEInput <- reactive({
        if (is.null(input$gene.de)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        data.table(read.table(input$gene.de$datapath, sep="\t", header=T))
    })
    
    output$gene.de.table <- renderTable({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(head(data[order(pval)]))
    })
    
    
    metDEInput <- reactive({
        if (is.null(input$met.de)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        data.table(read.table(input$met.de$datapath, sep="\t", header=T))
    })
    
    output$met.de.table <- renderTable({
        data <- metDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(head(data[order(logFC)]))
    })
    
    esInput <- reactive({
        input$preprocess
        network <- networks[[isolate(input$network)]]
        gene.de <- isolate(geneDEInput())
        gene.ids <- gene.ids[isolate(input$gene.ids)]
        met.de <- isolate(metDEInput())
        met.ids <- met.ids[isolate(input$met.ids)]
        if (is.null(gene.de) || is.null(met.de)) {
            return(NULL)
        }
        
        makeExperimentSet(network=network, met.de=met.de, gene.de=gene.de, met.ids=met.ids, gene.ids=gene.ids, plot=F)
    })
    
    output$network.summary <- renderPrint({
        es <- esInput()
        es$subnet
    })
    
    moduleInput <- reactive({
        input$find
        met.fdr <- isolate(input$met.fdr)
        gene.fdr <- isolate(input$gene.fdr)
        absent.met.score=isolate(input$absent.met.score)
        
        es <- isolate(esInput())
        
        if (is.null(es)) {
            return(NULL)
        }

        findModules(es,
                    met.fdr=met.fdr,
                    gene.fdr=gene.fdr,
                    absent.met.score=absent.met.score,
                    heinz.py=heinz.py)
    })
    
    output$module.summary <- renderPrint({
        moduleInput()
    })
    
})