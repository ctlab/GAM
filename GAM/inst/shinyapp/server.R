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
        if (is.null(input$gene_de)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        data.table(read.table(input$gene_de$datapath, sep="\t", header=T))
    })
    
    output$gene_de_table <- renderTable({
        data <- geneDEInput()
        if (is.null(data)) {
            return(NULL)
        }
        format(head(data[order(pval)]))
    })
    
    
    metDEInput <- reactive({
        if (is.null(input$met_de)) {
            # User has not uploaded a file yet
            return(NULL)
        }
        
        data.table(read.table(input$met_de$datapath, sep="\t", header=T))
    })
    
    output$met_de_table <- renderTable({
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
        gene.ids <- gene.ids[isolate(input$gene_ids)]
        met.de <- isolate(metDEInput())
        met.ids <- met.ids[isolate(input$met_ids)]
        if (is.null(gene.de) || is.null(met.de)) {
            return(NULL)
        }
        
        makeExperimentSet(network=network, met.de=met.de, gene.de=gene.de, met.ids=met.ids, gene.ids=gene.ids, plot=F)
    })
    
    output$network_summary <- renderPrint({
        es <- esInput()
        es$subnet
    })
    
    moduleInput <- reactive({
        input$find
        met.fdr <- isolate(input$met_fdr)
        gene.fdr <- isolate(input$gene_fdr)
        absent.met.score=isolate(input$absent_met_score)
        
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
    
    output$module_summary <- renderPrint({
        moduleInput()
    })
    
    output$module <- renderGraph({
        moduleInput()
    })
    
})