library(shiny)

# Define UI for application that plots random distributions 
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Shiny GAM"),
    
    # Sidebar with a slider input for number of observations
    sidebarPanel(
        selectInput("network", "Select a network:",
                    choices=names(networks)),
        fileInput("gene.de", label="Select file with gene p-values"),
        selectInput("gene.ids", "Select type of gene IDs:",
                    choices=names(gene.ids)),
        fileInput("met.de", label="Select file with metabolite p-values"),
        selectInput("met.ids", "Select type of metabolite IDs:",
                    choices=names(met.ids)),
        br(),
        actionButton("preprocess", "Preprocess"),
        numericInput("gene.fdr", "FDR for genes", 1e-6, min=1e-100, max=1, step=1e-200),
        numericInput("met.fdr", "FDR for metabolites", 1e-6, min=1e-100, max=1, step=1e-200),
        numericInput("absent.met.score", "Score for absent metabolites", -20),
        br(),
        actionButton("find", "Find module"),
        br()
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
        tableOutput("gene.de.table"),
        tableOutput("met.de.table"),
        textOutput("network.summary"),
        textOutput("module.summary"),
        htmlOutput("graph.container")
    )
))