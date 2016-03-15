#' Plot a module to PDF
#' @param module Module to save
#' @param file File to save plot to
#' @param width Width of a page (in inches)
#' @param height Height of a page (in inches)
#' @examples
#' data(examplesGAM)
#' \dontrun{
#' saveModuleToPdf(module.re, "module.re.pdf")
#' }
#' @export
saveModuleToPdf <- function(module, file, width=15, height=15) {
    pdf(file, width, height)
    plotNetwork(module, attr.label="label", vertex.size=2)
    dev.off()
}

saveModuleToJson <- function(module, outputFilePrefix) {
    graphString <- getModuleJsonString(module)
    f <- file(paste0(outputFilePrefix, ".json"))
    writeLines(graphString, f)
    close(f)
}

#' Get values of multiple atrributes of multiple vertices
#' @param graph Graph
#' @param index Vertices to select
#' @param attrs Vector of attribute names
#' @return List of vectors of attribute values, one element per attribute
#' @examples
#' data(examplesGAM)
#' vertex.attributes <- get.vertex.attributes(module.re)
#' @export
get.vertex.attributes <- function(graph, index=V(graph), attrs=list.vertex.attributes(graph)) {
    data.frame(
        sapply(attrs,
               function(attr) get.vertex.attribute(graph, attr, index),
               simplify=FALSE,
               USE.NAMES=TRUE),
        stringsAsFactors=F)
}

#' Get values of multiple atrributes of multiple edges
#' @param graph Graph
#' @param index Edges to select
#' @param attrs Vector of attribute names
#' @param include.ends Whether include edge ends into the table or not
#' @return List of vectors of attribute values, one element per attribute
#' data(examplesGAM)
#' edge.attributes <- get.edge.attributes(module.re)
#' @export
get.edge.attributes <- function(graph, index=E(graph), attrs=list.edge.attributes(graph), include.ends=FALSE) {
    res <- data.frame(sapply(
            attrs, function(attr) get.edge.attribute(graph, attr, index),
            simplify = FALSE, USE.NAMES = TRUE),
            stringsAsFactors=F)
        
    if (include.ends) {
        ends <- as.data.frame(get.edgelist(graph)[index,], stringsAsFactors=F)
        res <- cbind(from=ends[,1], to=ends[,2], res, stringsAsFactors=F)
    }     
    res
}

#' Converts a module from igraph to a list of nodes and links
#' @param module igraph module 
#' @return list of two elements: list of nodes and list of edges
#' data(examplesGAM)
#' graph.list <- module2list(module.re)
#' @export
module2list <- function(module) {
    getNodeObject <- function(i) {
        return(c(
            list(index = i - 1),
            get.vertex.attributes(module, i)))
    }
    
    getEdgeObject <- function(i) {
        es <- get.edges(module, i)
        return(c(
            list(source=es[1]-1, target=es[2]-1),
            get.edge.attributes(module, i)))
    }
    
    graphObject <- list(
        nodes=lapply(V(module), getNodeObject),
        edges=lapply(E(module), getEdgeObject)
        )
    graphObject
}

#' Get json string for a module
#' @param module Module to convert to JSONstring
#' data(examplesGAM)
#' graph.json <- getModuleJsonString2list(module.re)
#' @export
getModuleJsonString <- function(module) {
    if (!require(rjson)) {
        stop("getModuleJsonString needs rjson module to work")
    }
    graphObject <- module2list(module)    
    return(toJSON(graphObject))
}


#' Save network to an XGMML file
#' @param network Network to save
#' @param file File to save to
#' @param name Name of the network
#' @examples
#' data(examplesGAM)
#' \dontrun{
#' saveModuleToXgmml(module.re, "M1 vs M2", "module.re.xgmml")
#' }
#' @export
saveModuleToXgmml <- function(network, file, name=NULL) {
    if (is.null(name)) {
        name <- deparse(substitute(network))    
    }    
    s <- getGraphXmlString(network, name)
    write(s, file)
}

getGraphXmlString <- function(network, name) {
    res <- c()
    res <- c(res, '<?xml version="1.0"?>\n')
    res <- c(res, paste0('<graph label="', name, '" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">\n'))
    res <- c(res, '  <att name="documentVersion" value="1.1"/>\n')
    res <- c(res, '  <att name="networkMetadata">\n')
    res <- c(res, getRDFXmlString(network, name, indent="    "))
    res <- c(res, '  </att>\n')
    res <- c(res, getNodeXmlStrings(network, indent="  "))
    res <- c(res, getEdgeXmlStrings(network, indent="  "))
    res <- c(res, '</graph>\n')
    
    paste(res, collapse="")
}

xmlNodeString <- function(name, text) {
    paste0("<", name, ">", text, "</", name, ">")
}

getRDFXmlString <- function(network, name, indent="") {
    res <- c()
    res <- c(res, paste0("<rdf:RDF>\n"))
    res <- c(res, paste0("  ", "<rdf:Description rdf:about=\"http://www.cytoscape.org/\">\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:type", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:description", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:identifier", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:date", Sys.time()), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:title", name), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:format", "Cytoscape-XGMML"), "\n"))
    res <- c(res, paste0("  ", "</rdf:Description>\n"))
    res <- c(res, paste0("</rdf:RDF>\n"))
    
    res <- paste0(indent, res)
    paste(res, collapse="")
}


getAttrXmlStrings <- function(attr.values, indent="") {
    attr.xmlStrings <- list() 
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.rtype <- is(attr.values[[i]])[1]
        if (attr.rtype == "character") {
            type <- "string"
        } else if(attr.rtype == "integer") {
            type <- "integer"
        } else if(attr.rtype == "numeric") {
            type <- "real"
            attr.values[[i]] <- sub("Inf", "Infinity", as.character(attr.values[[i]]), fixed=TRUE)
        } else {
            type <- "string"
            attr.values[[i]] <- as.character(attr.values[[i]])
        }
        
        attr.xmlStrings[[i]] <- paste0(indent, "<att type=", "\"", type, "\"", " name=", "\"", attr.names[i], "\"", " value=", "\"", attr.values[[i]], "\"/>\n")
        attr.xmlStrings[[i]][is.na(attr.values[[i]])] <- NA
    }
    names(attr.xmlStrings) <- attr.names
    do.call("cbind", attr.xmlStrings)
}

getNodeXmlStrings <- function(network, indent="") {
    if (length(V(network)) == 0) {
        return(NULL)
    }
    attr.values <- get.vertex.attributes(network)
    attr.xmlStrings <- getAttrXmlStrings(attr.values, paste0(indent, "  "))
    
    if(is.null(V(network)$name))
    {
        V(network)$name <- as.character(V(network))
    }
    node.label <- V(network)$name
    node.id <- as.vector(V(network))
    xmlHeaders <- paste0(indent, 
                         "<node",
                         " label=", "\"", node.label, "\"", 
                         " id=", "\"", node.id, "\"", 
                         ">\n")
    xmlFooter <- paste0(indent, "</node>\n")
    xmlStrings <- paste0(xmlHeaders,
                              apply(attr.xmlStrings, 1, function(x) paste(na.omit(x), collapse="")),
                              xmlFooter)
                                    
    xmlStrings
}

getEdgeXmlStrings <- function(network, indent="") {
    if (length(E(network)) == 0) {
        return(NULL)
    }
    attr.values <- get.edge.attributes(network)
    attr.xmlStrings <- getAttrXmlStrings(attr.values, paste0(indent, "  "))
    
    edgelist.names <- get.edgelist(network, names=TRUE)
    edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
    edgelist.ids <- get.edgelist(network, names=FALSE)
    
    xmlHeaders <- paste0(indent, 
                         "<edge",
                         " label=", "\"", edgelist.names, "\"", 
                         " source=", "\"", edgelist.ids[,1], "\"", 
                         " target=", "\"", edgelist.ids[,2], "\"", 
                         ">\n")
    xmlFooter <- paste0(indent, "</edge>\n")
    xmlStrings <- paste0(xmlHeaders,
                              if (is.null(attr.xmlStrings)) {
                                  NULL
                              } else {
                                  apply(attr.xmlStrings, 1, function(x) paste(na.omit(x), collapse=""))
                              },
                              xmlFooter)
                                    
    xmlStrings
}

    
#' Save module to different formats
#' @param module Module to save
#' @param outputFilePrefix Path to save to (without extension)
#' @param types Vector of file types, "pdf" or one of the supported by BioNet::saveNetwork function
#' @examples
#' data(examplesGAM)
#' \dontrun{
#' saveModuleToPdf(module.re, "module.re", types=c("pdf", "XGMML", "dot"))
#' }
#' @export
saveModule <- function(module, outputFilePrefix, types=c("pdf", "XGMML", "dot")) {
    outdir <- dirname(outputFilePrefix)
    
    if (!file.exists(outdir)) {
        dir.create(outdir)
    }
    
    
    for (type in types) {
        if (type == "pdf") {
            saveModuleToPdf(module, paste0(outputFilePrefix, ".pdf"))
        } else if (type == "XGMML") {
            saveModuleToXgmml(module, file=paste0(outputFilePrefix, ".xgmml"), name=basename(outputFilePrefix))
        } else if (type == "dot") {
            saveModuleToDot(module, file=paste0(outputFilePrefix, ".dot"), name=basename(outputFilePrefix))
        } else {
            saveNetwork(module, file=outputFilePrefix,type=type)
        }
        
    }
}

#' Save network to a graphviz dot file
#' @param network Network to save
#' @param file File to save to
#' @param name Name of the network
#' @examples
#' data(examplesGAM)
#' \dontrun{
#' saveModuleToDot(module.re, "module.re.dot", name="M1 vs M2")
#' }
#' @export
saveModuleToDot <- function(network, file, name=NULL) {
    if (is.null(name)) {
        name <- deparse(substitute(network))    
    }    
    s <- getGraphDotString(network, name)
    write(s, file)
}

getGraphDotString <- function(network, name) {
    res <- c()
    res <- c(res, sprintf('graph "%s" {\n', name))    
    res <- c(res, "outputorder=edgesfirst\n")    
    res <- c(res, "bgcolor=transparent\n")    
    res <- c(res, "pad=\"2,0.25\"\n")   
    res <- c(res, getEdgeDotStrings(network, indent="  "))
    res <- c(res, getNodeDotStrings(network, indent="  "))
    res <- c(res, '}\n')
    
    paste(res, collapse="")
}


getAttrDotStrings <- function(attr.values) {
    attr.dotStrings <- list() 
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.dotStrings[[i]] <- sprintf("%s = \"%s\"", attr.names[i], attr.values[[i]])
        attr.dotStrings[[i]][is.na(attr.values[[i]])] <- NA
    }
    names(attr.dotStrings) <- attr.names
    do.call("cbind", attr.dotStrings)
}

upRamp <- colorRamp(c("#cccccc", "#ff3333"))
downRamp <- colorRamp(c("#cccccc", "green"))

getDotColor <- function(log2FC) {
    if (is.na(log2FC)) {
        return("#7777ff")
    }
    
    if (log2FC > 0) {
        col <- upRamp(min(1, abs(log2FC) / 2))
    } else {
        col <- downRamp(min(1, abs(log2FC) / 2))
    }
    rgb(col[,1], col[,2], col[,3], maxColorValue=255)
}

getDotSize <- function(logPval) {
    logPval[is.na(logPval)] <- 0
    return(pmin(0.2 - logPval/100/0.3, 0.5))
}


nodeShapeMap <- c(met="circle", rxn="square")
edgeStyleMap <- c(main="solid", trans="dashed")

getDotNodeStyleAttributes <- function(attrs) {            
    with(attrs, data.frame(
        label=if (!is.null(attrs$label)) label else "",
        shape=if (!is.null(attrs$nodeType)) nodeShapeMap[nodeType] else "circle",
        fixedsize="true",
        style="filled",            
        width=sapply(logPval, getDotSize),
        fontsize=sapply(logPval, getDotSize) * 45,
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "black",
        fillcolor=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "white"
    ))
}

getDotEdgeStyleAttributes <- function(attrs) {    
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(        
        label=if (!is.null(attrs$label)) label else "",
        style=if (!is.null(attrs$rptype)) edgeStyleMap[rptype] else "solid",
        penwidth=sapply(logPval, getDotSize) * 20,        
        fontsize=sapply(logPval, getDotSize) * 45,        
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "grey"        
    ))
}

getDotTooltip <- function(attr.values) {
    attr.strings <- list() 
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.strings[[i]] <- sprintf("%s: %s", attr.names[i], attr.values[[i]])
    }
    names(attr.strings) <- attr.names
    apply(do.call("cbind", attr.strings), 1, paste0, collapse="&#10;")
}

getNodeDotStrings <- function(network, indent="") {
    if (length(V(network)) == 0) {
        return(NULL)
    }
    attr.values <- get.vertex.attributes(network)
    style.attr.values <- getDotNodeStyleAttributes(attr.values)
    # ignoring technical nodeType and big pathway attributes
    tooltip <- getDotTooltip(attr.values[, !colnames(attr.values) %in% c("pathway", "nodeType")])
    
    attr.dotStrings <- getAttrDotStrings(cbind(style.attr.values, tooltip=tooltip))
    
    if(is.null(V(network)$name))
    {
        V(network)$name <- as.character(V(network))
    }
    
    
    node.label <- V(network)$name
    node.id <- as.vector(V(network))
    node.attrs <- apply(attr.dotStrings, 1, function(x) paste(na.omit(x), collapse=", "))
    nodeStrings <- sprintf("%sn%s [ %s ];\n", indent, node.id, node.attrs)
    nodeStrings
}

getEdgeDotStrings <- function(network, indent="") {
    if (length(E(network)) == 0) {
        return(NULL)
    }
    attr.values <- get.edge.attributes(network, include.ends = TRUE)
    style.attr.values <- getDotEdgeStyleAttributes(attr.values)
    
    # ignoring big pathway attribute
    tooltip.values <- attr.values[, !colnames(attr.values) %in% c("pathway")]
    # for readability replacing node IDs with labels
    tooltip.values$from <- V(network)[tooltip.values$from]$label
    tooltip.values$to <- V(network)[tooltip.values$to]$label
    tooltip <- getDotTooltip(tooltip.values)
    
    attr.dotStrings <- getAttrDotStrings(cbind(style.attr.values, tooltip=tooltip, labeltooltip=tooltip))
    
    edgelist.names <- get.edgelist(network, names=TRUE)
    edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
    edgelist.ids <- get.edgelist(network, names=FALSE)
    
    edge.attrs <- apply(attr.dotStrings, 1, function(x) paste(na.omit(x), collapse=", "))
    edgeStrings <- sprintf("%sn%s -- n%s [ %s ];\n",
                           indent,                           
                           edgelist.ids[,1],
                           edgelist.ids[,2],
                           edge.attrs)
    edgeStrings
}
