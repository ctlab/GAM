networkFromSif <- function(file.basename) {
    print("Loading network...")
    network <- newEmptyObject()
    class(network) <- "network"
    network$edges <- read.table(paste(file.basename, "sif", sep="."), header=FALSE, colClasses=c("character", "factor", "character"), comment.char="#", col.names=c("u", "e", "v"))
    
    for (f in list.files(dirname(file.basename), paste(basename(file.basename), "_\\w+.NA", sep=""), full.names=T)) {
        meta_type.cur <- sub(paste("^", file.basename, "_(\\w+).NA", sep=""), "\\1", f)        
        print(paste("Loading attribute", meta_type.cur))
        meta_temp <- unique(read.table(f, skip=1, comment.char="#", colClasses="character"))
        meta.cur <- data.frame(meta_temp[,1], meta_temp[,3], stringsAsFactors=F)
        colnames(meta.cur) <- c("v", meta_type.cur)
        if (is.null(network$meta)) {
            network$meta <- meta.cur
        } else {
            network$meta <- merge(network$meta, meta.cur, by="v", all=T)
        }        
    }
    
    return(network)
}

networkToSif <- function(network, file.basename) {
    write.table(network$edges[,c("u", "e", "v")], paste(file.basename, "sif", sep="."),
                col.names=F, row.names=F, quote=F, sep="\t")
    for (m in colnames(network$meta)) {
        if (m != "v") {
            x <- network$meta[,m, drop=F]
            rownames(x) <- network$meta$v
            x <- na.omit(x)
            f <- paste(file.basename, "_", colnames(x), ".NA", sep="")
            write.table(x, f, row.names=T, quote=T, col.names=T, sep=" = ")
        }
    }
}

networkToTab <- function(edges, file.name) {
    write.table(edges[,c("u", "v")], file.name, col.names=F, row.names=F, quote=F, sep="\t")
}


#' @export
saveModuleToPdf <- function(module, file) {
    pdf(file, width=15, height=15)
    plotNetwork(module, attr.label="label", vertex.size=2)
    dev.off()
}

saveModuleToJson <- function(module, outputFilePrefix) {
    graphString <- getModuleJsonString(module)
    f <- file(paste0(outputFilePrefix, ".json"))
    writeLines(graphString, f)
    close(f)
}

#' @export
get.vertex.attributes <- function(graph, index=V(graph), attrs=list.vertex.attributes(graph)) {
    sapply(attrs,
           function(attr) get.vertex.attribute(graph, attr, index),
           simplify=F,
           USE.NAMES=T)
}

#' @export
get.edge.attributes <- function(graph, index=E(graph), attrs=list.edge.attributes(graph)) {
    sapply(attrs, 
           function(attr) get.edge.attribute(graph, attr, index),
           simplify=F,
           USE.NAMES=T)
}

#' @export
module2list <- function(module) {
    getNodeObject <- function(i) {
        return(c(
            list(index = i - 1),
            get.vertex.attributes(imodule, i)))
    }
    
    getEdgeObject <- function(i) {
        es <- get.edges(imodule, i)
        return(c(
            list(source=es[1]-1, target=es[2]-1),
            get.edge.attributes(imodule, i)))
    }
    
    graphObject <- list(
        nodes=lapply(V(imodule), getNodeObject),
        links=lapply(E(imodule), getEdgeObject)
        )
    graphObject
}

#' @export
#' @importFrom rjson toJSON
getModuleJsonString <- function(module) {
    graphObject <- module2list(module)    
    return(toJSON(graphObject))
}

#' @export
#' @importFrom XML append.xmlNode addAttributes
saveModuleToXgmml <- function(network, name, file) {
    require(XML)
    top <- .XGMML.destription(name=name)
    print("...adding nodes")
    # append nodes
    nodes <- .XGMML.nodes(network=network)
    top <- append.xmlNode(top, nodes) 
    print("...adding edges")
    # append edges
    edges <- .XGMML.edges(network=network)
    top <- append.xmlNode(top, edges) 
    # save to file as xgmml
    print("...writing to file")
    saveXML(top, file, encoding="UTF-8")
    if("package:XML" %in% search()){detach("package:XML")}
}
# internal method to create the first part of the XGMML-file, description
.XGMML.destription <- function(name)
{
    # create top node
    # top part of xml
    top <- xmlNode("graph", attrs = c(label=name, "xmlns:dc"="http://purl.org/dc/elements/1.1/", "xmlns:xlink"="http://www.w3.org/1999/xlink", "xmlns:rdf"="http://www.w3.org/1999/02/22-rdf-syntax-ns#", "xmlns:cy"="http://www.cytoscape.org", xmlns="http://www.cs.rpi.edu/XGMML"))
    top <- append.xmlNode(top, xmlNode("att", attrs=c(name="documentVersion", value="1.1")))
    
    d <- xmlNode("rdf:Description", attrs=c("rdf:about"="http://www.cytoscape.org/"))
    d <- append.xmlNode(d,  xmlNode("dc:type", "Protein-Protein Interaction"))
    d <- append.xmlNode(d,  xmlNode("dc:description", "N/A"))
    d <- append.xmlNode(d,  xmlNode("dc:identifier", "N/A"))
    d <- append.xmlNode(d,  xmlNode("dc:date", Sys.time()))
    d <- append.xmlNode(d,  xmlNode("dc:title", name))
    d <- append.xmlNode(d,  xmlNode("dc:format", "BioNet-Cytoscape-XGMML"))
    
    c <- xmlNode("att", attrs=c(name="networkMetadata"), xmlNode("rdf:RDF", d))
    top <- append.xmlNode(top, c)
    return(top)
}

# internal method for the addition of nodes to xml
#' @importFrom igraph get.vertex.attribute
.XGMML.nodes <- function(network)
{
    # create node-nodes
    c.node <- rep("node", length(V(network)-1))
    nodes <- lapply(c.node, xmlNode)
    
    # create node attributes
    attrib <- list.vertex.attributes(network)
    node.attribs <- matrix(data=NA, nrow=length(attrib), ncol=length(V(network)))
    for(i in 1:length(attrib))
    {
        attrib.values <- get.vertex.attribute(network, attrib[i])
        if(is(attrib.values)[1] == "character") {
            type <- "string"
        } else if(is(attrib.values)[1] == "integer") {
            type <- "integer"
        } else if(is(attrib.values)[1] == "numeric") {
            type <- "real"
        } else {
            type <- "string"
            attrib.values <- as.character(attrib.values)
        }
        node.attribs[i,] <- paste("att type=", "\"", type, "\"", " name=", "\"", attrib[i], "\"", " value=", "\"", attrib.values, "\"", sep="")
        node.attribs[i, is.na(attrib.values)] <- NA
    }
    node.attribs.tmp <- matrix(lapply(node.attribs, xmlNode), nrow = length(attrib), ncol = length(V(network)))
    node.attribs.tmp[is.na(node.attribs)] <- NA
    node.attribs <- node.attribs.tmp
    if(is.null(V(network)$name))
    {
        V(network)$name <- as.character(V(network))
    }
    node.label <- V(network)$name
    node.id <- as.vector(V(network))
    
    # append node attributes
    for(i in 1:length(V(network)))
    {
        nodes[[i]] <- addAttributes(nodes[[i]], label = node.label[i], id=node.id[i])
        nodes[[i]] <- append.xmlNode(nodes[[i]], node.attribs[,i][!is.na(node.attribs[,i])])
    }
    
    return(nodes)
}

# internal method for the addition of edges to XGMML
#' @importFrom igraph get.edge.attribute get.edgelist list.edge.attributes
.XGMML.edges <- function(network)
{
    # create edge-nodes
    c.edge <- rep("edge", length(E(network)-1))
    edges <- lapply(c.edge, xmlNode)
    
    edgelist.names <- get.edgelist(network, names=TRUE)
    edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
    edgelist.ids <- get.edgelist(network, names=FALSE)
    
    # create edge attributes
    attrib <- list.edge.attributes(network)
    edge.attribs <- matrix(data=NA, nrow=length(attrib), ncol=length(E(network)))
    for(i in 1:length(attrib))
    {
        attrib.values <- get.edge.attribute(network, attrib[i])
        if(is(attrib.values)[1] == "character") {
            type <- "string"
        } else if(is(attrib.values)[1] == "integer") {
            type <- "integer"
        } else if(is(attrib.values)[1] == "numeric") {
            type <- "real"
        } else {
            type <- "string"
            attrib.values <- as.character(attrib.values)
        }
        edge.attribs[i,] <- paste("att type=", "\"", type, "\"", " name=", "\"", attrib[i], "\"", " value=", "\"", attrib.values, "\"", sep="")
        edge.attribs[i, is.na(attrib.values)] <- NA
    }
    edge.attribs.tmp <- matrix(lapply(edge.attribs, xmlNode), nrow = length(attrib), ncol = length(E(network)))
    edge.attribs.tmp[is.na(edge.attribs)] <- NA
    edge.attribs <- edge.attribs.tmp
    
    # append edge attributes
    for(i in 1:length(E(network)))
    {
        edges[[i]] <- addAttributes(edges[[i]], label=edgelist.names[i], source=edgelist.ids[i,1], target=edgelist.ids[i,2])
        edges[[i]] <- append.xmlNode(edges[[i]], edge.attribs[,i][!is.na(edge.attribs[,i])])
    }
    
    return(edges)
}

    
#' Save module to different formats
#' @param module Module to save
#' @param outputFilePrefix Path to save to (without extension)
#' @param types Vector of file types, "pdf" or one of the supported by BioNet::saveNetwork function
#' @export
saveModule <- function(module, outputFilePrefix, types=c("pdf", "XGMML")) {
    # :ToDO: fix saving to XGMML (NAs, escaping, trivial modules)
    outdir <- dirname(outputFilePrefix)
    
    if (!file.exists(outdir)) {
        dir.create(outdir)
    }
    
    
    #module <- subNetwork(nodes(module), network=subnet)
    module    

    
    for (type in types) {
        if (type == "pdf") {
            saveModuleToPdf(module, paste0(outputFilePrefix, ".pdf"))
        } else if (type == "XGMML") {
            saveModuleToXgmml(module, name=basename(outputFilePrefix), paste0(outputFilePrefix, ".xgmml"))
        } else {
            saveNetwork(module, file=outputFilePrefix,type=type)
        }
        
    }
    
    
}
