graphOutputBinding = new Shiny.OutputBinding();
graphOutputBinding.find = function(scope) {
    return $(scope).find(".graph-output");
};

graphOutputBinding.renderValue = function(el, data) {
    loadGraph(data);
}


Shiny.outputBindings.register(graphOutputBinding, "alserg.graphOutputBinding");



var svg;
var container = "#graph_container";
var force;

positiveFCScale = d3.scale.linear().clamp(true).domain([0,2]).range(["#cccccc","#ff0000"]);
negativeFCScale = d3.scale.linear().clamp(true).domain([0,-2]).range(["#cccccc","#00ff00"]);
function getColor(d) {
    if (d.logFC === parseFloat(d.logFC)) {
        if (d.logFC >= 0) {
            return positiveFCScale(d.logFC);
        } else {
            return negativeFCScale(d.logFC);
        }
    }
    return "#00acad";
}

pvalScale = d3.scale.linear().clamp(true).domain([0, -100]).range([8, 20]);
function getSize(d) {
    if (d.logPval === parseFloat(d.logPval)) {
        return pvalScale(d.logPval);
    }
    return pvalScale(0);
}

window.onload = function() {
    d3.select(window).on("resize", sizeChange);

    sizeChange();


    var width = $(container).width(),
        height = $(window).height();


    force = d3.layout.force()
        .charge(-300)
        .linkDistance(100)
        .size([width, height]);

    svg = d3.select(container).append("svg")
        .attr("width", width)
        .attr("height", height)
        .append("g");


    svg.append("rect")
        .attr("class", "overlay")
        .attr("width", width)
        .attr("height", height);

    zoom = d3.behavior.zoom().scaleExtent([1/8, 8]).on("zoom", zoomed);
    drag = force.drag()
        .on("dragstart", function(d) { d3.event.sourceEvent.stopPropagation(); } );

    svg = svg.call(zoom)
      .append("g");


}

function loadGraph(graph) {
    svg.text(null)
    force.nodes(graph.nodes)
        .links(graph.links)
        .start();

    var link = svg.selectAll(".link")
        .data(graph.links)
        .enter().append("line")
        .attr("class", "link")
        .style("stroke", getColor)
        .style("stroke-width", getSize);

    var node = svg.selectAll(".node")
        .data(graph.nodes)
        .enter()
        .append("g")
        .attr("class", "node")
        .call(drag);

    node.append("circle")
        .attr("cx", 0)
        .attr("cy", 0)
        .attr("r", getSize)
        .style("fill", getColor);

    node.append("text")
        .attr("x", 0)
        .attr("dy", ".35em")
        .attr("text-anchor", "middle")
        .style("font-size", getSize)
        .text(function(d) { return d.label });

    force.on("tick", function() {
        link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });
    });
}

function zoomed() {
  svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
}

function sizeChange() {
    var width = $(container).width(),
        height = $(window).height();
    d3.select("svg")
        .attr("width", width)
        .attr("height", height);
    d3.select("rect")
        .attr("width", width)
        .attr("height", height);
}
