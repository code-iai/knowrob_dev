function BarChart (options) {
  options = options || {};
  var w = 200;
  var h = 200;
  var data = options.data || [];
  var where = options.where;
  var label = options.label;
  //console.log(data);

    //setup the svg
    var svg = d3.select(where).append("svg:svg")
        .attr("width", w+100)
        .attr("height", h+100)
    svg.append("svg:rect")
        .attr("width", "100%")
        .attr("height", "100%")
        .attr("stroke", "#000")
        .attr("fill", "none");

    var vis = svg.append("svg:g")
        //.attr("id", "barchart")
        .attr("transform", "translate(50,50)");

  this.remove = function() {
    svg.remove();
  }

  this.update = function(rawdata) {

    var data = rawdata.value2;
    console.log(data);
    max = d3.max(data, function(d) {return parseInt(d)});
    console.log(max);

    //nice breakdown of d3 scales
    //http://www.jeromecukier.net/blog/2011/08/11/d3-scales-and-color/
    x = d3.scale.linear()
        .domain([0, max])
        .range([0, w]);

    y = d3.scale.ordinal()
        .domain(d3.range(data.length))
        .rangeBands([0, h], .2);


    //var vis = d3.select("#barchart")
    
    //a good written tutorial of d3 selections coming from protovis
    //http://www.jeromecukier.net/blog/2011/08/09/d3-adding-stuff-and-oh-understanding-selections/
    var bars = vis.selectAll("rect.bar")
        .data(data);

    //update
    bars
        .attr("fill", "#0a0")
        .attr("stroke", "#050");

    //enter
    bars.enter()
        .append("svg:rect")
        .attr("class", "bar")
        .attr("fill", "#0a0") // #800
        .attr("stroke", "#050"); // #800


    //exit 
    bars.exit()
    .transition()
    .duration(300)
    .ease("exp")
        .attr("width", 0)
        .remove();


    bars
        .attr("stroke-width", 4)
    .transition()
    .duration(300)
    .ease("quad")
        .attr("width", x)
        .attr("height", y.rangeBand())
        .attr("transform", function(d,i) {
            return "translate(" + [0, y(i)] + ")"
        });

  }

  //this.update(data);

}
