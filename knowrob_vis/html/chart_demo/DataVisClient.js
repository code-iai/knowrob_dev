
function DataVisClient(options) {
  var ros = options.ros;
  var containerId = options.containerId;
  var topic = options.topic;
  var width = options.width || 300;
  var height = options.height || 300;

  var chartHandle = [];

  var rosTopic = new ROSLIB.Topic({
    ros : ros,
    name : topic,
    messageType : 'data_vis_msgs/DataVis'
  });


  rosTopic.subscribe(function(message) {
    if (chartHandle.findIndex(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        }) == -1) {

      var options = {
        data: message.values[0],
        where: containerId,
        label: message.title,
        width: width,
        height: height,
        radius: height*3/10,
        innerRadius: height*4/30
      };

      if (message.type == 0) {
        chartHandle.push({
          id: message.id,
          handle: new DonutChart(options)
        });

        chartHandle.find(function (element, index, array) {
            if(element.id == message.id) {return true} else {return false}
          })
          .handle.update(message.values[0]);
      } else if (message.type == 1) {

        chartHandle.push({
          id: message.id,
          handle: new BarChart(options)
        });

        chartHandle.find(function (element, index, array) {
            if(element.id == message.id) {return true} else {return false}
          })
          .handle.update(message.values[0]);
      }

    } else if (message.values[0].value2.length == 0) {

      chartHandle.find(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        })
        .handle.remove();
      
      chartHandle.splice(chartHandle.findIndex(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        }), 1);

    } else {

      chartHandle.find(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        })
        .handle.update(message.values[0]);
    }
    console.log(chartHandle);
  });
}
