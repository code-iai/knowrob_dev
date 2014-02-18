
DataVisClient = function(options) {
  var ros = options.ros;
  var containerId = options.containerId;

  var chartHandle = [];

  var rosTopic = new ROSLIB.Topic({
    ros : ros,
    name : '/data_visualisation',
    messageType : 'data_vis_msgs/DataVis'
  });


  rosTopic.subscribe(function(message) {

    if (chartHandle.findIndex(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        }) == -1) {

      var options = {
        data: message.values[0],
        where: containerId,
        label: message.title
      };

      chartHandle.push({
        id: message.id,
        handle: new DonutChart(options)
      });

      chartHandle.find(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        })
        .handle.update(message.values[0]);

    } else {

      chartHandle.find(function (element, index, array) {
          if(element.id == message.id) {return true} else {return false}
        })
        .handle.update(message.values[0]);
    }
  });
}
