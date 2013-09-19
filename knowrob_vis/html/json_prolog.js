
function JsonProlog(){
  
  var ros = new ROSLIB.Ros({
    url : 'ws://localhost:9090'
  });

  this.jsonQuery = function(query) {

      var qid = this.makeid();

      var jsonPrologQueryClient = new ROSLIB.Service({
        ros : ros,
        name : '/json_prolog/simple_query',
        serviceType : 'json_prolog/PrologQuery'
      });

      
      // send query
      var request = new ROSLIB.ServiceRequest({
        mode : 0,
        id : qid,
        query : query
      });

      console.log(request);
    
      jsonPrologQueryClient.callService(request, function(result) {

        // collect results
        var jsonPrologNextResultClient = new ROSLIB.Service({
          ros : ros,
          name : '/json_prolog/next_solution',
          serviceType : 'json_prolog/PrologNextSolution'
        });

      console.log(result);
    
        var request2 = new ROSLIB.ServiceRequest({
          id : qid
        });

        jsonPrologNextResultClient.callService(request2, function(result) {

          console.log(result);
          
//           var res = JSON.parse(result.solution);
// 
//           for(var obj in res) {
//             console.log(obj);
//           }

          
          var jsonPrologFinishClient = new ROSLIB.Service({
            ros : ros,
            name : '/json_prolog/finish',
            serviceType : 'json_prolog/PrologFinish'
          });

          request3 = new ROSLIB.ServiceRequest({
            id : qid
          });
          jsonPrologFinishClient.callService(request3, function(e) { });

        });
      });

  };

  this.makeid = function() {

    var text = "";
    var possible = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

    for( var i=0; i < 8; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
  };
}