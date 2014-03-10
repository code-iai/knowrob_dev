function Control(options) {
  var that = this;
  var sliderLow = options.sliderLow || "";
  var sliderHigh = options.sliderHigh;
  var initButton = options.initButton;
  var startTime = 0;
  var endTime = 0;

  this.init = function() {

    var startEndQuery = "findall(Time, task_start(T, Time), List), sort(List, Sorted), nth0(0, Sorted, X), last(Sorted, Y), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Firststring, X), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Laststring, Y), atom_number(Firststring, First), atom_number(Laststring, Last)";

    var prolog = new JsonProlog({raw: true});
    prolog.jsonQuery(startEndQuery, function(result){
      prolog.finishClient();

      startTime = result.First;
      endTime = result.Last;
      console.log(startTime);
      console.log(endTime);

      var range = endTime - startTime;

      that.update(0);

      document.getElementById(sliderHigh).value = 0;  
      document.getElementById(sliderHigh).max = range;
      document.getElementById(sliderHigh).onchange = function(x) {
        console.log(x.explicitOriginalTarget.value);
        that.update(x.explicitOriginalTarget.value);
      };

    });

    //prolog.finishClient();
    //var startend =
    //that.getStartEnd();

    //var startTime = startend.start;
    //var endTime = startend.end;

    /*var range = endTime - startTime;

    that.update(startTime);
  
    document.getElementById(sliderHigh).max = range;
    document.getElementById(sliderHigh).onchange = function(x) {
      console.log(x.explicitOriginalTarget.value);
      that.update(x.explicitOriginalTarget.value);
    };*/

  }

  // get start and end time of loaded experiment
  /*this.getStartEnd = function() {
    var startEndQuery = "findall(Time, task_start(T, Time), List), sort(List, Sorted), nth0(0, Sorted, X), last(Sorted, Y), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Firststring, X), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Laststring, Y), atom_number(Firststring, First), atom_number(Laststring, Last)";


    //var usefulNameForThis = this;

    var prolog = new JsonProlog({raw: true});
    prolog.jsonQuery(startEndQuery, function(result){
      startTime = result.First;
      endTime = result.Last;
      console.log(startTime);
      console.log(endTime);
    });

    prolog.finishClient();
    //console.log(startTime);
    //console.log(endTime);
    //return {
    //  start: startTime,
    //  end: endTime
    //}
  }*/

  this.update = function(timePoint) {

    var queryErrorOverall = "Labels = ['ManipulationPoseUnreachable', 'ObjectNotFound', 'ManipulationFailed', 'LocationNotReached', 'ObjectLost', 'ManipulationPoseOccupied'],findall(N, (member(Error, Labels),findall(E, (failure_class(E, C), C=knowrob:Error), List), length(List, N)),Occurence),add_diagram(overallerrors, 'error distribution for whole experiment', barchart, xlabel, ylabel, '210', '210', '12px', [[Labels,Occurence]])";

    var oldfilequeryErrorOverall = "Labels = ['MANIPULATION-POSE-UNREACHABLE', 'OBJECT-NOT-FOUND', 'MANIPULATION-FAILED', 'LOCATION-NOT-REACHED', 'OBJECT-LOST', 'MANIPULATION-POSE-OCCUPIED'],findall(N,(member(Error, Labels),findall([E,B], (failure_task(E, C), failure_attribute(E, rdfs:'label', X), [X]=[literal(type(A,B))], B=Error), List), length(List, N)),Occurence),add_diagram(overallerrors, 'error distribution for whole experiment', barchart, xlabel, ylabel, '210', '210', '12px', [[Labels,Occurence]])";

    var actualTime = parseInt(timePoint) + parseInt(startTime);

    var oldfilequeryErrorsUntilNow = "Time = " + actualTime.toString() + ",Labels = ['MANIPULATION-POSE-UNREACHABLE', 'OBJECT-NOT-FOUND', 'MANIPULATION-FAILED', 'LOCATION-NOT-REACHED', 'OBJECT-LOST', 'MANIPULATION-POSE-OCCUPIED'],findall(N,(member(Error, Labels),findall([E,B], (failure_task(E, C), failure_attribute(E, rdfs:'label', X), [X]=[literal(type(A,B))], B=Error, failure_attribute(E, knowrob:startTime, Y), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Etimestring, Y), atom_number(Etimestring, Etime), Etime =< Time), List), length(List, N)),Occurence),add_diagram(errorsuntilnow, 'error distribution until current time', barchart, xlabel, ylabel, '210', '210', '12px', [[Labels,Occurence]])";

    var queryErrorsUntilNow = "Time = " + actualTime.toString() + ",Labels = ['ManipulationPoseUnreachable', 'ObjectNotFound', 'ManipulationFailed', 'LocationNotReached', 'ObjectLost', 'ManipulationPoseOccupied'],findall(N, (member(Error, Labels),findall(E, (failure_class(E, C), C=knowrob:Error, failure_attribute(E, knowrob:startTime, Y), string_concat('http://ias.cs.tum.edu/kb/cram_log.owl#timepoint_', Etimestring, Y), atom_number(Etimestring, Etime), Etime =< Time), List), length(List, N)),Occurence),add_diagram(errorsuntilnow, 'error distribution until current time', barchart, xlabel, ylabel, '210', '210', '12px', [[Labels,Occurence]])";
    

    console.log("update world!");

    var prolog = new JsonProlog({raw: true});
    prolog.jsonQuery(oldfilequeryErrorOverall, function(result){
      console.log("update all error diagram!");
    });

    prolog.finishClient();

    var prolog = new JsonProlog({raw: true});
    prolog.jsonQuery(oldfilequeryErrorsUntilNow, function(result){
      console.log("update errors until now!");
    });

    prolog.finishClient();

  }

  document.getElementById(initButton).onclick = function() {
    document.getElementById(initButton).disabled = true;
    that.init();
  };
}
