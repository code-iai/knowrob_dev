/** <module> knowrob_vis

  Description:
    Module providing visualisation capabilities


  Copyright (C) 2013 by Moritz Tenorth

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

@author Moritz Tenorth
@license GPL
*/

:- module(knowrob_vis,
    [
      visualisation_canvas/0,
      clear_canvas/0,
      add_object/1,
      add_object/2,
      add_object_with_children/1,
      add_object_with_children/2,
      remove_object/1,
      remove_object_with_children/1,
      highlight_object/1,
      highlight_object/2,
      highlight_object/3,
      highlight_object/6,
      highlight_object_with_children/1,
      highlight_object_with_children/2,
      reset_highlight/0,
      diagram_canvas/0,
      add_diagram/9,
      remove_diagram/1,
      clear_diagram_canvas/0,
      show_trajectory/3,
      remove_trajectory/0
    ]).

:- use_module(library('semweb/rdfs')).
:- use_module(library('semweb/rdf_db')).
:- use_module(library('semweb/rdfs_computable')).
:- use_module(library('jpl')).


:- rdf_meta add_object(r),
            add_object(r,r),
            add_object_with_children(r),
            add_object_with_children(r,r),
            remove_object(r),
            remove_object_with_children(r),
            highlight_object(r),
            highlight_object(r,?),
            highlight_object(r,?,?),
            highlight_object(r,?,?,?,?,?),
            highlight_object_with_children(r),
            highlight_object_with_children(r,?),
            add_diagram(+,+,+,+,+,+,+,+,+),
            remove_diagram(+), % TODO: not sure if this is correct
            show_trajectory(+,+,+).


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Visualization canvas management
%

%% visualisation_canvas is det.
%
% Launch the visualization canvas
%
:- assert(v_canvas(fail)).
visualisation_canvas :-
    v_canvas(fail),
    jpl_new('org.knowrob.vis.MarkerVisualization', [], Canvas),
    retract(v_canvas(fail)),
    assert(v_canvas(Canvas)),!.
visualisation_canvas(Canvas) :-
    v_canvas(Canvas).


%% clear_canvas is det.
%
% Completely clears the scene
%
clear_canvas :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'clear', [], _).



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Add / remove objects
%

%% add_object(+Identifier) is nondet.
%
% Add object to the scene
%
% @param Identifier Object identifier, eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
add_object(Identifier) :-
    get_timepoint(Time),
    add_object(Identifier, Time).


%% add_object(+Identifier, +Time) is nondet.
%
% Add object to the scene with its position at time 'Time'
%
% @param Identifier Object identifier, eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
add_object(Identifier, Time) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'addObject', [Identifier, Time], _).



%% add_object_with_children(+Identifier)
%
% Adds objects to the scene, including all items that are reachable via knowrob:properPhysicalPartTypes
% or via knowrob:describedInMap
%
% @param Identifier eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
add_object_with_children(Identifier) :-
    get_timepoint(Time),
    add_object_with_children(Identifier, Time).


%% add_object_with_children(+Identifier, +Time)
%
% Adds objects to the scene, including all items that are reachable via knowrob:properPhysicalPartTypes
% or via knowrob:describedInMap, with their positions at time 'Time'
%
% @param Identifier eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
add_object_with_children(Identifier, Time) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'addObjectWithChildren', [Identifier, Time], _).


%% remove_object(+Identifier) is det.
%
% Remove object from the scene
%
% @param Identifier Object identifier, eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
remove_object(Identifier) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'removeObject', [Identifier], _).


%% remove_object_with_children(+Identifier) is det.
%
% Removes objects from the scene, including all items that are reachable via knowrob:properPhysicalPartTypes
% or via knowrob:describedInMap
%
% @param Identifier eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
remove_object_with_children(Identifier) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'removeObjectWithChildren', [Identifier], _).



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Highlighting objects
%

%% highlight_object(+Identifier) is det.
%% highlight_object(Identifier, Highlight) is det.
%% highlight_object(Identifier, Highlight, Color) is det.
%% highlight_object(Identifier, Highlight, R, B, G, Alpha) is det.
%
% Different methods for highlighting objects. By default, objects are drawn in bright red
% if they are highlighted, but different colors can be specified using either one integer
% value (e.g. #00FF00) or separate values for red, green, and blue.
%
% The parameter Highlight specifies if the highlighting shall be activated or reset; if
% it is missing, a value of @(true) is assumed.
%
% If the object detection was uncertain, its probability can be visualized using the Prob
% parameter. This is done e.g. using the alpha channel or the hue value in HSV space
% (ignoring, in this case, the parameters R, B, G).
%
% @param Identifier eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
% @param Highlight  @(true) = highlight; @(false)=remove highlighting
% @param Color      Color value as integer, e.g. #AARRBBGG
% @param R          Red color value
% @param B          Blue color value
% @param G          Green color value
% @param Prob       Object existence probability
%
highlight_object(Identifier) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlight', [Identifier, @(true)], _).

highlight_object(Identifier, Highlight) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlight', [Identifier, Highlight], _).

highlight_object(Identifier, Highlight, Color) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlight', [Identifier, Highlight, Color], _).

highlight_object(Identifier, Highlight, R, B, G, Prob) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlight', [Identifier, Highlight, R, B, G, Prob], _).



%% highlight_object_with_children(+Identifier) is det.
%% highlight_object_with_children(+Identifier, +Highlight) is det.
%
% Highlights an object and everything that is reachable from it via knowrob:properPhysicalPartTypes
%
% The parameter Highlight specifies if the highlighting shall be activated or reset; if
% it is missing, a value of @(true) is assumed.
%
% @param Identifier eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
highlight_object_with_children(Identifier) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlightWithChildren', [Identifier, @(true)], _).

highlight_object_with_children(Identifier, Highlight) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'highlightWithChildren', [Identifier, Highlight], _).



%% reset_highlighting is det.
%
% Reset all highlighted objects in the canvas.
%
reset_highlight :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'clearHighlight', [], _).





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% Diagram canvas management
%



%% diagram_canvas is det.
%
% Launch the diagram data publisher
%
:- assert(d_canvas(fail)).
diagram_canvas :-
    d_canvas(fail),
    jpl_new('org.knowrob.vis.DiagramVisualization', [], Canvas),
    retract(d_canvas(fail)),
    assert(d_canvas(Canvas)),!.
diagram_canvas(Canvas) :-
    d_canvas(Canvas).


%% add_diagram(+Identifier, +Time) is nondet.
%
% Add object to the scene with its position at time 'Time'
%
% @param Identifier Object identifier, eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls"
%
add_diagram(Id, Title, Type, Xlabel, Ylabel, Width, Height, Fontsize, ValueList) :-
    d_canvas(Canvas),
    lists_to_arrays(ValueList, ValueArr),
    jpl_call(Canvas, 'addDiagram', [Id, Title, Type, Xlabel, Ylabel, Width, Height, Fontsize, ValueArr], _),!.

remove_diagram(Id) :-
    d_canvas(Canvas),
    jpl_call(Canvas, 'removeDiagram', [Id], _),!.

clear_diagram_canvas :-
    d_canvas(Canvas),
    jpl_call(Canvas, 'clear', [], _).

%% show_trajectory(+Starttime, +Endtime, +Interval)
% TODO description
show_trajectory(Starttime, Endtime, Interval) :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'showTrajectory', [Starttime, Endtime, Interval], _).

remove_trajectory :-
    v_canvas(Canvas),
    jpl_call(Canvas, 'removeTrajectory', [], _).

