/** <module> knowrob_mesh_reasoning

  Description:
    Module providing mesh reasoning capabilities


  Copyright (C) 2012 by Stefan Profanter

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

@author Stefan Profanter
@license GPL
*/

:- module(knowrob_mesh_reasoning,
    [
      mesh_reasoning/2,
      mesh_reasoning_path/2,
      mesh_element_types/2,
      mesh_find_annotations/3,
      mesh_find_supporting_planes/2,
      mesh_is_supporting_plane/1,
      mesh_is_supporting_plane/2,
      mesh_reasoning_highlight/2,
      mesh_reasoning_clear_highlight/1,
      mesh_find_handle/2,
      mesh_find_handle/4,
      listsplit/3,
      jpl_set_to_list/2,
      comp_physical_parts/2,
      annotation_area/2,
      annotation_area_coverage/2,
      annotation_pose_list/2,
      annotation_cone_radius_avg/2,
      annotation_cone_radius_max/2,
      annotation_cone_radius_min/2,
      annotation_cone_volume/2,
      annotation_cone_height/2,
      annotation_cone_direction/2,
      annotation_sphere_radius/2,
      annotation_sphere_is_concave/2,
      annotation_sphere_volume/2,
      annotation_plane_normal/2,
      annotation_plane_longside/2,
      annotation_plane_shortside/2,
      annotation_container_direction/2,
      annotation_container_volume/2
    ]).

:- use_module(library('semweb/rdfs')).
:- use_module(library('semweb/rdf_db')).
:- use_module(library('semweb/rdfs_computable')).
:- use_module(library('knowrob_objects')).
:- use_module(library('knowrob_coordinates')).
:- use_module(library('jpl')).


:- rdf_meta comp_physical_parts(r,r),
            mesh_annotator_handle(r,?),
            mesh_annotation_handle(r,?),
            mesh_reasoning(r,-),
            comp_physical_parts(r,r),
            annotation_area(r,?),
            annotation_area_coverage(r,?),
            annotation_pose_list(r,?),
            annotation_cone_radius_avg(r,?),
            annotation_cone_radius_max(r,?),
            annotation_cone_radius_min(r,?),
            annotation_cone_volume(r,?),
            annotation_cone_height(r,?),
            annotation_cone_direction(r,?),
            annotation_sphere_radius(r,?),
            annotation_sphere_is_concave(r,?),
            annotation_sphere_volume(r,?),
            annotation_plane_normal(r,?),
            annotation_plane_longside(r,?),
            annotation_plane_shortside(r,?),
            annotation_container_direction(r,?),
            annotation_container_volume(r,?).



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Main
%

%% mesh_reasoning(+Identifier, -MeshReasoning) is det.
%
% Do mesh reasoning on cad model with given identifier.
%
% @param Identifier 	   eg. "http://ias.cs.tum.edu/kb/ias_semantic_map.owl#F360-Containers-revised-walls" or "knowrob:'Spoon'"
% @param MeshReasoning     MeshReasoning object
%
mesh_reasoning(Identifier, MeshReasoning) :-
  get_model_path(Identifier,Path),
  mesh_reasoning_path(Path, MeshReasoning).

%% mesh_reasoning_path(+Path, -MeshReasoning) is det.
%
% Do mesh reasoning on cad model with given path (supported: local, package://, http://, ftp://).
%
% @param Path eg.   "/home/user/model.dae" or "http://example.com/model.kmz"
% @param MeshReasoning     MeshReasoning object
%
mesh_reasoning_path(Path, MeshReasoning) :-
  mesh_reasoning_init(MeshReasoning),
  jpl_call(MeshReasoning, 'analyseByPath', [Path], _).

    
%% mesh_reasoning_init(-MeshReasoning,+WithCanvas) is det.
%% mesh_reasoning_init(-MeshReasoning) is det.
%
% Create mesh reasoning object. WithCanvas indicates if you want to show canvas window.
% WithCanvas defaults to true if not indicated
%
mesh_reasoning_init(MeshReasoning, WithCanvas) :-
  jpl_call('edu.tum.cs.vis.model.MeshReasoning', 'initMeshReasoning', [WithCanvas], MeshReasoning).
mesh_reasoning_init(MeshReasoning) :-
  mesh_reasoning_init(MeshReasoning, @(true)).



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Annotations
%

%% mesh_reasoning_highlight(+MeshReasoning,+[AnnotationHead|AnnotationTail]) is det
%% mesh_reasoning_highlight(+MeshReasoning,+[]) is det
%% mesh_reasoning_highlight(+MeshReasoning,+Annotation) is det
%
% Highlight/select the specified annotation in GUI
%
% @param MeshReasoning		reasoning container
% @param Annotation			annotation to highlight
mesh_reasoning_highlight(MeshReasoning,[AnnotationHead|AnnotationTail]) :-
	mesh_reasoning_highlight(MeshReasoning, AnnotationHead),!,
	mesh_reasoning_highlight(MeshReasoning, AnnotationTail),!.
mesh_reasoning_highlight(_,[]).
mesh_reasoning_highlight(MeshReasoning, Annotation) :-
	jpl_call(MeshReasoning, 'highlightAnnotation', [Annotation], _).
	
%% mesh_reasoning_clear_hightlight(+MeshReasoning) is det
%
% Clear all annotation highlights in GUI
%
% @param MeshReasoning		reasoning container
%
mesh_reasoning_clear_highlight(MeshReasoning) :-
	jpl_call(MeshReasoning, 'clearHightlight', [], _).
	

%% mesh_element_types(+MeshReasoning,-TypeList) is det
%
% Get list of all found annotation types for current model in MeshReasoning
%
% @param MeshReasoning		reasoning container
% @param TypeList			List with annotation types eg: ['Plane','Sphere','Cone','Container']. Values can be directly used in mesh_find_annotations as Type
%
mesh_element_types(MeshReasoning, TypeList) :-
	jpl_call(MeshReasoning, 'getAnnotationTypes', [], TypeListFlat),
    jpl_call(TypeListFlat, toArray, [], TypeListArr),
	jpl_array_to_list(TypeListArr, TypeList).


%% mesh_find_annotations(+MeshReasoning,+Type, -AnnotationsList) is det.
%
% Get a list of all annotations with given type
%
% @param MeshReasoning		reasoning container
% @param Type		String indicating annotation type (Plane,Sphere,Cone,Container)
% @param AnnotationsList List with found annotations
%
mesh_find_annotations(MeshReasoning,Type,AnnotationsList) :-
	concat('findAnnotations', Type, Method),
	jpl_call(MeshReasoning, Method, [], AnnotationsSet),
  jpl_set_to_list(AnnotationsSet, AnnotationsList).
	
	
%% mesh_find_supporting_planes(+MeshReasoning, -PlaneList) is det.
%
% Get list of all supporting planes
%
% @param MeshReasoning		reasoning container
% @param PlaneList			returning list which contains all supporting planes
%
mesh_find_supporting_planes(MeshReasoning, PlaneList) :-
	mesh_find_annotations(MeshReasoning,'Plane',AnnList),
	findall(P,(member(P, AnnList),mesh_is_supporting_plane(P)),PlaneList).


%% mesh_is_supporting_plane(+PlaneAnnotation) is det.
%% mesh_is_supporting_plane(+PlaneAnnotation, +Identifier) is det.
%
% Check if plane annotation surface normal is upwards.
% Means checking if abs(acos(n.z)) < 10 degree = PI/18 rad 
%
% @param PlaneAnnotation	A Java Reference Object for a plane annotation
% @param Identifier			If Identifier is instantiated the current pose of the object is also considered when checking if normal is upwards
%
mesh_is_supporting_plane(PlaneAnnotation, Identifier) :-
	jpl_call(PlaneAnnotation,'getPlaneNormal',[],Norm),
	jpl_get(Norm,'z',NormZ),
	(nonvar(Identifier) -> 
		% we only need M20,M21,M22 because '(rot) * (NormX,NormY,NormZ,0) = (_, _, NewNormZ, _)' for calculating angle
		jpl_get(Norm,'x',NormX),
		jpl_get(Norm,'y',NormY),
		current_object_pose(Identifier,[_, _, _, _, _, _, _, _, M20, M21, M22, _, _, _, _, _]),
		NormZ is M20 * NormX+M21 * NormY+M22 * NormZ,
		true
	; 
		true
	),
	abs(acos(NormZ)) < pi / 18.
mesh_is_supporting_plane(PlaneAnnotation) :-
	mesh_is_supporting_plane(PlaneAnnotation,_).

:- dynamic
        mesh_min_radius/1, 
        mesh_max_radius/1.

mesh_handle_comparator(Comp, W1, W2) :-
	jpl_call(W1,'getAreaCoverage',[],Cov1),
	jpl_call(W2,'getAreaCoverage',[],Cov2),
	jpl_call(W1,'getRadiusAvgUnscaled',[],Rad1),
	jpl_call(W2,'getRadiusAvgUnscaled',[],Rad2),
	( mesh_min_radius(MinRadius),mesh_max_radius(MaxRadius) -> 
		(
			(Rad1 > MinRadius , Rad1 < MaxRadius) -> 
			Rad1Ok = true
			; Rad1Ok = false
		)
		; Rad1Ok = false
	),
	( mesh_min_radius(MinRadius),mesh_max_radius(MaxRadius) -> 
		(
			(Rad2 > MinRadius , Rad2 < MaxRadius) -> 
			Rad2Ok = true
			; Rad2Ok = false
		)
		; Rad2Ok = false
	),
	(	Cov1 < 0.6 , Cov2 >= 0.6 -> Comp = '>'
	;	Cov1 >= 0.6, Cov2 < 0.6 -> Comp = '<'
	;	not(Rad1Ok), Rad2Ok -> Comp = '>'
	;	Rad1Ok, not(Rad2Ok) -> Comp = '<'
	;	jpl_call(W1,'getHeightUnscaled',[],H1),
		jpl_call(W2,'getHeightUnscaled',[],H2),
		(   H1 < H2 -> Comp = '>'
		;   H1 >= H2 -> Comp = '<')
	).

listsplit([H|T], H, T).

%% jpl_set_to_list(+Set,-List) is det.
%
%  Extracts all elements from a JPL set to a prolog list by calling jpl_set_element
%
% @param Set	The jpl set
% @param List	returning list
%
jpl_set_to_list(Set,List) :-
  findall(P,jpl_set_element(Set,P),List).

%% mesh_find_handle(+MeshReasoning, -HandleAnnotations) is det.
%% mesh_find_handle(+MeshReasoning, -HandleAnnotations, +MinRadius, +MaxRadius) is det.
%
% Returns a list which contains annotations sorted by its probability that they are the object handle.
% Sorting is archeived by calling mesh_handle_comparator which compares two annotations by its probability.
%
% @param MeshReasoning			reasoning container
% @param HandleAnnotations		the resulting sorted annotation list
% @param MinRadius				minimum radius which the handle should have
% @param MaxRadius				maximum radius which the handle should have
%
mesh_find_handle(MeshReasoning, HandleAnnotations) :-
  mesh_find_annotations(MeshReasoning,'Cone',AnnList),
  predsort(mesh_handle_comparator, AnnList, HandleAnnotations).

mesh_find_handle(MeshReasoning, HandleAnnotations, MinRadius, MaxRadius) :-
  assert(mesh_min_radius(MinRadius)),
  assert(mesh_max_radius(MaxRadius)),
  mesh_find_handle(MeshReasoning, HandleAnnotations),
  retractall(mesh_min_radius(_)),
  retractall(mesh_max_radius(_)).




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% Computable definitions
% 

%
% check if mesh annotator has already been initialized for this instance
% if not, initialize and assert for this object
% 

:- dynamic 
   mesh_annotator_handle/2,
   mesh_annotation_handle/2.

mesh_annotator(Obj, MeshAnnotator) :-
  ((mesh_annotator_handle(Obj, MeshAnnotator),!) ;
   (get_model_path(Obj, Path),
    mesh_reasoning_init(MeshAnnotator, @(false)),
    jpl_call(MeshAnnotator, 'analyseByPath', [Path], _),
    assert(mesh_annotator_handle(Obj, MeshAnnotator)))).



comp_physical_parts(Obj, PartInst) :-

  % avoid re-creation of object parts
  (\+ mesh_annotator_handle(Obj,_)),

  % locates mesh for the object
  % determines object components
  % assert java reference to mesh reasoning object
  mesh_annotator(Obj, MeshAnnotator),
  mesh_element_types(MeshAnnotator, ContainedTypes),
  member(Type, ContainedTypes),

  mesh_find_annotations(MeshAnnotator,Type,AnnotationsList),
  member(Annotation, AnnotationsList),

  % TODO: check if component already exists; otherwise:

  annotation_pose_list(Annotation, PoseList),
  annotation_to_knowrob_class(Type, KnowRobClass),
  create_object_perception(KnowRobClass, PoseList, ['MeshSegmentation'], PartInst),

  % assert sub-parts
  rdf_assert(Obj, knowrob:properPhysicalParts, PartInst),

  % assert Java annotation ID for each plane/sphere/... in order to retrieve their properties
  assert(mesh_annotation_handle(PartInst, Annotation)).


annotation_to_knowrob_class('Cone', 'Cone').
annotation_to_knowrob_class('Sphere', 'Sphere').
annotation_to_knowrob_class('Plane', 'FlatPhysicalSurface').
annotation_to_knowrob_class('Container', 'ContainerArtifact').


annotation_pose_list(PrimitiveAnnotation, PoseList) :-
  jpl_datum_to_type(PrimitiveAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],[_])),
  jpl_call(PrimitiveAnnotation,'getPoseMatrix',[],PoseMatrix),
  knowrob_coordinates:matrix4d_to_list(PoseMatrix, PoseList).





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Computables for reading information about annotation properties
% 

annotation_area(PartInst, Area) :-
  mesh_annotation_handle(PartInst, PrimitiveAnnotation),
  jpl_datum_to_type(PrimitiveAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],[_])),
  jpl_call(PrimitiveAnnotation,'getPrimitiveArea',[],Area).

annotation_area_coverage(PartInst, AreaCoverage) :-
  mesh_annotation_handle(PartInst, PrimitiveAnnotation),
  jpl_datum_to_type(PrimitiveAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],[_])),
  jpl_call(PrimitiveAnnotation,'getAreaCoverage',[],AreaCoverage).




% % % % % % % % % % % % % % % % % % % % % % % 
% CONES

annotation_cone_radius_avg(PartInst, RadiusAvg) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getRadiusAvg',[],RadiusAvg).

annotation_cone_radius_max(PartInst, RadiusMax) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getRadiusLarge',[],RadiusMax).

annotation_cone_radius_min(PartInst, RadiusMin) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getRadiusSmall',[],RadiusMin).

annotation_cone_volume(PartInst, Volume) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getVolume',[],Volume).

annotation_cone_height(PartInst, Height) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getHeight',[],Height).

annotation_cone_direction(PartInst, Direction) :-
  mesh_annotation_handle(PartInst, ConeAnnotation),
  jpl_datum_to_type(ConeAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['ConeAnnotation'])),
  jpl_call(ConeAnnotation,'getDirection',[],DirVec),
  knowrob_coordinates:vector3d_to_list(DirVec, Direction).




% % % % % % % % % % % % % % % % % % % % % % % 
% SPHERES

annotation_sphere_radius(PartInst, RadiusAvg) :-
  mesh_annotation_handle(PartInst, SphereAnnotation),
  jpl_datum_to_type(SphereAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['SphereAnnotation'])),
  jpl_call(SphereAnnotation,'getRadius',[],RadiusAvg).

annotation_sphere_is_concave(PartInst, Concave) :-
  mesh_annotation_handle(PartInst, SphereAnnotation),
  jpl_datum_to_type(SphereAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['SphereAnnotation'])),
  jpl_call(SphereAnnotation,'isConcav',[],Concave).

annotation_sphere_volume(PartInst, Volume) :-
  mesh_annotation_handle(PartInst, SphereAnnotation),
  jpl_datum_to_type(SphereAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['SphereAnnotation'])),
  jpl_call(SphereAnnotation,'getVolume',[],Volume).




% % % % % % % % % % % % % % % % % % % % % % % 
% Planes

annotation_plane_normal(PartInst, PlaneNormal) :-
  mesh_annotation_handle(PartInst, PlaneAnnotation),
  jpl_datum_to_type(PlaneAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['PlaneAnnotation'])),
  jpl_call(PlaneAnnotation,'getPlaneNormal',[],NormalVec),
  knowrob_coordinates:vector3d_to_list(NormalVec, PlaneNormal).

annotation_plane_longside(PartInst, LongSide) :-
  mesh_annotation_handle(PartInst, PlaneAnnotation),
  jpl_datum_to_type(PlaneAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['PlaneAnnotation'])),
  jpl_call(PlaneAnnotation,'getLongSide',[],LongSideVec),
  knowrob_coordinates:vector3d_to_list(LongSideVec, LongSide).

annotation_plane_shortside(PartInst, ShortSide) :-
  mesh_annotation_handle(PartInst, PlaneAnnotation),
  jpl_datum_to_type(PlaneAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation,primitive],['PlaneAnnotation'])),
  jpl_call(PlaneAnnotation,'getShortSide',[],ShortSideVec),
  knowrob_coordinates:vector3d_to_list(ShortSideVec, ShortSide).




% % % % % % % % % % % % % % % % % % % % % % % 
% Containers

annotation_container_direction(PartInst, OpeningDirection) :-
  mesh_annotation_handle(PartInst, ContainerAnnotation),
  jpl_datum_to_type(ContainerAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation],['ContainerAnnotation'])),
  jpl_call(ContainerAnnotation,'getDirection',[],DirVec),
  knowrob_coordinates:vector3d_to_list(DirVec, OpeningDirection).

annotation_container_volume(PartInst, Volume) :-
  mesh_annotation_handle(PartInst, ContainerAnnotation),
  jpl_datum_to_type(ContainerAnnotation, 
      class([edu,tum,cs,vis,model,uima,annotation],['ContainerAnnotation'])),
  jpl_call(ContainerAnnotation,'getVolume',[],Volume).

