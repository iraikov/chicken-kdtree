[[tags:egg]]

== kd-tree

K-D tree implementation.

[[toc:]]

== Usage

(require-extension kd-tree)

== Documentation

This library implements a [[http://en.wikipedia.org/wiki/K-d_tree|k-d
tree]], a data structure for organizing and searching points in
k-dimensional space. 

The k-d tree is a binary search tree in which every branching node
contains a k-dimensional point, and every leaf node contains a set of
points. Every branching node represents a splitting hyperplane that
divides the space into two parts, known as half-spaces.

Points to the left of the splitting hyperplane are contained in the
left subtree of the node and points right of the hyperplane are
contained in the right subtree. The splitting hyperplane is chosen so
as to be perpendicular to one of the axes in the k-dimensional
space. The axis at each branching level is chosen in a round-robin
fashion. For instance, in 3-D space, at level 0, the chosen axis is X,
so points are divided according to their X-coordinates; at level 1,
the chosen axis is Y, so the points are divided according to their
Y-coordinates; at the next branch level the chosen axis is Z, and so
on.


=== Point

This library currently only supports points in 3D space.

<procedure>make-point3d:: DOUBLE * DOUBLE * DOUBLE -> POINT3D</procedure>

3D point constructor.

<procedure>point3d? :: POINT3D -> BOOL</procedure>

3D point predicate.

<procedure>point3d-x :: POINT3D -> DOUBLE</procedure>
<procedure>point3d-y :: POINT3D -> DOUBLE</procedure>
<procedure>point3d-z :: POINT3D -> DOUBLE</procedure>

Accessors for the x,y,z coordinates of a 3D point.

=== KD-tree

A K-D tree (short for k-dimensional tree) is a space-partitioning data
structure for organizing points in a k-dimensional space.

==== Constructors
   
<procedure>list->kd-tree:: POINT3D LIST  -> KD-TREE</procedure>

Given a list of points, constructs and returns a K-D tree object.

==== Predicates

<procedure>kd-tree? :: KD-TREE -> BOOL </procedure>

Returns {{#t}} if the given object is a K-D tree, {{#f}} otherwise.

<procedure>kd-tree-empty? :: KD-TREE -> BOOL  </procedure>

Returns {{#t}} if the given K-D tree object is empty, {{#f}} otherwise.

<procedure>kd-tree-is-valid? :: KD-TREE -> BOOL  </procedure>

Checks whether the K-D tree property holds for the given tree.
Specifically, it tests that all points in the left subtree lie to the
left of the plane, p is on the plane, and all points in the right
subtree lie to the right.

<procedure>kd-tree-all-subtrees-are-valid? :: KD-TREE -> BOOL </procedure>

Checks whether the K-D tree property holds for the given tree and all
subtrees.

==== Accessors

<procedure>kd-tree->list :: KD-TREE -> POINT3D LIST</procedure>

Returns a list with the points contained in the tree.

<procedure>kd-tree->list* :: KD-TREE -> (INT . POINT3D) LIST </procedure>

Returns a list where every element has the form {{(i . p)}}, where i
is the relative index of this point, and {{p}} is the point.

<procedure>kd-tree-subtrees :: KD-TREE -> KD-TREE LIST</procedure>

<procedure>kd-tree-point :: KD-TREE -> POINT3D  </procedure>

==== Combinators

<procedure>kd-tree-map </procedure>

<procedure>kd-tree-for-each </procedure>

<procedure>kd-tree-for-each* </procedure>

<procedure>kd-tree-fold-right </procedure>

<procedure>kd-tree-fold-right* </procedure>

==== Query procedures

<procedure>kd-tree-nearest-neighbor </procedure>

<procedure>kd-tree-near-neighbors </procedure>

<procedure>kd-tree-near-neighbors* </procedure>

<procedure>kd-tree-k-nearest-neighbors </procedure>

<procedure>kd-tree-slice </procedure>

<procedure>kd-tree-slice* </procedure>

==== Modifiers

<procedure>kd-tree-remove </procedure>


== Examples


== About this egg


=== Author

[[/users/ivan-raikov|Ivan Raikov]]

=== Version history

; 4.1-4.8 : Using f64vector for internal point representation
; 4.0-4.1 : Added with-distance? flag to kd-tree-near-neighbors
; 3.2 : Bug fix in kd-tree-near-neighbors
; 2.0 : Improvements to internal representation
; 1.0 : Initial release

=== License


 Copyright 2012-2015 Ivan Raikov
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 A full copy of the GPL license can be found at
 <http://www.gnu.org/licenses/>.

