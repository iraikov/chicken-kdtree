# kd-tree

K-D tree implementation in Chicken Scheme.

## Documentation

This library implements a K-D tree
(http://en.wikipedia.org/wiki/K-d_tree), which is a data structure for
organizing and searching points in k-dimensional space.

The K-D tree is a binary search tree in which every branching node
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


### K-dimensional point space

Module `kspace` provides facilities for managament of K-dimensional point spaces.

<procedure>make-space:: COORDS -> SPACE</procedure>

Given a list of coordinate collections of length K, constructs a yasos  K-dimensional point space object. The coordinate collections can be SRFI-4 f32vectors, or collection objects as defined in the yasos collections module.

<procedure>space? :: OBJECT -> BOOL</procedure>

K-dimensional point space predicate.

<procedure>dimension :: OBJECT -> INT</procedure>
Returns the dimensionality of the point space.

<procedure>point :: OBJECT * INT -> REAL LIST</procedure>
Returns the coordinates of the point at the given index.

<procedure>coord :: OBJECT * INT * INT -> FLOAT</procedure>
Returns the k'th coordinate of i'th point, starting from 0.


<procedure>compare-coord :: SPACE * INT * INT * INT -> INT</procedure>

Given the indices of two points and a coordinate index, compares the
respective coordinates of the two points and returns -1, 0, or 1,
depending on whether the coordinates are less than, equal, or greater
than each other.

<procedure>squared-distance :: SPACE * INT * INT -> FLOAT</procedure>

Returns the square of the Euclidean distance between the points at the given indices.

<procedure>compare-distance :: SPACE * INT * INT -> INT</procedure>

Compares the square of the Euclidean distance between the points at the given indices.

### K-D tree interface

#### Constructors
   
<procedure>make-kd-tree:: SPACE -> OBJECT</procedure>

Given a `kspace` object, constructs and returns a yasos spatial map object.

#### Predicates

<procedure>spatial-map? :: OBJECT -> BOOL </procedure>

Returns `#t` if the given object is a spatial map, `#f` otherwise.

<procedure>empty? :: OBJECT -> BOOL  </procedure>

Returns `#t` if the given spatial map object is empty, `#f` otherwise.

<procedure>kd-tree-is-valid? :: OBJECT -> BOOL  </procedure>

Checks whether the K-D tree property holds for the given spatial map.
Specifically, it tests that all points in the left subtree lie to the
left of the plane, p is on the plane, and all points in the right
subtree lie to the right.

<procedure>kd-tree-all-subtrees-are-valid? :: KD-TREE -> BOOL </procedure>

Checks whether the K-D tree property holds for the given spatial map
and all subtrees.

#### Accessors

<procedure>get-kspace :: OBJECT -> KSPACE</procedure>

Returns the underlying `kspace` object of the map.

<procedure>spatial-map->list :: KD-TREE -> POINT LIST</procedure>

Returns a list with the points contained in the spatial map.

#### Query procedures

<procedure>nearest-neighbor :: SPATIAL-MAP * POINT -> POINT</procedure>

Nearest neighbor of a point.

<procedure>near-neighbors :: SPATIAL-MAP * POINT * RADIUS -> POINT LIST </procedure>

Neighbors of a point within radius r.

<procedure>k-nearest-neighbors :: SPATIAL-MAP * POINT * INT -> POINT LIST </procedure>

K nearest neighbors of a point.

<procedure>slice :: SPATIAL-MAP * AXIS * COORD * COORD -> POINT LIST </procedure>

Returns all points between two planes.

#### Combinators

<procedure>spatial-map-for-each :: SPATIAL-MAP * F -> VOID </procedure>

Point iterator.

<procedure>spatial-map-fold :: SPATIAL-MAP * F -> AX </procedure>

<procedure>spatial-map-fold-right :: SPATIAL-MAP * F * INIT -> INIT </procedure>

Fold on points.

<procedure>spatial-map-fold-right* :: SPATIAL-MAP * F * INIT -> INIT </procedure>

Fold on points and their indices.

#### Modifiers

<procedure>kd-tree-remove :: SPATIAL-MAP * POINT -> SPATIAL-MAP </procedure>


## Examples


```scheme

(import scheme (chicken base) 
        yasos random-mtzig kspace kd-tree)

(let* (
       (n (inexact->exact 1e5)) (k 40) (r 1.0) (randst (init 9))
       
       ;; generate random coordinates
       (xs (randn/f32! n randst))
       (ys (randn/f32! n randst))
       (zs (randn/f32! n randst))
       ;; create a kspace
       (pts (list xs ys zs))
       (kspace3d (make-space pts))
       ;; create the spatial map
       (kdt (make-kd-tree kspace3d))
       ;; choose a random point index
       (xi (inexact->exact (modulo (random! randst) n)))
       ;; retrieve the coordinates of the chosen point
       (x  (point kspace3d xi))
       )

    (print "nearest-neighbor of " x ": " (nearest-neighbor kdt x))
    (print k " nearest neighbors of " x ": " (k-nearest-neighbors kdt x k))
    (print "near neighbors of " x " within " r ": " (near-neighbors kdt k r))

)


```

## Version history

- 6.0 : refactored to use yasos, ported to CHICKEN 5
- 5.0 : added list->kd-tree* procedure to KdTree class
- 4.1-4.8 : Using f64vector for internal point representation
- 4.0-4.1 : Added with-distance? flag to kd-tree-near-neighbors
- 3.2 : Bug fix in kd-tree-near-neighbors
- 2.0 : Improvements to internal representation
- 1.0 : Initial release

## License

>
> Copyright 2012-2019 Ivan Raikov
> 
>  This program is free software: you can redistribute it and/or modify
>  it under the terms of the GNU General Public License as published by
>  the Free Software Foundation, either version 3 of the License, or (at
>  your option) any later version.
>  
>  This program is distributed in the hope that it will be useful, but
>  WITHOUT ANY WARRANTY; without even the implied warranty of
>  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
>  General Public License for more details.
> 
>  A full copy of the GPL license can be found at
>  <http://www.gnu.org/licenses/>.

