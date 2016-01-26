;;
;; An implementation of the K-d tree spatial indexing data structure.
;; 
;; http://en.wikipedia.org/wiki/K-d_tree
;;
;; The k-d tree is a binary search tree in which every branching node
;; contains a k-dimensional point, and every leaf node contains a set
;; of points. Every branching node represents a splitting hyperplane
;; that divides the space into two parts, known as half-spaces.
;;
;; Points to the left of the splitting hyperplane are contained in the
;; left subtree of the node and points right of the hyperplane are
;; contained in the right subtree. The splitting hyperplane is chosen
;; so as to be perpendicular to one of the axes in the k-dimensional
;; space. The axis at each branching level is chosen in a round-robin
;; fashion. For instance, in 3-D space, at level 0, the chosen axis is
;; X, so points are divided according to their X-coordinates; at level
;; 1, the chosen axis is Y, so the points are divided according to
;; their Y-coordinates; at the next branch level the chosen axis is Z,
;; and so on.
;;
;;
;; This code is based on the Haskell kd-tree library implementation of
;; K-D trees.
;;
;; Copyright 2012-2016 Ivan Raikov
;;
;; This program is free software: you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation, either version 3 of the
;; License, or (at your option) any later version.
;;
;; This program is distributed in the hope that it will be useful, but
;; WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;; General Public License for more details.
;;
;; A full copy of the GPL license can be found at
;; <http://www.gnu.org/licenses/>.
;;

(module kd-tree
	
  (
   <Point> default-<Point> Point3d Point2d
   point? make-point

   <KdTree> default-<KdTree> KdTree3d KdTree2d

   kd-tree? 
   kd-tree-empty?
   kd-tree->list
   kd-tree->list*
   kd-tree-map
   kd-tree-for-each
   kd-tree-for-each*
   kd-tree-fold-right
   kd-tree-fold-right*
   kd-tree-subtrees
   kd-tree-node-points
   kd-tree-node-indices
   kd-tree-min-index
   kd-tree-max-index
   kd-tree-size
  )

  (import scheme chicken data-structures foreign)
  
  (require-library srfi-1 srfi-4 extras cis)
  (require-extension typeclass datatype)

  (import (only srfi-1 fold list-tabulate split-at span every fold-right take filter filter-map remove zip)
	   (only srfi-4 f64vector? f64vector f64vector-ref)
	   (only extras fprintf pp)
	   (prefix cis cis:))

  (define log2 (foreign-lambda double "log2" double))

  (define cdist2
    (foreign-lambda* double ((int dimension) (f64vector a) (f64vector b))
#<<EOF
     int i; double sum, delta;

     sum = 0.0;

     for (i=0; i<dimension; i++)
     {
        delta = (a[i] - b[i]);
        sum = sum + delta*delta;
     }

     C_return(sum);
EOF
))

  (define-class <Point> 
    ;; dimension has the number of coordinates of a point.
    dimension ;; Int

    ;; gets the k'th coordinate, starting from 0.
    coord ;; Int * Point -> Double

    ;; compares the given coordinates
    compare-coord ;; Int * Point * Point -> Bool

    ;; returns the squared distance between two points.
    dist2 ;; Point * Point -> Double

    ;; returns 0, negative or positive number depending on the
    ;; distance between two points
    compare-distance

    )

  (define (minimum-by lst less? . rest)
    (if (null? lst) #f
	(if (null? rest)
	    (let recur ((lst (cdr lst)) (m (car lst)))
	      (if (null? lst) m
		  (if (less? (car lst) m)
		      (recur (cdr lst) (car lst))
		      (recur (cdr lst) m)
		      ))
	      )
	    (let recur ((lst (cdr lst)) 
			(rest (map cdr rest))
			(m (map car (cons lst rest))))
	      (if (null? lst) m
		  (if (less? (car lst) (car m))
		      (recur (cdr lst) (map cdr rest) (map car (cons lst rest)))
		      (recur (cdr lst) (map cdr rest) m)
		      ))
	      )
	    )))


  (define (kd-elt-index x) (let ((xp (car x))) (if (number? xp) xp (car xp))))
  (define (kd-elt-value x) (and (pair? (car x)) (cadar x)))
  (define (kd-elt-point x) (cadr x))


  (define (default-<Point> dimension coord)

    (let* (
	   (dist2 (lambda (a b) (cdist2 dimension a b)))

	   (compare-distance
	    (lambda (p a b . reltol)
	      (let ((delta (- (dist2 p a) (dist2 p b))))
		(if (null? reltol) 
		    delta 
		    (if (<= delta (car reltol)) 0 delta)))))

	   (compare-coord 
	    (lambda (c a b)
	      (< (coord c a) (coord c b))))
	      
	   )
	    
    (make-<Point> 
     dimension coord compare-coord dist2 compare-distance)
    ))

  (define point? f64vector?)
  (define make-point f64vector)


  (define Point3d
    (default-<Point> 3 (lambda (i p) (f64vector-ref p i)) ))

  (define Point2d
    (default-<Point> 2  (lambda (i p) (f64vector-ref p i))))


  (define-class <KdTree> 

    ;; constructs a kd-tree from a list of points
    list->kd-tree
    ;; constructs a kd-tree from a list of points with indices and values
    list->kd-tree*
    ;; nearest neighbor of a point
    kd-tree-nearest-neighbor
    ;; the index of the nearest neighbor of a point
    kd-tree-nearest-neighbor*
    ;; neighbors of a point within radius r
    kd-tree-near-neighbors
    ;; neighbors of a point within radius r (using point indices)
    kd-tree-near-neighbors*
    ;; k nearest neighbors of a point
    kd-tree-k-nearest-neighbors
    ;; removes a point from the tree
    kd-tree-remove
    ;; retrieves all points between two planes
    kd-tree-slice
    ;; retrieves all points between two planes (using point indices)
    kd-tree-slice*
    ;; checks that the kd-tree properties are preserved
    kd-tree-is-valid?
    kd-tree-all-subtrees-are-valid?

    )


  (define (positive-or-zero-integer? x)
    (and (integer? x) (or (zero? x) (positive? x))))


  (define (positive-integer? x)
    (and (integer? x) (positive? x)))


  (define-datatype kd-tree kd-tree?
    (KdNode (left  kd-tree?)
	    (p     point?)
	    (i     positive-or-zero-integer?)
            (v     (lambda (v) (or (not v) v)))
	    (right kd-tree?)
	    (axis  positive-or-zero-integer?)
	    (ci    cis:cis?)
	    )
    (KdLeaf (ii cis:cis?) 
	    (pp (lambda (lst) (every point? lst))) 
	    (vv (lambda (v) (or (list? v) (not v))))
	    (axis positive-or-zero-integer?) )
    )


  (define (kd-tree-empty? t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) (cis:empty? ii))
	   (else #f)))

  
  (define (kd-tree-map f t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) 
		   (KdLeaf ii (map f pp) vv axis))
	   (KdNode (l x i v r axis ci)
		   (KdNode (kd-tree-map f l)
			   (f x)
			   i v
			   (kd-tree-map f r)
			   axis ci))
	   ))
  
  (define (kd-tree-for-each f t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) (for-each f pp))
	   (KdNode (l x i v r axis ci)
		   (begin
		     (kd-tree-for-each f l)
		     (f x)
		     (kd-tree-for-each f r)
		     ))
	   ))


  (define (kd-tree-for-each* f t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis)
                   (if vv (for-each f (zip (reverse (cis:elements ii)) vv) pp)
                       (for-each f (reverse (cis:elements ii)) pp)))
	   (KdNode (l x i v r axis ci)
		   (begin
		     (kd-tree-for-each* f l)
		     (if v (f (list i v) x) (f i x))
		     (kd-tree-for-each* f r)
		     ))
	   ))

  
  (define (kd-tree-fold-right f init t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) 
		   (fold-right f init pp))
	   (KdNode (l p i v r axis ci)
		   (let* ((init2 (kd-tree-fold-right f init r))
			  (init3 (f p init2)))
		     (kd-tree-fold-right f init3 l)))
	   ))


  (define (kd-tree-fold-right* f init t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) 
		   (if vv
		       (fold-right f init (zip (reverse (cis:elements ii)) vv) pp)
		       (fold-right f init (reverse (cis:elements ii)) pp)))
	   (KdNode (l x i v r axis ci)
		   (let* ((init2 (kd-tree-fold-right* f init r))
			  (init3 (if v (f (list i v) x init2) (f i x init2))))
		     (kd-tree-fold-right* f init3 l)))
	   ))
  
  

  
  (define (kd-tree->list t)
    (kd-tree-fold-right cons '() t))

  
  (define (kd-tree->list* t #!key (fn list))
    (kd-tree-fold-right* 
     (lambda (iv p ax) (cons (fn iv p) ax))
     '() t))
  
  
  ;; Returns a list containing t and all its subtrees, including the
  ;; leaf nodes.
  
  (define (kd-tree-subtrees t)
    (cases kd-tree t
		  (KdLeaf (ii pp vv axis)  (list t))
		  (KdNode (l x i v r axis ci)
			  (append (kd-tree-subtrees l) 
				  (list t) 
				  (kd-tree-subtrees r)))
		  ))

  
  (define (kd-tree-node-points t)
    (cases kd-tree t
		  (KdLeaf (ii pp vv axis)  pp)
		  (KdNode (l x i v r axis ci) (list x))
		  ))
  
  (define (kd-tree-node-indices t)
    (cases kd-tree t
		  (KdLeaf (ii pp vv axis) ii)
		  (KdNode (l x i v r axis ci) ci)
		  ))

  (define (kd-tree-node-values t)
    (cases kd-tree t
		  (KdLeaf (ii pp vv axis) vv)
		  (KdNode (l x i v r axis ci) (list v))
		  ))

  (define (kd-tree-size t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) (cis:cardinal ii))
	   (KdNode (l x i v r axis ci) (cis:cardinal ci))))

  (define (kd-tree-min-index t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) (cis:get-min ii))
	   (KdNode (l x i v r axis ci) (cis:get-min ci))))

  (define (kd-tree-max-index t)
    (cases kd-tree t
	   (KdLeaf (ii pp vv axis) (cis:get-max ii))
	   (KdNode (l x i v r axis ci) (cis:get-max ci))))



  ;; construct a kd-tree from a list of points
  (define=> (make-list->kd-tree/depth <Point>)

    (lambda (make-value)

      (letrec (
	       (split
		(lambda (m n points depth)

		  (let* ((axis   (modulo depth dimension))
			 (cmpfn  (lambda (p0 p1) (compare-coord axis (car p0) (car p1))))
			 (sorted (sort points cmpfn))
			 (median-index (quotient (- n m) 2)))

		    (let-values (((lte gte) (split-at sorted median-index)))

		      (let ((median (car (car gte))))

			(let-values (((lt xeq) (span (lambda (x) (< (coord axis (car x)) (coord axis median))) lte)))
			  
			  (if (null? xeq)
			      (values (car gte) median-index lt (cdr gte))
			      (let ((split-index (length lt)))
				(values (car xeq) split-index lt (append (cdr xeq) gte))))
			))
		      ))
		  ))

		(list->kd-tree/depth
		 (lambda (m n points depth #!key (bucket-size (* 10 (max (log2 (- n m)) 1))) (offset 0))

		   (let ((k (- n m)))

		     (cond
		      ((null? points) (KdLeaf cis:empty '() '() depth))

		      ((<= k bucket-size)

		       (let* ((es (take points k))
			      (ps (map car es)) 
			      (ii (cis:shift offset (cis:interval m (- n 1))))
			      (vs (and make-value
				       (map (lambda (i e) (make-value i (cdr e)))
					    (reverse (cis:elements ii)) 
					    es)))
			      )
                        
			(KdLeaf ii ps vs (modulo depth dimension))
			))
		    
		      ((null? (cdr points))
		       (let* ((e (car points))
			      (ps (list (car e)) )
			      (vs (and make-value (list (make-value m (cdr e))))))
			 
			 (KdLeaf (cis:shift offset (cis:singleton m) )
				 ps vs
				 (modulo depth dimension))
			 ))
		      
		      (else
		       (let-values (((median median-index lt gte)
				     (split m n points depth)))
			 
			 
			 (let* ((depth1 (+ 1 depth))
				(i (+ m median-index offset))
				(p (car median))
				(v (and make-value (make-value i (cdr median))))
				(axis (modulo depth dimension)))
			   
			   (KdNode (list->kd-tree/depth m (+ m median-index) lt depth1 
							bucket-size: bucket-size)
				   p i v
				   (list->kd-tree/depth (+ m median-index 1) n gte depth1 
							bucket-size: bucket-size)
				   axis 
				   (cis:shift offset (cis:interval m (- n 1))))))
		       ))
		     ))
		 ))
	list->kd-tree/depth
	))
      )

  
  ;; A variant of list->kd-tree specialized for
  ;; kd-tree-remove. Specifically, it preserves the
  ;; original point indices.  This variant of
  ;; list->kd-tree assumes that the input list of points
  ;; has the form ( (INDEX POINT VALUE) ... )

  (define=> (make-list->kd-tree/depth* <Point>)

      (letrec (
               
	       (split
		(lambda (points depth)
                  
		  (let* ((axis   (modulo depth dimension))
			 (cmpfn  (lambda (p0 p1) (compare-coord axis (kd-elt-point p0) (kd-elt-point p1))))
			 (sorted (sort points cmpfn))
			 (median-index (quotient (length sorted) 2))
			 )

		    (let-values (((lte gte) (split-at sorted median-index)))

		      (let ((median (kd-elt-point (car gte))))

			(let-values (((lt xeq) (span (lambda (x) (< (coord axis (kd-elt-point x)) (coord axis median))) lte)))
			  
			  (if (null? xeq)
			      (values (car gte) lt (cdr gte))
			      (let ((split-index (length lt)))
				(values (car xeq) lt (append (cdr xeq) gte))))
			))
		      ))
		  ))

		(list->kd-tree/depth
		 (lambda (points depth bucket-size)

		   (cond

		    ((null? points) (KdLeaf cis:empty '() '() depth))
		    
		    ((<= (length points) bucket-size)

		       (let* (
			      (ps (map kd-elt-point points)) 
			      (ii (fold (lambda (x ax) (cis:add x ax)) cis:empty (map kd-elt-index points)))
			      (vs (let ((vs (map kd-elt-value points))) (and (every identity vs) vs)))
			     )
                        
			(KdLeaf ii ps vs (modulo depth dimension))
			)
		       )
		    
		    ((null? (cdr points))

		     (let* ((ps (map kd-elt-point points))
			    (ii (fold (lambda (x ax) (cis:add x ax)) cis:empty (map kd-elt-index points)))
			    (vs (map kd-elt-value points)))
		       
		       (KdLeaf ii ps vs (modulo depth dimension))
		       )
		     )
		    
		   (else
		    (let-values (((median lt gte)
				  (split points depth)))
		      
		      (let* ((depth1 (+ 1 depth))
			     (i (kd-elt-index median))
			     (p (kd-elt-point median))
			     (v (kd-elt-value median))
			     (axis (modulo depth dimension))
			     (l (list->kd-tree/depth lt depth1 bucket-size))
			     (r (list->kd-tree/depth gte depth1 bucket-size))
			     )

			(KdNode l p i v r axis 
				(cis:add i (cis:union (kd-tree-node-indices l) (kd-tree-node-indices r))))
			))
		    ))
                   ))
                )
        list->kd-tree/depth
        ))
  
  ;; Returns the nearest neighbor of p in tree t.
  
  (define=> (make-kd-tree-nearest-neighbor <Point>)
    (define (tree-empty? t) (cases kd-tree t (KdLeaf (ii pp vv axis) (cis:empty? ii)) (else #f)))
    (letrec ((find-nearest
	      (lambda (t1 t2 p probe xp x-probe)

		(let* ((candidates1 
			(let ((best1 (nearest-neighbor t1 probe)))
			  (or (and best1 (list best1 p)) (list p))))
		       
		       (sphere-intersects-plane? 
			(let ((v (- x-probe xp)))
			  (< (* v v) (dist2 probe (car candidates1)))))

		       (candidates2
			(if sphere-intersects-plane?
			    (let ((nn (nearest-neighbor t2 probe)))
			      (if nn (append candidates1 (list nn)) candidates1))
			    candidates1)))

		  (minimum-by candidates2 (lambda (a b) (negative? (compare-distance probe a b))))
		  )))

	     
	     (nearest-neighbor
	      (lambda (t probe)
		(cases kd-tree t
		       (KdLeaf (ii pp vv axis) 
			       (minimum-by pp (lambda (a b) (negative? (compare-distance probe a b)))))

		       (KdNode (l p i v r axis ci)
			       (if (and (tree-empty? l)
					(tree-empty? r)) p
					(let ((x-probe (coord axis probe))
					      (xp (coord axis p)))

					  (if (< x-probe xp) 
					      (find-nearest l r p probe xp x-probe) 
					      (find-nearest r l p probe xp x-probe))
					  ))
			       ))
		)))
      
      nearest-neighbor
      ))
  
  ;; Applies a user supplied function to the nearest neighbor of p in tree t.
  
  (define=> (make-kd-tree-nearest-neighbor* <Point>)

    (define (tree-empty? t) (cases kd-tree t (KdLeaf (ii pp vv axis) (cis:empty? ii)) (else #f)))

    (letrec ((find-nearest
	      (lambda (t1 t2 i v p probe xp x-probe fn)

		(let* ((candidates1 
			(let ((best1 (nearest-neighbor t1 probe fn)))
			  (or (and best1 (list best1 (if v (list (list i v) p) (list i p))))
			      (list (if v (list (list i v) p) (list i p))))
                          ))
		       
		       (sphere-intersects-plane? 
			(let ((v (- x-probe xp)))
			  (< (* v v) (dist2 probe (cadar candidates1)))))

		       (candidates2
			(if sphere-intersects-plane?
			    (let ((nn (nearest-neighbor t2 probe fn)))
			      (if nn (append candidates1 (list nn)) candidates1))
			    candidates1)))

		  (let ((v (minimum-by candidates2 (lambda (a b) (negative? (compare-distance probe (cadr a) (cadr b)))))))
		    v)
		  )))

	     
	     (nearest-neighbor
	      (lambda (t probe #!key (fn list))
		(cases kd-tree t
		       (KdLeaf (ii pp vv axis) 
			       (let ((res
				      (if vv
					  (minimum-by pp (lambda (a b) (negative? (compare-distance probe a b))) 
						      (zip (reverse (cis:elements ii)) vv ))
					  (minimum-by pp (lambda (a b) (negative? (compare-distance probe a b))) 
						      (reverse (cis:elements ii))))))
				 (and res (apply fn (reverse res)))))

		       (KdNode (l p i v r axis ci)

			       (if (and (tree-empty? l)
					(tree-empty? r))

                                   (if v (fn (list i v) p) (fn i p))

				   (let ((x-probe (coord axis probe))
					 (xp (coord axis p))
					 (xi i))
				     (if (< x-probe xp) 
					 (find-nearest l r i v p probe xp x-probe fn) 
					 (find-nearest r l i v p probe xp x-probe fn))
				     ))
			       ))
		)))
      
      nearest-neighbor
      ))
  
    
  ;; near-neighbors t r p returns all neighbors within distance r from p in tree t.

  (define=> (make-kd-tree-near-neighbors <Point>)
    (define-inline (tree-empty? t) (cases kd-tree t (KdLeaf (ii pp vv axis) (cis:empty? ii)) (else #f)))
    (define-inline (filter-fn probe pp d2)
      (filter-map (lambda (p) 
		    (let ((pd (dist2 probe p)))
		      (and (<= pd d2) (list p (sqrt pd) )))) pp))
    (letrec ((near-neighbors
	      (lambda (t radius probe)
		(cases kd-tree t
		       (KdLeaf (ii pp vv axis)  
			       (let ((r2 (* radius radius)))
				 (filter-fn probe pp r2)))

		       (KdNode (l p i v r axis ci)
			       (let ((maybe-pivot (filter-fn probe (list p) (* radius radius))))
				 
				 (if (and (tree-empty? l)
					  (tree-empty? r))

				     maybe-pivot

				     (let ((x-probe (coord axis probe))
					   (xp (coord axis p)))

				       (if (<= x-probe xp)

					   (let ((nearest (append maybe-pivot (near-neighbors l radius probe))))
					     (if (> (+ x-probe (abs radius)) xp)
						 (append (near-neighbors r radius probe) nearest)
						 nearest))

					   (let ((nearest (append maybe-pivot (near-neighbors r radius probe))))
					     (if (< (- x-probe (abs radius)) xp)
						 (append (near-neighbors l radius probe) nearest)
						 nearest)))
				       ))))
		       ))
	      ))
      
	  near-neighbors
	  ))


  (define=> (make-kd-tree-near-neighbors* <Point>)
    (define-inline (tree-empty? t) (cases kd-tree t (KdLeaf (ii pp vv axis) (cis:empty? ii)) (else #f)))
    (define-inline (filter-fn probe pp ii vv r2 fn)
      (if vv 
	  (filter-map (lambda (i v p) 
			(let ((pd (dist2 probe p))) 
			  (and (<= pd r2) (fn (list i v) p (sqrt pd)))))
		      (reverse (cis:elements ii)) vv pp)
	  (filter-map (lambda (i p)
			(let ((pd (dist2 probe p))) 
			  (and (<= pd r2) (fn i p (sqrt pd)))))
		      (reverse (cis:elements ii)) pp)
	  ))

    (letrec ((near-neighbors
	      (lambda (t radius probe fn)
		(cases kd-tree t

		       (KdLeaf (ii pp vv axis)  
			       (let ((r2 (* radius radius)))
				 (let ((res (filter-fn probe pp ii vv r2 fn)))
				   res)))

		       (KdNode (l p i v r axis ci)
			       (let ((maybe-pivot (filter-fn probe (list p) (cis:singleton i) (and v (list v)) (* radius radius) fn)))

				 (if (and (tree-empty? l)
					  (tree-empty? r))

				     maybe-pivot

				     (let ((x-probe (coord axis probe))
					   (xp (coord axis p)))

				       (if (<= x-probe xp)

					   (let ((nearest (append maybe-pivot (near-neighbors l radius probe fn))))
					     (if (> (+ x-probe (abs radius)) xp)
						 (append (near-neighbors r radius probe fn) nearest)
						 nearest))

					   (let ((nearest (append maybe-pivot (near-neighbors r radius probe fn))))
					     (if (< (- x-probe (abs radius)) xp)
						 (append (near-neighbors l radius probe fn) nearest)
						 nearest)))


				       ))
				 ))
		       ))
	      ))
      (lambda (t radius probe #!key (fn list))
	(near-neighbors t radius probe fn)
	))
    )
  


  
  ;; Returns the k nearest points to p within tree.
  (define=> (make-kd-tree-k-nearest-neighbors <Point>)
    (lambda (kd-tree-remove kd-tree-nearest-neighbor)
      (letrec ((k-nearest-neighbors
		(lambda (t k probe)
		  (cases kd-tree t

		       (KdLeaf (ii pp vv axis) 
			       (let recur ((res '()) (pp pp) (k k))
				 (if (or (<= k 0) (null? pp)) 
				     res
				     (let ((nearest (minimum-by pp (lambda (a b) (negative? (compare-distance probe a b))))))
				       (recur (cons nearest res)
					      (remove (lambda (p) (equal? p nearest)) pp)
					      (- k 1))
				       ))
				 ))

		       (else
			(if (<= k 0) '()
			    (let* ((nearest (kd-tree-nearest-neighbor t probe))
				   (tree1 (kd-tree-remove t nearest)))
			      (cons nearest (k-nearest-neighbors tree1 (- k 1) probe)))
			    ))
		       ))
		))
	k-nearest-neighbors)))

  
  ;; removes the point p from t.
  (define=> (make-kd-tree-remove <Point>)
    
    (lambda (list->kd-tree/depth)

      (letrec (

	       (tree-remove
		(lambda (t p-kill #!key (bucket-size (* 10 (max (log2 (kd-tree-size t)) 1))) (tol 1e-9))
                  
                  (let ((tol^2 (* tol tol)))
                  
;                  (fprintf (current-error-port) "kd-tree-remove: t = ~A~%" t)
;                  (fprintf (current-error-port) "kd-tree-remove: p-kill = ~A~%" p-kill)

		  (cases kd-tree t
			 (KdLeaf (ii pp vv axis) 

;                                 (fprintf (current-error-port) "kd-tree-remove (KdLeaf): ii = ~A~%" (cis:elements ii))
;                                 (fprintf (current-error-port) "kd-tree-remove (KdLeaf): vv = ~A~%" vv)
;				 (fprintf (current-error-port) "kd-tree-remove (KdLeaf): pp = ~A~%" pp)
					 
				 (if vv
				     (let ((ipvs (filter-map
						  (lambda (i p v) (and (< (dist2 p p-kill) tol^2) (list i p v)))
						  (reverse (cis:elements ii)) pp vv)))

;					 (fprintf (current-error-port) "kd-tree-remove (KdLeaf): ipvs = ~A~%" ipvs)

				       (and (pair? ipvs)
					    (let ((ii1 (fold (lambda (i ax) (cis:remove i ax)) 
							     ii (map car ipvs)))
						  (pp1 (fold (lambda (x ax) (remove (lambda (p) (equal? x p)) ax))
							     pp (map cadr ipvs)))
						  (vv1 (fold (lambda (x ax)
							       (remove (lambda (v) (equal? x v)) ax))
							     vv (map caddr ipvs)))
						  )
;                                              (fprintf (current-error-port) "kd-tree-remove (KdLeaf): ii1 = ~A~%" (cis:elements ii1))
;                                              (fprintf (current-error-port) "kd-tree-remove (KdLeaf): pp1 = ~A~%" pp1)
;                                              (fprintf (current-error-port) "kd-tree-remove (KdLeaf): vv1 = ~A~%" vv1)

					      (KdLeaf ii1 pp1 vv1 axis))
					    ))

				     (let ((ips (filter-map (lambda (i p) (and (< (dist2 p p-kill) tol^2) (list i p)))
							    (reverse (cis:elements ii)) pp)))
				       
;				       (fprintf (current-error-port) "kd-tree-remove (KdLeaf): ips = ~A~%" ips)

				       (and (pair? ips)
					    (let ((ii1 (fold (lambda (i ax) (cis:remove i ax)) 
							     ii (map car ips)))
						  (pp1 (fold (lambda (x ax) (remove (lambda (p) (equal? x p)) ax))
							     pp (map cadr ips)))
						  )

;                                              (fprintf (current-error-port) "kd-tree-remove (KdLeaf): ii1 = ~A~%" (cis:elements ii1))
;                                              (fprintf (current-error-port) "kd-tree-remove (KdLeaf): pp1 = ~A~%" pp1)

                                              (KdLeaf ii1 pp1 vv axis))
					    ))
				     ))

			 (KdNode (l p i v r axis ci)

				 (cond ((< (dist2 p p-kill) tol^2)
					(let ((pts1 (append (kd-tree->list* l) (kd-tree->list* r))))
					  (list->kd-tree/depth pts1 axis bucket-size)))

				       (else

					(if (< (coord axis p-kill) (coord axis p))

					 (let* ((l1   (tree-remove l p-kill))
						(l1-is (and l1 (kd-tree-node-indices l1)))
						(r-is  (kd-tree-node-indices r))
						(ci1   (cis:add i (cis:union l1-is r-is))))

					   (and l1 (KdNode l1 p i v r axis ci1))
					   )
					 
					 (let* ((r1   (tree-remove r p-kill))
						(r1-is (and r1 (kd-tree-node-indices r1)))
						(l-is  (kd-tree-node-indices l))
						(ci1   (cis:add i (cis:union r1-is l-is))))

					   (and r1 (KdNode l p i v r1 axis ci1))

					   )))
				     
				     ))
			 ))
		)))
	tree-remove))
    )

  ;; Checks whether the K-D tree property holds for a given tree.
  ;;
  ;; Specifically, it tests that all points in the left subtree lie to
  ;; the left of the plane, p is on the plane, and all points in the
  ;; right subtree lie to the right.
  
  (define=> (make-kd-tree-is-valid? <Point>)
    (lambda (t)
      (cases kd-tree t
	     (KdLeaf (ii pp vv axis)  #t)

	     (KdNode (l p i v r axis ci)
		     (let ((x (coord axis p)))
		       (and (every (lambda (y) (< (coord axis y) x ))
				   (kd-tree->list l))
			    (every (lambda (y) (>= (coord axis y) x))
				   (kd-tree->list r)))))
	     )))
  
  
  ;; Checks whether the K-D tree property holds for the given tree and
  ;; all subtrees.
  
  (define (make-kd-tree-all-subtrees-are-valid? kd-tree-is-valid?)
    (lambda (t) (every kd-tree-is-valid? (kd-tree-subtrees t))))
  

  (define=> (make-kd-tree-slice <Point>)
    (lambda (x-axis x1 x2 t)
      (let recur ((t t)  (pts '()))
	(cases kd-tree t

	       (KdLeaf (ii pp vv axis) 
		       (append (filter (lambda (p) 
					 (and (<= x1 (coord x-axis p))
					      (<= (coord x-axis p) x2)))
				       pp)
			       pts))
			   

	       (KdNode (l p i v r axis ci)
		       (if (= axis x-axis)
			   
			   (cond ((and (<= x1 (coord axis p))
				       (<= (coord axis p) x2))
				   (recur l (cons p (recur r pts))))
				 
				 ((< (coord axis p) x1)
				  (recur r pts))
			       
				 ((< x2 (coord axis p))
				  (recur l pts)))
			   
			   (if (and (<= x1 (coord x-axis p))
				    (<= (coord x-axis p) x2))
			       (recur l (cons p (recur r pts)))
			       (recur l (recur r pts)))
			   ))
	       ))
      ))
  
  
  (define=> (make-kd-tree-slice* <Point>)
    (lambda (x-axis x1 x2 t #!key (fn list))
      (let recur ((t t)  (pts '()))
	(cases kd-tree t
	       (KdLeaf (ii pp vv axis) 
		       (append
                        (if vv
                            (filter-map (lambda (i v p) 
                                          (and (<= x1 (coord x-axis p))
                                               (<= (coord x-axis p) x2)
                                               (fn (list i v) p)))
                                        (reverse (cis:elements ii)) vv pp)
                            (filter-map (lambda (i p) 
                                          (and (<= x1 (coord x-axis p))
                                               (<= (coord x-axis p) x2)
                                               (fn i p)))
                                        (reverse (cis:elements ii)) pp))
                        pts))

	       (KdNode (l p i v r axis ci)
		       (if (= axis x-axis)
			   
			   (cond ((and (<= x1 (coord axis p))
				       (<= (coord axis p) x2))
				   (recur l (cons (if v (fn (list i v) p) (fn i p)) (recur r pts))))
				 
				 ((< (coord axis p) x1)
				  (recur r pts))
			       
				 ((< x2 (coord axis p))
				  (recur l pts)))
			   
			   (if (and (<= x1 (coord x-axis p))
				    (<= (coord x-axis p) x2))
			       (recur l (cons (if v (fn (list i v) p) (fn i p)) (recur r pts)))
			       (recur l (recur r pts)))
			   ))
	       ))
      ))

  (define (default-<KdTree> point-class)
    (let* ((list->kd-tree/depth
	    (make-list->kd-tree/depth point-class))
           (list->kd-tree/depth*
	    (make-list->kd-tree/depth* point-class))
	   (kd-tree-remove
	    ((make-kd-tree-remove point-class) list->kd-tree/depth*))
	   (kd-tree-nearest-neighbor
	    (make-kd-tree-nearest-neighbor point-class)))

      (make-<KdTree> 

       (lambda (points #!key
                       (make-point identity) 
                       (make-value #f)
                       (offset 0)
                       )
	 ((list->kd-tree/depth make-value) 
          0 (length points) 
	  (map (lambda (p) (cons (make-point p) p)) points) 0 
          offset: offset))

       (lambda (points) (list->kd-tree/depth* points 0))

       (make-kd-tree-nearest-neighbor point-class)
       (make-kd-tree-nearest-neighbor* point-class)
       (make-kd-tree-near-neighbors point-class)
       (make-kd-tree-near-neighbors* point-class)
       ((make-kd-tree-k-nearest-neighbors point-class)
	kd-tree-remove kd-tree-nearest-neighbor)
       kd-tree-remove
       (make-kd-tree-slice point-class)
       (make-kd-tree-slice* point-class)
       (make-kd-tree-is-valid? point-class)
       (make-kd-tree-all-subtrees-are-valid? 
	(make-kd-tree-is-valid? point-class)) 
       )))

  (define KdTree3d 
    (default-<KdTree> Point3d))

  (define KdTree2d 
    (default-<KdTree> Point2d))



)

