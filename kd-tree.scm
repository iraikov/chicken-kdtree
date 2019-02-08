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
;; Copyright 2012-2019 Ivan Raikov
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
   make-kd-tree
   spatial-map?
   empty?
   size
   dimension
   spatial-map-for-each
   spatial-map-fold-right
   spatial-map-fold-right*
   nearest-neighbor
   near-neighbors
   k-nearest-neighbors
   remove
   slice
   get-kspace
   spatial-map->list
   is-valid?
   all-subtrees-are-valid?
   )

  (import scheme (chicken base) (chicken foreign) (prefix (chicken sort) list.)
          datatype yasos yasos-collections 
          (only srfi-1 fold list-tabulate split-at span every fold-right take filter filter-map zip)
	   (only srfi-4 f32vector? f32vector f32vector-ref)
	   (only (chicken format) fprintf) (only (chicken pretty-print) pp)
           kspace)


  (define log2 (foreign-lambda double "log2" double))

  (define-predicate  spatial-map?)

  ;; nearest neighbor of a point
  (define-operation (nearest-neighbor smap p))
  ;; neighbors of a point within radius r
  (define-operation (near-neighbors smap p r))
  ;; k nearest neighbors of a point
  (define-operation (k-nearest-neighbors smap p k))
  ;; removes a point from the tree
  (define-operation (remove smap p))
  ;; retrieves all points between two planes
  (define-operation (slice smap l u))

  (define (range m n)
    (if (< n m) (range n m)
        (list-tabulate (- n m) (lambda (i) (+ i m)))))

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
            ))
    )

  
  (define (split kspace sorted axis)

    (let* ((median-index (quotient (length sorted) 2)))

      (let-values (((lte gte) (split-at sorted median-index)))

        
        (let ((median (car gte)))
          
          (let-values (((lt xeq) (span (lambda (x) (< (coord kspace x axis)
                                                      (coord kspace median axis)))
                                       lte)))
            
            (if (null? xeq)
                (values median lt (cdr gte))
                (let ((split-index (length lt)))
                  (values (car xeq) lt (append (cdr xeq) gte))))
            ))
        ))
    )




  (define (positive-or-zero-integer? x)
    (and (integer? x) (or (zero? x) (positive? x))))


  (define (positive-integer? x)
    (and (integer? x) (positive? x)))


  (define-datatype kd-tree kd-tree?
    (KdNode (left  kd-tree?)
            (i     positive-or-zero-integer?)
            (right kd-tree?)
            (axis  positive-or-zero-integer?)
            )
    (KdLeaf (ii list?)
            (axis positive-or-zero-integer?) )
    )


  (define (kd-tree-empty? t)
    (cases kd-tree t
           (KdLeaf (ii axis) (null? ii))
           (else #f)))
  
  (define (kd-tree-map f kspace t)
    (cases kd-tree t
           (KdLeaf (ii axis) 
                   (KdLeaf (map (lambda (i) (f (point kspace i))) ii) axis))
           (KdNode (l i r axis)
                   (KdNode (kd-tree-map f kspace l)
                           (f (point kspace i))
                           (kd-tree-map f kspace r)
                           axis))
           ))
  
  (define (kd-tree-for-each f kspace t)
    (cases kd-tree t
           (KdLeaf (ii axis) (for-each (lambda (i) (f (point kspace i))) ii))
           (KdNode (l i r axis)
                   (begin
                     (kd-tree-for-each f kspace l)
                     (f (point kspace i))
                     (kd-tree-for-each f kspace r)
                     ))
           ))


  
  (define (kd-tree-fold-right f init kspace t)
    (cases kd-tree t
           (KdLeaf (ii axis) 
                   (fold-right (lambda (i ax) (f (point kspace i) ax)) init ii))
           (KdNode (l i r axis)
                   (let* ((init2 (kd-tree-fold-right f init kspace r))
                          (init3 (f (point kspace i) init2)))
                     (kd-tree-fold-right f init3 kspace l)))
           ))

  
  (define (kd-tree-fold-right* f init kspace t)
    (cases kd-tree t
           (KdLeaf (ii axis) 
                   (fold-right (lambda (i ax) (f i (point kspace i) ax)) init ii))
           (KdNode (l i r axis)
                   (let* ((init2 (kd-tree-fold-right* f init kspace r))
                          (init3 (f i (point kspace i) init2)))
                     (kd-tree-fold-right* f init3 kspace l)))
           ))

  (define (kd-tree->list kspace t)
    (kd-tree-fold-right* (lambda (i p ax) (cons i ax)) '() kspace t))

  
  ;; Returns a list containing t and all its subtrees, including the
  ;; leaf nodes.
  
  (define (kd-tree-subtrees t)
    (cases kd-tree t
           (KdLeaf (ii axis)  (list t))
           (KdNode (l i r axis)
                   (append (kd-tree-subtrees l) 
                           (list t)
                           (kd-tree-subtrees r)))
           ))

  

  (define (kd-tree-size t)
    (cases kd-tree t
           (KdLeaf (ii axis) (length ii))
           (KdNode (l i r axis) (+ 1 (length l) (length r)))))
  
  
  ;; (define (insert pseq tree k i #!key (leaf-size (* 4 (max (log2 (kd-tree-size tree)) 1))))
  ;;   (let ((point (elt-ref pseq i)))
  ;;     (cases kd-tree t
  ;;            (KdLeaf (ii axis) 
  ;;                    (if (null? ii)
  ;;                        (KdLeaf (list i) i axis)
  ;;                        (let ((ii1 (merge (list i) ii (lambda (x y) (< (coord x axis) (coord y axis ))))))
  ;;                          (if (> (length ii1) leaf-size)
  ;;                              (let
  ;;                                  ((axis1   (modulo (+ 1 axis) k))
  ;;                                   (median  (list-ref ii1 (quotient (length ii1) 2)))
  ;;                                   (medianc (coord median axis)))
  ;;                                (let-values (((ii-left ii-right)
  ;;                                              (partition (lambda (x) (< (coord x axis) medianc))
  ;;                                                         ii1)))
  ;;                                  (let ((left (KdLeaf ii-left axis1))
  ;;                                        (right (KdLeaf (remove median ii-right) axis1)))
  ;;                                    (KdNode left median right axis))))
  ;;                              (KdLeaf ii1 pp1 axis)))
  ;;                        ))
  ;;          (KdNode (left i right axis)::ts => 
  ;;                 if n > nodesize
  ;;                 then (case lst of 
  ;;                           u::rest => addPoint' (nodesize,leafsize) (P,[KdLeaf{ii=[],axis=0}],j,n-1,joinTrees P (t,u)::rest)
  ;;                         | [] => addPoint' (nodesize,leafsize) (P,[KdLeaf{ii=[],axis=0}],j,n,[t]))
  ;;                 else addPoint' (nodesize,leafsize) (P,ts,j,n+1,t::lst)
  ;;               | [] => (KdLeaf {ii=[j],axis=0}) :: lst

  
  
  
  
  ;; Construct a kd-tree from an kspace
  (define (make kspace #!key (axis 0) (points (range 0 (size kspace))))

    (letrec (
             (k (dimension kspace))
             (axial-compare (lambda (axis) (lambda (p0 p1) (compare-coord kspace p0 p1 axis))))
             (make/depth
              (lambda (points depth #!key (bucket-size (* 10 (max (log2 (length points)) 1))))

                (let* ((axis (modulo depth k))
                       (ii (list.sort points (axial-compare axis)))
                       (extent (length ii)))
                  
                  (if (<= extent bucket-size)

                      (KdLeaf ii axis)
                      
                      (let-values (((median lt gte) (split kspace ii axis)))
                        
                        (KdNode (make/depth lt (add1 depth) bucket-size: bucket-size)
                                median
                                (make/depth gte (add1 depth) bucket-size: bucket-size)
                                axis )
                        ))
                  ))
              ))
             
      (make/depth points axis)
      )
    )
    
    ;; Returns the nearest neighbor of p in tree t.
    
    (define (kd-tree-nearest-neighbor kspace t probe)

      (let ((find-nearest
             (lambda (t1 t2 p probe xp x-probe)
               
               (let* ((candidates1 
                       (let ((best1 (kd-tree-nearest-neighbor kspace t1 probe)))
                         (or (and best1 (list best1 p)) (list p))))
                      (sphere-intersects-plane? 
                       (let ((v (- x-probe xp)))
                         (< (* v v) (squared-distance kspace probe (car candidates1)))))
                      
                      (candidates2
                       (if sphere-intersects-plane?
                           (let ((nn (kd-tree-nearest-neighbor kspace t2 probe)))
                             (if nn (append candidates1 (list nn)) candidates1))
                           candidates1)))
                 
                 (minimum-by
                  candidates2
                  (lambda (a b) (negative? (compare-distance kspace probe a b))))
                 ))
             ))

        (cases kd-tree t
               (KdLeaf (ii axis)
                       (let ((res (minimum-by
                                   ii
                                   (lambda (a b) (negative? (compare-distance kspace probe a b))))))
                         (and res (point kspace res))))

               (KdNode (l i r axis)
                       (let ((x-probe (list-ref probe axis))
                             (xp (coord kspace i axis)))
                         (if (< x-probe xp)
                             (find-nearest l r i probe xp x-probe) 
                             (find-nearest r l i probe xp x-probe)
                             ))
                       ))
        ))

    
    ;; near-neighbors t r p returns all neighbors within distance r from p in tree t.

    (define (kd-tree-near-neighbors kspace t radius probe)

      (define (filter-fn probe pp d2)
        (filter-map (lambda (p) 
                      (let ((pd (squared-distance kspace probe p)))
                        (and (<= pd d2) p )))
                    pp))

      (define (get-point p) (point kspace p))
      
      (cases kd-tree t
             (KdLeaf (ii axis)  
                     (let ((r2 (* radius radius)))
                       (map get-point (filter-fn probe ii r2))))
             
             (KdNode (l p r axis)
                     (let ((maybe-pivot (filter-fn probe (list p) (* radius radius))))
                       
                       (if (and (kd-tree-empty? l)
                                (kd-tree-empty? r))
                           
                           (map get-point maybe-pivot)
                           
                           (let ((x-probe (coord kspace probe axis))
                                 (xp (coord kspace p axis)))
                             
                             (if (<= x-probe xp)
                                 
                                 (let ((nearest
                                        (append (map get-point maybe-pivot)
                                                (kd-tree-near-neighbors kspace l radius probe))))
                                   (if (> (+ x-probe (abs radius)) xp)
                                       (append (kd-tree-near-neighbors kspace r radius probe) nearest)
                                       nearest))
                                 
                                 (let ((nearest
                                        (append (map get-point maybe-pivot)
                                                (kd-tree-near-neighbors kspace r radius probe))))
                                   (if (< (- x-probe (abs radius)) xp)
                                       (append (kd-tree-near-neighbors kspace l radius probe) nearest)
                                       nearest)))
                             ))
                       ))
             ))




    
    ;; Returns the k nearest points to p within tree.
    (define (kd-tree-k-nearest-neighbors kspace t k probe)

      (define (get-point p) (point kspace p))

      (cases kd-tree t
             
             (KdLeaf (ii axis) 
                     (let recur ((res '()) (pp pp) (k k))
                       (if (or (<= k 0) (null? pp))
                           (map get-point res)
                           (let ((nearest
                                  (minimum-by pp
                                    (lambda (a b) (negative? (compare-distance kspace probe a b))))))
                             (recur (cons nearest res)
                                    (filter (lambda (p) (not (equal? p nearest))) ii)
                                    (- k 1))
                             ))
                       ))
             
             (else
              (if (<= k 0) '()
                  (let* ((nearest (kd-tree-nearest-neighbor kspace t probe))
                         (tree1 (kd-tree-remove kspace t nearest 1e-3)))
                    (cons nearest (kd-tree-k-nearest-neighbors kspace tree1 (- k 1) probe)))
                  ))
             ))


    
    ;; removes the point p from t.
    (define (kd-tree-remove kspace t p-kill tol)
      
      (let ((tol^2 (* tol tol)))

        (cases kd-tree t
               (KdLeaf (ii axis)
                       (let ((ii1
                              (filter
                               (lambda (p)
                                 (> (squared-distance kspace p p-kill) tol^2))
                               ii)))
                         
                         (KdLeaf ii1 axis)))

               (KdNode (l p r axis)
                       (cond ((< (squared-distance kspace p p-kill) tol^2)
                              (let ((pts1 (append (kd-tree->list kspace l)
                                                  (kd-tree->list kspace r))))
                                (make kspace points: pts1 axis: axis)))
                             (else
                              
                              (if (< (if (list? p-kill) (list-ref p-kill axis)
                                         (coord kspace p-kill axis))
                                     (coord kspace p axis))
                                  
                                  (let* ((l1 (kd-tree-remove kspace l p-kill tol)))
                                    (and l1 (KdNode l1 p r axis))
                                    )
                                  
                                  (let* ((r1   (kd-tree-remove kspace r p-kill tol)))
                                    (and r1 (KdNode l p r1 axis))
                                    ))
                             
                              ))
               ))
      ))

    
    ;; Checks whether the K-D tree property holds for a given tree.
    ;;
    ;; Specifically, it tests that all points in the left subtree lie to
    ;; the left of the plane, p is on the plane, and all points in the
    ;; right subtree lie to the right.
    
    (define (kd-tree-is-valid? kspace t)
      (cases kd-tree t
             (KdLeaf (ii axis)  #t)
             
             (KdNode (l p r axis)
                     (let ((x (coord kspace p axis)))
                       (and (every (lambda (y) (< (coord kspace y axis) x ))
                                   (kd-tree->list kspace l))
                            (every (lambda (y) (>= (coord kspace y axis) x))
                                   (kd-tree->list kspace r)))))
             ))
    
    
    ;; Checks whether the K-D tree property holds for the given tree and
    ;; all subtrees.
    
    (define (kd-tree-all-subtrees-are-valid? kspace t)
      (every (lambda (t) (kd-tree-is-valid? kspace t))
             (kd-tree-subtrees t)))
    

    (define (kd-tree-slice kspace x-axis x1 x2 t)
      (define (get-point p) (point kspace p))
      (let recur ((t t)  (pts '()))
        (cases kd-tree t

               (KdLeaf (ii axis) 
                       (append
                        (map get-point
                             (filter (lambda (p) 
                                       (and (<= x1 (coord kspace p x-axis))
                                            (<= (coord kspace p x-axis) x2)))
                                     ii))
                        pts))
               

               (KdNode (l p r axis)
                       (if (= axis x-axis)
                           
                           (cond ((and (<= x1 (coord kspace p axis))
                                       (<= (coord kspace p axis) x2))
                                  (recur l (cons (get-point p) (recur r pts))))
                                 
                                 ((< (coord kspace p axis) x1)
                                  (recur r pts))
                                 
                                 ((< x2 (coord kspace p axis))
                                  (recur l pts)))
                           
                           (if (and (<= x1 (coord kspace p x-axis))
                                    (<= (coord kspace p x-axis) x2))
                               (recur l (cons (get-point p) (recur r pts)))
                               (recur l (recur r pts)))
                           ))
               ))
      )
    
  ;; kspace of the spatial map
  (define-operation (get-kspace smap))
  ;; nearest neighbor of a point
  (define-operation (nearest-neighbor smap p))
  ;; neighbors of a point within radius r
  (define-operation (near-neighbors smap p r))
  ;; k nearest neighbors of a point
  (define-operation (k-nearest-neighbors smap p k))
  ;; removes a point from the tree
  (define-operation (remove smap p tol))
  ;; retrieves all points between two planes
  (define-operation (slice smap axis l u))
  (define-operation (is-valid? smap))
  (define-operation (all-subtrees-are-valid? smap))

  ;; standard iterators
  (define-operation (spatial-map->list smap))
  (define-operation (spatial-map-fold-right smap f init))
  (define-operation (spatial-map-fold-right* smap f init))
  (define-operation (spatial-map-for-each smap f))

    
  (define (make-kd-tree kspace)

    (let ((kdt (make kspace)))
      
      (object
       ((spatial-map? self) #t)
       ((empty? self) (kd-tree-empty? kdt))
       ((get-kspace self) kspace)
       ((size self) (size kspace))
       ((dimension self) (dimension kspace))
       ((nearest-neighbor self p)
        (kd-tree-nearest-neighbor kspace kdt p))
       ((near-neighbors self p r)
        (kd-tree-near-neighbors kspace kdt r p))
       ((k-nearest-neighbors self p k)
        (kd-tree-k-nearest-neighbors kspace kdt k p))
       ((remove self p tol)
        (kd-tree-remove kspace kdt p tol))
       ((slice self axis l u)
        (kd-tree-slice kspace axis l u kdt))
       ((is-valid? self)
        (kd-tree-is-valid? kspace kdt))
       ((all-subtrees-are-valid? self)
        (kd-tree-all-subtrees-are-valid? kspace kdt))
       ((spatial-map-for-each self f)
        (kd-tree-for-each f kspace kdt))
       ((spatial-map-fold-right self f init)
        (kd-tree-fold-right f init kspace kdt))
       ((spatial-map-fold-right* self f init)
        (kd-tree-fold-right* f init kspace kdt))
       ((spatial-map->list self)
        (kd-tree->list kspace kdt))
       ))
    )





)
