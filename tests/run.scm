;; Unit tests for kd-tree


(import scheme (chicken base) (chicken sort)
        test srfi-1 srfi-4 yasos random-mtzig kspace kd-tree)

(define (consistency-tests kdt pts x k r)
  (define kspace (get-kspace kdt))

  (define (sort-points x pts) 
    (sort pts (lambda (a b)
                (negative? (compare-distance kspace x a b)))))
  
;    (define sorted-points (sort-points pts))
;    (define sliced-points (sort-points (filter (lambda (p) (and (<= 0. (coord kspace 0 p))
;                                                                (<= (coord kspace 0 p) 0.5)))
;                                               pts)))
    (define tree-points (spatial-map->list kdt))

;    (define nn (car sorted-points))

;    (define knn (if (> (length sorted-points) k)
;		    (take sorted-points k) 
;		    sorted-points))

;    (define rnn (let ((r2 (* r r)))
;		  (sort-points-rnn
;		   (filter-map
;		    (lambda (y) (let ((d (dist2 x y))) (and (<= d r2) (list y (sqrt d)))))
;		    sorted-points)
;		  )))

;    (define knn-removed 
;      (remove (lambda (x) (member x knn)) sorted-points))
    
    (test-group
     "KD tree consistency"
     
     (test-assert "tree-is-valid?"
	   (is-valid? kdt))

     (test-assert "tree-all-subtrees-are-valid?"
	   (all-subtrees-are-valid? kdt))

     (test-assert "tree-size"
	   (= (size kdt) (f32vector-length (car pts))))
     )
    
    (print "nearest-neighbor of " x ": " (nearest-neighbor kdt x))
    (print k " nearest neighbors of " x ": " (k-nearest-neighbors kdt x k))
    (print "near neighbors of " x " within " r ": " (near-neighbors kdt k r))
)

 
;;      (test "kd-tree->list"
;; 	  sorted-points tree-points)

;;      (test "nearest-neighbor"
;; 	   nn
;; 	   (nearest-neighbor kdt x))

;;      (test "k-nearest-neighbors"
;; 	   knn
;; 	   (sort-points (k-nearest-neighbors kdt k x)))
     
;;      (test "near-neighbors"
;; 	   rnn
;; 	   (sort-points-rnn (near-neighbors kdt r x )))


;;      (test "slice"
;; 	   sliced-points
;; 	   (sort-points (kd-tree-slice 0 0. 0.5 t)))

;;      (test "remove"
;; 	   knn-removed
;; 	   (sort-points (kd-tree->list (fold (lambda (x ax) (kd-tree-remove ax x)) t knn))))

;;      ))

;; (define elt< (lambda (i ai j bj) (< ai bj)))


(let ((n (inexact->exact 1e5)) (k 40) (r 1.0) (randst (init 9)))
  
  (let recur ((ntrials 1))

    (let* ((xs (randn/f32! n randst))
	   (ys (randn/f32! n randst))
	   (zs (randn/f32! n randst)))

      (print "random coordinates generated...")
      
      (let* (
             (pts1 (list xs ys))
             (pts2 (list xs ys zs))
             (kspace2d (make-space pts1))
             (kspace3d (make-space pts2))
             (dd    (print "tree construction..."))
             (t1    (time (make-kd-tree kspace2d)))
             (dd    (print "tree 1 constructed!"))
             (t2    (time (make-kd-tree kspace3d)))
             (dd    (print "tree 2 constructed!"))
             (x1    (inexact->exact (modulo (random! randst) n)))
             (x2    (inexact->exact (modulo (random! randst) n)))
             (xx1    (list (+ 0.1 (coord kspace2d x1 0)) 
                           (- (coord kspace2d x1 1) 0.1)))
             (xx2    (list (+ 0.1 (coord kspace2d x2 0)) 
                           (- (coord kspace3d x2 1) 0.1) 
                           (+ 0.1 (coord kspace3d x2 2))))
             )

        (time (consistency-tests t1 pts1 xx1 k r))
        (time (consistency-tests t2 pts2 xx2 k r))
        
        (if (positive? ntrials)
            (recur (- ntrials 1)))
        ))
    ))
  


(test-exit)

