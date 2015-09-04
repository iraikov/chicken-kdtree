;; Unit tests for kd-tree


(use typeclass kd-tree test srfi-1 random-mtzig)

(define=> (make-consistency-tests <Point> <KdTree>)
  (lambda (pts t x k r)

    (define (sort-points pts) 
      (sort pts (lambda (a b) (negative? (compare-distance x a b)))))
    (define (sort-points-rnn pts) 
      (sort pts (lambda (a b) (negative? (compare-distance x (car a) (car b))))))
    (define (sort-points* pts) 
      (sort pts (lambda (a b) (negative? (compare-distance x (cadr a) (cadr b))))))

    (define sorted-points (sort-points pts))
    (define sliced-points (sort-points (filter (lambda (p) (and (<= 0. (coord 0 p)) (<= (coord 0 p) 0.5))) pts)))
    (define tree-points (sort-points (kd-tree->list t)))
    (define tree-points-indices (kd-tree->list* t))

    (define nn (car sorted-points))

    (define knn (if (> (length sorted-points) k)
		    (take sorted-points k) 
		    sorted-points))

    (define rnn (let ((r2 (* r r)))
		  (sort-points-rnn
		   (filter-map
		    (lambda (y) (let ((d (dist2 x y))) (and (<= d r2) (list y (sqrt d)))))
		    sorted-points)
		  )))

    (define rnn* (map cons (map (lambda (nn) (list-index (lambda (x) (equal? (cadr x) (car nn))) tree-points-indices)) rnn) rnn))

    (define knn-removed 
      (remove (lambda (x) (member x knn)) sorted-points))
    
    (test-group
     "KD tree consistency"
     
     (test-assert "monotonically increasing indices"
	   (sorted? tree-points-indices (lambda (x y) (< (car x) (car y)))))

     (test-assert "tree-is-valid?"
	   (kd-tree-is-valid? t))

     (test-assert "tree-all-subtrees-are-valid?"
	   (kd-tree-all-subtrees-are-valid? t))

     (test-assert "tree-size"
	   (= (kd-tree-size t) (length pts)))

     (test "kd-tree->list"
	  sorted-points tree-points)

     (test "nearest-neighbor"
	   nn
	   (kd-tree-nearest-neighbor t x))

     (test "nearest-neighbor*"
	   (list (list-index (lambda (x) (equal? (cadr x) nn)) tree-points-indices) nn)
	   (kd-tree-nearest-neighbor* t x))

     (test "k-nearest-neighbors"
	   knn
	   (sort-points (kd-tree-k-nearest-neighbors t k x)))
     
     (test "near-neighbors"
	   rnn
	   (sort-points-rnn (kd-tree-near-neighbors t r x )))

     (test "near-neighbors*"
	   rnn*
	   (sort-points* (kd-tree-near-neighbors* t r x)))

     (test "slice"
	   sliced-points
	   (sort-points (kd-tree-slice 0 0. 0.5 t)))

     (test "remove"
	   knn-removed
	   (sort-points (kd-tree->list (fold (lambda (x ax) (kd-tree-remove ax x)) t knn))))

     )))

(define consistency-tests/3d (make-consistency-tests Point3d KdTree3d))
(define consistency-tests/2d (make-consistency-tests Point2d KdTree2d))
(define elt< (lambda (i ai j bj) (< ai bj)))


(define (test1)
  (let ((n (inexact->exact 1e5)) (k 40) (r 0.2) (randst (init)))
  
  (let recur ((ntrials 1))

    (let* ((xs (f64vector-randn! n randst))
	   (ys (f64vector-randn! n randst))
	   (zs (f64vector-randn! n randst)))

      (print "random coordinates generated...")
      
      (let* (
             (pts1 (let precur ((i (- n 1)) (ax '()))
                     (let ((p (make-point (f64vector-ref xs i)
                                          (f64vector-ref ys i)
                                          (f64vector-ref zs i))))
                       (if (positive? i)
                           (precur (- i 1) (cons p ax))
                           (cons p ax)))))
             (pts2 (let precur ((i (- n 1)) (ax '()))
                     (let ((p (make-point (f64vector-ref xs i)
                                          1.0
                                          (f64vector-ref zs i))))
                       (if (positive? i)
                           (precur (- i 1) (cons p ax))
                           (cons p ax)))))
             (dd    (print "tree construction..."))
             (t1    (time (with-instance ((<KdTree> KdTree3d)) (list->kd-tree pts1))))
             (dd    (print "tree 1 constructed!"))
             (t2    (time (with-instance ((<KdTree> KdTree3d)) (list->kd-tree pts2))))
             (dd    (print "tree 2 constructed!"))
             (xi1    (inexact->exact (modulo (random! randst) n)))
             (xi2    (inexact->exact (modulo (random! randst) n)))
             (x1     (list-ref pts1 xi1))
             (x2     (list-ref pts2 xi2))
             (xx1    (with-instance ((<Point> Point3d)) 
                                    (make-point (+ 0.1 (coord 0 x1)) 
                                                (- (coord 1 x1) 0.1) 
                                                (+ 0.1 (coord 2 x1)))))
             (xx2    (with-instance ((<Point> Point3d)) 
                                    (make-point (+ 0.1 (coord 0 x2)) 
                                                (- (coord 1 x2) 0.1) 
                                                (+ 0.1 (coord 2 x2)))))
             )

        (time (consistency-tests/3d pts1 t1 xx1 k r))
        (time (consistency-tests/3d pts2 t2 xx2 k r))
        ;;      (consistency-tests/2d pts t xx k r)
        
        (if (positive? ntrials)
            (recur (- ntrials 1)))
        ))
    ))
  )


(test1)

