;;
;; K-dimensional point spaces.
;;
;; Copyright 2019 Ivan Raikov
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

(module kspace
	
  (
   
   space? make-space
          dimension
          size
          point
          coord
          compare-coord
          squared-distance
          compare-distance
  )

  (import scheme (chicken base) (chicken foreign)
          yasos yasos-collections
          (only srfi-1 fold)
          (only srfi-4 f32vector? make-f32vector
                f32vector-ref f32vector-set!
                f32vector-length))


  (define-syntax map/iter
    (syntax-rules ()
      ((map/iter index start end body ...)
       (let loop ((index start) (ax '()))
         (if (< index end)
             (let ((v (begin body ...)))
               (loop (add1 index) (cons v ax)))
             (reverse ax))
         ))
      ))

  (define (collection->f32vector c)
    (if (f32vector? c) c
        (let ((res (make-f32vector (size c)))
              (i (make-parameter 0)))
          (begin
            (for-each-elt (lambda (v)
                            (f32vector-set! res (i) v) (i (add1 (i)))) c)
            res))
        ))

  (define-predicate space?)

  ;; dimension of the space
  (define-operation (dimension space))
    
  ;; gets the k'th coordinate of i'th point, starting from 0.
  (define-operation (coord space i k))
    
  ;; gets the i'th point in the space
  (define-operation (point space i))

  ;; compares the given coordinates
  (define-operation (compare-coord space point1 point2 i))

  ;; returns the squared distance between two points.
  (define-operation (squared-distance space point1 point2)) ;; Point * Point -> Double

  ;; returns 0, negative or positive number depending on the
  ;; distance between two points relative to origin
  (define-operation (compare-distance space origin point1 point2))

  (define (make-space coords)

    (let* ((k (length coords))

           (points (map collection->f32vector coords))

           (n (f32vector-length (car points)))

           (point-ref (lambda (p) (map/iter i 0 k (f32vector-ref (list-ref points i) p))))

           (coord-ref (lambda (p i) (f32vector-ref (list-ref points i) p)))
           
           (squared-coord-distance
            (lambda (x y)
              (let ((delta (- x y)))
                (* delta delta))))
           
           (fn-squared-distance
            (lambda (a b)
              (fold (lambda (x y ax)
                      (+ ax (squared-coord-distance x y)))
                    0.0 a b)))
           
           (fn-compare-distance
            (lambda (p a b . reltol)
              (let ((delta (- (fn-squared-distance p a) (fn-squared-distance p b))))
                (if (null? reltol) 
                    delta 
                    (if (<= delta (car reltol)) 0 delta)))))
           
           )
      (object
       ((space? self) #t)
       ((dimension self) k)
       ((size self) n)
       ((coord self p i)
        (coord-ref p i))
       ((point self p)
        (point-ref p))
       ((compare-coord self p1 p2 i)
        (< (coord-ref p1 i) (coord-ref p2 i)))
       ((squared-distance self p1 p2)
        (fn-squared-distance (if (list? p1) p1 (point-ref p1))
                             (if (list? p2) p2 (point-ref p2))))
       ((compare-distance self origin p1 p2)
        (fn-compare-distance origin
                             (if (list? p1) p1 (point-ref p1))
                             (if (list? p2) p2 (point-ref p2))))
       ))
    )

  )

