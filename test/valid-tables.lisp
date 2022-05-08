(in-package :gauss-quad/test)

(defun legendre-table (n)
  "Tabulated nodes and weights for Legendre polynomials quadrature
   on n points. Sorted by ascending value of node."
  (let ((tj) (w))
    (multiple-value-setq (tj w)
      (case n
        (1 (values (list 0d0) (list 2d0)))
        (2 (values (list (/ -1d0 (sqrt 3d0)) (/ 1d0 (sqrt 3d0)))
                   (list 1d0 1d0)))
        (3 (values (list (- (sqrt (/ 3d0 5))) 0d0 (sqrt (/ 3d0 5)))
                   (list (/ 5d0 9) (/ 8d0 9) (/ 5d0 9))))
        (4 (values (list (- (sqrt (+ (/ 3d0 7) (* (/ 2 7) (sqrt (/ 6d0 5))))))
                         (- (sqrt (- (/ 3d0 7) (* (/ 2 7) (sqrt (/ 6d0 5))))))
                         (sqrt (- (/ 3d0 7) (* (/ 2 7) (sqrt (/ 6d0 5)))))
                         (sqrt (+ (/ 3d0 7) (* (/ 2 7) (sqrt (/ 6d0 5))))))
                   (list (/ (- 18d0 (sqrt 30d0)) 36)
                         (/ (+ 18d0 (sqrt 30d0)) 36)
                         (/ (+ 18d0 (sqrt 30d0)) 36)
                         (/ (- 18d0 (sqrt 30d0)) 36))))
        (5 (values (list (* (/ -1d0 3) (sqrt (+ 5 (* 2 (sqrt (/ 10d0 7))))))
                         (* (/ -1d0 3) (sqrt (- 5 (* 2 (sqrt (/ 10d0 7))))))
                         0d0
                         (* (/ 1d0 3) (sqrt (- 5 (* 2 (sqrt (/ 10d0 7))))))
                         (* (/ 1d0 3) (sqrt (+ 5 (* 2 (sqrt (/ 10d0 7)))))))
                   (list (/ (- 322 (* 13 (sqrt 70d0))) 900)
                         (/ (+ 322 (* 13 (sqrt 70d0))) 900)
                         (/ 128d0 225)
                         (/ (+ 322 (* 13 (sqrt 70d0))) 900)
                         (/ (- 322 (* 13 (sqrt 70d0))) 900))))
        (otherwise
         (error "legendre-table: n = ~A outside of available range." n))))
    (values (make-array n :element-type 'double-float :initial-contents tj)
            (make-array n :element-type 'double-float :initial-contents w))))

(defun legendre-abc (n)
  "Three-term recurrence coefficients for Legendre polynomials"
  (let ((a (make-array n :element-type 'double-float))
        (b (make-array n :element-type 'double-float))
        (c (make-array n :element-type 'double-float)))
    (loop for i from 0 upto (1- n) do
      (psetf (aref a i) (/ (1+ (* 2d0 i)) (1+ i))
             (aref b i) 0d0
             (aref c i) (/ i (+ 1d0 i))))
    (values a b c)))
