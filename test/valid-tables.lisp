;;;; Handwritten validation tables for some quadrature rules
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

(defun laguerre-alpha075-n10 ()
  "Nodes and weights of the Gauss-Laguerre quadrature with n=10, alpha=-0.75.
   From GW paper, citing Concus et al."
  (let ((tj) (w))
    (setf tj (list 2.76665586707972d-2
                   4.54784422605949d-1
                   1.38242576115899d0
                   2.833980012092697d0
                   4.850971448764914d0
                   7.500010942642825d0
                   1.0888408023834404d1
                   1.5199478044237603d1
                   2.0789214621070107d1
                   2.8573060164922106d1)
          w (list 2.566765557790772d0
                  7.73347970344341d-1
                  2.33132834973219d-1
                  4.64367470895670d-2
                  5.54912350203625d-3
                  3.65646662677638d-4
                  1.18687985710245d-5
                  1.58441094205678d-7
                  6.19326672679684d-10
                  3.03775992651750d-13)
          tj (make-array 10 :element-type 'double-float :initial-contents tj)
          w (make-array 10 :element-type 'double-float :initial-contents w))
    (values tj w)))

(defun laguerre-abc (n alpha)
  "Compute three-term recurrence coefficients for Laguerre polynomial with
   parameters alpha, up to order n."
  (let ((a (make-array n :element-type 'double-float))
        (b (make-array n :element-type 'double-float))
        (c (make-array n :element-type 'double-float)))
    (loop for i from 0 upto (1- n) do
      (psetf (aref a i) (/ -1d0 (1+ i))
             (aref b i) (/ (+ (* 2d0 i) alpha 1) (1+ i))
             (aref c i) (/ (+ i alpha) (+ i 1d0))))
    (values a b c)))
