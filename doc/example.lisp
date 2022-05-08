(asdf:load-system "gauss-quad")

;;; Gauss-Legendre quadrature example
;; Recurrence coefficients
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

;; Computing the quadrature
(let ((n 5) (muzero 2d0)
      (a) (b) (c) (asymm) (bsymm)
      (tj) (w))
  (multiple-value-setq (a b c) (legendre-abc n))
  (multiple-value-setq (asymm bsymm) (gauss-quad:abc-to-symm a b c))
  (multiple-value-setq (tj w) (gauss-quad:gw asymm bsymm n muzero))
  (format t "~&Nodes: ~A~%Weights: ~A~%" tj w))
