;;;; Golub-Welsch algorithm, transposed from the original ALGOL.

(defun gw (a-coeffs b-coeffs n muzero)
  "a, b are the coefficients in the symmetric matrix J.
   a is of size N+1, the index 0 is unused.
   b if of size N, the index 0 is unused.
   n is the order of the quadrature scheme.
   Return tj and w, the nodes and weights of the computed quadradure scheme."
  (let ((w (make-array (1+ n) :element-type 'double-float))  ; Weights
        (tj (make-array (1+ n) :element-type 'double-float)) ; Nodes
        (a (make-array (1+ n) :element-type 'double-float))
        (b (make-array n :element-type 'double-float))
        (eps 0d0)                       ; Relative zero tolerance
        (lambd 0d0) (lambd1 0d0) (lambd2 0d0) (rho 0d0))
    ;; Copy and index shifting by 1 of coefficients a b
    (setf (subseq a 1) a-coeffs (subseq b 1) b-coeffs)
    ;; Find maximum row sum norm: compute eps
    (let ((norm 0d0))
      (setf (aref b 0) 0d0)
      (loop for i from 1 upto (1- n) do
        (setf norm (max norm (+ (abs (aref b (1- i)))
                                (abs (aref a i))
                                (abs (aref b i))))))
      (setf norm (max norm (+ (abs (aref a n)) (abs (aref b (1- n)))))
            eps (* norm double-float-epsilon))
      (psetf lambd norm lambd1 norm lambd2 norm rho norm))
    (setf (aref w 1) 1d0)
    ;; INSPECT: look for convergence of lower diagonal element
    (let ((k 0) (i 0) (m1 0) (m n))
      (loop until (= m 0) do
        (psetf i (1- m) k (1- m) m1 (1- m))
        (block inspect-iter
          (if (>= m1 1)
              (when (<= (abs (aref b m1)) eps)
                (setf (aref tj m) (aref a m)
                      (aref w m) (* muzero (expt (aref w m) 2))
                      rho (if (< lambd1 lambd2) lambd1 lambd2)
                      m m1)
                (return-from inspect-iter))
              (progn
                (setf (aref tj 1) (aref a 1)
                      (aref w 1) (* muzero (expt (aref w 1) 2))
                      m 0)
                (return-from inspect-iter)))
          ;; Small off diagonal element means matrix can be split
          (if (<= (abs (aref b i)) eps)
              (setf k i)
              (setf k 1))
          ;; Find eigenvalues of lower 2x2 and select accelerating shift
          (let ((b2 0d0) (det 0d0) (aa 0d0) (eigmax 0d0))
            (setf b2 (expt (aref b m1) 2)
                  det (sqrt (+ (expt (- (aref a m1) (aref a m)) 2)
                               (* 4 b2)))
                  aa (+ (aref a m1) (aref a m))
                  lambd2 (/ (if (>= aa 0) (+ aa det) (- aa det)) 2)
                  lambd1 (/ (- (* (aref a m1) (aref a m)) b2) lambd2)
                  eigmax (max lambd1 lambd2))
            (if (<= (abs (- eigmax rho)) (/ (abs eigmax) 8))
                (setf lambd eigmax rho eigmax)
                (setf rho eigmax)))
          ;; Transform block from K to M
          (let ((cj 0d0) (r 0d0) (st 0d0) (ct 0d0) (aj 0d0)
                (f 0d0) (q 0d0) (wj 0d0))
            (setf cj (aref b k)
                  (aref b (1- k)) (- (aref a k) lambd))
            (loop for j from k upto m1 do
              (setf r (sqrt (+ (expt cj 2) (expt (aref b (1- j)) 2)))
                    st (/ cj r) ct (/ (aref b (1- j)) r) aj (aref a j))
              (when (> j 1) (setf (aref b (1- j)) r))
              (when (< j (1- n))
                (setf cj (* (aref b (1+ j)) st)
                      (aref b (1+ j)) (* (aref b (1+ j)) (- ct))))
              (setf
               f (+ (* aj ct) (* (aref b j) st))
               q (+ (* (aref b j) ct) (* (aref a (1+ j)) st))
               (aref a j) (+ (* f ct) (* q st))
               (aref b j) (- (* f st) (* q ct))
               wj (aref w j)
               (aref a (1+ j)) (+ aj (aref a (1+ j)) (- (aref a j)))
               (aref w j) (+ (* wj ct) (* (aref w (1+ j)) st))
               (aref w (1+ j)) (- (* wj st) (* (aref w (1+ j)) ct))))))))
    ;; Sort the nodes in ascending order (along with weights)
    (let ((zipped-results (map 'list (lambda (x y) (list x y))
                               (subseq tj 1) (subseq w 1))))
      (setf zipped-results (sort zipped-results #'< :key #'car))
      (setf tj (map 'vector #'car zipped-results)
            w (map 'vector #'cadr zipped-results)))
    (values tj w)))

;;; Manual test with Legendre polynomials
(let* ((muzero 2d0) (n 5)
       (acoeffs (make-array n :element-type 'double-float))
       (bcoeffs (make-array (1- n) :element-type 'double-float))
       (tj) (w))
  (loop for i from 1 upto (1- n) do
    (setf (aref bcoeffs (1- i)) (/ i (sqrt (1- (* 4d0 (expt i 2)))))))
  (format t "~&Legendre polynomials, n = ~A~%" n)
  (multiple-value-setq (tj w) (gw acoeffs bcoeffs n muzero))
  (format t "~&Nodes: ~A~%Weights: ~A~%" tj w))
