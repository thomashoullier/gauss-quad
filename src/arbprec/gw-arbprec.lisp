;;;; Golub-Welsch algorithm in arbitrary precision.
(in-package :gauss-quad/arbprec)

(defun coerce-vec-ratio (vec)
  "Coerce a vector into a vector with elements ratio."
  (map 'vector #'rational vec))

(defun abc-to-symm-creal (acoef bcoef ccoef)
  "Return symmetrized coefficients a, b from the three-term recurrence
   coefficients arrays of size n.
   Pn(x) = (An.x + Bn).Pn-1(x) - Cn.Pn-2(x).
   creal computations."
  (let ((a (coerce-vec-ratio acoef))
        (b (coerce-vec-ratio bcoef))
        (n (length acoef)))
    (if (equalp acoef ccoef)
        ;; Already symmetric
        (values a (subseq b 0 (1- n)))
        ;; Symmetrize
        (let ((ai 0) (c (coerce-vec-ratio ccoef)))
          (loop for i from 0 below (1- n) do
            (setf ai (aref a i)
                  (aref a i) (/ (- (aref b i)) ai)
                  (aref b i)
                  (sqrt-r (/ (aref c (1+ i)) (* ai (aref a (1+ i)))))))
          (setf (aref a (1- n)) (/r (-r (aref b (1- n))) (aref a (1- n))))
          (values a (subseq b 0 (1- n)))))))

(defun gw-creal (a-coeffs b-coeffs n muzero sigbits)
  "a, b are the coefficients in the symmetric matrix J.
   a is of size N
   b if of size N-1
   n is the order of the quadrature scheme.
   sigbits is the number of significant bits used to compute
   the epsilon in comparisons (eg. 53 for double-float).
   muzero must be a creal.
   Return tj and w, the nodes and weights of the computed quadradure scheme.
   creal computations."
  (let ((w (make-array (1+ n) :element-type 'creal))  ; Weights
        (tj (make-array (1+ n) :element-type 'creal)) ; Nodes
        (a (make-array (1+ n) :element-type 'creal))
        (b (make-array n :element-type 'creal))
        (eps 0)                       ; Relative zero tolerance
        (numeps (/ 1 (expt 2 sigbits)))
        (lambd 0) (lambd1 0) (lambd2 0) (rho 0))
    ;; Copy and index shifting by 1 of coefficients a b
    (psetf (subseq a 1) a-coeffs
           (subseq b 1) b-coeffs)
    ;; Find maximum row sum norm: compute eps
    (let ((norm 0))
      (setf (aref b 0) 0)
      (loop for i from 1 upto (1- n) do
        (setf norm (max-r sigbits norm (+r (abs-r sigbits (aref b (1- i)))
                                           (abs-r sigbits (aref a i))
                                           (abs-r sigbits (aref b i))))))
      (setf norm (max-r sigbits norm (+r (abs-r sigbits (aref a n))
                                         (abs-r sigbits (aref b (1- n)))))
            eps (* (rational-approx-r norm sigbits) numeps))
      (psetf lambd norm lambd1 norm lambd2 norm rho norm))
    (setf (aref w 1) 1)
    ;; INSPECT: look for convergence of lower diagonal element
    (let ((k 0) (i 0) (m1 0) (m n))
      (loop until (= m 0) do
        (psetf i (1- m) k (1- m) m1 (1- m))
        (block inspect-iter
          (if (>= m1 1)
              (when (<= (abs (rational-approx-r (aref b m1) sigbits)) eps)
                (setf (aref tj m) (aref a m)
                      (aref w m) (*r muzero (sqr-r (aref w m)))
                      rho (min-r sigbits lambd1 lambd2)
                      m m1)
                (return-from inspect-iter))
              (progn
                (setf (aref tj 1) (aref a 1)
                      (aref w 1) (*r muzero (sqr-r (aref w 1)))
                      m 0)
                (return-from inspect-iter)))
          ;; Small off diagonal element means matrix can be split
          (if (<= (abs (rational-approx-r (aref b i) sigbits)) eps)
              (setf k i)
              (setf k 1))
          ;; Find eigenvalues of lower 2x2 and select accelerating shift
          (let ((b2 0) (det 0) (aa 0) (eigmax 0))
            (setf b2 (sqr-r (aref b m1))
                  det (sqrt-r (+r (sqr-r (-r (aref a m1) (aref a m)))
                                  (*r 4 b2)))
                  aa (+r (aref a m1) (aref a m))
                  lambd2 (/r (if (>= (rational-approx-r aa sigbits) 0)
                                 (+r aa det) (-r aa det)) 2)
                  lambd1 (/r (-r (*r (aref a m1) (aref a m)) b2) lambd2)
                  eigmax (max-r sigbits lambd1 lambd2))
            (if (<= (abs (rational-approx-r (-r eigmax rho) sigbits))
                    (/ (abs (rational-approx-r eigmax sigbits)) 8))
                (setf lambd eigmax rho eigmax)
                (setf rho eigmax)))
          ;; Transform block from K to M
          (let ((cj 0) (r 0) (st 0) (ct 0) (aj 0) (f 0) (q 0) (wj 0))
            (setf cj (aref b k)
                  (aref b (1- k)) (-r (aref a k) lambd))
            (loop for j from k upto m1 do
              (setf r (sqrt-r (+r (sqr-r cj) (sqr-r (aref b (1- j)))))
                    st (/r cj r) ct (/r (aref b (1- j)) r) aj (aref a j))
              ;; Intermediate approximation for performance.
              (setf r (rational-approx-r r sigbits)
                    st (rational-approx-r st sigbits)
                    ct (rational-approx-r ct sigbits))
              (when (> j 1) (setf (aref b (1- j)) r))
              (when (< j (1- n))
                (setf cj (*r (aref b (1+ j)) st)
                      (aref b (1+ j)) (*r (aref b (1+ j)) (-r ct))))
              (setf
               f (+r (*r aj ct) (*r (aref b j) st))
               q (+r (*r (aref b j) ct) (*r (aref a (1+ j)) st))
               (aref a j) (+r (*r f ct) (*r q st))
               (aref b j) (-r (*r f st) (*r q ct))
               wj (aref w j)
               (aref a (1+ j)) (+r aj (aref a (1+ j)) (-r (aref a j)))
               (aref w j) (+r (*r wj ct) (*r (aref w (1+ j)) st))
               (aref w (1+ j)) (-r (*r wj st) (*r (aref w (1+ j)) ct))))))))
    ;; Converting nodes and weights to ratio with sigbits precision.
    (psetf tj (map 'vector (lambda (x) (rational-approx-r x sigbits)) tj)
           w (map 'vector (lambda (x) (rational-approx-r x sigbits)) w))
    ;; Sort the nodes in ascending order (along with weights)
    (let ((zipped-results (map 'list (lambda (x y) (list x y))
                               (subseq tj 1) (subseq w 1))))
      (setf zipped-results (sort zipped-results #'< :key #'car))
      (setf tj (map 'vector #'car zipped-results)
            w (map 'vector #'cadr zipped-results)))
    (values tj w)))
