;;;; Shortcuts for common quadrature schemes
(in-package :gauss-quad)

;;; Common helpers
(defun common-abc (n abc-funs)
  "Generate tables of recurrence relation coefficients of size n.
   The functions provided take only i as a parameter from 0 upto n-1."
  (let ((a (make-array n))
        (b (make-array n))
        (c (make-array n)))
    (loop for i from 0 upto (1- n) do
      (loop for coefs in (list a b c)
            for coef-fun in abc-funs do
              (setf (aref coefs i) (funcall coef-fun i))))
    (values a b c)))

(defun common-quad (n abc-funs muzero)
  "Compute a quadrature scheme defined by the polynomials' recurrence relation
   coefficients and muzero."
  (let ((a) (b) (c) (asymm) (bsymm)
        (tj) (w))
    (multiple-value-setq (a b c) (common-abc n abc-funs))
    (multiple-value-setq (asymm bsymm) (abc-to-symm a b c))
    (multiple-value-setq (tj w) (gw asymm bsymm n muzero))
    (values tj w)))

;;; Legendre
(defvar *legendre-muzero* 2)

;; Recurrence relation coefficients
(defun legendre-an (n) (/ (1+ (* 2 n)) (1+ n)))
(defun legendre-bn (n) (declare (ignore n)) 0)
(defun legendre-cn (n) (/ n (+ 1 n)))
(defvar *legendre-abc-funs* (list #'legendre-an #'legendre-bn #'legendre-cn))

(defun legendre (n)
  "Compute the Legendre polynomials quadrature scheme with n points."
  (common-quad n *legendre-abc-funs* *legendre-muzero*))
