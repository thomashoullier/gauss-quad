(in-package :gauss-quad/arbprec)

(defun common-quad-creal (n abc-funs muzero sigbits)
  "common-quad in arbitrary precision with creal numbers.
   sigbits is the number of significant bits for computations."
  (let ((a) (b) (c) (asymm) (bsymm)
        (tj) (w))
    (multiple-value-setq (a b c) (common-abc n abc-funs))
    (multiple-value-setq (asymm bsymm) (abc-to-symm-creal a b c))
    (multiple-value-setq (tj w) (gw-creal asymm bsymm n muzero sigbits))
    (values tj w)))

(defun legendre-creal (n sigbits)
  "Gauss-Legendre quadrature in arbitrary precision.
   sigbits is the number of significant bits for computations
   (eg. would be 53 for double-float)."
  (common-quad-creal n *legendre-abc-funs* *legendre-muzero* sigbits))
