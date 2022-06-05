(defpackage :gauss-quad
  (:use :cl)
  (:export #:gw
           #:abc-to-symm
           ;; Common quadrature schemes
           #:legendre
           ;; For other packages
           *legendre-abc-funs*
           *legendre-muzero*
           #:common-abc))
