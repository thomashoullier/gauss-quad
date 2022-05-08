;;;; Rove test suite for gauss-quad
(in-package :gauss-quad/test)

(defparameter *valid-eps-legendre* 1d-14)
(defparameter *valid-eps-laguerre* 1d-12)

;; Computation of zero moment for Laguerre weight function requires
;; a gamma function routine.
(defparameter *laguerre-alpha075-muzero* 3.625609908221908d0)
(defparameter *laguerre-alpha05-muzero* 1.772453850905516d0)

(deftest node-weight-tables
  (testing "Legendre"
    (let ((nl (list 1 2 3 4 5))
          (muzero 2d0)
          (a) (b) (c) (asymm) (bsymm)
          (tj) (w) (tjval) (wval) (errtj) (errw))
      (loop for n in nl do
        (multiple-value-setq (a b c) (legendre-abc n))
        (multiple-value-setq (asymm bsymm) (gauss-quad:abc-to-symm a b c))
        (multiple-value-setq (tj w) (gauss-quad:gw asymm bsymm n muzero))
        (multiple-value-setq (tjval wval) (legendre-table n))
        (psetf errtj (max-error-arr tj tjval)
               errw (max-error-arr w wval))
        (ok (and (<= (max-error-arr tj tjval) *valid-eps-legendre*)
                 (<= (max-error-arr w wval) *valid-eps-legendre*))
            (format nil "n = ~A~%tj error: ~A~%w error: ~A" n errtj errw)))))

  (testing "Laguerre alpha=-0.75"
    (let ((n 10)
          (alpha -0.75)
          (muzero *laguerre-alpha075-muzero*)
          (a) (b) (c) (asymm) (bsymm)
          (tj) (w) (tjval) (wval) (errtj) (errw))
      (multiple-value-setq (a b c) (laguerre-abc n alpha))
      (multiple-value-setq (asymm bsymm) (gauss-quad:abc-to-symm a b c))
      (multiple-value-setq (tj w) (gauss-quad:gw asymm bsymm n muzero))
      (multiple-value-setq (tjval wval) (laguerre-alpha075-n10))
      (psetf errtj (max-error-arr tj tjval)
             errw (max-error-arr w wval))
      (ok (and (<= (max-error-arr tj tjval) *valid-eps-laguerre*)
               (<= (max-error-arr w wval) *valid-eps-laguerre*))
          (format nil "n = ~A~%tj error: ~A~%w error: ~A" n errtj errw)))))
