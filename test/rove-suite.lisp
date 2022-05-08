;;;; Rove test suite for gauss-quad
(in-package :gauss-quad/test)

(defparameter *valid-eps* 1d-14)

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
        (ok (and (<= (max-error-arr tj tjval) *valid-eps*)
                 (<= (max-error-arr w wval) *valid-eps*))
            (format nil "n = ~A~%tj error: ~A~%w error: ~A" n errtj errw))))))
