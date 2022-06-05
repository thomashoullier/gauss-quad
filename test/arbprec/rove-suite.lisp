;;;; Rove test suite for arbitrary precision gauss-quad.
(in-package :gauss-quad/arbprec/test)

(defvar *valid-eps-legendre-creal* 1d-14)
(defvar *valid-legendre-sigbits* 53)

(deftest short
  (testing "Legendre vs table"
    (let ((nl (list 1 2 3 4 5))
          (tj) (w) (tjval) (wval) (errtj) (errw))
      (loop for n in nl do
        (multiple-value-setq (tj w) (legendre-creal n *valid-legendre-sigbits*))
        (multiple-value-setq (tjval wval) (legendre-table n))
        (psetf errtj (max-error-arr tj tjval)
               errw (max-error-arr w wval))
        (ok (and (<= (max-error-arr tj tjval) *valid-eps-legendre-creal*)
                 (<= (max-error-arr w wval) *valid-eps-legendre-creal*))
            (format nil "n = ~A~%tj error: ~A~%w error: ~A" n errtj errw))))))
