;;;; Test utilities
(in-package :gauss-quad/test)

(defun max-error-arr (arr1 arr2)
  "Compute the maximum error of arr1 - arr2"
  (loop for el1 across arr1 for el2 across arr2
        maximizing (abs (- el1 el2))))
