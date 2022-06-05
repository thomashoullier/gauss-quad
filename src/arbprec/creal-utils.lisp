(in-package :gauss-quad/arbprec)

(defun vec-creal-to-ratio (seq sigbits)
  "Convert a vector of creal to a vector of ratio with sigbits
  significant bits."
  (map 'vector (lambda (x) (rational-approx-r x sigbits)) seq))

(defun <-r (sigbits num1 num2)
  "< comparison for creal approximated to sigbits significant bits."
  (< (rational-approx-r num1 sigbits)
     (rational-approx-r num2 sigbits)))

(defun max-r (sigbits number &rest more-numbers)
  "max function for creal. Numbers are first converted to
   ratio with sigbits significant bits and a max is done."
  (let ((numlist (cons number more-numbers))
        (numbers)
        (imax 0) (num-max 0))
    (setf numbers (vec-creal-to-ratio numlist sigbits)
          num-max (aref numbers 0))
    (loop for i from 0
          for el across numbers do
            (when (> el num-max) (setf imax i num-max el)))
    (elt numlist imax)))

(defun min-r (sigbits number &rest more-numbers)
  "min function for creal. Numbers are first converted to
   ratio with sigbits significant bits and a min is done."
  (let ((numlist (cons number more-numbers))
        (numbers)
        (imin 0) (num-min 0))
    (setf numbers (vec-creal-to-ratio numlist sigbits)
          num-min (aref numbers 0))
    (loop for i from 0
          for el across numbers do
            (when (< el num-min) (setf imin i num-min el)))
    (elt numlist imin)))

(defun abs-r (sigbits number)
  "abs function for creal. Number is converted to ratio with sigbits
   in order to determine sign."
  (let ((num (rational-approx-r number sigbits)))
    (if (< num 0) (-r number) number)))

(defun sqr-r (number)
  "Square a creal."
  (*r number number))
