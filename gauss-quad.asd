(defsystem gauss-quad
  :name "gauss-quad"
  :author "Thomas HOULLIER"
  :depends-on ("alexandria")
  :components
  ((:module "src"
    :components ((:file "package")
                 (:file "gw" :depends-on ("package"))
                 (:file "short" :depends-on ("gw")))))
  :in-order-to ((test-op (test-op "gauss-quad/test"))))

(defsystem gauss-quad/test
  :name "gauss-quad/test"
  :depends-on ("rove" "gauss-quad")
  :components
  ((:module "test"
    :components ((:file "package")
                 (:file "utils" :depends-on ("package"))
                 (:file "valid-tables" :depends-on ("package"))
                 (:file "rove-suite" :depends-on ("valid-tables" "utils")))))
  :perform (test-op (o c) (symbol-call :rove '#:run c)))

;;; Arbitrary precision
(defsystem gauss-quad/arbprec
  :name "gauss-quad/arbprec"
  :author "Thomas HOULLIER"
  :depends-on ("gauss-quad" "computable-reals")
  :components
  ((:module "src/arbprec"
    :components ((:file "package")
                 (:file "creal-utils" :depends-on ("package"))
                 (:file "gw-arbprec" :depends-on ("creal-utils"))
                 (:file "short" :depends-on ("gw-arbprec")))))
  :in-order-to ((test-op (test-op "gauss-quad/arbprec/test"))))

(defsystem gauss-quad/arbprec/test
  :name "gauss-quad/arbprec/test"
  :author "Thomas HOULLIER"
  :depends-on ("gauss-quad/arbprec" "gauss-quad/test")
  :components
  ((:module "test/arbprec"
    :components ((:file "package")
                 (:file "rove-suite" :depends-on ("package")))))
  :perform (test-op (o c) (symbol-call :rove '#:run c)))
