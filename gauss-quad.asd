(defsystem gauss-quad
  :name "gauss-quad"
  :author "Thomas HOULLIER"
  :depends-on ("alexandria")
  :components
  ((:module "src"
    :components ((:file "package")
                 (:file "gw" :depends-on ("package")))))
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
