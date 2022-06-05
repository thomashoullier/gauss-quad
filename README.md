# Gaussian quadrature scheme generation
The system `gauss-quad` is aimed at computing nodes and weights for Gaussian
quadrature schemes, with respect to an arbitrary weight function.

The system `gauss-quad/arbprec` implements the routines with arbitrary precision
computations.

Please refer to [7] for context on Gaussian quadrature.

## Implementation
We use the original Golub-Welsch algorithm [1], since it is general to many
quadratures and allows computing for low orders n.

The original ALGOL60 routines can be found in [2]. They were also
transposed into Matlab by Meurant [5] as part of the book [4]. This
latest version [5] is the most useful for transposing the algorithm into
other languages.

We observed that many modern implementations run the full eigenvectors
computation, instead of the specialized QR algorithm used by Golub
and Welsch [1], which only computes the first component of each eigenvector.
This could be because the algorithm is very succinctly expressed when
leveraging an existing eigenvectors routine. However this is a bit of a shame
since it turns an O(n^2) algorithm into an O(n^3) one.

## Usage
We present the usage of `gauss-quad`. The computations in this system are
performed with `double-float` numbers. The routines for `gauss-quad/arbprec`
are a transposition of the routines for `gauss-quad` with an additional
argument for the number of significant bits.

### Common polynomials
We define a few shortcuts for computing quadrature schemes for common
common polynomials. These functions are detailed in a section below.

```common-lisp
(legendre n)
```

### Low level routines for the general case
Say we want to compute the Gaussian quadrature nodes (roots) `tj` and
weights `w` with respect to some weight function, with `n` quadrature
points.

We need to know two things:
* The coefficients `a`, `b`, `c` of the three-term recurrence for the
  polynomials which correspond to our weight function. Tables of such
  coefficients are found in [6], we use the same conventions. `a`, `b`, and `c`
  are of length `n`.
* The zero moment of the weight function (eg. 2 for Legendre).

For instance, let's say we want to compute the Gauss-Legendre quadrature with
5 points (see [example](doc/example.lisp)).

Start by defining the arrays of recurrence coefficients:

```common-lisp
(defun legendre-abc (n)
  "Three-term recurrence coefficients for Legendre polynomials"
  (let ((a (make-array n :element-type 'double-float))
        (b (make-array n :element-type 'double-float))
        (c (make-array n :element-type 'double-float)))
    (loop for i from 0 upto (1- n) do
      (psetf (aref a i) (/ (1+ (* 2d0 i)) (1+ i))
             (aref b i) 0d0
             (aref c i) (/ i (+ 1d0 i))))
    (values a b c)))
```

The moment zero for the Legendre weight function (1 over [-1, 1]) is 2.
We compute the nodes of weights of the quadrature with:

```common-lisp
(let ((n 5) (muzero 2d0)
      (a) (b) (c) (asymm) (bsymm)
      (tj) (w))
  (multiple-value-setq (a b c) (legendre-abc n))
  (multiple-value-setq (asymm bsymm) (gauss-quad:abc-to-symm a b c))
  (multiple-value-setq (tj w) (gauss-quad:gw asymm bsymm n muzero))
  (format t "~&Nodes: ~A~%Weights: ~A~%" tj w))
```

The nodes are sorted in ascending order, along with corresponding weights.
The values for Gauss-Legendre, n=5, can be checked at [8].

```text
Nodes: #(-0.9061798459386642d0 -0.5384693101056827d0 2.220446049250313d-16
         0.5384693101056829d0 0.906179845938664d0)
Weights: #(0.23692688505618867d0 0.47862867049936786d0 0.5688888888888896d0
           0.4786286704993665d0 0.2369268850561889d0)
```

## Functions
### `gauss-quad`
#### Legendre
**legendre** n => tj w

Compute the Gauss-Legendre quadrature nodes and weights with n points.

#### abc-to-symm
**abc-to-symm** a b c => asymm bsymm

Perform a symmetrization on the three-term recurrence coefficients a, b, c.
a, b, c are vectors of size n.
Returns asymm and bsymm, the components of the symmetric tridiagonal matrix
constructed in [1].

#### gw
**gw** asymm bsymm n muzero => tj w

Compute the quadrature rule defined by asymm and bsymm the symmetrized
recurrence coefficients, muzero the zero moment of the weight function,
n the number of quadrature points.
Returns the vectors of nodes (roots) tj, sorted in ascending order, along
with corresponding weights w.

### `gauss-quad/arbprec`
The argument sigbits is the number of significant bits used to carry out
the computations. Note that this is not the guaranteed accuracy of the
results for the routines below as intermediate approximations are done
for the sake of performance. Setting sigbits to 53 is roughly equivalent to
computing with the `gauss-quad` system in double-float precision.

The nodes and weights computed are given as a `ratio`.

#### Legendre
**legendre-creal** n sigbits => tj w

#### abc-to-symm-creal
**abc-to-symm-creal** a b c => asymm bsymm

In the general case, asymm and bsymm are `creal` numbers.

#### gw-creal
**gw-creal** asymm bsymm n muzero sigbits => tj w

muzero must be of type `creal`.

## Tests
To launch tests, run:

```common-lisp
(asdf:test-system "gauss-quad")
(asdf:test-system "gauss-quad/arbprec")
```

The included validations are:
* `gauss-quad`:
  * Gauss-Legendre (n from 1 to 5)
  * Gauss-Laguerre example from [1] (alpha = -0.75, n=10).
* `gauss-quad/arbprec`:
  * Gauss-Legendre (n from 1 to 5)

## Dependencies
* `gauss-quad`
  * [alexandria](https://github.com/keithj/alexandria)
* `gauss-quad/test`
  * [rove](https://github.com/fukamachi/rove)
* `gauss-quad/arbprec`
  * [computable-reals](https://github.com/stylewarning/computable-reals)

## References
1. G. H. Golub and J. H. Welsch, “Calculation of Gauss quadrature
   rules,” Mathematics of Computation 23, 221–230 (1969).
   doi.org/10.1090/s0025-5718-69-99647-1
2. R. Chan, C. Greif, and D. O’Leary, Milestones in Matrix
   Computation: The Selected Works of Gene H. Golub with Commentaries
   (OUP Oxford, 2007)
3. G. Meurant and A. Sommariva, “Fast variants of the Golub and Welsch
   algorithm for symmetric weight functions in Matlab,” Numerical
   Algorithms 67, 491–506 (2013).
   doi.org/10.1007/s11075-013-9804-x
4. G. H. Golub and G. Meurant, Matrices, Moments and Quadrature with
   Applications (Princeton University Press, 2009).
   doi.org/10.1515/9781400833887
5. [mmq_toolbox](https://github.com/gegemeu/mmq_toolbox)
6. [NIST Digital Library of Mathematical Functions - 18.9](https://dlmf.nist.gov/18.9)
7. [Wikipedia - Gaussian Quadrature](https://en.wikipedia.org/w/index.php?title=Gaussian_quadrature&oldid=1083985268)
8. [Wikipedia - Gauss-Legendre quadrature](https://en.wikipedia.org/w/index.php?title=Gauss–Legendre_quadrature&oldid=1012677738)

## See also
* [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
* [GaussQuadrature.jl](https://github.com/billmclean/GaussQuadrature.jl)
