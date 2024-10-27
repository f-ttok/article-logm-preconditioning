#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#set page(
  paper: "a4", margin: 1.5cm,
)
#set text(font: ("Libertinus Serif", "Noto Serif JP"), size: 12pt)
#show math.equation: set text(font: "Libertinus Math")
#set math.equation(numbering: "(1)")
// #show ref: it => {
//   let eq = math.equation
//   let el = it.element
//   if el != none and el.func() == eq {
//     // Override equation references.
//     link(el.location(),numbering(el.numbering, ..counter(eq).at(el.location())
//     ))
//   } else {
//     // Other references as usual.
//     it
//   }
// }

#let title = [Imprementation details]
#let author = [Fuminori Tatsuoka]

#align(center, text(1.5em)[*#title*])
#align(center, text(1.2em)[
  #author#footnote([email: #link("mailto:f-tatsuoka@na.nuap.nagoya-u.ac.jp", "f-tatsuoka@na.nuap.nagoya-u.ac.jp"), github: #link("https://github.com/f-ttok", "@f-ttok")])
])

#let argmin = math.op("argmin", limits: true)
#let asinh = math.op("asinh", limits: false)
#let atanh = math.op("atanh", limits: false)
#let tildeA = {$tilde(A)$}
#let tildeP1 = {$tilde(P)_1$}
#let cprime = {$c prime$}


= Computing $log(A)$ without preconditioning
We consider the computation of $log(A)$ for HPD matrices using the Gauss--Legendre (GL) quadrature to
$
  log(A) = integral_(-1)^(1) F(t; A) dif t
  #h(2em)
  (F(t; A) = (A-I)[(1-t)I + (1+t)A]^(-1)).
$ <eq.logm.ir>

Initially, we scale the matrix $A$ to $tildeA := c A$ with $c = 1 slash sqrt(lambda_max lambda_min)$ to optimize the convergence rate.
The scaling changes the extreme eigenvalues of $tildeA$ to $sqrt(kappa)$ and $1 slash sqrt(kappa)$, where $kappa$ is the condition number of $A$.

Subsequently, we select the number of abscissas $m$ according to the given error tolerance $epsilon$.
Here, the error of the GL quadrature for $log(tildeA)$ can be estimated as that for scalar logarithms:
$
  norm(log(tildeA) - sum_(k=1)^m F(t_k; tildeA))_2
  = max_(lambda in Lambda)
    abs(log(c lambda) - sum_(k=1)^m F(t_k; c lambda)),
$
where $Lambda$ is the spectrum of $A$.
When $m$ is sufficiently large, the error associated with the extreme eigenvalues will dominate.
Thus, we presume the error of the $m$-point GL quadrature for $log(tildeA)$ as
$
  norm(log(tildeA) - sum_(k=1)^m F(t_k; tildeA))_2
  approx 
    abs(log(sqrt(kappa)) - sum_(k=1)^m F(t_k; sqrt(kappa))).
$
In our implementation, for the error tolerance $epsilon$, we find the minimum $m$ satisfying
$
  abs(log(sqrt(kappa)) - sum_(k=1)^m F(t_k; sqrt(kappa))) <= epsilon.
$
Then we compute $log(A)$ with the $m$-point GL quadrature.


= Computing $log(A)$ with the preconditioning
When we use the preconditioning, $log(A)$ is computed via
$
  log(A) = log(cprime tildeA tildeP1) - log(cprime tildeP1) - log(c)I,
$
where $cprime = sqrt((kappa^(1 slash 2)+1)(kappa^(-1 slash 2) + 1))$.
In practice, it is not necessary to explicitly compute the inverse in $tildeP1$.
For the computation of $log(cprime tildeA tildeP1)$, the integrand can be rewritten as
$
  F(t; cprime tildeA tildeP1)
  & = [cprime tildeA tildeP1 - I][(1-t)I + (1+t)cprime tildeA tildeP1]^(-1)\
  & = [cprime tildeA (tildeA + I)^(-1) - I][(1-t)I + (1+t)cprime tildeA (tildeA + I)^(-1)]^(-1) \
  & = [cprime tildeA - (tildeA + I)][(1-t)(tildeA + I) + (1+t)cprime tildeA]^(-1) \
  & = [(cprime - 1)tildeA - I][((cprime-1)t + 1+cprime)tildeA + (1-t)I]^(-1),
$
and the integrand $F(t; cprime tildeP1)$ can be rewritten as 
$
  F(t; cprime tildeP1)
  & = [cprime tildeP1 - I][(1-t)I + (1+t)cprime tildeP1]^(-1) \
  & = [cprime (tildeA+I)^(-1) - I][(1-t)I + (1+t)cprime (tildeA+I)^(-1)]^(-1) \
  & = [cprime I - (tildeA+I)][(1-t)(tildeA+I) + (1+t)cprime I]^(-1) \
  & = [-tildeA + (cprime-1)I][(1-t)tildeA + ((cprime-1)t + cprime+1)I]^(-1).
$

= Imprementation of the DE formula
The DE formula exploits the trapezoidal rule after applying a specialized change of variables.
Because the convergence analysis in the non-reviewed report @tatsuoka2020convergence is conducted using the variable transformation $tanh(pi sinh(x) slash 2)$, we also use $tanh(pi sinh(x) slash 2)$ while the original paper @tatsuoka2020algorithms uses the change of variables $t = t(x) = tanh(sinh(x))$.

In the experiments, in order to use the different change of variables and bound the absolute error, we modifed the $m$-point DE formula algorithm @tatsuoka2020algorithms[Alg. 1] as follows:

- $F_(upright(D E))(x) = cosh(x) sech^2(sinh(x))) [(1 + tanh(sinh(x)))(A-I) + 2I]^(-1)$\
  $-> quad F_(upright(D E))(x) = (pi cosh(x) sech^2(pi sinh(x) slash 2))/(2) [(1 + tanh(pi sinh(x)slash 2))(A-I) + 2I]^(-1)$
- $theta = |log(rho(A))| quad -> quad theta = 1$
- $l = asinh(atanh(2 a - 1)) quad -> quad l = asinh(2 atanh(2 a - 1) slash pi)$
- $r = asinh(atanh(2 b - 1)) quad -> quad r = asinh(2 atanh(2 b - 1) slash pi)$

#bibliography("reference.bib")