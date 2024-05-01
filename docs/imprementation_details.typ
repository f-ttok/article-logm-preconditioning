#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm

#set page(
  paper: "a4", margin: 1.5cm,
)
#set text(font: "Libertinus", size: 12pt)
#show math.equation: set text(font: "Libertinus Math")

#let title = [Imprementation details]
#let author = [Fuminori Tatsuoka]

#align(center, text(1.5em)[*#title*])
#align(center, text(1.2em)[
  #author#footnote("f-tatsuoka@na.nuap.nagoya-u.ac.jp")
])

= Computing $log(A)$ for 
For given $A$, whose maximum (minimum) eigenvalue is $lambda_max$ ($lambda_min$), let we difine
$ tilde(A) = tilde(A) $

$ log(A) = (A - I) integral_(-1)^1 [(1+u)A + (1-u)I]^(-1) dif u $

If we set $P = (A + sqrt(kappa(A))I)$
