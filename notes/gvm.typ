#import "@preview/bananote:0.1.2": *
#import "@preview/physica:0.9.8": *

#show: note.with(
  title: [Notas sobre Gauss Von Mises],
)

= Distribución Gauss Von Mises (GVM)

Consideremos $n$ variables euclídeas, distribuidas conjuntamente en una normal multivariable, y una variable angular $L$. La distribución GVM es el producto de la normal multivariable y una distribución de Von Mises. La distribución normal tiene función de densidad de probabilidad @horwoodGaussMisesDistribution2014

$
  p^N (va(x)) = 1 / sqrt(det(2 pi P)) exp[-1/2 (va(x) - va(mu))^TT P^(-1) (va(x) - va(mu))],
$

mientras que la Von Mises es @horwoodGaussMisesDistribution2014

$
  p^V (va(x), theta) = 1 / (2 pi exp(-kappa) I_0 (kappa)) exp[kappa cos(theta - alpha - va(beta)^TT A^(-1) (va(x) - va(mu)))] ,
$

dónde $A$ es una matriz transformación de blanqueamiento, típicamente la descomposición triangular inferior de $P$, tal que $P = A A^TT$ @horwoodGaussMisesDistribution2014.

La distribución GVM entonces es

$
  p^"GVM" (va(x), theta) = p^N (va(x)) med p^V (va(x), theta)
$

== Caso con $kappa -> infinity$

Cuando $kappa >> 1$, la distribución de Von Mises tiende a una normal @horwoodGaussMisesDistribution2014

$
  lim_(kappa -> infinity) p^V (va(x), theta) =
  1 / sqrt((2pi) / kappa) exp[-1/2 kappa (theta - alpha - beta^TT A^(-1) (va(x) - va(mu)))^2],
$

que equivale a una distribución normal con media $alpha - beta^TT A^(-1) (va(x) - va(mu))$ y varianza $1 / kappa$.

Por lo tanto, si consideramos toda la distribución de Gauss Von Mises como una normal multivariable, tendríamos

$
  p(va(y)) = 1 / sqrt(det(2 pi Q)) exp[-1/2 (va(y) - va(nu))^TT Q^(-1) (va(x) - va(nu))]
$

dónde @horwoodGaussMisesDistribution2014

$
  va(y) = mat(va(x); alpha) \
  Q = mat(A A^TT, A beta; beta^TT A^T, beta^TT beta + 1 / kappa). \
$

Por lo tanto, sí $kappa$ es grande, la distribución Von Mises es equivalente a considerar la variable angular como otra variable Gaussiana, así como su correlación con el resto de variables y su propia varianza.


= La distribución marginal en el ángulo

Investigaremos en este apartado la distribución marginal de una Gauss Von Mises para el ángulo.

Para ello, utilizamos la función característica para la distribución (con $Gamma$ nulo), que es @horwoodGaussMisesDistribution2014

$
  phi^"GVM" (va(xi), m) = (I_(abs(m))(kappa)) / (I_0(kappa)) exp[i (va(mu)^TT va(xi) + m alpha) - 1/2 (A^TT va(xi) + m beta)^TT (A^TT va(xi) + m va(beta))],
$

dónde $va(xi) in RR^n$ y $m in ZZ$. Recordando la definición de la función característica,

$
  phi^"GVM" (va(xi), m) = EE[e^i(va(xi)^TT va(x) + m theta)]
$

se hace evidente porqué $m$ queda limitado a $ZZ$, esto se debe a que, para $m$ no entero, $e^(i m theta)$ no sería cíclico en $theta$, y la expresión carecería de sentido @mardiaDirectionalStatistics1999.

La distribución marginal a partir de la función característica se obtiene poniendo el parámetro sobre el que marginalizamos a $0$. Por lo tanto, para $va(xi) = 0$, tenemos

$
  phi^"GVM" (0, m) = (I_abs(m) (kappa)) / (I_0 (kappa)) exp[
    i m alpha - m^2/2 va(beta)^TT va(beta)
  ].
$

Podemos entender esta expresión como el producto de una distribución Von Mises (con media nula, y valor $kappa$) y una distribución normal "enrollada" con media $alpha$ y varianza $sigma^2 = va(beta)^TT va(beta)$:

$
  phi^"GVM" (0, m) = underbrace((I_abs(m) (kappa)) / (I_0 (kappa)), "Von Mises") underbrace(exp[i m alpha - m^2 / 2 sigma^2], "Gaussiana en círculo").
$

== Aproximación de la función característica por otra Von Mises

La anteriormente presentada función característica marginal no es una Von Mises. Ya que la función característica se puede convertir en la función distribución de probabilidad mediante la transformada inversa de Fourier, podemos entender, mediante el teorema de la convolución, este marginal como la convolución de una Von Mises y una distribución normal (enrollada en el círculo, debido al argumento $m$ siendo entero en vez de un número real).

Para aproximar esta distribución con una Von Mises, y notando que la distribución Von Mises tiene dos parámetros, parece razonable utilizar la dirección media $alpha$ y la longitud media resultante $rho$, definida para una distribución en el círculo @mardiaDirectionalStatistics1999

$
  phi^"GVM" (0, 1) = rho exp[i alpha].
$

Si igualamos este "momento" del marginal con el de una Von Mises (con parámetros $hat(alpha)$ y $hat(kappa)$), obtenemos

$
  (I_1 (kappa)) / (I_0 (kappa)) exp [i alpha - 1 / 2 sigma^2] =
  (I_1 (hat(kappa))) / (I_0 (hat(kappa))) exp [i hat(alpha)].
$

Ya que la ecuación es una igualdad de numeros complejos, podemos igualar su norma y argumento por separado.

El argumento es inmediato,

$
  exp[i alpha] = exp[i hat(alpha)] => alpha = hat(alpha),
$

pero la norma requiere resolver

$
  (I_1 (kappa)) / (I_0 (kappa)) exp[-1/2 sigma^2] = (I_1 (hat(kappa))) / (I_0 (hat(kappa))).
$

Es posible que esta ecuación tenga una solución analítica, pero consideraremos su solución numérica. Si escribimos $A(kappa) = (I_1 (kappa)) / (I_0 (kappa))$, la ecuación equivale a

$
  A(hat(kappa)) - A(kappa) exp[-1 / 2 sigma^2] = 0,
$

que afortunadamente tiene derivada analítica en cuanto a $hat(kappa)$ (calculada con Wolfram Mathematica)

$
  dv(A, kappa) = (I_0 (kappa) + I_2 (kappa)) / (2 I_0 (kappa)) -A(kappa)^2,
$

lo que permite escribir una expresión iterativa para Newton-Rhapson.

= Planteamiento de un nuevo método de propagación de incertidumbre GVM

Con lo anteriormente desarollado, se plantea el siguiente esquema de propagación de incertidumbre:

- Propagar

= Efecto de $Gamma != 0$

Tomando de @horwoodGaussMisesDistribution2014 la función característica completa con $Gamma != 0$ es algo más compleja





#bibliography("../writeup/refs.bib")
