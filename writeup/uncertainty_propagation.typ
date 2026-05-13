#import "@preview/physica:0.9.8": *
#import "utils.typ": *


== Conceptos de estadística en el lenguaje geométrico

Consideremos un punto aleatorio $x in X tilde.equiv RR^n$, cuya función de densidad de probabilidad es $p : X -> RR$, cumpliendo que $p(x) >= 0$ y que #footnote[La integral es n-dimensional, no unidimensional, pero usamos el símbolo de integral simple por simplicidad notacional.]

$
  integral p(x) dd(x) = 1.
$

La primera forma de extraer información sobre $x$ es la media $mu in X$. Esta se define como una suma ponderada por $p(x)$ de gran cantidad de todas las realizaciones de $x$:

$
  mu = EE[x] = integral x p(x) dd(x).
$

Consideremos ahora el vector aleatorio $(x - mu) in T X$. Su valor esperado es nulo:

$
  EE[x - mu] = integral (x - mu) p(x) dd(x) = integral x p(x) dd(x) - mu integral p(x) dd(x) = mu - mu = 0.
$

Consideremos también una 1-forma no nula $alpha in T^* X$, actuando sobre $(x - mu)$. Su valor esperado es también nulo, ya que por linealidad

$
  EE[alpha(x - mu)] = integral alpha(x - mu) p(x) dd(x) = alpha (integral (x - mu) p(x) dd(x)) = alpha(0) = 0
$

Es aparente que es necesario una operación no lineal para obtener datos útiles sobre la distribución de $x$ más allá de la media. En particular, consideremos lo que sucede si calculamos el valor esperado de la anterior 1-forma, pero elevada al cuadrado:

$
  EE[alpha(x - mu)^2] = integral alpha(x - mu)^2 p(x) dd(x) > 0.
$

Generalizando, consideremos dos 1-formas $alpha, beta : T^* X$, y construyamos el producto de su evaluación en $(x - mu)$:

$
  EE[alpha(x - mu) beta(x - mu)] = integral alpha(x - mu) beta(x - mu) p(x) dd(x).
$

Consideremos ahora el anterior objecto como una función de dos 1-forma arbitrarias, $P : T^* X times T^* X -> RR$, tal que

$
  P(alpha, beta) = EE[alpha(x - mu) beta(x - mu)].
$

Este objeto es conocido como el (2,0)-tensor de covarianza, ya que consume un par de 1-formas. #footnote[Este mismo concepto se puede generalizar para mayor cantidad de 1-formas. Considerando 3, obtenemos el (3, 0)-tensor de asimetría, y utilizando 4, el (4, 0)-tensor de curtosis. Por desgracia, para distribuciones de varias variables, estos no se pueden escribir fácilmente sobre el papel; el tensor de asimetría requiere un "cubo" de números, y el de curtosis un hipercubo.] Es fácil demostrar, utilizando la linealidad de las 1-formas y el operador $EE$ que $P$ es bilineal, es decir,

$
  P(lambda alpha_1 + alpha_2, beta) = lambda P(alpha_1, beta) + P(alpha_2, beta) \
  P(alpha, lambda beta_1 + beta_2 ) = lambda P(alpha, beta_1) + P(alpha, beta_2).
$


Podemos obtener la expresión matricial del tensor, $vb(P)$ evaluándolo contra todos los pares de una base de 1-formas, por ejemplo $dd(x_1), dd(x_2), ... dd(x_n)$:

$
  vb(P) = mat(
    P(dd(x_1), dd(x_1)), ..., P(dd(x_1), dd(x_n));
    dots.v, dots.v, dots.v;
    P(dd(x_n), dd(x_1)), ..., P(dd(x_n), dd(x_n))
  )
$

Es interesante que, si escribimos la 1-forma $alpha$ como $alpha_1 dd(x_1) + ...+ alpha_n dd(x_n)$, y $beta$ como $beta_1 dd(x_1) + ... + beta_n dd(x_n)$, podemos expresar

$
  P(alpha, beta) = mat(alpha_1, ..., alpha_n) vb(P) mat(beta_1, ..., beta_n)^TT = vb(alpha) vb(P) vb(beta)^TT,
$

#example[Expresión en coordenadas de $P(alpha, beta)$][
  Para hacer esto evidente, consideremos un ejemplo bidimensional,

  $
    vb(alpha) vb(P) vb(beta)^TT & = mat(alpha_1, alpha_2) med
                                  mat(
                                    P(dd(x_1), dd(x_1)), P(dd(x_1), dd(x_2));
                                    P(dd(x_2), dd(x_1)), P(dd(x_2), dd(x_2))
                                  ) med
                                  mat(beta_1; beta_2) = \
                                & = mat(alpha_1, alpha_2) med
                                  mat(
                                    P(dd(x_1), dd(x_1)) beta_1 + P(dd(x_1), dd(x_2)) beta_2;
                                    P(dd(x_2), dd(x_1)) beta_1 + P(dd(x_2), dd(x_2)) beta_2
                                  )
  $

  Por la bilinealidad de $P$, tenemos que

  $
    P(dd(x_1), dd(x_1)) beta_1 + P(dd(x_1), dd(x_2)) beta_2 =
    P(dd(x_1), beta_1 dd(x_1) + beta_2 dd(x_2)) & = \
                                                & = P(dd(x_1), beta),
  $

  y de forma similar en la segunda fila,

  $
    P(dd(x_2), dd(x_1)) beta_1 + P(dd(x_2), dd(x_2)) beta_2 =
    P(dd(x_2), dd(x_1) beta_1 + dd(x_2) beta_2) & = \
                                                & = P(dd(x_2), beta).
  $

  Multiplicando las dos matrices remanentes

  $
    mat(alpha_1, alpha_2) med
    mat(
      P(dd(x_1), beta);
      P(dd(x_2), beta)
    ) =
    alpha_1 P(dd(x_1), beta) + alpha_2 P(dd(x_2), beta) =
  $

  que, de nuevo, por bilinealidad de $P$,

  $
    = P(alpha_1 dd(x_1) + alpha_2 dd(x_2), beta) = P(alpha, beta).
  $

]


En adelante, utilizaremos esta notación, donde el símbolo sin negrita representa un expresión libre de coordenadas, y el símbolo con negrita un objecto con coordenadas en una base implícita.


=== Mapas afines sobre una distribución

La ventaja de introducir la nomenclatura de la geometría diferencial y evitar asumir coordenadas es la claridad notacional que se puede conseguir. Este apartado busca recordar los pushforward y pullbacks previamente definidos, y esclarecer su utilidad.

Supongamos que tenemos un mapa afín $cal(A): X -> X$, y consideremos $tilde(x) = cal(A)(x)$. La realización con coordenadas más general de este mapa es $L(vb(x)) = vb(A) vb(x) + vb(b)$, siendo $vb(A)$ una matriz y $vb(b)$ un vector.

La media de $tilde(x)$ se obtiene directamente gracias a la afinidad de este, $tilde(mu) = E[cal(A)(x)] = cal(A)(EE[x]) = cal(A)(mu)$. En coordenadas, $vb(tilde(mu)) = vb(A) vb(mu) + vb(b)$.

Por su parte, consideremos un vector $v in T X$, y consideremos su transformación bajo el mapa. Al tratarse de un vector, definimos que este se transforma a través del pushforward $dd(cal(A)) : T X -> T X$, es decir, $tilde(v)= dd(cal(A)) (v)$. Con coordenadas escribiríamos, $vb(tilde(v)) = vb(A) vb(v)$, siempre que escribamos los vectores como vectores columna.

De forma similar, consideremos una 1-forma $alpha in T^* X$, e impongamos que exista una 1-forma $tilde(alpha)$ tal que $tilde(alpha)(tilde(v)) = alpha(v)$. Ya que $tilde(v) = dd(cal(A))(v)$, la anterior igualdad implica que

$
  tilde(alpha) comp dd(cal(A)) = alpha
$

Este mapa que nos lleva de $tilde(alpha)$ a $alpha$ lo denominamos el "pullback" y lo escribimos $dd(cal(A))^*: T^* X -> T^* X$, tal que $dd(cal(A))^*(tilde(alpha)) = alpha$. En coordenadas, buscamos que $vb(tilde(alpha)) tilde(vb(v)) = vb(tilde(alpha)) vb(A) vb(v) = vb(alpha) vb(v)$. Identificando términos, escribimos $vb(tilde(alpha)) vb(A) = vb(alpha)$, siempre que escribamos las 1-formas como vectores fila.

Finalmente, con estas construcciones, estamos listos para transformar el tensor de covarianza. Intuitivamente, definimos que

$
  tilde(P)(tilde(alpha), tilde(beta)) = EE[tilde(alpha)(tilde(x) - tilde(mu)) tilde(beta)(tilde(x) - tilde(mu))],
$

pero estamos interesados en evaluar $tilde(P)$ en función de las 1-formas originales (en esencia, para obtener su nueva forma matricial). Ya que $cal(A)$ es afín, podemos afirmar que $tilde(x) - tilde(mu) = dd(cal(A))(x - mu)$. Por lo tanto, podemos escribir

$
  tilde(P)(tilde(alpha), tilde(beta)) = EE[alpha(x - mu), beta(x - mu)] = P(alpha, beta).
$

A su vez, recordamos que $dd(cal(A))^*(tilde(alpha)) = alpha$, y análogamente para $beta$, lo que nos permite escribir

$
  tilde(P)(tilde(alpha), tilde(beta)) = P(alpha, beta) = P(dd(cal(A))^*(tilde(alpha)), dd(cal(A))^*(tilde(beta))).
$

Equivalentemente a las 1-formas, definimos la actuación del pullback en este como una propiedad distributiva del pullback dentro del tensor, es decir,

$
  tilde(P)(tilde(alpha), tilde(beta)) = dd(cal(A))^* P(alpha, beta) = P(dd(cal(A))^* alpha, dd(cal(A))^* beta).
$

Expresando las anteriores igualdades en coordenadas,

$
  vb(tilde(alpha)) vb(tilde(P)) vb(tilde(beta))^TT
  = vb(alpha) vb(P) vb(beta)^TT
  = (vb(tilde(alpha)) vb(A)) vb(P) (vb(tilde(beta)) vb(A))^TT = vb(tilde(alpha)) ( vb(A) vb(P) vb(A^TT)) vb(tilde(beta))^TT,
$

identificando términos, $vb(tilde(P)) = vb(A) vb(P) vb(A)^TT$, una expresión que usaremos abundantemente en la implementación numérica.

== Métodos basados en muestreo

=== Monte Carlo convencional (MC)

Quizás una de las maneras más conocidas y sencillas para la propagación de incertidumbres es el método de Monte Carlo (MC por sus siglas). Este método consiste en tomar una gran cantidad de puntos de partida, muestreando a cierta distribución, y propagar las ecuaciones del sistema dinámico para cada punto. La posición final de cada punto se puede estudiar estadísticamente.

Si bien el método es prácticamente trivial de implementar, no es óptimo computacionalmente, y la nube final de puntos no es sencilla de estudiar y depende en gran medida del muestreo original.

=== Quasi Monte Carlo (QMC)

== Métodos Gaussianos

#let distN(x, mu, P) = $cal(N)(#x\; med med #mu, #P)$

Consideremos un mapa $Psi: X^n -> X^n$, tal que $tilde(x) = Psi(x)$, y asumamos que $x$ se distribuye según $p(x) = distN(x, mu, P)$. El problema de la propagación de incertidumbres consiste en encontrar una distribución $distN(tilde(x), tilde(mu), tilde(P))$ que aproxime satisfactoriamente la distribución real de $tilde(x)$.

=== Distribución Gaussiana multivariable


=== Cuadratura Gaussiana

=== Distancia de Mahalanobis

=== Matriz de transición de estado (STM)

El método de la matriz de transición de estado consiste en notar que la distribución normal es cerrada en cuanto a mapas afines, es decir, si

$
  cal(A)_b : X -> X, quad cal(A)_b (vb(x)) = tilde(vb(x)) = vb(A) vb(x) + vb(b),
$

y $p(vb(x)) = distN(vb(x), vb(mu), vb(P))$, entonces

$
  p(vb(tilde(x))) = distN(vb(tilde(x)), vb(tilde(mu)), vb(tilde(P))),
$

donde $vb(tilde(mu)) = vb(A) vb(mu) + vb(b)$ y $vb(tilde(P)) = vb(A) vb(P) vb(A)^TT$ como previamente hemos demostrado.





=== Transformación Unscented (UT)


== Método de Gauss-Von Mises

#let cite_gvm(x) = cite(
  <horwoodGaussMisesDistribution2014>,
  supplement: [pág. #x],
)


El método de Gauss-Von Mises es similar a la Transformación Unscented. En este caso, en vez de utilizar una distribución Gaussiana para la variable aleatoria que es el estado del satélite, se utiliza la distribución Gauss-Von Mises, que desarrollaremos a continuación. Usaremos como base teórica #cite(<horwoodGaussMisesDistribution2014>), si bien adaptaremos la nomenclatura y utilizaremos un planteamiento basado en geometría diferencial.

=== Distribución de Von Mises

La distribución de Von Mises queda definida por su función de densidad de probabilidad #cite_gvm(283)

#let distVM(theta, alpha, kappa) = $cal(V M)(#theta\; med med #alpha, #kappa)$

$
  p(theta) = distVM(theta, alpha, kappa) = (e^(kappa cos(theta - alpha))) / (2 pi I_0 (kappa))
$

para una variable aleatoria $theta in RR$. Notamos que la distribución es periódica con periodo $2 pi$, ya que $p(theta) = p(theta + 2 pi)$.

$I_n (x)$ es la función de Bessel modificada de orden $n in NN$, definida según @weissteinModifiedBesselFunction

$
  I_n (kappa) = 1 / pi integral_0^pi e^(kappa cos theta) cos (n theta) dd(theta).
$

Para el caso en el que $n = 0$, tenemos

$
  I_0 (kappa) = 1 / pi integral_0^pi e^(kappa cos theta) dd(theta).
$


A su vez, por simetría del coseno

$
  I_0 (kappa) = 1 / (2 pi) integral_(-pi)^pi e^(kappa cos theta) dd(theta).
$

Notamos que

$
  1 / (2 pi) integral_(-pi)^pi e^(kappa cos (theta - alpha)) dd(theta) = integral_(-pi)^pi e^(kappa cos theta cos alpha + kappa sin theta sin alpha)dd(theta) =
$

Esta última integral es conocida @ListIntegralsExponential2026, y toma por valor

$
  = 2 pi I_0(sqrt(kappa^2 cos^2 alpha + kappa^2 sin^2 alpha)) = 2 pi I_0(kappa).
$

Por lo tanto, podemos comprobar que la distribución de probabilidad integra a la unidad,

$
  integral_(-pi)^(pi) p(theta) dd(theta) = (2 pi I_0 (kappa)) / (2 pi I_0 (kappa)) = 1.
$

Esto, junto con $p(theta) >= 0 med thick forall med theta$ y la periodicidad, convierten a la distribución de Von Mises en una distribución definida para el círculo #cite_gvm(283).

=== Distribución Gauss-Von Mises canónica

Consideremos las variables aleatorias $(z, phi)$ donde $z in Z subset.eq RR^n$, $phi in S^1$. Denominaremos a ambas variables como distribuidas conjuntamente bajo una distribución Gauss-Von Mises canónica sí cumplen que

$
  p(z, phi) = distN(z, 0, I) med distVM(phi, 0, kappa),
$

donde $distN(x, 0, I)$ es una distribución normal canónica y $distVM(theta, 0, kappa)$ es una distribución Von Mises con media angular nula y "afilamiento" $kappa$.

Notamos que, bajo esta distribución canónica, ya que no existe acoplamiento entre $z$ y $phi$, el marginal en $z$ de la distribución es la propia distribución normal, y el marginal en $phi$ es la propia distribución de Von Mises.

=== Distribución Gauss-Von Mises general

Consideremos el mapa afín $cal(C)_mu : Z -> X subset.eq RR^n$, tal que $cal(C)_mu (z) = x = C z + mu$, donde $C$ es una matriz triangular inferior. Observamos que, si $z$ se distribuye según $distN(z, 0, I)$, $x$ a su vez se distribuye siguiendo $distN(x, mu, P)$, con $P = C C^TT$. El mapa inverso es $z = C^(-1) (x - mu)$ y lo denominaremos $cal(C)_mu^(-1): X -> Z$.

Por otra parte, consideremos la familia de funciones $cal(Theta)_(alpha, beta, Gamma) : Z -> RR$ generadas por

$
  Theta_(alpha, beta, Gamma)(z) = alpha + beta(z) + Gamma(z, z),
$

donde $alpha in RR$, $beta: Z -> RR$ es una 1-forma y $Gamma: Z -> RR$ es una forma cuadrática, y construyamos la siguiente distribución

#let distGVM(x, theta, mu, p, alpha, beta, Gamma, kappa) = {
  $cal(G V M)(#x, #theta\; med med #mu, #p, #alpha, #beta, #Gamma, #kappa)$
}

$
  p(x, theta) = distGVM(x, theta, mu, P, alpha, beta, Gamma, kappa) = distN(x, mu, P) med distVM(theta, Theta_(alpha, beta, Gamma)(cal(C)_mu^(-1) (x)), kappa).
$

En este caso, y al contrario de la distribución canónica, hemos logrado unir íntimamente la variable $x$ con $theta$ (siempre que $beta$ y $Gamma$ no sean ambos trivialmente nulos).

La distribución marginal en $x$ se mantiene normal,

$
  p_x (x) = integral_(-pi)^(pi) distN(x, mu, P) e^(kappa cos(theta - Theta^*(x))) / (2 pi I_0 (kappa)) dd(theta) = distN(x, mu, P)
$

sin importar $Theta^*(x) = Theta(A^(-1) (x - mu))$, que es constante en cuanto a la integración. Pero, hemos perdido la simplicidad del marginal en $theta$, ya que

$
  p_theta (theta) = integral_(-infinity)^(infinity) dots integral_(-infinity)^(infinity) distN(x, mu, P) distVM(theta, Theta(x), kappa) dd(x_1) dots dd(x_n)
$

no tiene una expresión analítica (o por lo menos, una expresión razonable de desarrollar), y se trata de una mezcla suave, ponderada por distribuciones normales, de una infinidad de distribuciones Von Mises con diferentes medias angulares.

Notamos que, por como hemos construido esta distribución, siempre es posible convertirla de vuelta a su forma canónica $(z, phi)$ mediante los mapas $cal(C)_mu^(-1)$ y $Theta_(alpha, beta, Gamma)^(-1)$.

=== Técnica de propagación de incertidumbre

Consideremos un mapa $Psi : X^n times S^1 -> X^n times S^1$, y denominemos

$
  (tilde(x), tilde(theta)) = Psi(x, theta),
$

donde $(x, theta)$ están conjuntamente distribuidos siguiendo Gauss-Von Mises.

Ya que el mapa $Psi$ es arbitrario, $(tilde(x), tilde(theta))$, en general, tan solo estarán aproximadamente distribuidos según Gauss Von-Mises. Nuestro objetivo es entonces encontrar una distribución Gauss Von-Mises justificadamente correcta.

==== Distancia de Mahalanobis-Von Mises

==== Cuadratura de Gauss-Von Mises

==== Algoritmo de propagación


== Método por expansión en caós polinómico
