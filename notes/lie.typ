#import "@preview/physica:0.9.8": *

#import "@preview/ilm:1.4.2": *

#set text(lang: "es")
// #set text(font: "Noto Sans")
// #show math.equation: set text(font: "Fira Math", fallback: true)
#set quote(block: true)

#show: ilm.with(
  title: [Notas sobre propagación de incertidumbre],
  author: "José Antonio Mayo García",
  bibliography: bibliography("../writeup/refs.bib"),
  figure-index: (enabled: true),
  table-index: (enabled: false),
  listing-index: (enabled: false),
  abstract: [
    Se parte del paper @markovicWrappingKalmanFilter2017 y se particulariza a las necesidades de la propagación orbital, considerando un grupo de Lie "cilíndrico" $"SO"(2) times RR^5$ como el espacio en el que realizar la propagación orbital. Posteriormente, se considera el vector de estado como una variable aleatoria con media sobre este espacio cilíndrico, y incertidumbre sobre el álgebra de Lie siguiendo el tratamiento estadístico detallado en @bourmaudContinuousDiscreteExtendedKalman2015. Se deducen las características de la propagación bajo la hipótesis de normal concentrada, y se introducen varias distancias de Mahalanobis apropiadas para el problema.],
  paper-size: "a4",
)

#set math.equation(numbering: "(1)")


#outline()

#pagebreak()
= Notación

Se utilizará la misma notación que en @markovicWrappingKalmanFilter2017.

Denominaremos con $G$ a un grupo de Lie. Su álgebra de Lie asociada es $frak(g)$, y tenemos el mapa exponencial y logarítmico

$
  exp_G: frak(g) -> G \
  log_G: G -> frak(g).
$

El álgebra de Lie es el espacio tangente al grupo de Lie en su origen. Como espacio tangente, es euclídeo. Denominaremos a la conversión de elementos del álgebra a "vectores convencionales" y su inversa

$
  [.]_G^or: frak(g) -> RR^p \
  [.]_G^and: RR^p -> frak(g).
$

== Particularización e intuición para $"SO"(2)$

$G = "SO"(2)$ es isomórfico al círculo, es decir, un elemento de $"SO"(2)$ se puede entender como un punto del círculo. Equivalentemente, podemos imaginar los elementos de $"SO"(2)$ como las matrices de rotación que nos llevan de un vector identidad a todos los puntos del círculo. Por lo tanto, podemos expresar los elementos $X_k = R_(theta_k) in "SO"(2)$ como

$
  X_k = mat(cos(theta_k), -sin(theta_k); sin(theta_k), cos(theta_k)).
$

El grupo $"SO"(2)$ tiene operación multiplicativa mediante la multiplicación de estas matrices, que es una operación suave, y cerrada en el grupo (el producto de varias matrices de rotación es otra matriz de rotación).

A su vez, las matrices de rotación son siempre invertibles, ya que $X_k^TT = X_k^(-1)$ para una matriz de rotación. Esta operación es también suave, lo que nos confirma que las matrices de rotación efectivamente forman un grupo de Lie.

== El álgebra de Lie para $"SO"(2)$

Como se ha introducido previamente, el álgebra de Lie es el espacio tangente al grupo en su origen. Podemos calcular un vector tangente cualquiera al grupo mediante la derivada,

$
  lim_(h -> 0) (R_(theta + h) - R_(theta)) / (h) =
  mat(-sin(theta), -cos(theta); cos(theta), -sin(theta)).
$

El álgebra de Lie es el espacio tangente en el origen. Sus elementos son por tanto proporcionales a la anterior derivada evaluada en $theta = 0$,

$
  [theta_k]_G^and = theta_k mat(0, -1; 1, 0) = mat(0, -theta_k; theta_k, 0).
$

Notamos que $theta_k in RR$, y se aprecia como el operador $[.]_G^and$ nos lleva del espacio isomórfico a $frak("so")(2)$, los reales, al propio algebra de Lie.

El mapa inverso es por tanto,

$
  [mat(0, -theta_k; theta_k, 0)]_G^or = theta_k.
$

== El mapa exponencial y logarítmico

Bajo esta notación, el mapa exponencial y logarítmico son simplemente

$
  exp_G [[theta_k]_G^and] = R_(theta_k), \
  log_G [R_(theta_k)] = [theta_k]_G^and.
$

El mapa exponencial es, literalmente, el exponencial matricial

$
  exp_G [[theta_k]_G^and] = exp[mat(0, -theta_k; theta_k, 0)] =
  sum_(n = 0)^infinity mat(0, -theta_k; theta_k, 0)^n / (n!) = \
  mat(1, 0; 0, 1)
  + mat(0, -theta_k; theta_k, 0)
  + 1 / 2 mat(-theta_k^2, 0; 0, -theta_k^2)
  + 1 / 6 mat(0, theta_k^3; -theta_k^3, 0)
  + ...
$

que identificamos con el desarrollo en serie del seno y coseno en cada posición de la matriz, para obtener

$
  exp_G [[theta_k]_G^and] = mat(cos(theta_k), -sin(theta_k); sin(theta_k), cos(theta_k)) = R_(theta_k)
$

El logaritmo es más complicado de desarrollar, ya que presenta varias ramas. Para el desarrollo que realizaremos, no es necesaria intuición sobre el logaritmo matricial, así que no entraremos en detalle.


#pagebreak()
= Grupo de Lie cilíndrico $"SO"(2) times RR^n$

En @markovicWrappingKalmanFilter2017 se considera que el grupo cilíndrico incluye tan solo dos coordenadas cartesianas. Estas son la velocidad y aceleración angular, por lo tanto íntimamente relacionadas con la coordenada angular. Para la propagación orbital, esta relación no existe, las coordenadas cartesianas no están tan claramente relacionadas con el parámetro angular.

== Grupo de Lie en expresión matricial

De todas formas, podemos expandir $"SO"(2)$ con las coordenadas cartesianas $mat(a_k, b_k, ...) in RR^n$ si consideramos la matriz bloque

$
  X_k = mat(
    R_(theta_k), , , ;
    , mat(1, a_k; 0, 1), , ;
    , , mat(1, b_k; 0, 1), ;
    , , , dots.down
  ).
$

Entonces, por las propiedades de las matrices bloque diagonales,

$
  X_1 X_2 =
  mat(
    R_(theta_1) R_(theta_2), , , ;
    , mat(1, a_1 + a_2; 0, 1), , ;
    , , mat(1, b_1 + b_2; 0, 1), ;
    , , , dots.down
  ), \
  X_k^(-1) =
  mat(
    R_(theta_k)^(-1), , , ;
    , mat(1, -a_k; 0, 1), , ;
    , , mat(1, -b_k; 0, 1), ;
    , , , dots.down
  )
$

Podemos apreciar que estas matrices forman un grupo de Lie ya que su operación producto y inversa son suaves y cerradas.

== Algebra de Lie en expresión matricial

La construcción del algebra de Lie es muy similar a la que realizamos previamente. Consideremos la derivada de la matriz respecto a una perturbación en $theta_k$, $a_k$, $b_k$, etc... Abusando notación, escribimos

$
  delta X_k = mat(
    mat(-sin(theta_k), -cos(theta_k); cos(theta_k), sin(theta_k)), , , ;
    , mat(0, 1; 0, 0), , ;
    , , mat(0, 1; 0, 0), ;
    , , , dots.down
  ),
$

evaluando en el origen, obtenemos un elemento del álgebra de Lie, asociado a su vector de $RR^(n+1)$, $va(x_k) = mat(theta_k, a_k, b_k, ...)^TT$

$
  [va(x_k)]_G^and = mat(
    mat(0, -theta_k; theta_k, 0), , , ;
    , mat(0, a_k; 0, 0), , ;
    , , mat(0, b_k; 0, 0), ;
    , , , dots.down
  ).
$

Podemos verificar que el mapa exponencial sigue siendo el exponencial matricial, notando que el exponencial de una matriz de bloques con los elementos fuera de la diagonal vacíos no es más que la matriz de los bloques exponenciados:

$
  exp[[va(x_k)]_G^and] =
  mat(
    exp[mat(0, -theta_k; theta_k, 0)], , , ;
    , exp[mat(0, a_k; 0, 0)], , ;
    , , exp[mat(0, b_k; 0, 0)], ;
    , , , dots.down
  ) = \
  = mat(
    R_(theta_k), , , ;
    , mat(1, a_k; 0, 1), , ;
    , , mat(1, b_k; 0, 1), ;
    , , , dots.down
  ) = X_k. \
$

De nuevo, no demostraremos que la expresión inversa es cierta para el logaritmo matricial por la complicación que esto supone.

#pagebreak()
= Interpretación de la dinámica orbital en $"SO"(2) times RR^n$

El modelo de un sistema en un filtro de Kalman en un grupo de Lie se escribe @markovicWrappingKalmanFilter2017

$
  X_(k + 1) = X_k exp_G [[Omega(X_k, u_k) + n_k]_G^and].
$

Para nuestro caso, no necesitamos toda la mecánica del filtro de Kalman. Consideraremos ruido aditivo nulo, $n_k = 0$. La función $Omega$ es el desplazamiento del estado, y asumiremos que no existe control, es decir, $u_k = 0$. Esto simplifica la expresión a

$
  X_(k + 1) = X_k exp_G [[Omega(X_k)]_G^and].
$

En nuestro caso particular de la dinámica orbital,

$
  X_k in G = "SO"(2) times RR^5 \
  Omega: G -> RR^6 \
  [.]_G^and: RR^6 -> frak(g) = frak("so")(2) times RR^5 \
  exp_G [.]: frak(g) -> G.
$

Intuitivamente:

- $Omega$ convierte de un elemento del grupo de Lie a un desplazamiento en su espacio tangente, expresado este como vector euclídeo convencional.
- $[.]_G^and$ convierte este vector a un elemento del álgebra de Lie.
- $exp_G [.]$ "integra" el elemento del álgebra de Lie para convertirlo en otro elemento del grupo de Lie.
- El producto del elemento del grupo de Lie del estado anterior con este nuevo, por como está definido el grupo de Lie, nos ofrece la nueva versión del estado.

== Tratamiento de la transición de estado como un elemento del grupo de Lie

No es immediatamente obvio que relación tiene $Omega$ con nuestro propagador orbital. Denominemos #box[$va(x_k) = (theta_k, a_k, b_k, c_k, d_k, e_k) in RR^6$] a nuestro vector de estado orbital (ignorando la realidad angular de $theta_k$, $[va(x_k)]_G^and in frak(g)$ es el verdadero elemento del álgebra de Lie que trata el ángulo correctamente).

Consideremos, heredando de la notación del filtro de Kalman utilizada en @markovicWrappingKalmanFilter2017, que el propagador orbital se puede escribir como una actualización de estado

$
  va(x_(k + 1)) = f(va(x_k)) = va(x_k) + hat(f_k)(va(x_k)),
$

notando que no se pierde generalidad por poder tener el término $-va(x_k)$ dentro de $hat(f_k)(x_k)$ @markovicWrappingKalmanFilter2017.

Ahora, consideremos que significado tiene la adición de cada término vectorial:

- La suma de ángulos $theta_k$ se translada a una multiplicación de matrices de rotación en el grupo de Lie
- La suma de elementos euclídeos $a_k$ se translada al producto de matrices $mat(1, a_1; 0, 1) mat(1, a_2; 0, 1) = mat(1, a_1 + a_2; 0, 1)$

Por lo tanto, siempre será posible escribir $Omega$ a partir de $hat(f_k)$. Con esto, garantizamos que el desarrollo que vamos a realizar dentro de la maquinaria del grupo de Lie sea aplicable en el caso numérico.

#pagebreak()
= Estadística en $"SO"(2) times RR^n$

== Media extrínseca (circular) en grupo $"SO"(2)$

Afortunadamente, no estamos tratando con cualquier grupo de Lie, si no que $G = "SO"(2) times RR^5$. Para este caso particular, podemos tomar prestado de la estadística direccional el concepto de media circular. En este apartado, demostraremos que la media de las matrices de rotación, bajo una interpretación típica de la acción de estas sobre un vector, es exactamente equivalente a la media circular.

Se define, para un conjunto de $N$ ángulos $theta_n$, su media circular como el número complejo @mardiaDirectionalStatistics1999

$
  z = 1/N sum_(n = 1)^N exp [i theta_n],
$

pudiendose obtener la dirección media como el argumento de este número (si existe). De forma similar, podríamos considerar la media "extrínsica" sobre un grupo de Lie a través de una construcción similar

$
  va(x) = 1 / N sum_(n = 1)^N X_n va(1)
$

dónde $va(1)$ es un vector unitario apropiado. Por ejemplo, para nuestro caso, podría ser

$
  va(1) = mat(1, 0, bar, 0, 1, bar, 0, 1, bar ...)^TT.
$

Definimos que la acción de un miembro de $G$ en este vector sea la multiplicación matricial

$
  X_k va(1) = mat(
    R_(theta_k), , , ;
    , mat(1, a_k; 0, 1), , ;
    , , mat(1, b_k; 0, 1), ;
    , , , dots.down
  ) mat(mat(1; 0); mat(0; 1); mat(0; 1); dots.v) =
  mat(R_(theta_k) mat(1; 0); mat(a_k; 1); mat(b_k; 1); dots.v).
$

Entonces, entendemos intuitivamente que el vector resultante tiene
- Por primera submatriz, el vector $mat(1, 0)^TT$ rotado $theta_k$ (equivalente a la exponenciación de $e^(i theta_k)$ si consideramos $1$ como el vector $mat(1, 0)^TT$ y la magnitud imaginaria como el vector $mat(0, 1)^TT$)
- Por siguientes submatrices, un vector afín "transladado y escalado", ya que la escala es $1$, simplemente extrae el valor de la matriz.

Por otra parte, por álgebra matricial, podemos escribir equivalentemente

$
  va(x) = (1 / N sum_(n = 1)^N X_n) va(1),
$

dónde denominaremos por $mu_E$ a la media de matrices, que de forma general no será miembro de $"SO"(2) times RR^n$. Ahora, podemos definir un operador proyección que convierte esta matriz en un miembro del grupo de Lie, $mu = cal(P)(mu_E)$. En el caso de $"SO"(2)$, la siguiente proyección es apropiada @khanMeansRandomVariables2025:

$
  cal(P)(A) = A (A^TT A)^(-1/2).
$

Para entender porqué esta transformación tiene sentido, consideremos una matriz fruto de sumar matrices de rotación escaladas,

$
  A = mat(x, -y; y, x),
$

y aplicamos la transformación. Primero de todo,

$
  A^TT A = mat(x^2 + y^2, 0; 0, x^2 + y^2),
$

que tiene raiz cuadrada inversa siempre por ser diagonal. Podemos escribir

$
  (A^TT A)^(-1/2) = I_(2 times 2) / sqrt(x^2 + y^2),
$

y finalmente

$
  cal(P)(A) = A / sqrt(x^2 + y^2).
$

¿Que sentido tiene esta transformación? Consideremos ahora la acción de $cal(P)(A)$ en el vector unitario

$
  cal(P)(A) mat(1; 0) = 1 / sqrt(x^2 + y^2) mat(x; y),
$

que es unitario y forma un ángulo con el eje $x$ de $"atan2"(y, x)$. Notamos que, por ser $A$ una media de matrices de rotación

$
  x = 1 / N sum_(n=1)^N cos(theta_n) quad quad y = 1 / N sum_(n=1)^N sin(theta_n),
$

que podemos escribir como un número complejo $z = x + i y$

$
  z = 1 / N sum_(n=1)^N e^(i theta_n),
$

es decir, la media de matrices de rotación es exactamente equivalente a la media circular definida en @mardiaDirectionalStatistics1999.

$cal(P)$ solo actúa en la parte de rotación, la media de las magnitudes euclídianas es directa (no requiere proyección). Por lo tanto, se podría construir un $cal(P)$ para $"SO"(2) times RR^n$.

== Dispersión en $"SO"(2) times RR^n$

Del anterior apartado concluimos que existe una definición práctica de la media en $"SO"(2)$, análoga a la media circular de la estadística direccional. Consideremos ahora el caso de las covarianzas y varianzas.

En su forma euclídea, las covarianzas se definen, para un conjunto de $N$ muestras en $RR^n$,

$
  "Cov"(i, j) = 1 / N sum_(n=1)^N ((X_n)_i - mu_i) ((X_n)_j - mu_j),
$

siendo $mu_i$ y $mu_j$ las medias en cada variable. En un grupo de Lie, esta noción de sustración (o distancia a la media) se representa mediante el producto con la inversa de $mu$,
$
  "Cov"_G (i, j) = 1 / N sum_(n = 1)^N ([log_G [mu^(-1) X_n]]_G^or)_i ([log_G [mu^(-1) X_n]]_G^or)_j.
$

Notamos la aparición del logaritmo para permitirnos extraer coordenadas. Podemos escribir la matriz de covarianza entera a través de

$
  P = 1 / N sum_(n = 1)^N xi_n xi_n^TT quad quad xi_n = [log_G [mu^(-1) X_n]]_G^or.
$

La aparición del logaritmo es problemática, sí las muestras se encuentran aproximadamente diametralmente opuestas a $mu$, el operator $[ log_G[.] ]_G^or$ causará problemas por la discontinuidad. Por lo tanto, esta expresión de la covarianza es apropiada solo para distribuciones concentradas.

== Dispersión en el espacio tangente

Como hemos visto, la matriz de covarianza no está claramente definida en el grupo de Lie, debido a la aparición del logaritmo. Por otra parte, supongamos que las muestras obedecen

$
  X_n = mu exp_G [ [epsilon_n]_G^and ],
$<map_samples>

es decir, las muestras se construyen mediante la adición a la media de una desviación euclídea. Este cambio de punto de vista permite definir la covarianza como

$
  P = 1 / N sum_(n = 1) epsilon_n epsilon_n^(TT),
$

que está definido únicamente al ser $epsilon_n in RR^n$. Este planteamiento es, esencialmente, asumir que la distribución de los puntos no vive en el grupo de Lie, sino en el álgebra de Lie (intuitivamente, el espacio tangente a $mu$).


== Propagación orbital bajo normal concentrada

En este caso asumimos que los puntos $X_n$ se distribuyen siguiendo una normal concentrada, es decir, sus valores en el álgebra de Lie están suficientemente concentrados como para que el mapa exponencial resulte en una distribución aproximadamente normal sobre el grupo de Lie.

Este es el caso que se considera típicamente en robótica @markovicWrappingKalmanFilter2017, ya que permite un tratamiento directo de la covarianza en el álgebra de Lie, y por tanto permite aplicar directamente el filtro de Kalman. Presentamos a continuación los resultados de @markovicWrappingKalmanFilter2017 adaptados al caso de la propagación orbital.

Consideremos que el estado $X_1 in G$ esté formado por dos componentes

$
  X_1 = mu_1 exp_G [[epsilon_1]_G^and],
$<locura1>

dónde $mu_1 in G$ y $epsilon_1 in RR^6$ representa una incertidumbre aproximada por una distribución normal multivariable en $RR^6$ de media nula y matriz de covarianza $P_1$, es decir, $epsilon_1 tilde cal(N)(0, P_1)$. Notamos que estamos asumiendo que la distribución estadística vive en el álgebra de Lie, no en el grupo.

Consideremos ahora el efecto en esta incertidumbre de nuestra función desplazamiento $Omega$, para obtener $X_2 in G$. Asumiremos que tras la propagación, el término de incertidumbre $epsilon_2$ sigue siendo normal, de media nula y matriz de covarianza $P_2$.

El procedimiento más convencional considera que la media se propaga ignorando el término $epsilon$ @markovicWrappingKalmanFilter2017

$
  mu_2 = mu_1 exp_G [[Omega(X_1)]_G^and].
$

Por otra parte, la propagación de la matriz de covarianza, siguiendo la notación de @markovicWrappingKalmanFilter2017, y particularizando para el caso sin ruido aditivo, se escribe

$
  P_2 = cal(F)_1 P_1 cal(F)_1^TT
$

dónde, para un grupo abeliano (la multiplicación de nuestras matrices es commutativa)

$
  cal(F)_1 = 1 + cal(C),
$

dónde @markovicWrappingKalmanFilter2017

$
  cal(C) = I_(n times n) + evaluated(pdv(, epsilon) Omega(mu_1 exp_G [ [epsilon]^and_G ]))_(epsilon = 0)
$

se trata del Jacobiano respecto a una perturbación en el espacio de $epsilon$.

Notamos que la aparición de la unidad aditiva no es más que un artefacto de cómo hemos relacionado $f(x)$ con $Omega(X)$, recordando la ecuación @locura1 y el párrafo que sigue.

Ahora, ¿que sucede si los puntos dan una vuelta completa al círculo? La respuesta es, absolutamente nada. La distribución estadística vive en el espacio tangente, plano, construido alrededor de la media. Con esta construcción, la evolución de la matriz de covarianza es totalmente independiente del grupo de Lie en el que vive la dinámica real #footnote[Nuestro caso de estudio es commutativo. Un contraejemplo sería $"SO"(3)$ que no es abeliano, y dónde todo es mucho más complicado.] y por lo tanto, la propagación de incertidumbre sobre un espacio euclídeo no requiere un tratamiento especial, siempre que se cumpla la hipótesis de Gaussiana concentrada.


== Distancias de Mahalanobis en grupos de Lie

La idea de propagar la incertidumbre en el espacio tangente es muy práctica, ya que permite utilizar un propagador convencional y permite la utilización del filtro de Kalman (extendido, unscented, etc...).

Consideremos ahora la computación de una distancia entre una distribución, $X_1 ~ (mu, P)$, con #box[$mu in G$] y #box[$P in RR^n$], y un punto #box[$X_2 in G$].

=== Distancia de Mahalanobis "naive"

Una primera idea podría ser construir la distribución en $RR^n$ a través de una normalización de la media, y comparar esta con la proyección del punto. Podemos escribir

$
  va(mu) = [log_G [mu]]_G^or, quad "tal que" quad va(x_1) ~ cal(N)(va(mu), P) \
  va(x_2) = [log_G [X]]_G^or.
$

Entonces, se define la distancia de Mahalanobis entre el punto y la distribución como

$
  d = sqrt((va(x_2) - va(mu))^TT P^(-1) (va(x_2) - va(mu))).
$

El problema de esta definición es que es sensible a los errores de enrollamiento, como se aprecia en @fig:discontinuity.

#figure(
  image("../writeup/img/statistics/naive_mahalanobis.svg", width: 90%),
  caption: [Distancia de Mahalanobis, se aprecia que, al no tener en cuenta la geometría real del espacio, surge una discontinuidad que afecta gravemente a la evaluación de la distancia al punto naranja.],
) <fig:discontinuity>

=== Distancia de Mahalanobis enrollada

Como alternativa, consideremos una operación muy similar, pero esta vez proyectando el punto al espacio tangente a $mu$,

$
  va(x_1) ~ cal(N)(0, P) \
  va(x_2) = [log_G [mu^(-1) X]]_G^or,
$

y se define distancia de Mahalanobis al origen

$
  d = sqrt(va(x_2)^TT P^(-1) va(x_2)).
$

Notamos que la acción del logaritmo es, esencialmente, tomar la distancia a la media enrollada al intervalo $[-pi, pi)$, y por lo tanto este método es equivalente a la distancia de Mahalanobis enrollada utilizada en robótica.

El comportamiento de esta transformación es mejor, ya que se aleja la discontinuidad del logaritmo de la zona de interés de la distribución (asumiendo una normal concentrada), como se aprecia en la @fig:discontinuity2.

#figure(
  image("../writeup/img/statistics/wrapped_mahalanobis.svg", width: 90%),
  caption: [Distancia de Mahalanobis alrededor de la media. Se reduce el efecto de la discontinuidad, siempre que la distribución esté suficientemente concentrada.],
) <fig:discontinuity2>


Por supuesto, si la distribución normal no es concentrada, surge el mismo problema que antes.

== Distancia de Mahalanobis Von Mises

La distribución Gauss Von Mises, introducida en @horwoodGaussMisesDistribution2014, consiste en una distribución conjunta Gaussiana y Von Mises, dónde la media angular de la distribución Von Mises depende lineal y cuadráticamente en las variables Gaussianas. En nuestro caso de estudio, asumiremos una dependencia únicamente lineal. La distribución es entonces

$
  cal("GVM")(theta, va(x)) = cal(N)(va(mu), A A^TT) cal("VM")(Theta(va(x)), kappa),
$

dónde

$
  Theta(va(x)) = alpha + va(beta)^TT A^(-1) (va(x) - va(mu))
$

y $va(mu) in RR^n$ es la media de las variables "euclídeas", $A in RR^(n times n)$ es la descomposición inferior de Cholesky de la matriz de covarianza $P in RR^(n times n)$, $alpha in RR$ es la "media angular", $va(beta) in RR^n$ representa el acomplamiento lineal entre las distribuciones, y $kappa in RR$ es la dispersión de Von Mises.

En @horwoodGaussMisesDistribution2014 se define la distancia Mahalanobis Von Mises de un punto como la suma de la distancia de Mahalanobis en el espacio euclídeo, más una distancia angular. Para un punto $(theta, va(x))$, escribimos

$
  d^2 = (va(x) - va(mu))^TT (A A^(-1))^(-1) (va(x) - va(mu)) + 4 kappa sin^2 (1/2 (theta - alpha - beta^TT A^(-1) (va(x) - va(mu)))).
$

== Obtención de la distribución Gauss Von Mises a partir de una Gaussiana

Tenemos la expresión para calcular la distancia de Mahalanobis Von Mises, pero carecemos de un método para convertir nuestra Gaussiana en el álgebra de Lie en una distribución Gauss Von Mises.

La función característica para esta distribución es @horwoodGaussMisesDistribution2014

$
  phi^"GVM" (va(xi), m) = (I_(abs(m))(kappa)) / (I_0(kappa)) exp[i (va(mu)^TT va(xi) + m alpha) - 1/2 (A^TT va(xi) + m va(beta))^TT (A^TT va(xi) + m va(beta))],
$

dónde $va(xi) in RR^n$ y $m in ZZ$. Recordando la definición de la función característica,

$
  phi^"GVM" (va(xi), m) = EE[e^i(va(xi)^TT va(x) + m theta)]
$

se hace evidente porqué $m$ queda limitado a $ZZ$, esto se debe a que, para $m$ no entero, $e^(i m theta)$ no sería cíclico en $theta$, y la expresión carecería de sentido @mardiaDirectionalStatistics1999.

Ahora, consideremos la siguiente descomposición del espacio $RR^n$ característico del álgebra de Lie. Sea nuestra variable aleatoria $X = mu exp_G [ [va(epsilon)]_G^and ]$, con media proyectada $va(mu) = [log_G [mu]]_G^or$ y $va(epsilon) ~ cal(N)(0, P)$. Consideremos además el muestreo $va(x) = va(mu) + va(epsilon)$, que notamos diferente a la proyección de la variable aleatoria al álgebra de Lie, ya que esta $va(x)$ no sufre discontinuidades del logaritmo.

Podemos descomponer en bloques

$
  va(mu) = mat(hat(alpha); va(mu)_e), \
  hat(P) = mat(
    sigma^2, va(gamma)^TT;
    va(gamma), hat(P_e)
  ),\
  va(x) = mat(theta; va(x_e))
$

tal que podemos escribir la función característica de la distribución normal multivariable (ver anexo `expand_normal_characteristic.py`)

$
  phi^N (va(eta))
  = exp[i va(mu)^TT va(eta) - 1/2 va(eta)^TT P va(eta)] = \ =
  phi^N (va(eta_e), n) =
  exp[i (va(mu_e)^TT va(eta_e) + n hat(alpha)) - 1/2 sigma^2 n^2 - n va(eta_e)^TT va(gamma)- 1/2 va(eta_e)^TT hat(P_e) va(eta_e)]
$

dónde $va(eta) in RR^n$, $va(eta_e) in RR^(n-1)$ y $n in RR$. Si realizamos la factorización de Cholesky de $hat(P_e)$, y la denominamos $hat(A)_e$, podemos asimilar esta expresión a la de Gauss Von Mises.

Podemos para ello realizar la multiplicación (ver anexo `expand_normal_characteristic.py`)

$
  1/2 (A^TT va(xi) + m va(beta))^TT (A^TT va(xi) + m va(beta)) = \
  =
  1/2 [va(xi)^TT A A^TT va(xi) + va(xi)^TT A m va(beta) +
    va(beta)^TT m A^TT va(xi) + m^2 va(beta)^TT va(beta)].
$

Notamos que

$
  (va(beta)^TT m A^TT va(eta)_e)^TT = va(eta)_e^TT A m va(beta),
$

y esta última magnitud es un escalar (es decir, se cumple $x^TT = x$ trivialmente), por lo que podemos agrupar los dos términos.

Por otra parte, recordamos que $A A^TT = hat(P_e)$ por ser su descomposición de Cholesky. Por lo tanto

$
  va(xi)^TT A A^TT va(xi) = va(xi)^TT hat(P)_e va(xi)
$

Tenemos entonces

$
  1/2 (A^TT va(xi) + m va(beta))^TT (A^TT va(xi) + m va(beta)) = \
  m va(xi)^TT A va(beta) + 1/2 m^2 va(beta)^TT va(beta) + 1/2 va(xi)^TT P va(xi).
$

=== Ajuste por moment-matching si $kappa -> infinity$

Vamos a comparar ahora las expresiones de $phi^N$ y $phi^"GVM"$, bajo la anterior expansión

$
  phi^"GVM" & = (I_(abs(m))(kappa)) / (I_0(kappa)) &exp[ & i (va(mu)^TT va(xi) + m alpha) & - m va(xi)^TT A va(beta) &- 1/2 m^2 va(beta)^TT va(beta) &- 1/2 va(xi)^TT P va(xi).
  ], \
  phi^N & = &exp[&i (va(mu_e)^TT va(eta_e) + n hat(alpha)) &- n va(eta_e)^TT va(gamma) &- 1/2 sigma^2 n^2 &- 1/2 va(eta_e)^TT hat(P_e) va(eta_e)].
$

Observamos una marcada semejanza, lo que justifica una aproximación entre ambas. Comenzemos por el caso $kappa -> infinity$, ya que este es presentado en @horwoodGaussMisesDistribution2014 (si bien mediante un desarrollo alternativo). En este caso, $lim_(kappa -> infinity) (I_(abs(m))(kappa)) / (I_0(kappa)) = 1$.

Dentro del exponencial, asemejamos, para $va(mu) = va(mu_e)$, $va(xi) = va(eta_e)$ y $m = n$,

$
  hat(alpha) = alpha \
  hat(P_e) = P, \
  hat(P_e) = hat(A)_e hat(A)_e^TT,
$

dónde $P = A A^TT$ de la Gauss Von Mises. Es importante notar que es necesario realizar de nuevo la descomposición de Cholesky para la submatriz "euclídea". Por otra parte, igualando términos según su orden en $n = m$

$
  va(beta)^TT va(beta) = sigma^2,
$

y la ecuación (sustituyendo $va(xi) = va(eta_e)$)

$
  va(xi)^TT A va(beta) = va(xi^TT) va(gamma) \
$

de donde concluimos por independencia de la igualdad en $va(xi)$
$
  va(gamma) = A va(beta)
$

lo que nos permite obtener

$
  hat(P) = mat(
    va(beta)^TT va(beta), va(beta)^TT A^TT;
    A va(beta), P
  ),
$

expresión idéntica hasta el orden de filas con la presentada en @horwoodGaussMisesDistribution2014 cuando $kappa -> infinity$.

=== Ajuste por moment-matching con $kappa$ finito

En el caso de $kappa$ finito, asumamos igual que antes $va(mu) = va(mu_e)$, $va(xi) = va(eta_e)$ y $m = n$ para realizar el moment-matching. Ahora, no es posible igualar término a término debido a la aparición del factor de escala fuera de la exponencial, por lo tanto la distribución normal obtenida no será exacta (no podría serlo de ninguna forma). En particular, plantearemos que la distribución normal tenga idénticos momentos de primer orden angulares (es decir, media angular y dispersión) que la Gauss Von Mises.

Consideremos el caso $va(xi) = 0$ y $m = 1$, tenemos entonces

$
  phi^"GVM" & = (I_1(kappa)) / (I_0(kappa)) &exp[ & i alpha & - 1/2 va(beta)^TT va(beta) ], \
  phi^N & = &exp[ & i hat(alpha) & - 1/2 sigma^2 ]. \
$

Igualando las fases, garantizamos que ambas distribuciones tengan la misma media angular @mardiaDirectionalStatistics1999, lo que permite obtener $alpha = hat(alpha)$. Igualando las normas, garantizamos que ambas tengan la misma dispersión, y obtenemos la siguiente ecuación para $kappa$ y $beta$,

$
  (I_1(kappa)) / (I_0(kappa)) = exp[(va(beta)^TT va(beta) - sigma^2) / 2],
$<numericaBessel>

que no tiene solución analítica, pero se comporta bien numéricamente (más adelante lo demostramos).

Ahora, consideremos el caso con $va(xi) != va(0)$ y $m = 0$, tenemos entonces

$
  phi^"GVM" & = & exp[ & i va(mu)^TT va(xi)   &         - 1/2 va(xi)^TT P va(xi).
                                                        ], \
      phi^N & = & exp[ & i va(mu_e)^TT va(xi) & - 1/2 va(xi)^TT hat(P_e) va(xi)]. \
$

De nuevo, igualando fases podemos relacionar las medias y obtener

$
  va(mu_e) = va(mu),
$

y igualando las magnitudes, relacionamos las matrices de covarianza y obtenemos

$
  P = hat(P_e) \
  P = A A^TT,
$

notando que es necesario realizar la descomposición de Cholesky a la submatriz. Por último, para obtener $va(beta)$, estudiaremos el caso en el que $va(xi) -> 0$ a través del gradiente de la función característica en el origen. Para ello, sustituimos las igualdades anteriores,

$
  phi^"GVM" & = (I_(abs(m))(kappa)) / (I_0(kappa)) &exp[ & i (va(mu)^TT va(xi) + m alpha) & - m va(xi)^TT A va(beta) &- 1/2 m^2 va(beta)^TT va(beta) &- 1/2 va(xi)^TT P va(xi).
  ], \
  phi^N & = &exp[&i (va(mu)^TT va(xi) + m alpha) &- m va(xi)^TT va(gamma) &- 1/2 sigma^2 m^2 &- 1/2 va(xi)^TT P va(xi)]. \
$

y totamos el gradiente respecto a $xi$, evaluado en $xi = 0$,


$
  nabla_xi phi^"GVM" & = (I_(abs(m))(kappa)) / (I_0(kappa)) &[&i va(mu)^TT &- m A va(beta) &- P va(xi)] &exp[ & i m alpha & - 1/2 m^2 va(beta)^TT va(beta)], \
  nabla_xi phi^N & = &[&i va(mu)^TT &- m va(gamma) &- P va(xi)]&exp[ & i m alpha & - 1/2 sigma^2 m^2 ]. \
$

Igualando ambas expresiones, y sustituyendo la @numericaBessel, tenemos

$
  A va(beta) = va(gamma) => va(beta) = A^(-1) va(gamma).
$

=== Evaluación del moment-matching

Para evaluar las anteriores experiones, planteamos un método Monte Carlo en el que se muestreará una normal multivariable, su Gauss Von Mises correspondiente obtenida por el método previamente introducido, y se evaluará la similitud de las nubes de puntos.

Para ello se ha planteado un método "Classifier two-sample test (C2ST)", basado en entrenar un clasificador con muestras de ambas distribuciones @lopez-pazRevisitingClassifierTwoSample2018.

Para ello, se parte de los conjunto de puntos $S_cal(N) ~ cal(N)(va(mu), P)$ y $S_cal("GVM") ~ cal("GVM")(va(mu), P, va(beta), kappa)$, y se marcan las muestras de cada uno con un indicator positivo o negativo, respectivamente. Posteriormente, se entrena un clasificador que, dado un conjunto de muestras, prediga si el indicador es positivo o negativo. Si el clasificador es exitoso (clasifica mejor que al azar), consideraremos ambas muestras tomadas de distribuciones similares. Si por el contrario, el clasificador no logra superar el azar, se considerarán que las muestras son de distribuciones diferentes. Para el entrenamiento simplemente se muestrean subconjuntos de ambos conjuntos de puntos, junto con sus indicadores.

La ventaja de estos métodos es que escalan linealmente con el número de puntos (al contrario de pruebas estadísticas más rigurosas, que escalan con el cuadrado) y son muy fáciles de implementar partiendo de las librerías de aprendizaje por máquina disponibles en Python. Todo el código se encuentra en el anexo `moment_matching_evaluation.py`.

