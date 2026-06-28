#import "@preview/bananote:0.1.2": *
#import "@preview/physica:0.9.8": *

#show: note.with(
  title: [Filtro de Kalman enrollado],
)

#abstract[Se parte del paper @markovicWrappingKalmanFilter2017 y se particulariza a las necesidades de la propagación orbital.]

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

= Particularización e intuición para $"SO"(2)$

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
  exp_G[[theta_k]_G^and] = exp[mat(0, -theta_k; theta_k, 0)] =
  sum_(n = 0)^infinity mat(0, -theta_k; theta_k, 0)^n / (n!) = \
  mat(1, 0; 0, 1)
  + mat(0, -theta_k; theta_k, 0)
  + 1 / 2 mat(-theta_k^2, 0; 0, -theta_k^2)
  + 1 / 6 mat(0, theta_k^3; -theta_k^3, 0)
  + ...
$

que identificamos con el desarrollo en serie del seno y coseno en cada posición de la matriz, para obtener

$
  exp_G[[theta_k]_G^and] = mat(cos(theta_k), -sin(theta_k); sin(theta_k), cos(theta_k)) = R_(theta_k)
$

El logaritmo es más complicado de desarrollar, ya que presenta varias ramas. Para el desarrollo que realizaremos, no es necesaria intuición sobre el logaritmo matricial, así que no entraremos en detalle.


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

evaluando en el origen, obtenemos un elemento del álgebra de Lie, asociado a su vector de $RR^(n+1)$ $x_k = mat(theta_k, a_k, b_k, ...)^TT$

$
  [x_k]_G^and = mat(
    mat(0, -theta_k; theta_k, 0), , , ;
    , mat(0, a_k; 0, 0), , ;
    , , mat(0, b_k; 0, 0), ;
    , , , dots.down
  ).
$

Podemos verificar que el mapa exponencial sigue siendo el exponencial matricial, notando que el exponencial de una matriz de bloques con los elementos fuera de la diagonal vacíos no es más que la matriz de los bloques exponenciados:

$
  exp[[x_k]_G^and] =
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

= Interpretación de la dinámica orbital en $"SO"(2) times RR^n$

El modelo de un sistema en un filtro de Kalman en un grupo de Lie se escribe @markovicWrappingKalmanFilter2017

$
  X_(k + 1) = X_k exp_G [[Omega(X_k, u_k) + n_k]_G^and].
$

Para nuestro caso, no necesitamos toda la mecánica del filtro de Kalman. Consideraremos ruido aditivo nulo, $n_k =$. La función $Omega$ es el desplazamiento del estado, y asumiremos que no existe control, es decir, $u_k = 0$. Esto simplifica la expresión a

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

No es immediatamente obvio que relación tiene $Omega$ con nuestro propagador orbital. Denominemos #box[$x_k = (theta_k, a_k, b_k, c_k, d_k, e_k) in RR^6$] a nuestro vector de estado orbital (ignorando la realidad angular de $theta_k$, $[x_k]_G^and in frak(g)$ es el verdadero elemento del álgebra de Lie que trata el ángulo correctamente).

Consideremos, heredando de la notación del filtro de Kalman utilizada en @markovicWrappingKalmanFilter2017, que el propagador orbital se puede escribir como una actualización de estado

$
  x_(k + 1) = f(x_k) = x_k + hat(f_k)(x_k),
$

notando que no se pierde generalidad por poder tener el término $-x_k$ dentro de $hat(f_k)(x_k)$ @markovicWrappingKalmanFilter2017.

Ahora, consideremos que significado tiene la adición de cada término vectorial:

- La suma de ángulos $theta_k$ se translada a una multiplicación de matrices de rotación en el grupo de Lie
- La suma de elementos euclídeos $a_k$ se translada al producto de matrices $mat(1, a_1; 0, 1) mat(1, a_2; 0, 1) = mat(1, a_1 + a_2; 0, 1)$

Por lo tanto, siempre será posible escribir $Omega$ a partir de $hat(f_k)$. Con esto, garantizamos que el desarrollo que vamos a realizar dentro de la maquinaria del grupo de Lie sea aplicabe en el caso numérico.

= Tratamiento estadístico del grupo de Lie





#bibliography("../writeup/refs.bib")

