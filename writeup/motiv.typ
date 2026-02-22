#import "@preview/physica:0.9.8": *
#import "utils.typ": *

#let cite_mech(x) = cite(<arnold_mathematical_1989>, supplement: [pág. #x])

== Conceptos básicos de topología

Comenzamos la sección con un repaso breve de varios conceptos que serán necesarios en el resto del trabajo. Si bien las definiciones no tienen un rigor matemático completo, son suficientes para su uso en mecánica clásica. Todas ellas han sido tomadas y adaptadas de @arnold_mathematical_1989.

#definition[Una *variedad topológica* se define como un conjunto de puntos, junto con una cantidad contable de cartas, tal que todo punto del conjunto esté representado por, al menos, una carta #cite_mech(77)] <def:manifold>

#definition[Una *carta* se puede definir como un subconjunto abierto $U$ de puntos de un espacio de coordenadas euclídeo, junto con una función biyectiva $f$ (uno a uno en ambas direcciones) que permite asignar a cada punto de $U$ un punto de la variedad #cite_mech(77). Representaremos la carta con la notación $(U, f)$.] <def:chart>

Equivalentemente, podemos hablar de una carta como el mapa de un cierto conjunto de puntos de la variedad al espacio euclídeo. Intuitivamente, cada carta nos permite asignar coordenadas a un cierto subconjunto de los puntos de la variedad. Esta última definición es más típica en los textos de matemáticas, mientras que la anterior es más apropiada en mecánica clásica.

#example[Carta trivial para el espacio euclídeo][
  El espacio euclídeo $RR^n = {(q_1, ..., q_n) | q_i in RR}$ es una variedad topológica trivialmente, tomando $U = RR^n$, y el mapa identidad, formando la carta $(RR^n, "Id")$.
]

#example[Carta trivial para la esfera][
  Un ejemplo menos trivial es la esfera $S^2 = {(x_1, x_2, x_3) in RR^3 | x_1^2 + x_2^2 + x_3^2 = 1}$, ya que $S^2$ no es un espacio euclídeo, no podemos utilizar la carta trivial que utilizamos en el ejemplo anterior.
]

#example[Las coordenadas esféricas NO son una carta para la esfera][
  Retomando el ejemplo anterior, consideremos $U = {(theta, phi) in RR^2 | theta in [0, pi), phi in [0, 2 pi)}$, y el mapa
  $
    x(theta, phi) = mat(x_1; x_2; x_3) = mat(
      sin theta cos phi;
      sin theta sin phi;
      cos theta
    ),
  $

  Si bien $x$ hace posible asignar a cada punto del rectángulo en $RR^2$ un punto de la esfera, el sistema no es una carta. Por una parte, incumplimos el requisito de que $U$ sea un subconjunto abierto de $RR^2$. Además, el mapa no es biyectivo. Por ejemplo, para $mat(x_1, x_2, x_3)^TT = (0, 0, 1)^TT$ tenemos $theta = 0$ y $phi in [0, 2 pi)$. Gráficamente:

  #text(red)[Figura]
]

#definition[Una *variedad diferenciable* es una variedad topológica que puede ser totalmente cubierta por un atlas #cite_mech(77) #footnote[Convencionalmente, se define la variedad diferencial a partir del atlas, con cierta pérdida del valor intuitivo.]] <def:diff_manifold>

#definition[Un *atlas* es la unión de varias cartas compatibles #cite_mech(78).]
#definition[Dos cartas $(U, f)$ y $(V, g)$ de la misma variedad $M$ son *compatibles* sí, para los subconjuntos $X subset U$ y $Y subset V$, ambos con la misma imagen en $M$, la función que lleva desde cada punto de $X$ a $Y$ es diferenciable #cite_mech(78).]

En nuestro caso de estudio particular, podremos asumir siempre que $X$ y $Y$ son espacios de las mismas dimensiones, digamos $n$, y por lo tanto el mapa de $X$ a $Y$ serán $n$ funciones de $n$ variables. Por lo tanto, que el mapa sea diferenciable no es más que estas funciones sean independientemente diferenciables en cuanto a cada variable.

Además, podemos obtener el mapa notando que las cartas son biyectivas, por lo tanto, para un punto $x in X$, tenemos que $f(x) = m$, siendo $m$ un punto de la variedad, y $g^(-1)(m) = y$, siendo $y in Y$. Por lo tanto, el mapa de $X$ a $Y$ no es más que $g^(-1)(f(x))$, expresado de forma general como $g^(-1) compose f$.

#text(red)[Figura]

#example[Atlas para la esfera][
  El mapa estereográfico de la esfera respecto a un punto $p in S^2$ se genera proyectando cada punto de esta, excepto $p$, sobre el plano tangente a la esfera en el punto $q$ diametralmente opuesto a $p$:

  #text(red)[Figura]

  Reaizando el desarrollo con coordenadas, consideremos el sistema de coordenadas de la figura. Comenzamos con una carta $(U, f)$, formada por la proyección estereográfica de la esfera desde el punto $n = (0, 0, 1)^TT in S^2$ (el polo norte). $f^(-1): S^2 \\{n} -> U$ se puede obtener por trigonometría según

  $
    mat(u_1; u_2) = mat(
      x_1 / (1 - x_3);
      x_2 / (1 - x_3)
    ),
  $

  aplicable para toda $S^2$ menos $n$.

  De forma análoga, podemos obtener una segunda carta $(V, g)$ de igual forma, pero proyectando desde $s = (0, 0, -1)^TT in S^2$ (el polo sur), y obteniendo #box[$g^(-1): S^2 \\ {s} -> V$] según

  $
    mat(v_1; v_2) = mat(
      x_1 / (1 + x_3);
      x_2 / (1 + x_3)
    ).
  $

  Ahora, a excepción de los polos norte (punto $n$) y sur (punto $s$), todo el resto de puntos de $S^2$ son compartidos entre ambas cartas, y por lo tanto, para que estas formen un atlas, el mapa de una carta a otra debe ser diferenciable en todo punto compartido. Previamente hemos deducido que el mapa debe ser $g^(-1) compose f$. Tras algo de manipulación algebraica #footnote[Disponible en el anexo "Sphere.nb"] <annex:sphere> , este mapa queda como

  $
    mat(v_1; v_2) = mat(
      u_1 / (u_1^2 + u_2^2);
      u_2 / (u_1^2 + u_2^2);
    ),
  $

  cuyo Jacobiano @annex:sphere existe en todo punto menos $mat(0, 0)^T$.

  Concluimos por tanto que la esfera es una variedad diferenciable, ya que hemos encontrado un atlas que la cubre por completo.
] <ej:sphereAtlas>

La introducción del concepto de diferenciabilidad de una variedad es vital, ya que nos permite definir el espacio tangente a un punto de la variedad.

#definition[Consideremos una variedad diferenciable $M$, y dos curvas denominadas $phi: RR -> M$ y $psi: RR -> M$, así como una carta en la que ambas tomen coordenadas (en la carta) $vb(phi)(t)$ y $vb(psi)(t)$ respectivamente. Ambas se denominan *curvas equivalentes en un punto* sí $vb(phi)(0) = vb(psi)(0)$ y $lim_(t -> 0) (vb(phi)(t) - vb(psi)(t)) / t = 0$ #cite_mech(81).]


#example[Dos curvas equivalentes en la esfera][
  Como ejemplo, consideremos las curvas cuya proyección en la carta estereográfica desde el polo norte $(U, f)$ del @ej:sphereAtlas, en particular

  $
    vb(phi)(t) = mat(phi_1(t); phi_2(t)) = mat(
      1 + t;
      1 + t^2;
    ) \
    vb(psi)(t) = mat(psi_1(t); psi_2(t)) = mat(
      1 + sin(t);
      2 - cos(t);
    ).
  $

  Notamos que $vb(phi)(0) = mat(1, 1)^TT$ y $vb(psi)(0) = mat(1, 1)^TT$, y a su vez

  $
    lim_(t -> 0) (vb(phi)(t) - vb(psi)(t)) / t =
    lim_(t -> 0) 1 / t mat(
      t - sin(t);
      t^2 + cos(t) - 1
    ) =
    lim_(t -> 0) mat(
      1 - sin(t)/t;
      t + (cos(t) - 1)/t
    ) =
    mat(0; 0).
  $

  Por lo tanto, ambas curvas son equivalentes para $t = 0$. En la @fig:sphere_equivalent_curves se pueden ver representadas ambas curvas en las dos cartas de la esfera y en la propia variedad @annex:sphere.

  #figure(
    image("img/motiv/sphere_equivalent_curves.svg", width: 80%),
    caption: [Dos curvas equivalentes en la carta estereográfica desde el polo norte (izquierda), desde el polo sur (derecha) y en la esfera representada dentro de $RR^3$ (derecha).],
  ) <fig:sphere_equivalent_curves>

]

== Mecánica Lagrangiana para una partícula en órbita

Consideremos una partícula orbitando la Tierra en una órbita elíptica. Con el fin de evitar definir espacios de coordenadas excesivamente temprano, denominemos a la variedad diferenciable de posibles posiciones del satélite $Q$, y la posición real del satélite es un punto de esta variedad diferenciable, $p in Q$. En términos de mecánica clásica, $Q$ es el espacio de configuración del satélite #cite_mech(53).

#example[Espacio de configuración de un péndulo esférico][

  Como ejemplo, consideremos un péndulo simple en tres dimensiones, es decir, una partícula unida rígidamente a un eje que puede girar en todas las direcciones. Consideremos que sobre la partícula actúa la gravedad, asumida unidireccional.

  Por intuición física, podemos deducir que la partícula puede situarse en una esfera dentro del volumen euclídeo. Esto nos lleva a pensar que el espacio de configuración $Q$ del péndulo se corresponde a la esfera. De forma más rigurosa, afirmamos que $Q$ es isomórfico a la esfera $S^2$, expresado como $Q tilde.equiv S^2$.

  Esta definición anterior es totalmente libre de coordenadas. En particular, podemos ofrecer una visión de la configuración del péndulo en cuanto a su vector posición euclídeo $vb(x) = x_1 vb(e_1) + x_2 vb(e_2) + x_3 vb(e_3)$ tal que $x_1^2 + x_2^2 + x_3^2 = 1$ definiendo el centro del péndulo como el origen del sistema de coordenadas, y la base vectorial ${vb(e_1), vb(e_2), vb(e_3)}$ como vectores de longitud igual a la longitud de la barra rígida del péndulo, uno de ellos apuntando en la dirección de la gravedad, y los otros dos completando el sistema de coordenadas por regla de la mano derecha:

  #text(red)[FIGURA]

  Por otra parte, podríamos ofrecer una definición del sistema en cuanto a tan solo dos coordenadas angulares, $theta in [0, pi)$ y $phi in [0, 2 pi]$, eliminando la ecuación de restricción:

  #text(red)[FIGURA]

] <ej:conf_pendulo>

El espacio de configuración no es suficiente si consideramos la evolución temporal del sistema, ya que no representa las velocidades. Asumiendo que $Q$ es una variedad diferenciable, el espacio de estados se define entonces como el fibrado de los espacios tangentes a $Q$, $T Q$ #cite_mech(53).

#example[Espacio de estados de un péndulo esférico][
  Retomando el @ej:conf_pendulo, el espacio tangente en cada punto $T_p Q$ es el espacio de velocidades de la partícula en cada uno de estos puntos. Debido a que la partícula está a una distancia fija del centro del péndulo, este espacio tangente es isomórfico a un plano, $T_p Q tilde.equiv RR^2$:

  #text(red)[FIGURA]

  Para realizar el fibrado de estos espacios tangentes,









] <ej:fase_pendulo>
