#import "@preview/physica:0.9.8": *
#import "utils.typ": *

#let cite_mech(x) = cite(<arnold_mathematical_1989>, supplement: [pág. #x])

== Conceptos básicos de topología

#definition[Una *variedad topológica* se define como un conjunto de puntos, junto con una cantidad contable de cartas, tal que todo punto del conjunto esté representado por, al menos, una carta #cite_mech(77)] <def:manifold>

#definition[Una *carta* se puede definir como un subconjunto abierto $U$ de puntos de un espacio de coordenadas euclídeo, junto con una función biyectiva $f$ (uno a uno en ambas direcciones) que permite asignar a cada punto de $U$ un punto de la variedad #cite_mech(77). Representaremos la carta con la notación $(U, f)$.] <def:chart>

Equivalentemente, podemos hablar de una carta como el mapa de un cierto conjunto de puntos de la variedad al espacio euclídeo. Intuitivamente, cada carta nos permite asignar coordenadas a un cierto subconjunto de los puntos de la variedad. Esta última definición es más típica en los textos de matemáticas, mientras que la anterior es más apropiada en mecánica clásica.

#example[Carta para el espacio euclídeo][
  El espacio euclídeo $RR^n = {(q_1, ..., q_n) | q_i in RR}$ es una variedad topológica trivialmente, tomando $U = RR^n$, y el mapa identidad, formando la carta $(RR^n, "Id")$.
]

#example[Carta trivial para la esfera][
  Un ejemplo menos trivial es la esfera $S^2 = {(x_1, x_2, x_3) in RR^3 | x_1^2 + x_2^2 + x_3^2 = 1}$, ya que $S^2$ no es un espacio euclídeo, no podemos utilizar la carta trivial que utilizamos en el ejemplo anterior.
]

#example[Las coordenadas esféricas NO son una carta para la esfera][
  Retomando el ejemplo anterior, consideremos $U = {(theta, phi) in RR^2 | theta in [0, pi), phi in [0, 2 pi)}$, y el mapa
  $
    mat(x_1; x_2; x_3) = mat(
      sin theta cos phi;
      sin theta sin phi;
      cos theta
    ),
  $

  por una parte, incumplimos el requisito de que $U$ sea un subconjunto abierto de $RR^2$. Además, el mapa no es biyectivo. Por ejemplo, para $mat(x_1, x_2, x_3)^TT = (0, 0, 1)^TT$ tenemos $theta = 0$ y $phi in [0, 2 pi)$. Gráficamente:

  #text(red)[Figura]
]

#definition[Una *variedad diferenciable* es una variedad topológica que es totalmente cubierta por un atlas #cite_mech(77).] <def:diff_manifold>

#definition[Un *atlas* es la unión de varias cartas compatibles #cite_mech(78).]
#definition[Dos cartas $(U, f)$ y $(U', f')$ de la misma variedad $M$ son *compatibles* sí, para los subconjuntos $V subset U$ y $V' subset U'$, ambos con la misma imagen en $M$, la función que lleva desde cada punto de $V$ a $V'$ es diferenciable #cite_mech(78).]

#example[Atlas para la esfera][

]




== Mecánica clásica para una partícula orbitando la Tierra

Consideremos una partícula orbitando la Tierra en una órbita elíptica. Con el fin de evitar definir espacios de coordenadas excesivamente temprano, denominemos a la variedad diferenciable de posibles posiciones del satélite $Q$, y la posición real del satélite es un punto de esta variedad diferenciable, $p in Q$. En términos de mecánica clásica, $Q$ es el espacio de configuración del satélite #cite_mech(53).

#example[Espacio de configuración de un péndulo esférico][

  Como ejemplo, consideremos un péndulo simple en tres dimensiones, es decir, una partícula unida rígidamente a un eje que puede girar en todas las direcciones. Consideremos que sobre la partícula actúa la gravedad, asumida unidireccional.

  Por intuición física, podemos deducir que la partícula puede situarse en una esfera dentro del volumen euclídeo. Esto nos lleva a pensar que el espacio de configuración $Q$ del péndulo se corresponde a la esfera. De forma más rigurosa, afirmamos que $Q$ es isomórfico a la esfera $S^2$, expresado como $Q tilde.equiv S^2$.

  Esta definición anterior es totalmente libre de coordenadas. En particular, podemos ofrecer una visión de la configuración del péndulo en cuanto a su vector posición euclídeo $vb(x) = x_1 vb(e_1) + x_2 vb(e_2) + x_3 vb(e_3)$ tal que $x_1^2 + x_2^2 + x_3^2 = 1$ definiendo el centro del péndulo como el origen del sistema de coordenadas, y la base vectorial ${vb(e_1), vb(e_2), vb(e_3)}$ como vectores de longitud igual a la longitud de la barra rígida del péndulo, uno de ellos apuntando en la dirección de la gravedad, y los otros dos completando el sistema de coordenadas por regla de la mano derecha:

  #text(red)[FIGURA]

  Por otra parte, podríamos ofrecer una definición del sistema en cuanto a tan solo dos coordenadas angulares, $theta in [0, pi)$ y $phi in [0, 2 pi]$, eliminando la ecuación de restricción:

  #text(red)[FIGURA]

  Más rigurosamente,

] <ej:conf_pendulo>

El espacio de configuración no es suficiente si consideramos la evolución temporal del sistema, ya que no representa las velocidades. Asumiendo que $Q$ es una variedad diferenciable, el espacio de estados se define entonces como el fibrado de los espacios tangentes a $Q$, $T Q$ #cite_mech(53).

#example[Espacio de estados de un péndulo esférico][
  Retomando el @ej:conf_pendulo, el espacio tangente en cada punto $T_p Q$ es el espacio de velocidades de la partícula en cada uno de estos puntos. Debido a que la partícula está a una distancia fija del centro del péndulo, este espacio tangente es isomórfico a un plano, $T_p Q tilde.equiv RR^2$:

  #text(red)[FIGURA]

  Para realizar el fibrado de estos espacios tangentes,









] <ej:fase_pendulo>
