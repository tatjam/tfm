#import "@preview/physica:0.9.8": *
#import "utils.typ": *

#let cite_mech(x) = cite(
  <arnoldMathematicalMethodsClassical1989>,
  supplement: [pág. #x],
)


== Mecánica orbital Lagrangiana

Tras esta breve introducción a las variedades diferenciales, vamos a explorar una aplicación directa de estas, alejándonos de la nomenclatura matemática y acercándonos a la física.

Consideremos una partícula orbitando la Tierra en una órbita elíptica. Denominemos a la variedad diferenciable de posibles posiciones de la partícula $Q$, y la posición real del satélite es un punto de esta variedad diferenciable, $p in Q$. En términos de mecánica clásica, $Q$ es el espacio de configuración del satélite #cite_mech(53).

#example[Espacio de configuración de un péndulo esférico][

  Como ejemplo, consideremos un péndulo simple en tres dimensiones, es decir, una partícula unida rígidamente a un eje que puede girar en todas las direcciones. Consideremos que sobre la partícula actúa la gravedad, asumida unidireccional.

  Por intuición física, podemos deducir que la partícula puede situarse en una esfera dentro del volumen euclídeo. Esto nos lleva a pensar que el espacio de configuración $Q$ del péndulo se corresponde a la esfera. Es decir, que $Q$ es isomórfico a la esfera $S^2$, expresado como $Q tilde.equiv S^2$.

  Esta definición anterior es totalmente libre de coordenadas. En particular, podemos ofrecer una visión de la configuración del péndulo en cuanto a su vector posición euclídeo $vb(x) = x_1 vb(e_1) + x_2 vb(e_2) + x_3 vb(e_3)$ tal que $x_1^2 + x_2^2 + x_3^2 = 1$ definiendo el centro del péndulo como el origen del sistema de coordenadas, y la base vectorial ${vb(e_1), vb(e_2), vb(e_3)}$ como vectores de longitud igual a la longitud de la barra rígida del péndulo, uno de ellos apuntando en la dirección de la gravedad, y los otros dos completando el sistema de coordenadas por regla de la mano derecha.

  Por otra parte, podríamos ofrecer una definición del sistema en cuanto a tan solo dos coordenadas angulares, $theta in [0, pi)$ y $phi in [0, 2 pi]$, eliminando la ecuación de restricción.

] <ej:conf_pendulo>

El espacio de configuración no es suficiente si consideramos la evolución temporal del sistema, ya que no representa las velocidades. Asumiendo que $Q$ es una variedad diferenciable, el espacio de estados se define entonces como el fibrado de los espacios tangentes a $Q$, $T Q$ #cite_mech(53).

#example[Espacio de estados de un péndulo esférico][

] <ej:estados_pendulo>

Con estos dos ejemplos para el péndulo, estamos listos para plantear el espacio de estados de la partícula en órbita. En particular, este espacio estará formado por

Una primera

== Mecánica orbital Hamiltoniana
