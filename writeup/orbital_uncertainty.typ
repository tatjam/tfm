#import "@preview/physica:0.9.8": *
#import "utils.typ": *

== Casos de estudio

Con el fin de evitar sobre-ajustar el método a un caso particular, estudiaremos tres casos con cierto realismo:

- Una órbita LEO prácticamente circular.
- Una órbita Molniya, en su perigeo y apogeo.

Tenemos los siguientes elementos MEE:

#align(
  center,
  table(
    columns: (auto, auto, auto, auto),
    inset: 10pt,
    align: horizon,
    table.header([*Coord.*], [*LEO*], [*Perigeo*], [*Apogeo*]),
    [p (km)], [1], [2], [3],
    [f], [1], [2], [3],
    [g], [1], [2], [3],
    [h], [1], [2], [3],
    [k], [1], [2], [3],
    [L (rad)], [1], [2], [3],
  ),
)

En cuanto a las incertidumbres, consideraremos a su vez tres casos. Con el fin de evitar repetirnos excesivamente, escribiremos estos casos en coordenadas RTN escaladas por el radio y velocidad orbital del punto de estudio, tal que

$
  Sigma = R S med overline(Sigma) med S^TT R^TT
$

siendo $R$ la matriz que lleva de ejes RTN a inerciales, $overline(Sigma)$ las matrices de la tabla que sigue y $S = diag(r, r, r, sqrt(mu / r), sqrt(mu / r), sqrt(mu / r))$ es la matriz de escala de cada caso de estudio.

Tenemos entonces

$
  overline(Sigma)_1 = mat(
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
  )
$
$
  overline(Sigma)_2 = mat(
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
  )
$
$
  overline(Sigma)_3 = mat(
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
    , , , , , ;
  )
$

Los casos buscan representar, respectivamente:
- Un caso de medición on-board con GPS, isotrópico y preciso
- Un caso de medición desde tierra con mayor error along-track
- Un caso con acoplamiento de incertidumbres y error más elevado


#text(red)[TODO: Visualizción]

Un estudio más detallado debería considerar, o bien un mayor número de matrices, o bien un estudio de sensibilidad en cada uno de sus grados de libertad. Ya que los métodos aplicados son costosos computacionalmente, nos limitamos a estudiar este total de 9 casos.

En cualquier caso, el cuerpo central es la Tierra, asumiendo la constante gravitacional $mu =$ y con única perturbación $J_2 =$. Un estudio más detallado debería considerar otras fuentes de no linealidad, como son los harmónicos gravitacionales de orden superior, o fuerzas no conservativas como el rozamiento con la atmósfera.

== Métodos en coordenadas euclídeas

Los métodos Monte Carlo, STM, UT, y PCE son directamente aplicables al problema orbital en coordenadas euclídeas. El método GVM carece de sentido, al no disponer de una coordenada angular.

=== Parámetros para el método UT

=== Parámetros para el método PCE

== Distribución de partida en MEE

Uno de los primeros problemas que surgen al estudiar la incertidumbre en coordenadas MEE es que la incertidumbre sobre la posición del objeto, generalmente, se da como una Gaussiana en coordenadas euclídeas. Por lo tanto, esta misma incertidumbre no es gaussiana en coordenadas MEE, y deberemos realizar algún tipo de transformación.

En general, es deseable que un muestro de la distribución de partida en MEE sea máximamente similar, tras realizar el mapa MEE a euclídeas, a la distribución original.

Los métodos basados en distribuciones gaussianas necesitan una gaussiana de partida. Para este primer caso, una opción bastante razonable es aplicar los métodos STM y UT al mapa que nos lleva de coordenadas euclídeas a MEE, obteniendo así una distribución gaussiana en MEE.

Por otra parte, el método GVM parte de una distribución Gauss-Von Mises. En este caso, procederemos a partir de la gaussiana previamente obtenida, y obtendremos por ajuste de momentos la Gauss-Von Mises más cercana.

El método PCE es el más directo, ya que no necesita una distribución de partida: el mapa euclídeo a MEE simplemente se abosorbe en la función de estudio y su no linealidad será capturada por la proyección de Galerkin.

En cualquier caso, utilizaremos el previamente descrito método de las energías para cuantificar la similitud entre la distribución de partida en coordenadas euclídeas y MEE. Para ello, para cada punto del caso de estudio, indicaremos el p-valor obtenido.

=== Distribución de partida usando STM

Presentamos a continuación los p-valores de la prueba de la energía para la distribución tras la transformación STM en cada caso de estudio:

#align(
  center,
  table(
    columns: (auto, auto, auto, auto),
    inset: 10pt,
    align: horizon,
    table.header([*Caso*], [*Matriz 1*], [*Matriz 2*], [*Matriz 3*]),
    [LEO], [1], [1], [1],
    [Perigeo], [1], [1], [1],
    [Apogeo], [1], [1], [1],
  ),
)


=== Distribución de partida usando UT

Al contrario del método anterior, utilizar el método UT implica la elección de los valores $alpha$, $kappa$ y $beta$. En el caso de $beta$, ya que la distribución de partida es gaussiana, tomamos $beta = 2$, que es el valor óptimo #text(red)[TODO: fuente].

Estudiaremos la variación de $alpha$ y $kappa$ en los rangos $[0.0001, 1]$ y $[-3, 3]$ respectivamente, para cada caso.

== UT

== PCE
