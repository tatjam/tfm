#import "@preview/frame-it:1.2.0": *
#import "@preview/ctheorems:1.1.3": *

#let (example) = frame(kind: "example", "Ejemplo", blue.lighten(60%))

#let theorem = thmbox("theorem", "Theorem", fill: rgb("#eeffee"))
#let corollary = thmplain(
  "corollary",
  "Corollary",
  base: "theorem",
  titlefmt: strong,
)

#let definition = thmbox(
  "definition",
  "Definición",
  inset: (
    x: 0em,
    top: 0em,
  ),
)

#let proof = thmproof("proof", "Proof")

