// Copyright (C) 2026 José Antonio Mayo García
// SPDX-License-Identifier: GPL-3.0-or-later

#import "@preview/ilm:1.4.2": *
#import "utils.typ": *
#import "@preview/ctheorems:1.1.3": *

#show: thmrules

#show: frame-style(kind: "example", styles.thmbox)
#show: frame-style(kind: "note", styles.thmbox)

#set text(lang: "es")
#set text(font: "Noto Sans")
#show math.equation: set text(font: "Fira Math", fallback: true)
#set quote(block: true)

#show figure: set block(breakable: true)

#show: ilm.with(
  title: [TFM (CAMBIAR PORTADA)],
  author: "José Antonio Mayo García",
  bibliography: bibliography("refs.bib"),
  figure-index: (enabled: true),
  table-index: (enabled: false),
  listing-index: (enabled: false),
)

= Fundamentos geométricos

#include "geometry.typ"

= Mecánica geométrica

#include "mechanics.typ"


= Propagación de incertidumbre

#include "uncertainty_propagation.typ"
