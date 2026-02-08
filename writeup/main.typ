// Copyright (C) 2026 José Antonio Mayo García
// SPDX-License-Identifier: GPL-3.0-or-later

#import "@preview/ilm:1.4.2": *
#import "utils.typ": *

#show: frame-style(kind: "example", styles.thmbox)
#show: frame-style(kind: "definition", styles.thmbox)
#show: thmrules

#set text(lang: "es")
#set text(font: "Noto Sans")
#show math.equation: set text(font: "Fira Math", fallback: true)
#set quote(block: true)



#show: ilm.with(
  title: [TFM (CAMBIAR PORTADA)],
  author: "José Antonio Mayo García",
  bibliography: bibliography("refs.bib"),
  figure-index: (enabled: true),
  table-index: (enabled: false),
  listing-index: (enabled: false),
)

= Motivación

#include "motiv.typ"




