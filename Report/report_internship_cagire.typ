// Test equation
// ------------------- config -------------------

#set page(margin:2.5cm)

#set text(font:"New Computer Modern", size:12pt)

#set par(justify: true)
#set text(hyphenate: false)

#show figure: set block(
  spacing: 2em,
)
// outline/table of content
#import "@preview/outrageous:0.4.0"
#show outline.entry: outrageous.show-entry


// ------------------- headings -------------------
#set heading(numbering: "1.")

#show heading: it => {
  if it.level == 1 {v(1em)} else {v(0.5em)}
  it
  if it.level == 1 {v(1em)} else {v(0.5em)}
}
// ------------------- equations -------------------

#set math.equation(numbering: "(1.)")
#set math.equation(supplement: none, numbering: it => {numbering("(1.1)", counter(heading).get().first(), it)})

//#set math.equation(numbering:none)
#show math.equation: set text(font: "New Computer Modern Math")


// ------------------- front matter -------------------

#page[#include "chapters/0-Cover.typ"]

#page[#include "chapters/0b-Abstract-ToC.typ"]

// ------------------- main content -------------------
#set page(numbering: "1")
#counter(page).update(1)
#counter(heading).update(0)
#show link: set text(fill: blue)

#include "chapters/1-Introduction.typ"
#include "chapters/2-Solver.typ"
#include "chapters/3-Numerical_Instablities.typ"

// ------------------- bibliography -------------------

#page[#bibliography("bibliography.bib", style:"ieee")]
#set page(numbering:none)
//#include "chapters/4-Appendix.typ"
