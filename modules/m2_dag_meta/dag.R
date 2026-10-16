# DAG identification with the genuine dagitty engine (jtextor/dagitty, R).
#
# Replaces modules/m2_dag_meta/dag.py (pydagitty has no verifiable Python
# API; dag.py is kept, deprecated, for the day one materializes).
# Two graphs, same nodes as dag.py:
#   G1: M/C/T confound the E->S effect (adjustment set expected).
#   G2: study inclusion I descends from E and S (selection collider).
# Reports: adjustmentSets() + d-connected open paths E-S per graph.
#
# Input:  none
# Output: results/m2/dag_adjustment.txt

suppressPackageStartupMessages({
  library(dagitty)
})

g1 <- dagitty("dag {
  E [exposure]
  S [outcome]
  M -> E
  M -> S
  C -> E
  C -> S
  T -> E
  T -> S
  E -> S
}")

g2 <- dagitty("dag {
  E [exposure]
  S [outcome]
  M -> E
  M -> S
  C -> E
  C -> S
  T -> E
  T -> S
  E -> S
  E -> I
  S -> I
}")

report <- function(g, label) {
  cat("---", label, "---\n")
  cat("Adjustment sets:\n")
  print(adjustmentSets(g))
  cat("E-S d-connected with empty set:", dconnected(g, "E", "S", c()), "\n")
  cat("Open paths E-S:\n")
  print(paths(g, "E", "S"))
  cat("\n")
}

dir.create("results/m2", showWarnings = FALSE, recursive = TRUE)
sink("results/m2/dag_adjustment.txt")
report(g1, "G1 confounding graph")
report(g2, "G2 selection-collider graph (I descends from E, S)")
sink()
cat("Wrote results/m2/dag_adjustment.txt\n")
