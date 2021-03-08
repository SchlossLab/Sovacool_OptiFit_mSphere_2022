library(DiagrammeR)

# Pat's MCC function 
# source: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017/blob/a8bc26855423bba85acc0b8e7cca075e5c94f533/submission/supplemental_text.Rmd#L26-L28
mcc <- function(tp, tn, fp, fn) {
   format(round((tp * tn - fp * fn) /	sqrt((tp + fp) * (tp + fn) * (tn + fp) *
                                              (tn + fn)), digits = 2),
          digits = 2,
          nsmall = 2L)
}

dat <- data.frame(`Seq. 1` = c("D", "F", "G", "H", "I", "I", "J", "J", "N", "O", "P", "P", "P", "Q", "Q"), 
                  `Seq. 2` = c("B", "E", "C", "A", "B", "D", "A", "H", "M", "L", "K", "L", "O", "E", "F"), 
                  Distance = c(0.024, 0.028, 0.028, 0.027, 0.016, 0.024, 0.028, 0.024, 0.024, 0.024, 0.016, 0.027, 0.027, 0.024, 0.028))

otus <- list(c('A', 'B'),
             c('C', 'D', 'E'),
             c('F'))

get_otus <- function(idx) {
  x <- otus[[idx]]
  paste0(x, collapse = ",")
}

grViz("
digraph algorithm {

graph [layout = dot,
       overlap = true]

{
   rank = same
   node [shape = circle]
   1 [label = '@@1'] 
   2 [label = '@@2'] 
   3 [label = '@@3'] 
}

1 -> 3 [dir = back; label = 0.7]
2 -> 3 [dir = back]
3 -> 3
}

[1]: get_otus(1)
[2]: get_otus(2)
[3]: get_otus(3)
", height = 100)
