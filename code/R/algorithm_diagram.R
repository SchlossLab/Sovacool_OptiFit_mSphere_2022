library(DiagrammeR)

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

2 -> 3
1 -> 3 [dir = back]
}

[1]: get_otus(1)
[2]: get_otus(2)
[3]: get_otus(3)
", height = 100)
