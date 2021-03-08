library(DiagrammeR)

dat <- data.frame(`Seq. 1` = c("D", "F", "G", "H", "I", "I", "J", "J", "N", "O", "P", "P", "P", "Q", "Q"), 
                  `Seq. 2` = c("B", "E", "C", "A", "B", "D", "A", "H", "M", "L", "K", "L", "O", "E", "F"), 
                  Distance = c(0.024, 0.028, 0.028, 0.027, 0.016, 0.024, 0.028, 0.024, 0.024, 0.024, 0.016, 0.027, 0.027, 0.024, 0.028))

otus <- list(1 = c('A','B'),
             2 = c('C', 'D', 'E'),
             3 = c('F'))

grViz("
digraph algorithm {

graph [layout = dot,
       overlap = true]

{
   rank = same
   node [shape = circle]
   1 [label = '@@1'] 
   2 
   3
}

2 -> 3
1 -> 3 [dir = back]

[1]: otus[[1]]
}
", height = 100)
