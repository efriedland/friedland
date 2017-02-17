UnitRootApply <-
function(data, k, ...){
  t(sapply(lapply(data, UnitRoot, k = k, ...), function(x) x[2:7]))
}
