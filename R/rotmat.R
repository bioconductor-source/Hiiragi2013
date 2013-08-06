rotmat <-
function(phi)
  cbind( c( cos(phi),-sin(phi)),
         c(+sin(phi), cos(phi)) )
