require(phaseR)

dy = function(t,y,parameters){
  k = parameters[1]
  k * y
}
expdecay  <- flowField(dy,xlim = c(0,2), ylim = c(0, 2),
                       parameters = -1,
                       points     = 21,
                       system     = "one.dim",
                       add        = FALSE, col = "blue")


