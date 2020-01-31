library(magrittr)
n = 10000
normDat = matrix(rnorm(n*2), ncol = 2) %>% scale()
expDat = matrix(sample(c(-1,1), 2*n, T)*rexp(n*2)^1.3, ncol = 2) %>% scale()
s = matrix(c(1, -2,-3,1), ncol = 2) %>% svd()
R = s$u
rexpDat = (expDat %*% R)
sc = 3



pdf(file = "rotInv.pdf", width = 10, height = 3.5)
par(mfrow =c(1,3), mar = c(.1,1,4,1), pty="s")
plot(normDat, pch = ".", xlim = sc*c(-1,1), ylim=sc*c(-1,1), axes=F, ylab = "", xlab="", 
     main = "Gaussian is\nrotationally invariant.", cex.main = 2)
arrows(x0 = c(-3,0),y0 = c(0,-3), x1 = c(3,0),y1 = c(0,3), code  =3, col = "grey43", lwd = 2)
# lines(c(-99,099), c(0,0), col = "grey43", lwd = 2)
# lines(c(0,0), c(-99,099), col = "grey43", lwd = 2)
points(normDat, pch = ".")
points(0,0, col = "red", cex=2, pch  =19)



plot(rexpDat, pch = ".", xlim = sc*c(-1,1), ylim=sc*c(-1,1), axes=F, ylab = "", xlab="", 
     main = "This non-Gaussian\nhas radial streaks.", cex.main = 2)
arrows(x0 = c(-3,0),y0 = c(0,-3), x1 = c(3,0),y1 = c(0,3), code  =3, col = "grey43", lwd = 2)
# lines(c(-99,099), c(0,0), col = "grey43", lwd = 2)
# lines(c(0,0), c(-99,099), col = "grey43", lwd = 2)
points(rexpDat, pch = ".")
points(0,0, col = "red", cex=2, pch  =19)


plot(rexpDat, pch = ".", xlim = sc*c(-1,1), ylim=sc*c(-1,1), axes=F, ylab = "", xlab="", 
     main = "Varimax correctly\nestimates the basis.", cex.main = 2)
v = varimax(rexpDat)
stdBasis = matrix(c(1,0,-1,0,0,1,0,-1), byrow=T, ncol = 2)
rb = stdBasis%*%R*3.5
arrows(x0 = rb[c(1,3),1],y0 = rb[c(1,3),2], x1 = rb[-c(1,3),1],y1 = rb[-c(1,3),2], code  =3, col = "grey43", lwd = 2)
# lines(x = rb[1:2,1], y = rb[1:2,2], col = "grey43", lwd = 2)
# lines(x = rb[3:4,1], y = rb[3:4,2], col = "grey43", lwd = 2)
points(rexpDat, pch = ".")
points(0,0, col = "red", cex=2, pch  =19)
dev.off()
