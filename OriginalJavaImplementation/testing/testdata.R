#
#  Clumpy
#
npts <- 150
R1 <- c(rnorm(npts), rnorm(npts), rnorm(npts))
R1 <- matrix(R1, nrow=npts, ncol=3)

R2 <- matrix(0, nrow=150, ncol=3)
R2[1:80,] <- matrix(20, nrow=80, ncol=3)
R2[1:80,2] <- 0
R2[81:150,] <- matrix(8, nrow=70, ncol=3)
C1 <- R1 + R2

R2[1:80,] <- matrix(12, nrow=80, ncol=3)
C2 <- R1 + R2

R2 <- matrix(0, nrow=150, ncol=3)
R2[1:50,] <- matrix(15, nrow=50, ncol=3)
R2[1:50,2] <- 0
R2[51:100,] <- matrix(8, nrow=50, ncol=3)
C3 <- R1 + R2

sc1 <- scagnostics3D(C1)
sc2 <- scagnostics3D(C2)
sc3 <- scagnostics3D(C3)
sc1
sc2
sc3
# 
# sc1
#
#    Outlying  7.223424e-02
#    Skewed  5.727241e-01
#    Clumpy  8.637553e-01
#    Sparse  6.251597e-05
#    Striated   5.664336e-01
#    Convex  3.069520e-01
#    Skinny  4.410324e-01 
#    Stringy  2.746250e-01
#    Monotonic  3.080593e-01
# 
 
# sc2
#    Outlying  1.190674e-01
#    Skewed  5.918202e-01 
#    Clumpy  7.956124e-01
#    Sparse  8.478095e-05 
#    Striated  4.068966e-01
#    Convex  8.013018e-01 
#    Skinny  2.617081e-01 
#    Stringy  2.226118e-01 
#    Monotonic  2.191561e-01 
#

# sc3 
#    Outlying  5.926901e-02 
#    Skewed  5.957959e-01 
#    Clumpy  5.711882e-01  
#    Sparse  6.755146e-05 
#    Striated  3.194444e-01 
#    Convex  6.078959e-01 
#    Skinny  4.166859e-01 
#    Stringy   2.962963e-01 
#    Monotonic   8.025560e-01 


# not-skinny
npts <- 100
CV <- c(rnorm(npts), rnorm(npts), rnorm(npts))
CV <- matrix(CV, nrow=npts, ncol=3)
s1 <- scagnostics3D(CV)
s1

#    Outlying  0.0741899087
#    Skewed    0.5313535237
#    Clumpy    0.0365686154
#    Sparse     0.0001685042
#    Striated  0.4421052632
#    Convex    0.7749889166
#    Skinny    0.1803681085
#    Stringy   0.2643510006
#    Monotonic 0.0468524826

#
# skinny - medium
#
d1 <- c(rnorm(100), rnorm(100), rnorm(100))
d1 <- matrix(d1, nrow=100, ncol=3)
svd <- svd(d1)
X2 <- svd$u %*% diag(c(100,10,1)) %*% svd$v

sk2 <- scagnostics3D(X2)
sk2

#    Outlying  7.611122e-02
#    Skewed    5.643034e-01
#    Clumpy   3.941309e-02
#    Sparse  6.643877e-05
#    Striated  3.263158e-01
#    Convex  8.214249e-01
#    Skinny 4.318609e-01
#    Stringy  3.348980e-01
#    Monotonic   9.111339e-01 

#
# skinny3, most
#

d1 <- c(rnorm(100), rnorm(100), rnorm(100))
d1 <- matrix(d1, nrow=100, ncol=3)
svd <- svd(d1)
X3 <- svd$u %*% diag(c(100,100,0.01)) %*% svd$v

sk3 <- scagnostics3D(X3)
sk3

#    Outlying  0.2508788891
#    Skewed  0.6840415712
#    Clumpy  0.0611772565
#    Sparse 0.0001212506
#    Striated 0.4555555556 
#    Convex  0.6358233687
#    Skinny  0.9379536765
#    Stringy 0.4814917182
#    Monotonic  0.8933471811

#
# convex - least 
#
pts = 33;
d <- matrix(0,nrow=3*pts+1, ncol=3)
d[1:pts,1] <- 0
d[1:pts,2] <- runif(pts)
d[1:pts,3] <- runif(pts)
b1 = pts+1
b2 = pts*2
d[b1:b2,1] <- runif(pts)
d[b1:b2,2] <- runif(pts)
d[b1:b2,3] <- 0
b3 = pts*2 +1
b4 = pts*3

d[b3:b4, 1] <- runif(pts)
d[b3:b4, 2] <- 0
d[b3:b4, 3] <- runif(pts)

lc<-scagnostics3D(d)
lc

#    Outlying 0.0000000000
#    Skewed   0.5522119633
#    Clumpy   0.0440931550
#    Sparse 0.0001879158
#    Striated  0.2857142857
#    Convex 0.1201119540
#    Skinny  0.5265370087
#    Stringy  0.4883730469
#    Monotonic  0.2291905846
 
#
# convex - medium
#
pts = 33;
d <- matrix(0,nrow=3*pts+1, ncol=3)
d[1:pts,1] <- 0
d[1:pts,2] <- runif(pts)
d[1:pts,3] <- runif(pts)
b1 = pts+1
b2 = pts*2
d[b1:b2,1] <- runif(pts)
d[b1:b2,2] <- runif(pts)
d[b1:b2,3] <- 0
b3 = pts*2 +1
b4 = pts*3

d[b3:b4, 1] <- runif(pts)
d[b3:b4, 2] <- runif(pts)
d[b3:b4, 3] <- runif(pts)

q<-scagnostics3D(d)
q

#    Outlying  0.0274204780 
#    Skewed  0.5772544032
#    Clumpy   0.0402523872
#    Sparse   0.0002222176 
#    Striated  0.4639175258 
#    Convex  0.4728369323
#    Skinny   0.4073486307
#    Stringy  0.3213227185
#    Monotonic  0.1366591646

 
#
# convex - most
#
npts <- 100
N1 <- c(rnorm(npts), rnorm(npts), rnorm(npts))
N1 <- matrix(N1, nrow=npts, ncol=3)
svd <- svd(N1)
N3 <- svd$u %*% diag(c(5,5,5)) %*% svd$v
sn3 <- scagnostics3D(N3)
sn3

#    Outlying 0.1044820190
#    Skewed   0.6220911825 
#    Clumpy  0.0287222207 
#    Sparse 0.0002126725
#    Striated  0.2795698925
#    Convex  0.8681868925 
#    Skinny   0.1732047307 
#    Stringy 0.3095613678
#    Monotonic  0.0204745811

 
#
# striated - many parallel planes
#
seed <- as.double(1)
RANDU <- function() {
    seed <<- ((2^16 + 3) * seed) %% (2^31)
    seed/(2^31)
}

npts <- 500

U <- matrix(nrow=npts, ncol=3)

for(i in 1:npts) {
    U[i,] <- c(RANDU(), RANDU(), RANDU())}

sa1 <- scagnostics3D(U)
sa1
    
#    Outlying  0.0042239445
#    Skewed    0.3050135154
#    Clumpy   0.0080589436
#    Sparse 0.0001200568
#    Striated  0.5010101010
#    Convex  0.9192159383
#    Skinny   0.2551027637
#    Stringy   0.2315152301 
#    Monotonic 0.0008410872
     
#
# striated - three parallel planes & 'lines'
#
d1 <- c(rnorm(80), rnorm(80), rnorm(80))
d1 <- matrix(d1, nrow=80, ncol=3)
svd <- svd(d1)
T1 <- svd$u %*% diag(c(10,1,1)) %*% svd$v
T2 <- svd$u %*% diag(c(100,1,100)) %*% svd$v
T3 <- svd$u %*% diag(c(100,10,1)) %*% svd$v

A <- matrix(0, nrow=240, ncol=3)
A[1:80,] <- T2
A[81:160,] <- T2+matrix(25, nrow=80,ncol=3)
A[161:240,] <- T2+matrix(10, nrow=80,ncol=3)

A1 <- A

A[1:80,] <- T2
A[81:160,] <- T2+matrix(10, nrow=80,ncol=3)
A2 <- A

st1 <- scagnostics3D(A1)
st2 <- scagnostics3D(A2)
st1
st2

#
##  sphere data
#
d <- matrix(rnorm(600),nrow=200, ncol=3)
dlen <- sqrt(diag(d%*%t(d)))
d <- d/dlen
sp <- scagnostics3D(d)
sp

#    Outlying  0.0000000000 
#    Skewed  0.4617564658
#    Clumpy  0.0236895500
#    Sparse  0.0001304427
#    Striated  0.2171717172
#    Convex  0.6026809852
#    Skinny   0.3916732924
#    Stringy   0.3744611561
#    Monotonic  0.0213962331
