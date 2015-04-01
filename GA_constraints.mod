# The metric is:
# G_c(h,k) = \int_{S^1} (...) ds, where
# (...) = pH0 * <h,h>
#       + pKappa * kappa^2 * <h,k>
#       + pH2 * < D^2_s h, D^2_s k >

# MODEL COEFFICIENTS
param pHor >=0 <=1;   # pHor=1 means only horizontal energy, pHor=0 means total energy
param pH0  >=0;       # the weight in front of the H^0 term of the metric
param pKappa >= 0;    # the weight in front of the kappa term of the metric
param pH2  =0;       # no weight in front of the H^2 term of the metric
param nKU integer;    # number of u-knots
param nKT integer;    # number of t-knots 
param nQU integer;    # number of u-quadrature points
param nQT integer;    # number of t-quadrature points

# KNOTS AND QUADRATURE POINTS
set KT := 1..nKT;          # knots in the t-dimension
set KU := 1..nKU circular; # knots in the u-dimension
set QT := 1..nQT;          # quadrature points in the t-dimension
set QU := 1..nQU circular; # quadrature points in the u-dimension

# COLLOCATION MATRICES
set iB        within   {QT,QU,KT,KU}; # the indices of non-zero entries of B 
set iBt       within   {QT,QU,KT,KU}; # the indices of non-zero entries of Bt
set iBtu      within   {QT,QU,KT,KU}; # the indices of non-zero entries of But 
set iBtuu     within   {QT,QU,KT,KU}; # the indices of non-zero entries of Buut 
set iBu       within   {QT,QU,KT,KU}; # the indices of non-zero entries of Bu 
set iBuu      within   {QT,QU,KT,KU}; # the indices of non-zero entries of Buu 
set iBuuu     within   {QT,QU,KT,KU}; # the indices of non-zero entries of Buuu 
set iBU       within   {QU,KU};       # the indices of non-zero entries of BU
set iBUu      within   {QU,KU};       # the indices of non-zero entries of BUu
set iBUuu     within   {QU,KU};       # the indices of non-zero entries of BUuu
set iBUuuu    within   {QU,KU};       # the indices of non-zero entries of BUuuu
param B       {iB}     default 0;     # basis of tensor product splines in the t and u dimension
param Bt      {iBt}    default 0;     # basis of tensor product splines in the t and u dimension
param Btu     {iBtu}   default 0;     # basis of tensor product splines in the t and u dimension
param Btuu    {iBtuu}  default 0;     # basis of tensor product splines in the t and u dimension
param Bu      {iBu}    default 0;     # basis of tensor product splines in the t and u dimension
param Buu     {iBuu}   default 0;     # basis of tensor product splines in the t and u dimension
param Buuu    {iBuuu}  default 0;     # basis of tensor product splines in the t and u dimension
param BU      {QU,KU}  default 0;     # basis of splines in the u-dimension used for h and m
param BUu     {QU,KU}  default 0;     # basis of splines in the u-dimension used for h and m
param BUuu    {QU,KU}  default 0;     # basis of splines in the u-dimension used for h and m
param BUuuu   {QU,KU}  default 0;     # basis of splines in the u-dimension used for h and m
param W       {QT,QU}  default 0;     # weight matrix for quadrature

# CONTROLS OF CURVE
param d0 {KU,1..2};
param d1 {KU,1..2};
param dmiddle {t in KT diff {1,nKT}, u in KU, i in 1..2} := 
  ((nKT-t)*d0[u,i]+(t-1)*d1[u,i])/(nKT-1);
param d {t in KT, u in KU, i in 1..2} := 
  (if t=1 then d0[u,i] else if t=nKT then d1[u,i] else dmiddle[t,u,i]);

# CURVE AND DERIVATIVES
param c {t in QT, u in QU, i in 1..2} :=
  sum {(t,u,tt,uu) in iB} B[t,u,tt,uu]*d[tt,uu,i];
param ct {t in QT, u in QU, i in 1..2} :=
  sum {(t,u,tt,uu) in iBt} Bt[t,u,tt,uu]*d[tt,uu,i];
param cu {t in QT, u in QU, i in 1..2} :=
  sum {(t,u,tt,uu) in iBu} Bu[t,u,tt,uu]*d[tt,uu,i];
param cuu {t in QT, u in QU, i in 1..2} :=
  sum {(t,u,tt,uu) in iBuu} Buu[t,u,tt,uu]*d[tt,uu,i];
param cuuu {t in QT, u in QU, i in 1..2} :=
  sum {(t,u,tt,uu) in iBuuu} Buuu[t,u,tt,uu]*d[tt,uu,i];
param norm_cu {t in QT, u in QU} := sqrt(sum {i in 1..2} cu[t,u,i]^2);

# CONTROLS OF HORIZONTAL PROJECTION 
var k {t in QT, u in KU, i in 1..2}; 

# HORIZONTAL PROJECTION
var h {t in QT, u in QU, i in 1..2} = 
  sum {(u,uu) in iBU} BU[u,uu]*k[t,uu,i];


# VERTICAL TEST FUNCTIONS
set iM := iBU union iBUu union iBUuu;
param m {t in QT,u in QU,i in 1..2,uu in KU} := 
  cu[t,u,i]*BU[u,uu];
param mu {t in QT,u in QU,i in 1..2,uu in KU} := 
  cuu[t,u,i]*BU[u,uu] + cu[t,u,i]*BUu[u,uu];
param muu {t in QT,u in QU,i in 1..2,uu in KU} := 
  cuuu[t,u,i]*BU[u,uu] + 2*cuu[t,u,i]*BUu[u,uu] + cu[t,u,i]*BUuu[u,uu];


var energy = sum {t in QT, u in QU} W[t,u] *
  (
    pH0 * sum {i in 1..2} h[t,u,i]^2
    +
    pKappa * (sum {i in 1..2} h[t,u,i]^2) 
      * (cu[t,u,1]*cuu[t,u,2]-cu[t,u,2]*cuu[t,u,1])^2/norm_cu[t,u]^6
  ) * norm_cu[t,u]
;

# normal part of velocity 
param a {t in QT, u in QU} := (ct[t,u,1]*cu[t,u,2]- ct[t,u,2]*cu[t,u,1])/norm_cu[t,u];

# difference of projection to true solution
#var error = sum {t in QT, u in QU} W[t,u] *(a[t,u]^2-(sum {i in 1..2} h[t,u,i]^2))^2*norm_cu[t,u];

# tangential part of h (should be zero)
#var error = sum {t in QT, u in QU} (W[t,u] *(sum {i in 1..2} h[t,u,i]*cu[t,u,i])^2)*norm_cu[t,u];


#Constraint value
var error = 
sum {t in QT, uu in KU}(
sum {(u,uu) in iBU} (W[t,u]*
  (
    (ct[t,u,1]-h[t,u,1])*m[t,u,2,uu]-(ct[t,u,2]-h[t,u,1])*m[t,u,1,uu]
  ) * norm_cu[t,u]))^2;







# MINIMIZATION
minimize constant: 
1;


# CONSTRAINTS
# h equals pHor times the horizontal projection of ct 

# force h to be normal by requiring G_c(h,m)=0, 
# for all vertical test functions m, at each t-quadrature point
subject to horizontal {t in QT, uu in KU}: 
  sum {(u,uu) in iM} (W[t,u] * (
    pH0 + pKappa  * (cu[t,u,1]*cuu[t,u,2]-cu[t,u,2]*cuu[t,u,1])^2/norm_cu[t,u]^6) * 
  ( sum {i in 1..2} h[t,u,i]*m[t,u,i,uu]
  ) * norm_cu[t,u]) = 0;

# force ct minus h to be vertical by requiring \int <ct-h,rot(m)> ds = 0,
# for all test functions m, at each t-quadrature point
subject to vertical {t in QT, uu in KU}:
  sum {(u,uu) in iBU} (W[t,u]*
  (
    (ct[t,u,1]-h[t,u,1])*m[t,u,2,uu]-(ct[t,u,2]-h[t,u,1])*m[t,u,1,uu]
  ) * norm_cu[t,u]) = 0;

