# The metric is:
# G_c(h,k) = \int_{S^1} (...) ds, where
# (...) = pH0 * <h,n>^2
#       + pKappa * kappa^2*<h,n>^2

# MODEL COEFFICIENTS
param pHor >=0 <=1;   # pHor=1 means only horizontal energy, pHor=0 means total energy
param pH0  >=0;       # the weight in front of the H^0 term of the metric
param pKappa >= 0;    # the weight in front of the kappa term of the metric
param pH2 = 0;    # no H2-term
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
var dmiddle {t in KT diff {1,nKT}, u in KU, i in 1..2} := 
  ((nKT-t)*d0[u,i]+(t-1)*d1[u,i])/(nKT-1);
var d {t in KT, u in KU, i in 1..2} = 
  (if t=1 then d0[u,i] else if t=nKT then d1[u,i] else dmiddle[t,u,i]);

# CURVE AND DERIVATIVES
var c {t in QT, u in QU, i in 1..2} =
  sum {(t,u,tt,uu) in iB} B[t,u,tt,uu]*d[tt,uu,i];
var ct {t in QT, u in QU, i in 1..2} =
  sum {(t,u,tt,uu) in iBt} Bt[t,u,tt,uu]*d[tt,uu,i];
var cu {t in QT, u in QU, i in 1..2} =
  sum {(t,u,tt,uu) in iBu} Bu[t,u,tt,uu]*d[tt,uu,i];
var cuu {t in QT, u in QU, i in 1..2} =
  sum {(t,u,tt,uu) in iBuu} Buu[t,u,tt,uu]*d[tt,uu,i];
var norm_cu {t in QT, u in QU} = sqrt(sum {i in 1..2} cu[t,u,i]^2);
var l {t in QT} = sum {u in QU} norm_cu[t,u];

#norm of velocity
var ctsquare {t in QT, u in QU} = (ct[t,u,1]*ct[t,u,1]+ ct[t,u,1]*cu[t,u,2]);

# normal part of velocity 
var a {t in QT, u in QU} = (ct[t,u,1]*cu[t,u,2]- ct[t,u,2]*cu[t,u,1])/norm_cu[t,u];

# penalty
var pen = sum {t in QT, u in QU} W[t,u]*(norm_cu[t,u]-l[t]/nQU)^2;

# energy
var energy = sum {t in QT, u in QU} (W[t,u] *
  (pH0 + pKappa * (cu[t,u,1]*cuu[t,u,2]-cu[t,u,2]*cuu[t,u,1])^2/norm_cu[t,u]^6)* a[t,u]^2*norm_cu[t,u])
;




# MINIMIZATION
minimize f:
 (energy);



