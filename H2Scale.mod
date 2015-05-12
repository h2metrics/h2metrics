# H2 Metric
# Letting $c_t = h$, the horizontal energy is: 
# $\int_0^1 \int_{S^1} ell^{-1}\langle h,h\rangle+ell \langle D^2_sh, D^2_sh \rangle +ell^3 \langle D^2_sh, D^2_sh \rangle ds dt$

# PARAMETERS 
param Pi default 3.141592653589793;
param nKP integer; # number of p-knots
param nKT integer; # number of t-knots
param nQP integer; # number of p-quadrature points
param nQT integer; # number of t-quadrature points
param L2coef; # L2-coefficient of the metric
param H1coef; # H1-coefficient of the metric
param H2coef; # H2-coefficient of the metric
param minTrans;

# KNOTS AND QUADRATURE POINTS
set KP := 1..nKP circular; # knots in the p-dimension
set KT := 1..nKT;          # knots in the t-dimension
set QP := 1..nQP circular; # quadrature points in the p-dimension
set QT := 1..nQT;          # quadrature points in the t-dimension

# COLLOCATION MATRICES
set iB        within  {QT,QP,KT,KP}; # the indices of non-zero entries of B 
set iBp       within  {QT,QP,KT,KP}; # the indices of non-zero entries of Bp 
set iBt       within  {QT,QP,KT,KP}; # the indices of non-zero entries of Bt
set iBpt      within  {QT,QP,KT,KP}; # the indices of non-zero entries of Bpt 
set iBpp      within  {QT,QP,KT,KP}; # the indices of non-zero entries of Bpp 
set iBppt     within  {QT,QP,KT,KP}; # the indices of non-zero entries of Bppt 
param B       {iB}    default 0;
param Bp      {iBp}   default 0;
param Bt      {iBt}   default 0;
param Bpt     {iBpt}  default 0;
param Bpp     {iBpp}  default 0;
param Bppt    {iBppt} default 0;
param W       {QT,QP} default 0;



#Rotations and Translations
var Trans {1..2} := 0;
var Scale:=1;

# CONTROLS
param d0 {KP,1..2};
param d1 {KP,1..2};
var dend  {p in KP,i in 1..2} = (if i=1 then Scale*d1[p,i]+minTrans*Trans[1] else Scale*d1[p,i]+minTrans*Trans[2]);
var dmiddle {t in KT diff {1,nKT}, p in KP, i in 1..2} := 
  ((nKT-t)*d0[p,i]+(t-1)*d1[p,i])/(nKT-1);
var d {t in KT, p in KP, i in 1..2} =  
  (if t=1 then d0[p,i] else if t=nKT then dend[p,i] else dmiddle[t,p,i]);

# CURVE AND DERIVATIVES
var c {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iB} B[t,p,tt,pp]*d[tt,pp,i];
var cp {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iBp} Bp[t,p,tt,pp]*d[tt,pp,i];
var ct {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iBt} Bt[t,p,tt,pp]*d[tt,pp,i];
var cpt {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iBpt} Bpt[t,p,tt,pp]*d[tt,pp,i];
var cpp {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iBpp} Bpp[t,p,tt,pp]*d[tt,pp,i];
var cppt {t in QT, p in QP, i in 1..2} =
  sum {(t,p,tt,pp) in iBppt} Bppt[t,p,tt,pp]*d[tt,pp,i];
var norm_cp {t in QT, p in QP} = sqrt(sum {i in 1..2} cp[t,p,i]^2);
var length {t in QT} = sum {p in QP} W[t,p] *norm_cp[t,p];

# MINIMIZATION

minimize energy: sum {t in QT, p in QP}
  W[t,p] *
  (
    L2coef/length[t]*(sum {i in 1..2} ct[t,p,i]^2) 
    + H1coef*length[t]*(
	sum {i in 1..2} cpt[t,p,i]^2/norm_cp[t,p]^2) 
    + H2coef*length[t]^3*(
	sum {i in 1..2} (cppt[t,p,i]/norm_cp[t,p]^2 
      - 
	(sum {j in 1..2} cp[t,p,j]*cpp[t,p,j])*cpt[t,p,i]/norm_cp[t,p]^4
			 )^2)
  ) * norm_cp[t,p];