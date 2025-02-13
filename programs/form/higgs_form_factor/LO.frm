#-
#: IncDir procedures
Off Statistics;
format 160;

#include declarations.h
.sort

l iM = -1/16/pi_^4*(-i_*mq/v*cOld(cOli1, cOli2))*i_*(g_(0, k) + g_(0, p1) + g_(0, p2) + mq)*DS(k + p1 + p2, mq)*(i_*g*g_(0, nu)*cOlT(a, cOli1, cOli3))
  *i_*(g_(0, k) + g_(0, p1) + mq)*DS(k + p1, mq)*(i_*g*g_(0, mu)*cOlT(b, cOli3, cOli2))*i_*(g_(0, k) + mq)*DS(k, mq)
  *(1 + replace_(p1, p2, p2, p1, mu, nu, nu, mu, a, b, b, a, k, -k - p1 - p2));

id DS(-k, mq?) = DS(k, mq);
id DS(-k - p1, mq?) = DS(k + p1, mq);
id DS(-k - p1 - p2, mq?) = DS(k + p1 + p2, mq);
id cOld(cOli1?, cOli2?)*cOlT(a?, cOli1?, cOli3?) = cOlT(a, cOli2, cOli3);
id cOlT(a?, cOli1?, cOli2?)*cOlT(b?, cOli2?, cOli1?) = TF*cOld(a, b);

tracen 0;
.sort
#call kinematics
.sort

l C = 4*pi_^2*v/i_/g^2*1/NA*cOld(a, b)*4/mH^4*den(d - 2)*(p2(mu)*p1(nu) - p1.p2*d_(mu, nu))*iM;
.sort
drop iM;
#call kinematics
id cOld(a?, b?)*cOld(a?, b?) = NA;

*partial fractioning
repeat;
  id once p2.k*DS(k + p1 + p2, mq?) = 1/2*(1 - DS(k + p1 + p2, mq)*(k.k + 2*p1.k + 2*p1.p2 - mq^2));
  id once p1.k*DS(k + p1, mq?) = 1/2*(1 - DS(k + p1, mq)*(k.k - mq^2));
  id once k.k*DS(k, mq?) = 1 + mq^2*DS(k, mq);
  #call kinematics
endrepeat;

if(match(DS(k, mq))==0);
  id DS(k + p1, mq) = DS(k + p1, mq)*replace_(k, k - p1);
endif;

if(count(DS, 1)==2);
* using the symmetry after k -> -k - p1 - p2 we find
  id p1.k*DS(k, mq)*DS(k + p1 + p2, mq) = -1/2*p1.p2*DS(k, mq)*DS(k + p1 + p2, mq);
endif;

#call kinematics;

*b k;
ab d;
print+s C;
.sort

if(count(k, 1)!=0);
  print "%t";
  exit "Found k in the numerator";
endif;

mul i_*pi_^2; * integral measure reads 1/(i_*pi_^{d/2}) [epsilon dimenisioality not tracked]
* Invariance with respect to on-shell momentum. The integral can only depend on mq^2
if(count(DS, 1)==2)
  id DS(k, mq)*DS(k + p2, mq) = DS(k, mq)*DS(k + p1, mq);

mul replace_(d, 4 - 2*ep);
if(count(DS, 1) == 2) id DS(k,mq)*DS(p1 + p2 + k,mq) = (1/ep + 2 - log(mq^2) - be*log((1 + be)*den(1 - be)));
if(count(DS, 1) == 3) id DS(k,mq)*DS(p1 + k,mq)*DS(p1 + p2 + k,mq) = (1/2/mH^2*log(-z*den(1 - z))^2*replace_(z, 1/2*(1 + be)));
#call epExpand
if(count(ep, 1) > 0)discard;
argument log;
  id den(1/2 - 1/2*be) = 2*den(1 - be);
  id be*den(1 - be) = -1 + den(1 - be);
endargument;

id mq^2 = mH^2/4/z;


format mathematica;
b DS, cOld, TF, g, v;
print+s ;
.end