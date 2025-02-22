#-
#: IncDir procedures
Off Statistics;
format 160;

#include declarations.h
.sort

*l iM0 = (i_*C0/v*(p1.p2*d_(mu, nu) - p1(nu)*p2(mu))*cOld(a, b));
l iM1 = (i_*(C0+C1)/v*(p1.p2*d_(mu, nu) - p1(nu)*p2(mu))*cOld(a, b))
  + 1/16/pi_^4*(DV(k, lori1, lori2)*VVV(k, p1, -p1 - k, lori2, mu, lori3, cOli1, a, cOli2)*DV(p1 + k, lori3, lori4)*VEffGGH(p1 + k, -k + p2, lori4, lori5)
    *DV(k - p2, lori5, lori6)*VVV(k - p2, p2, -k, lori6, nu, lori1, cOli2, b, cOli1))
  + 1/16/pi_^4*(DV(k + p1 + p2, lori1, lori2)*VEffGGH(k + p1 + p2, -k, lori2, lori3)*DV(k, lori3, lori4)*VVVV(k, p2, p1, -(k + p1 + p2), lori4, nu, mu, lori1, cOli1, b, a, cOli1))/2;

sum cOli1,...,cOli10;
sum lori1,...,lori10;

id DV(k?, lori1?, lori2?) = i_*(-d_(lori1, lori2))*DS(k);
id VVV(k?, p1?, p2?, lori1?, lori2?, lori3?, cOli1?, cOli2?, cOli3?)
  = g*cOlf(cOli1, cOli2, cOli3)*( (k(lori3) - p1(lori3))*d_(lori1, lori2) + (p1(lori1) - p2(lori1))*d_(lori2, lori3) + (p2(lori2) - k(lori2))*d_(lori1, lori3) );
id VVVV(p1?, p2?, p3?, p4?, mu?, nu?, rho?, sigma?, cOli1?, cOli2?, cOli3?, cOli4?)
  = (-i_)*g^2*(cOlf(cOli1, cOli2, cOli10)*cOlf(cOli3, cOli4, cOli10)*(d_(mu, rho)*d_(nu, sigma) - d_(mu, sigma)*d_(nu, rho))
              +cOlf(cOli1, cOli3, cOli10)*cOlf(cOli2, cOli4, cOli10)*(d_(mu, nu)*d_(rho, sigma) - d_(mu, sigma)*d_(nu, rho))
              +cOlf(cOli1, cOli4, cOli10)*cOlf(cOli2, cOli3, cOli10)*(d_(mu, nu)*d_(rho, sigma) - d_(mu, rho)*d_(nu, sigma)));
sum cOli10;
*id VEffGGH(p1?, p2?, lori1?, lori2?) = i_*(d_(lori1, lori2)*p1.p2 - p1(lori2)*p2(lori1) - p1(lori1)*p2(lori2))/v*C0;
id VEffGGH(k1?, k2?, lori1?, lori2?) = i_*k1.k2*(d_(lori1, lori2)*p1.p2 - p1(lori2)*p2(lori1) - p1(lori1)*p2(lori2))/(mH^2/2)/v*C0;

id DS(-k, mq?) = DS(k, mq);
id DS(-k - p1, mq?) = DS(k + p1, mq);
id DS(-k - p1 - p2, mq?) = DS(k + p1 + p2, mq);
id cOld(cOli1?, cOli2?)*cOlT(a?, cOli1?, cOli3?) = cOlT(a, cOli2, cOli3);
id cOlT(a?, cOli1?, cOli2?)*cOlT(b?, cOli2?, cOli1?) = TF*cOld(a, b);

.sort
#call kinematics
.sort

l C = pi_*v/i_/alphas*1/NA*cOld(a, b)*4/mH^4*den(d - 2)*(p2(mu)*p1(nu) - p1.p2*d_(mu, nu))*iM1;
.sort
drop iM1;
#call kinematics
id cOld(cOli1?, cOli2?)*cOlf(cOli2?, cOli3?, cOli4?) = cOlf(cOli1, cOli3, cOli4);
id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*cOld(cOli3, cOli4);
id cOld(a?, b?)*cOld(a?, b?) = NA;
id cOld(a?, a?) = NA;

if(match(DS(-p2 + k))) mul replace_(k, k + p2);

*partial fractioning
repeat;
  id once p2.k*DS(k + p1 + p2) = 1/2*(1 - DS(k + p1 + p2)*(k.k + 2*p1.k + 2*p1.p2));
  id once p1.k*DS(k + p1) = 1/2*(1 - DS(k + p1)*(k.k));
  id once p2.k*DS(k + p2) = 1/2*(1 - DS(k + p2)*(k.k));
  id once k.k*DS(k) = 1;
  if(match(DS(k))==0);
    if(count(DS, 1)==2) id DS(p2 + k)*DS(p1 + p2 + k) = DS(p2 + k)*DS(p1 + p2 + k)*replace_(k, -k - p1 - p2);
    if(count(DS, 1)==1);
      id DS(k + p1) = DS(k + p1)*replace_(k, k - p1);
      id DS(k + p2) = DS(k + p2)*replace_(k, k - p2);
      id DS(k + p1 + p2) = DS(k + p1 + p2)*replace_(k , k - p1 - p2);
    endif;
  endif;
  if(count(DS, 1)==3) id DS(k)*DS(p2 + k)*DS(p1 + p2 + k) = DS(k)*DS(p2 + k)*DS(p1 + p2 + k)*replace_(k, -k - p1 - p2);
  if(count(DS, 1)==2) id DS(k + p1)*DS(k + p1 + p2) = DS(k + p1)*DS(k + p1 + p2)*replace_(k, k - p1);
  id DS(-k) = DS(k);
  id DS(-k - p1) = DS(k + p1);
  id DS(-k - p1 - p2) = DS(k + p1 + p2);
endrepeat;

.sort
if(count(DS, 1)==1) id DS(k) = 0;
if(count(DS, 1)==2);
  id DS(k)*DS(k + p1) = 0;
  id DS(k)*DS(k + p2) = 0;
endif;

#call kinematics;

* IBP
if(count(k, 1)==0);
  if(count(DS, 1)==3)
    id DS(k)*DS(k + p1)*DS(k + p1 + p2) = DS(k)*DS(k + p1 + p2)*(-2*d+6)*den(d-4)/mH^2;
endif;

totensor k, kk;
if(count(DS, 1)==2);
  id kk(mu?)*DS(k)*DS(k + p1 + p2) = (p1(mu) + p2(mu))*(-1/2*DS(k)*DS(k + p1 + p2));
  id kk(mu?, nu?)*DS(k)*DS(k + p1 + p2) = (p1(mu) + p2(mu))*(p1(nu) + p2(nu))*d*den(d - 1)/mH^4*mH^4/4*DS(k)*DS(k + p1 + p2)
                                          + d_(mu, nu)*(-1)*den(d - 1)/mH^2*mH^4/4*DS(k)*DS(k + p1 + p2) ;
endif;

id d*den(-2 + d) = 1 + 2*den(-2 + d);
id d*den(-4 + d) = 1 + 4*den(-4 + d);

#call kinematics;

*b k;
b DS, k, kk;
print+s C;
.sort

if((count(k, 1)!=0)&&(count(kk, 1)!=0));
  print "%t";
  exit "Found k in the numerator";
endif;



if(occurs(DS)) mul i_*pi_^2; * integral measure reads 1/(i_*pi_^{d/2}) [epsilon dimenisioality not tracked]

mul replace_(d, 4 - 2*ep);
if(count(DS, 1) == 2) id DS(k)*DS(p1 + p2 + k) = 1/ep + (2 - log(-mH^2)) + ep*(4 - pi_^2/12 - 2*log(-mH^2) + 1/2*log(-mH^2)^2);

#call epExpand
if(count(ep, 1) > 0)discard;

id g^2 = alphas*4*pi_;

*id C0 = alphas/4/pi_*(-4/3);
*id C1 = (alphas/4/pi_)^2*(-44/3);

* For comparison with Dawson, 1991
*if(count(alphas, 1)==1) mul 2;
*mul 1/3*9*(1 - pi_^2/12*ep^2);

*id log(-mH^2) = log(mH^2) - i_*pi_;
*id log(-mH^2) = 0;
mul (1 + ep*log(-mH^2) + ep^2/2*log(-mH^2)^2);
*mul (1 + pi_^2/12*ep^2);
*mul replace_(CA, 3);
if(count(ep, 1) > 0)discard;
mul 3;

format mathematica;
b DS, cOld, TF, g, v, ep, log, alphas;
print+s ;
.end