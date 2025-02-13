#-
#: IncDir procedures
Off Statistics;
format 160;

#include declarations.h
.sort
#message "Computing Cross section for g(p1), g(p2) -> H(p3), g(p4)"
l iM = + E1(mu)*E2(nu)*VVV(p1, -p1 - p2, p2, mu, lori1, nu, c1, cOli1, c2)*DV(p1 + p2, lori1, lori2)*VEffGGH(p1 + p2, -p4, lori2, rho, cOli1, c4)*EStar4( rho)
  + E1(mu)*E2(nu)*VVV(p1, -p4, - (p1 - p4), mu, rho, lori1, c1, c4, cOli1)*DV(p1 - p4, lori1, lori2)*VEffGGH(p1 - p4, p2, lori2, nu, cOli1, c2)*EStar4( rho)
  + E1(mu)*E2(nu)*VVV(p2, -p4, - (p2 - p4), nu, rho, lori1, c2, c4, cOli1)*DV(p2 - p4, lori1, lori2)*VEffGGH(p2 - p4, p1, lori2, mu, cOli1, c1)*EStar4( rho)
  + E1(mu)*E2(nu)*VEffGGGH(p1, p2, -p4, mu, nu, rho, c1, c2, c4)*EStar4( rho);

l WardIdentity = iM;
sum cOli1,...,cOli10;
sum lori1,...,lori10;
.sort
if(expression(WardIdentity));
  id E4(mu?) = p4(mu);
  sum mu, nu, rho;
endif;
id VEffGGGH(p1?, p2?, p3?, mu?, nu?, rho?, c1?, c2?, c3?) = C0/v*VVV(p1, p2, p3, mu, nu, rho, c1, c2, c3);

id DV(k?, lori1?, lori2?) = i_*(-d_(lori1, lori2))*DS(k);
id VVV(k?, p1?, p2?, lori1?, lori2?, lori3?, cOli1?, cOli2?, cOli3?)
  = g*cOlf(cOli1, cOli2, cOli3)*( (k(lori3) - p1(lori3))*d_(lori1, lori2) + (p1(lori1) - p2(lori1))*d_(lori2, lori3) + (p2(lori2) - k(lori2))*d_(lori1, lori3) );

id VEffGGH(p1?, p2?, lori1?, lori2?, c1?, c2?) = i_*d_(c1, c2)*(d_(lori1, lori2)*p1.p2 - p1(lori2)*p2(lori1))/v*C0;


id DS(p?) = den(p.p);

.sort
#call mandelstams
argument den;
#call mandelstams
endargument;
id den(s?) = 1/s;
#do i=1,4
  id E`i'(p`i') = 0;
  id EStar`i'(p`i') = 0;
#enddo
b alphas, g, CA, NA, ep, den, C0, v, pi_;
print+s WardIdentity;
.sort
l iMConj = iM*replace_(N1_?, cOli1, N2_?, cOli3, E1, EStar1, E2, EStar2, EStar4, E4, i_, -i_, mu, alpha, nu, beta, rho, gamma);
.sort
l M2 = iM*iMConj;
l DawsonCheck = alphas^3/v^2*(32/3/pi_)*((mH^8 + s^4 + t^4 + u^4)/s/t/u*(1 - 2*ep)
  + ep/2*(mH^4 + s^2 + t^2 + u^2)^2/s/t/u) - M2;
.sort
drop iM, iMConj;


#do i=1,4
  id EStar`i'(mu?)*E`i'(nu?) = -d_(mu, nu);
#enddo

id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli1?, cOli2?, cOli4?) = CA*cOld(cOli3, cOli4);
id cOld(a?, b?)*cOld(a?, b?) = NA;
id cOld(a?, a?) = NA;

#call mandelstams
mul replace_(d, 4 - 2*ep);

id ep*den(1 - ep) = -1 + den(1 - ep);
id g^2 = alphas*4*pi_;
id C0 = alphas/4/pi_*(-4/3);
mul replace_(NA, 8, CA, 3);

id mH^2 = s + t + u;

format mathematica;
b alphas, g, CA, NA, ep, den, C0, v, pi_;
*ab marker;
print+s ;
.end