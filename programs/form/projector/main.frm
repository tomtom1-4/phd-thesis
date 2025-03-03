#-
#: IncDir procedures
Off Statistics;
format 160;

s d;
dimension d;

.sort
s s, t, u;
v p1, p2, p3;
auto i mu, nu, rho;
t met;
f V, U;
s F1, F2, F3, F4;



l t1 = g_(0, p3)*p1(rho1) - p1.p3*g_(0, rho1);
l t2 = g_(0, p3)*p2(rho1) - p2.p3*g_(0, rho1);
.sort
l t1Dagger = t1*replace_(rho1, rho2);
l t2Dagger = t2*replace_(rho1, rho2);
.sort
l t1t1 = t1Dagger*g_(0, p1)*met(p1, p3, rho1, rho2)*t1*g_(0, p2);
l t2t2 = t2Dagger*g_(0, p1)*met(p1, p3, rho1, rho2)*t2*g_(0, p2);
l t1t2 = t1Dagger*g_(0, p1)*met(p1, p3, rho1, rho2)*t2*g_(0, p2);
l t2t1 = t2Dagger*g_(0, p1)*met(p1, p3, rho1, rho2)*t1*g_(0, p2);
l test = (-d_(mu1, mu) + (p1(mu1)*p2(mu) + p1(mu)*p2(mu1))*2/s)
        *(-d_(nu1, nu) + (p2(nu1)*p3(nu) + p2(nu)*p3(nu1))*2/t)
        *(-d_(rho1, rho) + (p3(rho1)*p1(rho) + p3(rho)*p1(rho1))*2/u)
        *(d_(mu, nu)*p2(rho)*F1 + d_(mu, rho)*p1(nu)*F2 + d_(nu, rho)*p3(mu)*F3 + p3(mu)*p1(nu)*p2(rho)*F4);
.sort
sum rho;

id met(p1, p3, rho1?, rho2?) = -d_(rho1, rho2) + (p1(rho1)*p3(rho2) + p1(rho2)*p3(rho1))/t*2;
tracen 0;
*repeat;
*  id g_(0, rho?)*g_(0, p1) = -g_(0, p1)*g_(0, rho) + 2*p1(rho);
*  id g_(0, p2)*g_(0, rho?) = -g_(0, rho)*g_(0, p2) + 2*p2(rho);
*endrepeat;
*id g_(0, rho?, rho?) = d;
*id V*g_(0, p1) = 0;
*id g_(0, p2)*U = 0;
*id g_(0, p3?, p3?) = 0;

id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id p1.p2 = s/2;
id p1.p3 = t/2;
id p2.p3 = u/2;
id V*U = 1;
.sort
mul replace_(mu1, mu, nu1, nu, rho1, rho);
id p1(mu) = 0;
id p2(mu) = 0;
id p2(nu) = 0;
id p3(nu) = 0;
id p3(rho) = 0;
id p1(rho) = 0;


b p1, p2, p3, d_;
print+s t1t1, t1t2, t2t1, t2t2, test;
.end