#-
#: IncDir procedures
Off Statistics;
format 160;
s d;
dimension d;
.sort
v p;
f D;
s m, marker;
i mu, nu, rho, sig;

l num1 = ((2*i_*p(mu)*D(mu) + D(mu)*D(mu)*marker)*g_(0) - i_/2*i_/2*g_(0, mu, nu)*(D(mu)*D(nu) - D(nu)*D(mu)));
l num2 = num1*replace_(mu, rho, nu, sig);
.sort
l num = num1*num2;

sum mu, nu, rho, sig;
tracen 0;
.sort
id D(p)*D(p) = D(mu)*D(mu)*m^2/d;
sum mu;

print+s;
.sort

.end