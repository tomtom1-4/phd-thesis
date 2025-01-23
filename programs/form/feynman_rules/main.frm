#-
#: IncDir procedures
Off Statistics;
format 160;
#include declarations.h

g LG = -1/4*G(cOli1, mu, nu)*G(cOli1, mu, nu);
g LFP = -g*cOlf(cOli1, cOli2, cOli3)*D(CBar(cOli1), mu)*A(cOli3, mu)*C(cOli2);
g LM = g*PsiBar(cOli2)*g_(0, mu)*Psi(cOli3)*cOlT(cOli1, cOli2, cOli3)*A(cOli1, mu);
g LWeak = -1/4*sum_(a, 1, 3, W(a, mu, nu)*W(a, mu, nu)) - 1/4*B(mu, nu)*B(mu, nu);
g LH = DcovConj(Phi, mu)*Dcov(Phi, mu) -Lam*(Phi^4) + Mu^2*Phi^2;

sum cOli1, mu, nu;

#define Nper
#call definitions
id e_(1,2,3) = 1;
* Drop kinetic terms
if(count(D, 1) == 2) discard;
*if(count(v, 1) !=2) discard;
id sW = gY*sqrt_(1/(gY^2 + g2^2));
id cW = g2*sqrt_(1/(gY^2 + g2^2));
id Y = 1;
repeat id sqrt_(cW?)^2 = cW;
id gY^2*g2^2/(gY^2 + g2^2) = 1/2*((gY^2 + g2^2) - gY^4/(gY^2 + g2^2) - g2^4/(gY^2 + g2^2));

print+s LH;
.sort

*#call GFF.prc; * Gluon Fermion Fermion [psiBar(p1, c1), psi(p2, c2), A(p3, rho, c3)]
*#call GGG.prc; * Triple Gluon Vertex [A(p1, mu, c1), A(p2, nu, c2), A(p3, rho, c3)]
*#call GGGG.prc; * Quadruple Gluon Vertex [A(p1, mu, c1), A(p2, nu, c2), A(p3, rho, c3), A(p4, sigma, c4)]
*#call GCC.prc; * Gluon Ghost Ghost vertex [CBar(p1, c1), C(p2, c2), A(p3, rho, c3)]
#call ZZHH.prc;

.sort
l Vertex = i_*(L1+...+L`Nper'); *Factor i_ accounts for the factor in the exponent of the Gel-Mann--Low formula
renumber 1;
.sort
id ep1(cOli1?) = d_(mu, cOli1);
id ep2(cOli1?) = d_(nu, cOli1);
id ep3(cOli1?) = d_(rho, cOli1);
id ep4(cOli1?) = d_(sigma, cOli1);
b cOlf, g, g2, gY, sW, i_;
print+s ;
.end

