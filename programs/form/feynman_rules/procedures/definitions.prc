#procedure definitions

* Gluon
repeat;
  id once G(cOli1?, mu?, nu?) = D(A(cOli1, nu), mu) - D(A(cOli1, mu), nu) + g*cOlf(cOli1, cOli10, cOli11)*A(cOli10, mu)*A(cOli11, nu);
  sum cOli10, cOli11;
endrepeat;

* W-Boson
repeat;
  id once W(a?, mu?, nu?) = D(W(a, nu), mu) - D(W(a, mu), nu) + sum_(c, 1, 3, sum_(b, 1, 3, g2*e_(a, b, c)*W(b, mu)*W(c, nu)));
endrepeat;

* B-Boson
repeat;
  id once B(mu?, nu?) = D(B(nu), mu) - D(B(mu), nu);
  sum b, c;
endrepeat;

id Dcov(Phi, mu?) = D(Phi, mu) - i_*g2*sum_(a, 1, 3, W(a, mu)*I(a))*Phi + i_*gY*Y/2*B(mu)*Phi;
id DcovConj(Phi, mu?) = D(Phi, mu) + i_*g2*sum_(a, 1, 3, W(a, mu)*I(a))*Phi - i_*gY*Y/2*B(mu)*Phi;

* H-Boson
id Phi = sqrt_(1/2)*(v + H);
id Lam = Mu^2/v^2;

id I(a?)*I(b?) = 1/4*d_(a,b) + i_*sum_(c, 1, 3, e_(a,b,c)*I(c))/2;
* Lagrangian is sandwiched between (0, 1/sqrt(2)*(v + H)) and (0, 1/sqrt(2)*(v + H))^T. We can therefore kill all off-diagonals
if(count(I, 1) == 1);
  if(match(I(1))) discard;
  if(match(I(2))) discard;
  id I(3) = -1/2;
endif;

* Change basis
id W(1, mu?) = (Wp(mu) + Wm(mu))*sqrt_(1/2);
id W(2, mu?) = i_*(Wp(mu) - Wm(mu))*sqrt_(1/2);
id W(3, mu?) = cW*Z(mu) - sW*Agamma(mu);
id B(mu?) = cW*Agamma(mu) + sW*Z(mu);

id D(W(1, mu?), nu?) = (D(Wp(mu), nu) + D(Wm(mu), nu))*sqrt_(1/2);
id D(W(2, mu?), nu?) = i_*(D(Wp(mu), nu) - D(Wm(mu), nu))*sqrt_(1/2);
id D(W(3, mu?), nu?) = cW*D(Z(mu), nu) - sW*D(Agamma(mu), nu);
id D(B(mu?), nu?) = cW*D(Agamma(mu), nu) + sW*D(Z(mu), nu);
#endprocedure