#procedure GGH

#write "Effective Gluon-Gluon-Higgs vertex [A(p1, mu, c1), A(p2, nu, c2), H(p3)]"

if(count(g,1)!=0) discard;
if(count(D, 1)!=2) discard;
.store

#redefine Nper "2"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the gluons to be incoming

#$Ldummy = LHTL;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p1(nu)*ep1(mu)*d_(cOli1, c1);
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p2(nu)*ep2(mu)*d_(cOli1, c2);
#endif
#if `permutation'==2;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p2(nu)*ep2(mu)*d_(cOli1, c2);
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p1(nu)*ep1(mu)*d_(cOli1, c1);
#endif
id once H = 1;

.store
.sort
#enddo

#endprocedure