#procedure AZWW

#write "W-,W+,A,Z-Vertex [W-(p1, mu), W+(p2, nu), Agamma(p3, rho), A(p4, sigma)]"

if(count(g2,1)!=2) discard;
.sort
if(count(Agamma, 1)!=1) discard;
if(count(Z, 1)!=1) discard;
if(count(Wp, 1)!=1) discard;
if(count(Wm, 1)!=1) discard;
.store

#redefine Nper "1"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LWeak;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep2(mu);
  id once Agamma(mu?) = ep3(mu);
  id once Z(mu?) = ep4(mu);
#endif
.store
.sort
#enddo

#endprocedure