#procedure ZWW

#write "W-,W+,Z-Vertex [W-(p1, mu), W+(p2, nu), Z(p3, rho)]"

if(count(g2,1)!=1) discard;
.sort
s bosoncount;
id Z(mu?) = Z(mu)*bosoncount;
id Wp(mu?) = Wp(mu)*bosoncount;
id Wm(mu?) = Wm(mu)*bosoncount;
id D(Z(mu?), nu?) = D(Z(mu), nu)*bosoncount;
id D(Wp(mu?), nu?) = D(Wp(mu), nu)*bosoncount;
id D(Wm(mu?), nu?) = D(Wm(mu), nu)*bosoncount;
if(count(bosoncount, 1) != 3) discard;
id bosoncount = 1;
*if(count(D, 1) !=1) discard;
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
  id once Z(mu?) = ep3(mu);
  id once D(Wm(mu?), nu?) = (-i_)*p1(nu)*ep1(mu);
  id once D(Wp(mu?), nu?) = (-i_)*p2(nu)*ep2(mu);
  id once D(Z(mu?), nu?) = (-i_)*p3(nu)*ep3(mu);
#endif

.store
.sort
#enddo

#endprocedure