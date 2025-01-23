#procedure AAWW

#write "W-,W+,A,A-Vertex [W-(p1, mu), W+(p2, nu), Agamma(p3, rho), A(p4, sigma)]"

if(count(g2,1)!=2) discard;
.sort
if(count(Agamma, 1)!=2) discard;
id once Wp(mu?) = Wp(mu)*bosoncount;
id once Wm(mu?) = Wm(mu)*bosoncount;
if(count(bosoncount, 1) != 2) discard;
id bosoncount = 1;
*if(count(D, 1) !=1) discard;
.store

#redefine Nper "2"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LWeak;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep2(mu);
  id once Agamma(mu?) = ep3(mu);
  id once Agamma(mu?) = ep4(mu);
#endif
#if `permutation'==2;
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep2(mu);
  id once Agamma(mu?) = ep4(mu);
  id once Agamma(mu?) = ep3(mu);
#endif

.store
.sort
#enddo

#endprocedure