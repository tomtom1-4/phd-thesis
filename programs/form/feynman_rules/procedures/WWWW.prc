#procedure WWWW

#write "W-,W+,W-,W+-Vertex [W-(p1, mu), W+(p2, nu), W-(p3, rho), W+(p4, sigma)]"

if(count(g2,1)!=2) discard;
.sort
if(count(Wp, 1)!=2) discard;
if(count(Wm, 1)!=2) discard;
.store

#redefine Nper "4"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LWeak;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep2(mu);
  id once Wm(mu?) = ep3(mu);
  id once Wp(mu?) = ep4(mu);
#endif
#if `permutation'==2;
  id once Wm(mu?) = ep3(mu);
  id once Wp(mu?) = ep2(mu);
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep4(mu);
#endif
#if `permutation'==3;
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep4(mu);
  id once Wm(mu?) = ep3(mu);
  id once Wp(mu?) = ep2(mu);
#endif
#if `permutation'==4;
  id once Wm(mu?) = ep3(mu);
  id once Wp(mu?) = ep4(mu);
  id once Wm(mu?) = ep1(mu);
  id once Wp(mu?) = ep2(mu);
#endif

.store
.sort
#enddo

#endprocedure