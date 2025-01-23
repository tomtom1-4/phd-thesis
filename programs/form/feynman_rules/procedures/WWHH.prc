#procedure WWHH

#write "W-, W+, H, H-Vertex [W-(p1, mu), W+(p2, nu), H(p3), H(p4)]"

.sort
if(count(Wm, 1) != 1) discard;
if(count(Wp, 1) != 1) discard;
if(count(H, 1) != 2) discard;
.store

#redefine Nper "2"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LH;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
id once Wm(mu?) = ep1(mu);
id once Wp(mu?) = ep2(mu);
id once H = 1;
id once H = 1;
.store
.sort
#enddo

#endprocedure