#procedure ZZH

#write "Z,Z,H-Vertex [Z(p1, mu), Z(p2, nu), H(p3)]"

.sort
if(count(Z, 1) != 2) discard;
if(count(H, 1) != 1) discard;
.store

#redefine Nper "2"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LH;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once Z(mu?) = ep1(mu);
  id once Z(mu?) = ep2(mu);
  id once H = 1;
#endif
#if `permutation'==2;
  id once Z(mu?) = ep2(mu);
  id once Z(mu?) = ep1(mu);
  id once H = 1;
#endif


.store
.sort
#enddo

#endprocedure