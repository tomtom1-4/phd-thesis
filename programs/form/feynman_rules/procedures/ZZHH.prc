#procedure ZZHH

#write "Z, Z, H, H-Vertex [Z(p1, mu), Z(p2, nu), H(p3), H(p4)]"

.sort
if(count(Z, 1) != 2) discard;
if(count(H, 1) != 2) discard;
.store

#redefine Nper "4"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the bosons to be incoming

#$Ldummy = LH;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if((`permutation'==1) || (`permutation'==2))
  id once Z(mu?) = ep1(mu);
  id once Z(mu?) = ep2(mu);
  id once H = 1;
  id once H = 1;
#endif
#if((`permutation'==3) || (`permutation'==4))
  id once Z(mu?) = ep2(mu);
  id once Z(mu?) = ep1(mu);
  id once H = 1;
  id once H = 1;
#endif
.store
.sort
#enddo

#endprocedure