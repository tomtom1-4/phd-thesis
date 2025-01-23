#procedure HHH

#write "Triple Higgs Vertex [H(p1), H(p2), H(p3)]"

if(count(H, 1)!=3) discard;
.store

#redefine Nper "6"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the gluons to be incoming

#$Ldummy = LH;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#do j={1,2,3};
  id once H = 1;
#enddo

.store
.sort
#enddo

#endprocedure