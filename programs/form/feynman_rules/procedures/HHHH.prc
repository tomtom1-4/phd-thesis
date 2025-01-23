#procedure HHHH

#write "Quadruple Higgs Vertex [H(p1), H(p2), H(p3), H(p4)]"

if(count(H, 1)!=4) discard;
.store

#redefine Nper "24"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the gluons to be incoming

#$Ldummy = LH;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#do j={1,2,3,4};
  id once H = 1;
#enddo

.store
.sort
#enddo

#endprocedure