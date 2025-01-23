#procedure GCC

#write " Gluon Ghost Ghost vertex [CBar(p1, c1), C(p2, c2), A(p3, rho, c3)]"

if(count(g,1)!=1) discard;
if(count(A, 1)!=1) discard;
if(count(C, 1)!=1) discard;
.store

#redefine Nper "1"

* derivatives are replaced by i_*p in momentum space. The minus defines the momentum of the ghost to be outgoing

#$Ldummy = LFP;
#do permutation=1,1;
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once D(CBar(cOli1?), mu?) = i_*p1(mu)*d_(cOli1, c1);
  id once C(cOli1?) = d_(cOli1, c2);
  id once A(cOli1?, mu?) = ep3(mu)*d_(cOli1, c3);
#endif


.store
.sort
#enddo

#endprocedure