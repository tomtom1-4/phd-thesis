#procedure GGG

#write "Triple Gluon Vertex [A(p1, mu, c1), A(p2, nu, c2), A(p3, rho, c3)]"

if(count(g,1)!=1) discard;
if(count(A, 1)!=2) discard;
if(count(D, 1) !=1) discard;
.store

#redefine Nper "6"

* derivatives are replace by (-i_)*p in momentum space. The minus defines the gluons to be incoming

#$Ldummy = LG;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p1(nu)*ep1(mu)*d_(cOli1, c1);
  #do j={2,3};
    id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==2;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p1(nu)*ep1(mu)*d_(cOli1, c1);
  #do j={3,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==3;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p2(nu)*ep2(mu)*d_(cOli1, c2);
  #do j={1,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==4;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p2(nu)*ep2(mu)*d_(cOli1, c2);
  #do j={3,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==5;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p3(nu)*ep3(mu)*d_(cOli1, c3);
  #do j={1,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==6;
  id once D(A(cOli1?, mu?), nu?) = (-i_)*p3(nu)*ep3(mu)*d_(cOli1, c3);
  #do j={2,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif

.store
.sort
#enddo

#endprocedure