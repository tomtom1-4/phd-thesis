#procedure GGGGH

#write "Effectiv 4-Gluon-Higgs Vertex [A(p1, mu, c1), A(p2, nu, c2), A(p3, rho, c3), A(p4, sigma, c4), H(p5)]"

if(count(A,1)!=4) discard;
.store

#redefine Nper "24"

#$Ldummy = LHTL;
#do permutation=1,`Nper';
.sort
g L`permutation' = $Ldummy;
#if `permutation'==1;
  #do j={1,2,3,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==2;
  #do j={1,2,4,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==3;
  #do j={1,3,2,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==4;
  #do j={1,3,4,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==5;
  #do j={1,4,2,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==6;
  #do j={1,4,3,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==7;
  #do j={2,1,3,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==8;
  #do j={2,1,4,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==9;
  #do j={2,3,1,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==10;
  #do j={2,3,4,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==11;
  #do j={2,4,1,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==12;
  #do j={2,4,3,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==13;
  #do j={3,1,2,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==14;
  #do j={3,1,4,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==15;
  #do j={3,2,1,4};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==16;
  #do j={3,2,4,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==17;
  #do j={3,4,1,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==18;
  #do j={3,4,2,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==19;
  #do j={4,1,2,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==20;
  #do j={4,1,3,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==21;
  #do j={4,2,1,3};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==22;
  #do j={4,2,3,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==23;
  #do j={4,3,1,2};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
#if `permutation'==24;
  #do j={4,3,2,1};
  id once A(cOli1?, mu?) =   ep`j'(mu)*d_(cOli1, c`j');
  #enddo
#endif
id once H = 1;
.store
.sort
#enddo

#endprocedure