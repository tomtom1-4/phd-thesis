#procedure epExpand
mul replace_(d, 4 - 2*ep);
b den;
.sort
keep brackets;
#do i=1,1
  splitarg ((ep)) den;
  if(match(den(0,?args))) redefine i "0";
  repeat id den(mq?,mH?,be?,?args) = den(mq,mH+be,?args);
  id den(0,mq?) = 1/ep * den(mq/ep);
  b den;
  .sort
  keep brackets;
#enddo
.sort

#$max = 0;
if (count(ep,-1) > $max) $max = count_(ep,-1);

b den;
.sort
keep brackets;
#if `$max' > 0
  id den(mH?,be?) = den(mH)*sum_(mq,0,`$max',(-be*den(mH))^mq);
#else
  id den(mH?,be?) = den(mH);
#endif
.sort

if(count(ep,1)>0) discard;

factarg den;
repeat id den(mH?,be?,?args) = den(mH)*den(be)*den(?args);
id den = 1;

splitarg den;
id den(mq?) = 1/mq;
repeat id den(mH?,be?,?args) = den(mH+be,?args);

#endprocedure