#procedure mandelstams

id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = mH^2;
id p4.p4 = 0;

id p1.p2 = s/2;
id p3.p4 = s/2 - mH^2/2;
id p1.p3 = -t/2 + mH^2/2;
id p2.p4 = -t/2;
id p1.p4 = -u/2;
id p2.p3 = -u/2 + mH^2/2;

#endprocedure