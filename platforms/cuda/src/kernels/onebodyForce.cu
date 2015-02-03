real3 ROH1 = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
real3 ROH2 = make_real3(pos3.x-pos1.x, pos3.y-pos1.y, pos3.z-pos1.z);
real3 RHH  = make_real3(pos2.x-pos3.x, pos2.y-pos3.y, pos2.z-pos3.z);

ROH1 *= 10.;
ROH2 *= 10.;
RHH *= 10.;









energy += 1.;
real3 delta = make_real3(pos2.x-pos1.x, pos2.y-pos1.y, pos2.z-pos1.z);
real3 force1 = delta;
real3 force2 = delta;
real3 force3 = delta;
