v=100;
gamma=0.3;
r_wheel=0.02;


L=3;%lunghezza veicolo
d_r=0.3;
d_l=0.3;
d=L/tan(gamma);
w_cir=v/d;
v=[sqrt((d-d_l)^2+L^2),sqrt((d+d_r)^2+L^2),d-d_l,d+d_r]*w_cir;
w=v*r_wheel;
w_fl=w(1);
w_fr=w(2);
w_rl=w(3);
w_rr=w(4);