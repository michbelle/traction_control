function k=k(v,omega,alpha)
%longitudinal slip
% v : velocity
% omega : angular velocity
% alpha : side-slip angle
vx=v*cos(alpha);
k=-(vx-omega*re)/vx;
% re : radius of wheel