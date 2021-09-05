function Fy0=Fy0(alpha,mu_max)
%nominal force x

% alpha side slip angle
% mu attriction coefficent
Fy0=Fz*mu_max*sin(Cy*arctan(By*alpha-Ey*(By*alpha-arctan(By*alpha))));