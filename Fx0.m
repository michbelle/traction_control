function Fx0=Fx0(k,mu_max)
%nominal force x

% k longitudinal slip
% mu attriction coefficent
Fx0=Fz*mu_max*sin(Cx*arctan(Bx*k-Ex*(Bx*k-Ex*(Bx*k-arctan(Bx*k)))));