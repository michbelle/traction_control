function Gyk = Gyk(k)
%function used to calculate the factor to obtain y force from slip
Gyk=(cos(Cyk*arctan(Byk*(k+Shyk))))/(cos(Cyk*arctan(Byk*Shyk)));
