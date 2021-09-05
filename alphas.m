function [alphaFL,alphaFR,alphaRL,alphaRR]=alphas(vby,vbx,dphi,deltaFL,deltaFR);
%evaluation of side slip angles
% vbx,y = velocity of the body
% delta FL, FR = angular rotation along y axis of the tire
% lf,r =distances along y axis of the vehicle from center of gravity
alphaFL=(vby+dphi*lf)/vbx-deltaFL;
alphaFR=(vby+dphi*lf)/vbx-deltaFR;
alphaRL=(vby-dphi*lr)/vbx;
alphaRR=(vby-dphi*lr)/vbx;