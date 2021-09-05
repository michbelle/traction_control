function [Gxa] = Gxa(alpha)
%function used to calculate the factor to obtain x force from side-slip
Gxa=cos(Cxa*arctan(Bxa*alpha-Exa*(Bxa*alpha-arctan(alpha))));
