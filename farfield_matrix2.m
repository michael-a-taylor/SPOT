function [A_p,B_q] = farfield_matrix2(Nmax,theta,phi)
% farfield_matrix2.m
% Finds matrices which return the far field from VSWF expansion coefficients
%
% usage
% [A_p,B_q] = farfield_matrix(n,m,theta,phi)
%
% then, the electric field is given by
% E=reshape(A_p*(a+2*p2)+B_q*(b+2*q2),length(theta(:)),3);
%
% each row of E is the field (in spherical coordinates) in the
% (theta,phi) direction (assuming a distance scaling factor of kr)
%
%
% This file is part of the Quadrant Detection Package. It is a revised
% version of farfield.m, which now operates much faster and with improved
% numerical precision. Also, it skips radial field components (which are
% always zero) so it is less memory intensive.
%
% For details see:
%
% M. A. Taylor and W. P. Bowen, “A computational tool to characterize
% particle tracking measurements in optical tweezers,” J. Opt. 15, 085701
% (2013).  
%
% Copyright 2013-2016 The University of Queensland.
%
%
%
% The file is modified from "farfield.m" in the Optical tweezers toolbox 1.2,
% and calls on other functions from that software package.
%
% For details of the Optical tweezers toolbox, see:
% http://www.physics.uq.edu.au/people/nieminen/software.html
%(Copyright 2006-2013 The University of Queensland.)


% Assume m=-n:n

[theta,phi] = matchsize(theta,phi);

A_p = zeros(length(theta)*2,Nmax^2+2*Nmax);
B_q = zeros(length(theta)*2,Nmax^2+2*Nmax);

for n = 1:Nmax
   

[~,Ytheta,Yphi] = spharm(n,-n:n,theta,phi);


   Nn = 1/sqrt(n*(n+1));

   A_p(:,n^2:(n^2+2*n)) =  Nn * (-1i)^(n+1)*[Yphi;-Ytheta];
   B_q(:,n^2:(n^2+2*n)) =  Nn * (-1i)^n*[Ytheta;Yphi];   
end



return