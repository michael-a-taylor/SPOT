% Take an optimized solution and check performance, using parameters that
% can differ from those used in the original calculations. This file
% calculates lateral force Fx vs. x.

% This file is part of the SPOT package, as described in [publication].
% Please cite this if you use the code.

% This allows evaluation over different angular grids or with different
% Nmax (to assess numerical errors) or with different particle.

% This file should be run with the variables from SPOT optimization in the
% workspace. Specifically, it requires the variables:
% phase_opt, KX, KY.
%
% This file can also be used with different particle or optical parameters
% to those in the optimization, e.g. particle diameter or polarization

% Now what grid do we want to re-evaluate this over?
Nmax=26;
N_theta=100;
N_phi=120;


%Can comment this out if they are still in workspace, or set new values
%here
% radius=6/1.6; 
% NA=1.25;
% n_medium = 1.33; % Water
% n_particle = 1.46; %Silica
% polarisation = [1 1i];



% Particle displacement
x=linspace(-5,5,201);
dx=(x(2)-x(1));
zeq=0;



% Create the trapping field:

n_relative=n_particle/n_medium;
lg_mode=[0 0]; 
beam_offset = [ 0 0 0];
beam_angle = asin(NA/n_medium)*180/pi;
w0 = lg_mode_w0( lg_mode, beam_angle );

[n0,m0,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ lg_mode w0 1 polarisation 90 beam_offset ]);

[a,b,n,m] = make_beam_vector(a0,b0,n0,m0);

% Normalize total power of wave sum to 1.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
a=a/pwr;
b=b/pwr;



% Calculate the basis transform between (a,b) and E(theta,phi).
theta0=pi/(2*N_theta):pi/(N_theta):pi; %Use this instead of scaling from theta=0 -> theta=pi, since dA=0 for both endpoints. That is a waste of numerics.
phi0=pi/N_phi:2*pi/N_phi:2*pi;

[theta,phi]= meshgrid(theta0,phi0);
dA=sin(theta(:))*(theta0(2)-theta0(1))*(phi0(2)-phi0(1)); % The area of each grid point


[A_p,B_q]=farfield_matrix2(Nmax,theta(:),phi(:));

Aperture=theta(:)<asin(NA/n_medium);
Ap2=[Aperture;Aperture];


% Interpolate phase profile onto a grid.
[KXi,KYi]=meshgrid(2*pi*NA/n_medium*linspace(-1,1,501),2*pi*NA/n_medium*linspace(-1,1,501));
phase_grid = griddata(KX,KY,exp(1i*phase_opt),KXi,KYi);

phase_grid=phase_grid./abs(phase_grid); %Normalize amplitude
% imagesc(angle(phase_grid))


KX2=2*pi*sin(theta(Aperture)).*cos(phi(Aperture));
KY2=2*pi*sin(theta(Aperture)).*sin(phi(Aperture));


%Now interpolate back onto the new angular basis.
phase_test = angle(interp2(KXi,KYi,phase_grid,KX2,KY2));
phase_test(isnan(phase_test))=0;%Remove any NaNs that may appear at edge.



%Prepare coefficients for modified Gaussians
E_0=(A_p*(a)+B_q*(b)).*Ap2;

% Simple Gaussian, generated via basis transform. Should equal (a,b).
a1=A_p'*(E_0.*[dA;dA]);
b1=B_q'*(E_0.*[dA;dA]);

pwr = sqrt(sum( abs(a1).^2 + abs(b1).^2 ));
a1=a1/pwr;
b1=b1/pwr;

% Phase optimized
E2=E_0;
E2(Ap2)=E_0(Ap2).*exp(1i*[phase_test;phase_test]);

a_opt=A_p'*(E2.*[dA;dA]);
b_opt=B_q'*(E2.*[dA;dA]);

pwr = sqrt(sum( abs(a_opt).^2 + abs(b_opt).^2 ));
a_opt=a_opt/pwr;
b_opt=b_opt/pwr;


%Prepare for force calculations
T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,radius);


%calculate the x-axis coefficients for force calculation.
[rt,thetaOTT,phiOTT]=xyz2rtp(x,0,zeq);
Rx = z_rotation_matrix(pi/2,0);
Dx = wigner_rotation_matrix(Nmax,Rx);

Fx_0=zeros(size(x)); %Regular OTT with Gaussian trap
Fx_1=zeros(size(x)); %Convert (a,b) -> E(theta,phi) ->  (a,b). Allows estimation of the small numerical error associated with transforms. 
Fx_opt=zeros(size(x));% Phase optimized profile.

% Force calculation
for nr = 1:length(x)
    
    R = z_rotation_matrix(thetaOTT(nr),phiOTT(nr)); %calculates an appropriate axis rotation off z.
    D = wigner_rotation_matrix(Nmax,R);
    
    [A,B] = translate_z(Nmax,rt(nr));
    a2 = D'*(  A * D*a +  B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b +  B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_0(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
    
    a2 = D'*(  A * D*a1 +  B * D*b1 ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b1 +  B * D*a1 ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_1(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
    
    a2 = D'*(  A * D*a_opt +  B * D*b_opt ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b_opt +  B * D*a_opt ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_opt(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.    

end
Incl=abs(x)<=dx;

figure(1)
plot(x,Fx_0,x,Fx_1,x,Fx_opt)

% Gx=polyfit(x(Incl),Fx_0(Incl),1);
% kx_0=-Gx(1)*n_medium^2/(1064e-9*3e8)
% 
% Gx=polyfit(x(Incl),Fx_1(Incl),1);
% kx_1=-Gx(1)*n_medium^2/(1064e-9*3e8)

Gx=polyfit(x(Incl),Fx_opt(Incl),1);
kx_opt=-Gx(1)*n_medium^2/(1064e-9*3e8)*1e3

