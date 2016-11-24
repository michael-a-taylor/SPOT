tic
%This applies the SPOT algorithm to find a phase profile that optimizes the
%lateral trap stiffness when applied at the back-focal plane.  
%
% This file is part of the SPOT package, as described in [publication], and
% hosted online at:
% https://github.com/michael-a-taylor/SPOT
%
% Steps:
% 1) Define the parameters needed for the calculation.
% 2) Calculate the Mie scattering matrix in the E(theta,phi) basis.
% 3) Construct the matrix that predicts optical force for any incident field.
% 4) Calculate the 2D matrix that defines the spring constant.
% 5) Find optimized phase profile in the E(theta,phi) basis.
% 6) (Optional) Find optimized wavefront using Eigenmode method
% 7) Check results with Optical Tweezers Toolbox calculations


%% 1) Define the parameters needed for the calculation.

compute_Eigenmode=1;% A logical switch. If true, the Eigenmode solution will be computed and compared to phase-only optima.

% Particle size
radius =2/1.6; %In units of medium wavelength. 
% For 1064nm light in water, this is (Diameter in micron)/1.6


% Particle displacements at which we will assess optical forces 
x=linspace(-5,5,201);
dx=(x(2)-x(1));


%Numerical Aperture; this defines the input angles. 
NA=1.25;

% Polarisation. This parameter has quite a big influence.
% [ 1 1i ] and [1 -1i] are circularly polarised
% [ 1 0 ] is x-polarized, 
% [ 0 1 ] is y-polarized,
polarisation = [ 1 1i ];


%Refractive indices
n_medium = 1.33; % Water

% Specify refractive index of particle
n_particle = 1.46; %Silica
% n_particle = 1.58;% Polystyrene
% n_particle = 1.68;% Resin


wavelength=1064e-9; %Vacuum wavelength; this is just used to normalize to SI units

kx_norm=n_medium^2/(wavelength*3e8);% Normalize trap stiffness to mN/m/W. 


% Numerical expansion size. The formulae here calculate reasonable values
% based on the radius, so it is not strictly necessary to set these.

%To what order do we expand the vector spherical harmonics?
Nmax=max(ceil(7.5*radius+4),15); 
% Nmax=28;

% Angular basis: How many points along theta and phi do we use? 
N_theta=round(Nmax*2/3)+70;% This seems to work ok. 
N_phi=max(ceil((Nmax*2.35+10)/4)*4,60); 
% N_phi should be divisible by 4, so that each quadrant includes the same number of points.
% The VSH functions are aliased if N_phi<2*Nmax. Definitely use more than that



%% 2) Calculate the Mie scattering matrix in the E(theta,phi) basis.


% First step: calculate the basis transform between (a,b) and E(theta,phi).

% Define the angular grid
theta0=pi/(2*N_theta):pi/(N_theta):pi; %Use this instead of scaling theta from 0 to pi, since dA=0 for both endpoints. Those points are therefore excluded anyway.
phi0=pi/N_phi:2*pi/N_phi:2*pi;

[theta,phi]= meshgrid(theta0,phi0);

dA=sin(theta(:))*(theta0(2)-theta0(1))*(phi0(2)-phi0(1)); % The area of each grid point

% Compute basis transform matrices.
[A_p,B_q]=farfield_matrix2(Nmax,theta(:),phi(:)); 

% Define which input angles fall within the objective aperture
Aperture=theta(:)<asin(NA/n_medium);
Ap2=[Aperture;Aperture];

% Downsampled basis transform matrices, which only consider the input modes
% that are actually used
A_p_inc=A_p(Ap2,:);
B_q_inc=B_q(Ap2,:);


dA_ap=[dA(Aperture);dA(Aperture)];



% Next Step: Calculate T matrix and transform it to the angular basis.

% In principle it should be enough to use the transmission matrix, as
% described in the paper; but this can lead to excess numerical errors,
% particularly for small particles. It is better to treat the scattered and
% unscattered components separately. 

% Regular T matrix. This matrix is diagonal, so taking the diag() does not
% throw anything away.
T_ab = diag(tmatrix_mie(Nmax,2*pi,2*pi*n_particle/n_medium,radius));
% The T matrix from OTT defines the scattered light. The total transmission
% is given by (1+2*T)


% T matrix in the (theta, phi) basis.
T_E = A_p*diag(2*T_ab(1:(Nmax^2+2*Nmax)))*A_p_inc'+B_q*diag(2*T_ab((Nmax^2+2*Nmax+1):end))*B_q_inc';

% This applies the same basis transforms to the unscattered light. In
% principle it should be possible to remove this.
T0=A_p_inc*A_p_inc'+B_q_inc*B_q_inc';



%% 3) Construct the matrix that predicts optical force for any incident field.


% This is the transverse optical momentum. Here we work in OTT units of Q
% which already absorbs the missing n_medium factor   
PX=sin([theta(:);theta(:)]).*cos([phi(:);phi(:)]).*[dA;dA]; %

% This is the quadrant signal. Use this instead if optimizing SNR.
% PX=sign(cos([phi(:);phi(:)])).*Ap2.*[dA;dA];  

% Momentum after scattering is given by E'*(1+2*T)'*PX*(1+2*T)*E.  
% The force can be calculated from this matrix as E' * A_force_mat * E 
A_force_mat=T_E(Ap2,:)'*diag(PX(Ap2))*T0+T0'*diag(PX(Ap2))*T_E(Ap2,:)+T_E'*diag(PX)*T_E;


clear T_E T0; %These are very large matrices that we don't need anymore. 

%% 4) Calculate the 2D matrix that defines the spring constant.


% First create the trapping field amplitude.

% Spatial mode of the incident light in the Laguerre-Gauss basis
lg_mode=[0 0]; 

% Location of the focal point relative to the particle.
beam_offset = [ 0 0 0];

% Specify the beam width. 
beam_angle = asin(NA/n_medium)*180/pi;
w0 = lg_mode_w0( lg_mode, beam_angle );

% Use this for focused Gaussian field
[n0,m0,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ lg_mode w0 1 polarisation 90 beam_offset ]);

[a,b,n,m] = make_beam_vector(a0,b0,n0,m0);

% Normalize total power of wave sum to 1.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
a=a/pwr;
b=b/pwr;


% Convert to E(theta,phi) basis, remove light that is outside the aperture
E_0=(A_p*(a)+B_q*(b)).*Ap2;

% We sum over these, so we also need the area increment dA
E1=E_0.*[dA;dA];


% Transverse optical momentum
KX=2*pi*sin(theta(Aperture)).*cos(phi(Aperture));
KY=2*pi*sin(theta(Aperture)).*sin(phi(Aperture));


% Now, skip forces and go straight to calculation of spring constant k_x
D_MX=-1i*diag(E1(Ap2)')*(A_force_mat.*([KX;KX]*ones(1,sum(Ap2))-ones(sum(Ap2),1)*[KX;KX]'))*diag(E1(Ap2));

% Usage; spring constant kx=real(exp(1i*phase_XY)'*D_MX*exp(1i*phase_XY));
% Note, however, that the basis currently treats the two polarizations
% separately.

% Downsample the matrix so that the orthogonal polarization components are
% treated together, with the same phase always applied to both. 
D_MX_small=D_MX(1:length(KX),1:length(KX))+D_MX(length(KX)+1:end,length(KX)+1:end)+D_MX(1:length(KX),length(KX)+1:end)+D_MX(length(KX)+1:end,1:length(KX));


% Spring constants calculated using this matrix method. 
kx_mat=zeros(1,4);% 4 elements: Gaussian, phase-optimized, Eig-phase, Eig-opt
kx_mat(1)=real((ones(size(E1(Ap2)))'*D_MX*ones(size(E1(Ap2)))));


phase_opt=zeros(size(KX));

Spring_kx=kx_mat(1);% Track changes in spring constant with iteration number.


%% 5) Find optimized phase profile in the E(theta,phi) basis.

%Quadrants: use symmetries so that we only need to optimize part of the
%phase. Apply rotational symmetry to halve the space:
Q1=find(KX>0&KY>0);
% Q2=KX<0&KY>0;
% Q3=KX<0&KY<0;
% Q4=KX>0&KY<0;


Iteration_num=0;
run_optimization=1;

% Run optimization
while run_optimization
Iteration_num=Iteration_num+1;
for num=1:length(Q1)%length(KX)
    
    % 180 degree rotation symmetry
    N1=Q1(num);
    
    phase_opt(N1)=angle(D_MX_small(N1,:)*(exp(1i*[phase_opt])));
    phase_opt(N1+N_phi/2)=phase_opt(N1);
    
    N2=N1+N_phi/4;
    phase_opt(N2)=angle(D_MX_small(N2,:)*(exp(1i*[phase_opt])));
    phase_opt(N2+N_phi/2)=phase_opt(N2);

%     % x-y reflection symmetry
%     N1=Q1(num);
%     
%     NT=floor(N1/N_phi)*N_phi+1-rem(N1,N_phi);
%     
%     phase_opt(N1)=angle(D_MX_small(N1,:)*(exp(1i*[phase_opt])));
%     phase_opt(N_phi/2+NT)=phase_opt(N1);
%     
%     phase_opt(N1+N_phi/2)=phase_opt(N1);
%     phase_opt(N_phi+NT)=phase_opt(N1);

end

%Check new spring constant
Spring_kx(Iteration_num+1)=real((exp(1i*[phase_opt]))'*D_MX_small*(exp(1i*[phase_opt])));

if(Spring_kx(Iteration_num+1)/Spring_kx(Iteration_num)<1.01) %Ending criteria
    %End loop
    run_optimization=0;
end
end

kx_mat(2)=real((exp(1i*[phase_opt]))'*D_MX_small*(exp(1i*[phase_opt])));

%% 6) (Optional) Find optimized wavefront using Eigenmode method

if(compute_Eigenmode)
% For Eigenmode method the spring constant matrix D_MX must be normalized such
% that the input modes each contribute equally to the power. Since
% power=sum(abs(E).^2.*dA), the input modes must include normalization to
% the area, i.e. E.*dA^(1/2), not just E. 
%
% Recalculate the spring constant matrices with this normalization.


D_MX=-1i*diag(dA_ap.^(1/2))*(A_force_mat.*([KX;KX]*ones(1,sum(Ap2))-ones(sum(Ap2),1)*[KX;KX]'))*diag(dA_ap.^(1/2));

% Downsample this into fixed polarization basis. The polarization in (theta,phi) basis:
pol_0=[-polarisation(1)*cos(phi(Aperture)) - polarisation(2)*sin(phi(Aperture)),...
       polarisation(1)*sin(phi(Aperture)) - polarisation(2)*cos(phi(Aperture))];

D_MX_eig=diag(pol_0(:,1)')*D_MX(1:length(KX),1:length(KX))*diag(pol_0(:,1))...
    +diag(pol_0(:,2)')*D_MX(length(KX)+1:end,length(KX)+1:end)*diag(pol_0(:,2))...
    +diag(pol_0(:,1)')*D_MX(1:length(KX),length(KX)+1:end)*diag(pol_0(:,2))...
    +diag(pol_0(:,2)')*D_MX(length(KX)+1:end,1:length(KX))*diag(pol_0(:,1));


% Find Eigenmodes of the "spring constant" matrix
[Eigenvectors,Eigenvalues] = eig(D_MX_eig);

% Select the best solution
Opt_Eig=Eigenvectors(:,real(diag(Eigenvalues))==max(real(diag(Eigenvalues))));

% Extract the phase of this solution
Phase_EI=angle(Opt_Eig);

% Predict spring constants
kx_mat(3)=real((exp(1i*Phase_EI))'*D_MX_small*(exp(1i*Phase_EI)));% Phase-only
kx_mat(4)=real(Opt_Eig'*D_MX_eig*Opt_Eig);                        % Full wavefront

end

kx_mat=kx_mat*kx_norm; %Normalize kx into SI units N/m/W.




%% 7) Check results with Optical Tweezers Toolbox calculations


% Transform phase-optimized solution to (a,b) basis.
E2=E_0;
E2(Ap2)=E_0(Ap2).*exp(1i*[phase_opt;phase_opt]);

a_opt=A_p'*(E2.*[dA;dA]);
b_opt=B_q'*(E2.*[dA;dA]);

pwr = sqrt(sum( abs(a_opt).^2 + abs(b_opt).^2 ));
a_opt=a_opt/pwr;
b_opt=b_opt/pwr;


if(compute_Eigenmode)
    % Transform Eigenmode solution to (a,b) basis.
    E_eig=E_0;
    E_eig(Ap2)=[Opt_Eig;Opt_Eig].*pol_0(:).*dA_ap.^(-1/2);
    % We need dA_ap.^(-1/2) as Opt_Eig is defined in the basis E.*dA.^(1/2);
    % pol_0 sets the polarization.
    
    
    E_eig=E_eig/sqrt(sum(abs(E_eig).^2.*[dA;dA])); %Normalize power
    
    
    a_ei=A_p'*(E_eig.*[dA;dA]);
    b_ei=B_q'*(E_eig.*[dA;dA]);
    
    pwr = sqrt(sum( abs(a_ei).^2 + abs(b_ei).^2 ));
    a_ei=a_ei/pwr;
    b_ei=b_ei/pwr;
    
    
    % Transform phase-only part of Eigenmode solution to (a,b) basis.
    E2=E_0;
    E2(Ap2)=E_0(Ap2).*exp(1i*[Phase_EI;Phase_EI]);
    
    E2=E2/sqrt(sum(abs(E2).^2.*[dA;dA])); %Normalize power
    
    a_EIP=A_p'*(E2.*[dA;dA]);
    b_EIP=B_q'*(E2.*[dA;dA]);
    
    pwr = sqrt(sum( abs(a_EIP).^2 + abs(b_EIP).^2 ));
    a_EIP=a_EIP/pwr;
    b_EIP=b_EIP/pwr;
    
    % How closely can phase-only control approximate the optimized eigenmode?
    overlap=real(sum(conj(E_eig).*E2.*[dA;dA]));

end

% For good measure, transform the Gaussian field back into (a,b) basis.
% This should make no difference, but sometimes it does introduce small
% numerical changes - perhaps because OTT uses a true Gaussian profile,
% which has small but non-zero amplitudes incident at angles that are
% blocked by the objective aperture. Here these are set to zero.
a1=A_p'*(E_0.*[dA;dA]);
b1=B_q'*(E_0.*[dA;dA]);

pwr = sqrt(sum( abs(a1).^2 + abs(b1).^2 ));
a1=a1/pwr;
b1=b1/pwr;

n_relative=n_particle/n_medium;

% Calculate T matrix again. 
T = tmatrix_mie(Nmax,2*pi,2*pi*n_relative,radius);


% Pre-allocate the forces that are to be calculated
Fx_OTT = zeros(size(x));
Fx_OTT1 = zeros(size(x));
Fx_opt_OTT=zeros(size(x));
Fx_EIphase_OTT=zeros(size(x));
Fx_EI_OTT=zeros(size(x));




%calculate the x-axis coefficients for force calculation.
Rx = z_rotation_matrix(pi/2,0);
Dx = wigner_rotation_matrix(Nmax,Rx);

zeq=0;
[rt,thetaOTT,phiOTT]=xyz2rtp(x,0,zeq);

% Calculate all forces at once to save having to repeat calculations.
for nr = 1:length(x)
    
    % These matrices are used by all calculations below
    R = z_rotation_matrix(thetaOTT(nr),phiOTT(nr)); %calculates an appropriate axis rotation off z.
    D = wigner_rotation_matrix(Nmax,R);
    [A,B] = translate_z(Nmax,rt(nr));
    
    % The Gaussian as usual for OTT
    a2 = D'*(  A * D*a +  B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b +  B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_OTT(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
    
    
    % The Gaussian mode transformed to E(theta,phi) basis and then back again
    a2 = D'*(  A * D*a1 +  B * D*b1 ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b1 +  B * D*a1 ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_OTT1(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.

    
    % The phase optimized profile
    a2 = D'*(  A * D*a_opt +  B * D*b_opt ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b_opt +  B * D*a_opt ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    Fx_opt_OTT(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.    

    if(compute_Eigenmode)
        % The Eigenmode solution
        a2 = D'*(  A * D*a_ei +  B * D*b_ei ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
        b2 = D'*(  A * D*b_ei +  B * D*a_ei ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
        
        pq = T * [ a2; b2 ];
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
        
        Fx_EI_OTT(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.
        
        % The phase-only Eigenmode solution
        a2 = D'*(  A * D*a_EIP +  B * D*b_EIP ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
        b2 = D'*(  A * D*b_EIP +  B * D*a_EIP ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
        
        pq = T * [ a2; b2 ];
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
        
        Fx_EIphase_OTT(nr) = force_z(n,m,Dx*a2,Dx*b2,Dx*p,Dx*q); %Dx makes the z-force calculation the x-force calculation.

    end
end



Fx0_mat=zeros(size(x));
Fxopt_mat=zeros(size(x));

% Example force calculation using A_force_mat; use this to verify that it agrees
% with normal Optical Tweezers Toolbox calculations. Particle is displaced
% by applying a phase of KX * x to the light.
for num=1:length(x)
    
    %Gaussian trap
    phase_XY=KX*x(num); % don't include n_medium because x is already given in medium wavelengths 
    
    Fx0_mat(num)=-real((E1(Ap2).*exp(1i*[phase_XY;phase_XY]))'*A_force_mat*(E1(Ap2).*exp(1i*[phase_XY;phase_XY])));

    % Optimized trap
    phase_XY=phase_opt+KX*x(num); % applied phase + phase from displacement
    
    Fxopt_mat(num)=-real((E1(Ap2).*exp(1i*[phase_XY;phase_XY]))'*A_force_mat*(E1(Ap2).*exp(1i*[phase_XY;phase_XY])));

end



figure(1) % This is the "verify accuracy" figure
subplot(3,1,1)
semilogy(abs(T_ab))% Use this to check if N_max is large enough to describe the scatter. I like to truncate where T-matrix is <10^-8.
figure(1)
subplot(3,1,2)
plot(x,Fx0_mat,x,Fx_OTT,x,Fx_OTT1) % Gaussian trap, calculated 3 ways. Agreement indicates that A_force_mat is calculated correctly.
figure(1)
subplot(3,1,3)
plot(x,Fxopt_mat,x,Fx_opt_OTT)% Optimized trap, calculated 2 ways.


%Interpolate profiles onto a grid for a nicer viewing experience
[KXi,KYi]=meshgrid(2*pi*NA/n_medium*linspace(-1,1,501),2*pi*NA/n_medium*linspace(-1,1,501));
Z_opt = griddata(KX,KY,exp(1i*phase_opt),KXi,KYi);

figure(2) %Display phase-only profile
imagesc(angle(Z_opt))
title('Optimized phase')

if(compute_Eigenmode) %Display results of Eigenmode method
    
    figure(3)
    Z_eig = griddata(KX,KY,Opt_Eig.*dA(Aperture).^(-1/2),KXi,KYi);
    imagesc(angle(Z_eig))
    title('Eigenmode phase')

    figure(4)
    imagesc(abs(Z_eig))
    title('Eigenmode amplitude')

    figure(5)
    plot(x,Fx_OTT,x,Fx_opt_OTT,x,Fx_EIphase_OTT,x,Fx_EI_OTT)
    xlabel('x (\lambda/n_m)')
    ylabel('Force (Q)')
    legend('No shaping','phase-only','Eigenmode phase','full Eigenmode')
else
    figure(5)
    plot(x,Fx_OTT,x,Fx_opt_OTT)
    xlabel('x (\lambda/n_m)')
    ylabel('Force (Q)')
    legend('No shaping','phase-only')
end

% Now calculate kx from a linear fit to central part of the force-displacement curves.

kx_OTT=zeros(1,5);% Gaussian, phase-optimized, EI-phase, EI-opt, Gaussian OTT without basis transform

%Region to include in the fit
Incl=abs(x)<=dx;

Gx=polyfit(x(Incl),Fx_OTT1(Incl),1);
kx_OTT(1)=-Gx(1)*kx_norm;

Gx=polyfit(x(Incl),Fx_opt_OTT(Incl),1);
kx_OTT(2)=-Gx(1)*kx_norm;

Gx=polyfit(x(Incl),Fx_EIphase_OTT(Incl),1);
kx_OTT(3)=-Gx(1)*kx_norm;

Gx=polyfit(x(Incl),Fx_EI_OTT(Incl),1);
kx_OTT(4)=-Gx(1)*kx_norm;

Gx=polyfit(x(Incl),Fx_OTT(Incl),1);
kx_OTT(5)=-Gx(1)*kx_norm;


% kx_mat
% 
% kx_OTT

toc