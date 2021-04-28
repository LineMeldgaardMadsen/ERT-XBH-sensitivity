%% 2D ERT sensitivity distributions for surface and buried electrodes 
% By Line Meldgaard Madsen, Aarhus University, Febuary 2021
% Bugs and comments may be directed to linemeldgaard@geo.au.dk
%
% The sensitivity is computed analytically for a pole-pole configuration  
% and superposition then gives the quadrupole (dipole-dipole) sensitivity.
% 
% The pole-pole sensitivity is computed as the dot product between Ec and 
% Ep, which are the electric fields from the current electrodes 
% (modelled as a point source) and the  potential electrode (had it 
% functioned as a current electrode, e.i. the adjoint electrode).
% E = grad(V)
% 
% To account for the boundary conditions when using buried electrode, the
% potential is computed using image sources (the source mirrored in the 
% surface), so
% v_2D(x,z) = 1/(2*pi) * ( ln(r)+ln(r_im) )
% grad(v_2D(x,z)) = 1/(2*pi) * (d/dx,d/dy,d/dz)( ln(r)+ln(r_im) )
% where r and r_im are the distance fra the course and image source and to
% the point (x,z)
%
% x is zero at the surface and negative downwards
%

clear all; clc; close all;

%% Input
%Define the location of the electrode
A = [-1 -0.5]; %(x,z) of the current source
B = [-1 -2];   %(x,z) of the current sink 
M = [1 -0.5];  %(x,z) of the potential electrode
N = [1 -2];     %(x,z) of the potential electrode

%Define xmin, xmax, zmin and zmax for the computation domain
xmin = -2;  xmax = 2;
zmin = -4;  zmax = 0; %Note that zmax<=0

%Define numer of nodes in the x and z-direction
Nnx = 500;
Nnz = 500;

%Define calculation points of sensitivity (x,z)
xs = linspace(-2,2,500);
zs = linspace(-4,0,500);

%% Check input
if(A(2)>0); error('A(2) must be 0 or negative'); end
if(B(2)>0); error('B(2) must be 0 or negative'); end
if(M(2)>0); error('M(2) must be 0 or negative'); end
if(N(2)>0); error('N(2) must be 0 or negative'); end
if(zmax>0); zmax=0; disp('Warning: zmax has been put to 0'); end
if(A(1)>xmax||A(1)<xmin||A(2)>zmax||A(2)<zmin); disp('Warning: A is outside the domain'); end
if(B(1)>xmax||B(1)<xmin||B(2)>zmax||B(2)<zmin); disp('Warning: B is outside the domain'); end
if(M(1)>xmax||M(1)<xmin||M(2)>zmax||M(2)<zmin); disp('Warning: M is outside the domain'); end
if(N(1)>xmax||N(1)<xmin||N(2)>zmax||N(2)<zmin); disp('Warning: N is outside the domain'); end

%% Computations
%Define calculation points of sensitivity (x,z)
xs = linspace(xmin,xmax,Nnx);
zs = linspace(zmin,zmax,Nnz);

%Loop the quadrupoles
for pole = 1:4

if(pole==1) %AM
    r = A(1); t = A(2); %source in (r,t)
    l = M(1); n = M(2); %Potential in (l,n)
elseif(pole==2) %AN
    r = A(1); t = A(2); %source in (r,t)
    l = N(1); n = N(2); %Potential in (l,n)
elseif(pole==3) %BM
    r = B(1); t = B(2); %source in (r,t)
    l = M(1); n = M(2); %Potential in (l,n)
elseif(pole==4) %BN
    r = B(1); t = B(2); %source in (r,t)
    l = N(1); n = N(2); %Potential in (l,n)
end

for ii = 1:length(xs) %x
for jj = 1:length(zs) %z
    %point of sensitivity (x,y,z)
    x = xs(ii);
    z = zs(jj);

    %Electric field from source
    Ec_x = 1/(2*pi) * ( (r-x)/((r-x)^2+(t-z)^2) + (r-x)/((r-x)^2+(t+z)^2) );
    Ec_z = 1/(2*pi) * ( (t-z)/((r-x)^2+(t-z)^2) + (t+z)/((r-x)^2+(t+z)^2) );

    %Electroc field from potential
    Ep_x = 1/(2*pi) * ( (l-x)/(((l-x)^2+(n-z)^2)) + (l-x)/(((l-x)^2+(n+z)^2)) );        
    Ep_z = 1/(2*pi) * ( (n-z)/(((l-x)^2+(n-z)^2)) + (n+z)/(((l-x)^2+(n+z)^2)) );

    %Sensitivity
    S_2D(pole,ii,jj) = Ec_x*Ep_x + Ec_z*Ep_z;
end
end

end

% Superposition 
S_tot_2D(:,:) = S_2D(1,:,:)- S_2D(2,:,:)- S_2D(3,:,:)+ S_2D(4,:,:);

%% Plotting
%Discritize the domain
[X,Z] = meshgrid(xs,zs);

%Normalise the sensitivity
q = max(quantile(S_tot_2D,0.95),[],'all');
S_tot_2D = S_tot_2D/q;

%Plot the sensitivity distribution
figure
pcolor(X',Z',S_tot_2D); shading flat    
hold on;
contour(X',Z',S_tot_2D,[0.1, 0.1],'-','color',[0.6 0 0])
contour(X',Z',S_tot_2D,[-0.1, -0.1],'-','color',[0 0 0.6])
contour(X',Z',S_tot_2D,[0.01, 0.01],'--','color',[0.6 0 0])
contour(X',Z',S_tot_2D,[-0.01,-0.01],'--','color',[0 0 0.6])

%Change the colormap
colorbar
colormap(S_colorbar)
caxis([-1 1])

%Plot the labels for the electrode
plot(A(1),A(2),'kx','markersize',5,'linewidth',1)
t=text(A(1),A(2),'A','Color',[0 0 0],'fontsize',11,'backgroundcolor',[1 1 1 0.5]);
plot(B(1),B(2),'kx','markersize',5,'linewidth',1)
text(B(1),B(2),'B','Color',[0 0 0],'fontsize',11,'backgroundcolor',[1 1 1 0.5]);
plot(M(1),M(2),'kx','markersize',10,'linewidth',2)
text(M(1),M(2),'M','Color',[0 0 0],'fontsize',11,'backgroundcolor',[1 1 1 0.5]);
plot(N(1),N(2),'kx','markersize',10,'linewidth',2)
text(N(1),N(2),'N','Color',[0 0 0],'fontsize',11,'backgroundcolor',[1 1 1 0.5]);
title('Normalised 2D ERT sensitivity')

xlabel('x'); ylabel('z');
set(gca,'fontsize',11)
