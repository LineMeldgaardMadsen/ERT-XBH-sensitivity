%% 2D ERT sensitivity distributions for surface and buried electrodes 
% By Line Meldgaard Madsen, Aarhus University, Febuary 2021
% Bugs and comments may be directed to linemeldgaard@geo.au.dk
%
% The sensitivity is computed analytically for a pole-pole configuration  
% and superposition then gives the quadrupole (dipole-dipole) sensitivity.
% 
% The pole-pole sensitivity is computed as the dot product between Ec and 
% Ep, which are the electric fields from the current electrodes 
% (modelled as a current source) and the  potential electrode (had it 
% functioned as a current electrode, e.i. the adjoint electrode).
% E = grad(V)
% 
% To account for the boundary conditions when using buried electrode, the
% potential is computed using image sources (the source mirrored in the 
% surface), so
% v = 1/(4*pi) * (1/r + 1/r_im)
% grad(v) = 1/(4*pi) * (d/dx,d/dy,d/dz)((1/r + 1/r_im)) 
% where r and r_im are the distance fra the course and image source and to
% the point (x,y,z)
%
% Modifications
% April 2021,   Anders Kühl, Aarhus University
%               Added iso-surface plotting

clear all; clc; close all;

%% Input
%z is zero at the surface and negative in the ground.
%Define the location of the electrode
A = [ 1 -0.5 -1];    %(x,y,z) of the current source
B = [ 1 -0.5 -2];    %(x,y,z) of the current sink 
M = [-2 -2   -1];    %(x,y,z) of the potential electrode
N = [-2 -2   -2];    %(x,y,z) of the potential electrode

%Define xmin, xmax, zmin and zmax for the computation domain
xmin = -3;  xmax = 3;
ymin = -3;  ymax = 3;
zmin = -6;  zmax = 0; %Note that zmax<=0

%Define numer of nodes in the x and z-direction
Nnx = 80;
Nny = 80;
Nnz = 80;

%% Check input
if(A(3)>0); error('A(2) must be 0 or negative'); end
if(B(3)>0); error('B(2) must be 0 or negative'); end
if(M(3)>0); error('M(2) must be 0 or negative'); end
if(N(3)>0); error('N(2) must be 0 or negative'); end
if(zmax>0); zmax=0; disp('Warning: zmax has been put to 0'); end
if(A(1)>xmax||A(1)<xmin||A(2)>ymax||A(2)<ymin); disp('Warning: A is outside the domain'); end
if(B(1)>xmax||B(1)<xmin||B(2)>ymax||B(2)<ymin); disp('Warning: B is outside the domain'); end
if(M(1)>xmax||M(1)<xmin||M(2)>ymax||M(2)<ymin); disp('Warning: M is outside the domain'); end
if(N(1)>xmax||N(1)<xmin||N(2)>ymax||N(2)<ymin); disp('Warning: N is outside the domain'); end

%% Computations
%point of sensitivity (x,y,z)
xs = linspace(xmin,xmax,Nnx);
ys = linspace(ymin,ymax,Nny);
zs = linspace(zmin,zmax,Nnz);

%Initilise
S_3D = zeros(4,Nnx,Nny,Nnz);

%Loop the quadrupoles
for pole = 1:4 

if(pole==1) %AM
    r = A(1); s = A(2); t = A(3); %source in (r,s,t)
    l = M(1); m = M(2); n = M(3); %Potential in (l,m,n)
elseif(pole==2) %AN
    r = A(1); s = A(2); t = A(3); %source in (r,s,t)
    l = N(1); m = N(2); n = N(3); %Potential in (l,m,n)
elseif(pole==3) %BM
    r = B(1); s = B(2); t = B(3); %source in (r,s,t)
    l = M(1); m = M(2); n = M(3); %Potential in (l,m,n)
elseif(pole==4) %BN
    r = B(1); s = B(2); t = B(3); %source in (r,s,t)
    l = N(1); m = N(2); n = N(3); %Potential in (l,m,n)
end

%Compute electric fields
for ii = 1:Nnx %x
for jj = 1:Nny %y
for kk = 1:Nnz %z 
   
    %point of sensitivity (x,y,z)
    x = xs(ii);
    y = ys(jj);
    z = zs(kk);

    %Electric field from source
    Ec_x = -1/(4*pi)*( (r-x) / ((r-x)^2+(s-y)^2+(t-z)^2)^(3/2) + ...
                       (r-x) / ((r-x)^2+(s-y)^2+(t+z)^2)^(3/2) );
    Ec_y = -1/(4*pi)*( (s-y) / ((r-x)^2+(s-y)^2+(t-z)^2)^(3/2) + ...
                       (s-y) / ((r-x)^2+(s-y)^2+(t+z)^2)^(3/2) );
    Ec_z = -1/(4*pi)*( (t-z) / ((r-x)^2+(s-y)^2+(t-z)^2)^(3/2) + ...
                       (t+z) / ((r-x)^2+(s-y)^2+(t+z)^2)^(3/2) );

    %Electroc field from potential
    Ep_x = -1/(4*pi)*( (l-x) / ((l-x)^2+(m-y)^2+(n-z)^2)^(3/2) + ...
                       (l-x) / ((l-x)^2+(m-y)^2+(n+z)^2)^(3/2) );
    Ep_y = -1/(4*pi)*( (m-y) / ((l-x)^2+(m-y)^2+(n-z)^2)^(3/2) + ...
                       (m-y) / ((l-x)^2+(m-y)^2+(n+z)^2)^(3/2) );
    Ep_z = -1/(4*pi)*( (n-z) / ((l-x)^2+(m-y)^2+(n-z)^2)^(3/2) + ...
                       (n+z) / ((l-x)^2+(m-y)^2+(n+z)^2)^(3/2) );

    %Sensitivity
    S_3D(pole,jj,ii,kk) = Ec_x*Ep_x + Ec_y*Ep_y + Ec_z*Ep_z;
    
end
end
end

end

% Superposition 
S_tot(:,:,:) = S_3D(1,:,:,:)-S_3D(2,:,:,:)-S_3D(3,:,:,:)+S_3D(4,:,:,:);

%Normalise the sensitivity
q = max(quantile(S_tot,0.95),[],'all');
S = S_tot/q;

%% Plotting

%Mesh grid
[X,Y,Z]= meshgrid(xs,ys,zs);

%Plot sensitivity as slices
figure
h=slice(X,Y,Z,S,[A(1) B(1) M(1) N(1)],[],[A(3) B(3) M(3) N(3)]);
shading flat
colorbar
colormap(S_colorbar)
caxis([-1 1])
hold on
plot3(A(1),A(2),A(3),'kx','markersize',8)
plot3(B(1),B(2),B(3),'kx','markersize',8)
plot3(M(1),M(2),M(3),'kx','markersize',8)
plot3(N(1),N(2),N(3),'kx','markersize',8)
xlabel('x'); ylabel('y'); zlabel('z');

%Plot sensitivity as iso-surfaces
figure
colorbar
colormap(S_colorbar)
caxis([-1 1])
isosurface(X,Y,Z,S,0.9); %alpha(0.95)
isosurface(X,Y,Z,S,0.5); alpha(0.9)
isosurface(X,Y,Z,S,0.2); alpha(0.6)
isosurface(X,Y,Z,S,0.1); alpha(0.3)
isosurface(X,Y,Z,S,-0.9); %alpha(0.95)
isosurface(X,Y,Z,S,-0.5); alpha(0.9)
isosurface(X,Y,Z,S,-0.2); alpha(0.6)
isosurface(X,Y,Z,S,-0.1); alpha(0.3)
grid on; hold on
plot3(A(1),A(2),A(3),'kx','markersize',8)
plot3(B(1),B(2),B(3),'kx','markersize',8)
plot3(M(1),M(2),M(3),'kx','markersize',8)
plot3(N(1),N(2),N(3),'kx','markersize',8)
xlabel('x'); ylabel('y'); zlabel('z');
