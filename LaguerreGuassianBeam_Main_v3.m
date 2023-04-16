clear
% close all
%TODO: CHeck interval of angles on Laguerre, evaluation and ghow to
%increase in each step
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

mapgreen = AdvancedColormap('kbcw',256,[0 30 100 255]/255);  %color del haz

% Parameters for define all parameters of Laguerre
nu = 14;
mu = 10;
initialWaist = 100;
wavelength   = 0.6328;
initialZ     = 0;
units        = 'microns';

LPz0 = LaguerreParameters(initialZ, initialWaist, wavelength, nu, mu, units);
% LP give us all parameters

%paramater until last distance of propagation
k = LPz0.k;
RayleighDistance = LPz0.RayleighDistance;
finalZ = 0.5*LPz0.RayleighDistance;
LPzf = LaguerreParameters(finalZ, initialWaist, wavelength, nu, mu, units);

%% sampling of vectors 
%First, we estimate samplig in z-direction with propagation distance 
% z-direction
Dz = finalZ;  % z-window (propagation distance)
Nz = 2^8;     % number of points in z-direction
dz = Dz/Nz;   % Resolution in z
nz = 0:Nz-1;  % vector with N-points with resolution 1
z  = nz*dz;   % z-vector z of propagation 

% waist of Laguerre Gaussian Beam to distancec Dz
sigmaLZ = LPzf.laguerreWaist;

% y,x-direction
Nx  =  2^9;                 % Number of points in x,y axis
n   = -Nx/2+.05:Nx/2-1+.05; % vector with N-points with resolution 1
Dx  = (2*sigmaLZ)*1.37;     % Size of window 
dx  = Dx/Nx;                % Resolution
x   = n*dx;                 % Vector
y   = x;
[X] = meshgrid(x);

% frequency vectors
Du     = 1/dx;        
du     = 1/Dx;        
u      = n*du;         
[U]    = meshgrid(u);
% Vectores kx,ky
kx     = 2*pi*u;
[Kx]   = meshgrid(kx); 

[TH,R] = cart2pol(X,X');

Nth =  2^8;             % Number of points in x,y axis
nth = -Nth/2:Nth/2-1;   % vector with N-points with resolution 1
Dth = 2*pi;             % Size of window 
dth = Dth/Nth;          % Resolution
th  = nth*dth;          % Vector

% z-direction
Drho = Dx;        % z-window (propagation distance)
Nrho = 2^8;       % number of points in z-direction
drho = Drho/Nrho; % Resolution in z
nrho = 0:Nrho-1;  % vector with N-points with resolution 1
r  = nrho*drho;   % z-vector z of propagation 


difr = sqrt(drho^2);
% Normalization of vectors
scaleFactorX = 1/(sqrt(2)*initialWaist);
scaleFactorZ = 2/RayleighDistance;
xNormalized = scaleFactorX*x;
zNormalized = scaleFactorZ*z;

LGnm  =  LaguerreGaussBeam(nu,mu,initialWaist,wavelength,units,R,TH,0); 
XLGnm = XLaguerreGaussBeam(nu,mu,initialWaist,wavelength,units,R,TH,0); 

H1 = LGnm + 1i*XLGnm;
H2 = LGnm - 1i*XLGnm;

g = LGnm;
pxy=max(max(g));

sigmaLZ0  = LPz0.laguerreWaist; % cintura de Laguerre en z=0
radiusObs = sigmaLZ0/4.5;%%sigmaLZ0/4;%sigmaLo/4;        
%traslado
xt = 2.3*radiusObs;...1.1*radiusObs;%.15*sigmaLo;
yt = 0.01;   %xt=0;

[~,ro]      = cart2pol(X-xt,X'-yt);   
obstruction = double(ro<=radiusObs);
clear ro                         
g=g.*(1-obstruction);

total_rays = 50; 
raysH1 = get_init_rays_structure(xt,yt,radiusObs, total_rays);
raysH2 = raysH1;

prop = exp(1i*wavelength*dz*(Kx.^2+(Kx').^2)/(4*pi));

gx = zeros(Nx,length(z)); 
gy = zeros(Nx,length(z));
gg = zeros(Nx,Nz,Nx);
% Guardando campo transversal en z=0
gx(:,1)   = g(Nx/2+1,:);
gy(:,1)   = g(:,Nx/2+1);
gg(:,1,:) = g';

v = VideoWriter('LGmomentum.avi');
open(v);


for z_index = 2:length(z) 

    % ray of previous step on z
    rayH1 = raysH1(z_index-1);
    rayH2 = raysH2(z_index-1);
    zi    = z(z_index-1);
    
    fig = figure(6);
    get_plot_intensity(xNormalized, xNormalized, g, pxy);
    hold on
    plot_rays(rayH1, scaleFactorX, 'r');
    plot_rays(rayH2, scaleFactorX, 'y');
    hold off
    title([' z = ',num2str(scaleFactorZ*zi)])

    for ray_index = 1:total_rays
        
        % Ray traycing of H1
        %extract info of ray in previous step on z
        [ri, thi] = get_values_ray(rayH1, ray_index);
        
        hankeltype1 = 1;
        [H1r,H1th,H1z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                                 x, th, z, ...
                                                 ri, thi, zi, ...
                                                 units, ... 
                                                 hankeltype1);

        [grad_f] = CylindricalGradient(unwrap(angle(H1r )), ...
                                       unwrap(angle(H1th)), ...
                                       unwrap(angle(H1z )), ...
                                       x, th, z, ...
                                       ri, thi, zi, ...
                                       dx, dth, dz);

        [ri, zi, thi] = ray_tracing_cylindrical(ri, zi, thi, grad_f(1), grad_f(3), grad_f(2), drho, dz, dth);

        raysH1 = save_actual_values_on_ray(raysH1, z_index, ray_index, ...
                                           ri, thi,zi);

         % Ray traicing H2
        [ri, thi] = get_values_ray(rayH2, ray_index);
      
        [H2r,H2th,H2z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                                 x, th, z, ...
                                                 ri, thi, zi, ...
                                                 units, ... 
                                                 2);

        [grad_H2] = CylindricalGradient(unwrap(angle(H2r )), ...
                                       unwrap(angle(H2th)), ...
                                       unwrap(angle(H2z )), ...
                                       x, th, z, ...
                                       ri, thi, zi, ...
                                       dx, dth, dz);
        step = sqrt(dz^2+dx^2+dx^2);
  
        [ri, zi, thi] = ray_tracing_cylindrical(ri, zi, thi, grad_H2(1), grad_H2(3), grad_H2(2), drho, dz, dth);
        raysH2 = save_actual_values_on_ray(raysH2, z_index, ray_index, ...
                                           ri, thi,zi);

    end

    G = fftshift(fft2(g));
    g =ifft2(fftshift(G.*prop));
    gx(:,z_index)   = g(Nx/2+1,:);
    gy(:,z_index)   = g(:,Nx/2+1);
    gg(:,z_index,:) = g';

    pause(.01)
end

close(v);


function ray = save_actual_values_on_ray(ray, z_index, ray_index, r, th, z)
        ray(z_index).r(ray_index)    = r;        
        ray(z_index).th(ray_index)   = th;
        ray(z_index).x(ray_index)    = r*cos(th);
        ray(z_index).y(ray_index)    = r*sin(th);
        ray(z_index).z(ray_index)    = z;
end

function [r, th] = get_values_ray(ray,ray_index)
    r  = ray.r(ray_index);
    th = ray.th(ray_index);
end

function ray = get_default_ray(x,y,z)
    [th,r] = cart2pol(x,y);
    ray    = struct();
    ray.x  = x;
    ray.y  = y;
    ray.r  = r;
    ray.th = th;
    ray.z  = z;
end

function plot_rays(rays, scaleFactorX, color)
    for ray_index = 1:size(rays.x,2)
        xRay = scaleFactorX*rays.x(ray_index);
        yRay = scaleFactorX*rays.y(ray_index); 
        plot(xRay, yRay,'.','MarkerSize',10,'LineWidth',2,'color', color)   
    end
end

function plot_intensity(x,y,g)
    mapgreen = AdvancedColormap('kbcw',256,[0 30 100 255]/255);  %color del haz
    pcolor(x,y,abs(g).^2)
    axis square
    shading flat
    colormap(mapgreen)

end

function [resultStruct] = mergeStructure(mainStruct,struct2merge)
  if isempty(mainStruct)
    T1 = ([]);
  else
    T1=  struct2table(mainStruct);
  end

  if isempty(struct2merge)
    T2 = ([]);
  else
    T2=  struct2table(struct2merge);
  end
    T= [T1;T2];
    resultStruct = table2struct(T,"ToScalar",true);
end

function rays_structure = get_init_rays_structure(xt,yt,radiusObs,number_rays)
    rays_structure = ([]);

    for ray_index = 1:number_rays
        xi = xt + radiusObs*cos(ray_index*(2*pi)/(number_rays)); 
        yi = yt + radiusObs*sin(ray_index*(2*pi)/(number_rays));
        zi = 0;
    
        rayH1 = get_default_ray(xi, yi,zi);
        
        [rays_structure] = mergeStructure(rays_structure,rayH1);
    
    end
end

function get_plot_intensity(xNormalized, ...
                            yNormalized, ...
                            g, ...
                            pxy)

    mapgreen = AdvancedColormap('kbcw',256,[0 30 100 255]/255);  %color del haz

    pxyz   = g(1,1);
    g(1,1) = pxy;
    imagesc(xNormalized,yNormalized,abs(g))
    g(1,1) = pxyz;
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square


end