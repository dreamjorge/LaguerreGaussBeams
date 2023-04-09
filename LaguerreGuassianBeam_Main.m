clear
close all
%TODO: CHeck interval of angles on Laguerre, evaluation and ghow to
%increase in each step
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

mapgreen = AdvancedColormap('kbcw',256,[0 30 100 255]/255);  %color del haz

% Parameters for define all parameters of Laguerre
nu = 20;
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
Nz = 2^7;     % number of points in z-direction
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

% Coordenadas polares para el Laguerre
[TH,R] = cart2pol(X,X');

%
Nth =  2^10;            % Number of points in x,y axis
nth   = -Nth/2:Nth/2-1; % vector with N-points with resolution 1
Dth = 2*pi;             % Size of window 
dth = Dth/Nth;          % Resolution
th  = nth*dth;          % Vector

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
% Graficando Laguerre
figure(1)
plot_intensity(x,y,g)
pxy=max(max(g));

figure(2)
plot_intensity(x,y,H1)

%---------------Funcion a propagar con obstrucción en z=0-----------------%
sigmaLZ0  = LPz0.laguerreWaist; % cintura de Laguerre en z=0
radiusObs = sigmaLZ0/4;%%sigmaLZ0/4;%sigmaLo/4;        
%traslado
xt = 1.1*radiusObs;...1.1*radiusObs;%.15*sigmaLo;
yt = 0.01;   %xt=0;

[~,ro]      = cart2pol(X-xt,X'-yt);   
obstruction = double(ro<=radiusObs);
clear ro                         
% Función a propagar con obstrucción
g=g.*(1-obstruction);

figure(3)
plot_intensity(x,y,g)

pn=50; 
raysH1 = ([]);
raysH2 = ([]);

for ray_index = 1:pn
    xi = xt + radiusObs*cos(ray_index*(2*pi)/(pn)); 
    yi = yt + radiusObs*sin(ray_index*(2*pi)/(pn));
    zi = 0;
    [thi,ri] = cart2pol(xi,yi);

    rayH1 = get_default_ray(xi, yi);
    rayH2 = rayH1;

    hankeltype = 1;

    [H1r,H1th,H1z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                             x, th, z, ...
                                             ri, thi, zi, ...
                                             units, ... 
                                             hankeltype);

    % pendiente para el frente de onda para cada punto y guardandola
    [mzrH1,mzthH1,mrthH1] = gradientrthz(unwrap(angle(H1r )), ...
                                         unwrap(angle(H1th)), ...
                                         unwrap(angle(H1z )),k,dx,dth,dz,ri,thi,zi);
    rayH1.mzr  = mzrH1;
    rayH1.mzth = mzthH1;
    rayH1.mrth = mrthH1;

    hankeltype = 2;

    [H2r,H2th,H2z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                             x, th, z, ...
                                             ri, thi, zi, ...
                                             units, ... 
                                             hankeltype);

    [mzrH2,mzthH2,mrthH2] = gradientrthz(unwrap(angle(H2r)), ...
                                         unwrap(angle(H2th)), ...
                                         unwrap(angle(H2z)),k,dx,dth,dz,ri,thi,zi); 
    rayH2.mzr  = mzrH2; 
    rayH2.mzth = mzthH2; 
    rayH2.mrth = mrthH2; 
    
    [raysH1] = mergeStructure(raysH1,rayH1);
    [raysH2] = mergeStructure(raysH2,rayH2);

end

%-------------------------Propagagación física----------------------------%
% Propagador paraxial
prop = exp(1i*wavelength*dz*(Kx.^2+(Kx').^2)/(4*pi));
figure(5)
imagesc(u,u,(angle(prop)))
title('Propagador')
% Matrices de campos transversales para guardar los datos
gx = zeros(Nx,length(z)); 
gy = zeros(Nx,length(z));
gg = zeros(Nx,Nz,Nx);
% Guardando campo transversal en z=0
gx(:,1)   = g(Nx/2+1,:);
gy(:,1)   = g(:,Nx/2+1);
gg(:,1,:) = g';

z_check = 0.8*(1/5)*RayleighDistance;..., (1/3)*RayleighDistance

for z_index = 2:length(z) %corriendo todos los valores de sp
    zDistance = z(z_index);
    %graficando campo en z(ii-1)
    pxyz   = g(1,1);
    g(1,1) = pxy;
    fig = figure(6);
    set(gca,'un','n','pos',[0,0,1,1])
    imagesc(xNormalized,xNormalized,abs(g))
    axis off
    fig.MenuBar='none';
    fig.Position=[1,1,1354.0,733];
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square
    hold on
    %puntos propagados debido a Hankel 1 y 2
    plot_rays(raysH1, pn, scaleFactorX, z_index, 'r');
    plot_rays(raysH2, pn, scaleFactorX, z_index, 'y');
%     plotRays(rayo,pn,scaleFactorX,z_index)
    text(-6,6,[' z = ',num2str(scaleFactorZ*z(z_index-1))],'Color','y','FontSize',16)
    hold off
    g(1,1) = pxyz;

    % ray of previous step on z
    rayH1 = raysH1(z_index-1);
    rayH2 = raysH2(z_index-1);

    for ray_index = 1:pn
        %extract info of ray in previous step on z
        [rH1, thH1, mzrH1, mzthH1, mrthH1] = get_values_ray(rayH1, ray_index);

        % increase values
        rH1 = get_new_r(rH1,mzrH1,dz);
        thH1 = get_new_th(thH1,mzrH1,dth);
        
        hankeltype1 = 1;
        [H1r,H1th,H1z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, RayleighDistance, ...
                                                 x, th, z, ...
                                                 ri, thi, zi, ...
                                                 units, ... 
                                                 hankeltype1);

        [mzrH1,mzthH1,mrthH1] = gradientrthz(unwrap(angle(H1r )), ...
                                             unwrap(angle(H1th)), ...
                                             unwrap(angle(H1z )),k,dx,dth,dz,ri,thi,zDistance);
        %save values on actual z

        raysH1 = save_actual_values_on_ray(raysH1, z_index, ray_index, ...
                                           rH1, thH1, ...
                                           mzrH1, mzthH1, mrthH1);

        % Ray traicing H2
        [rH2, thH2, mzrH2, mzthH2, mrthH2] = get_values_ray(rayH2, ray_index);

        % increase r
        rH2  = get_new_r(rH2,mzrH2,dz);
        thH2 = get_new_th(thH2,mrthH2,dth);
        
        if  rH2<0
            hankeltype2 = 1;
            [H2r,H2th,H2z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, RayleighDistance, ...
                                                     x, th, z, ...
                                                     ri, thi, zi, ...
                                                     units, ... 
                                                     hankeltype2);

        else
            hankeltype2 = 2;
            
            [H2r,H2th,H2z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, RayleighDistance, ...
                                                     x, th, z, ...
                                                     ri, thi, zi, ...
                                                     units, ... 
                                                     hankeltype2);

        end

        [mzrH2,mzthH2,mrthH2] = gradientrthz(unwrap(angle(H2r)), ...
                                             unwrap(angle(H2th)), ...
                                             unwrap(angle(H2z)),k,dx,dth,dz,ri,thi,zDistance); 


        raysH2 = save_actual_values_on_ray(raysH2, z_index, ray_index, ...
                                           rH2, thH2, ...
                                           mzrH2, mzthH2, mrthH2);

    end
    %-----------------------Fin de Calculo de Rayos-----------------------%   
    %propagacion del campo
    G = fftshift(fft2(g));
    %obteniendo el campo propagado
    g =ifft2(fftshift(G.*prop));
    %guardando campo transversal
    gx(:,z_index)   = g(Nx/2+1,:);
    gy(:,z_index)   = g(:,Nx/2+1);
    gg(:,z_index,:) = g';

    pause(.01)
end


function r = get_new_r(r,m,dz)
        r = (1/m)*(dz)+r;
end

function th = get_new_th(th,m,dx)
    %TODO: check how to increase th
        th = -(1/m)*(dx)+th;

end

function ray = save_actual_values_on_ray(ray,z_index,ray_index,r,th,mzr,mzth,mrth)
        ray(z_index).r(ray_index)    = r;        
        ray(z_index).th(ray_index)   = th;
        ray(z_index).x(ray_index)    = r*cos(th);
        ray(z_index).y(ray_index)    = r*sin(th);
        ray(z_index).mzr(ray_index)  = mzr;
        ray(z_index).mzth(ray_index) = mzth;
        ray(z_index).mrth(ray_index) = mrth;
end

function [r, th, mzr, mzth, mrth] = get_values_ray(ray,ray_index)
        mzr  = ray.mzr(ray_index);
        mzth = ray.mzth(ray_index);
        mrth = ray.mrth(ray_index);
        r    = ray.r(ray_index);
        th   = ray.th(ray_index);

end

function ray = get_default_ray(x,y)
    [th,r] = cart2pol(x,y);

    ray = struct();
    
    ray.x  = x;
    ray.y  = y;
    ray.r  = r;
    ray.th = th;
end

function plot_rays(rays, total_rays, scaleFactorX, z_index, color)
    for ray_index = 1:total_rays
        xRay = scaleFactorX*rays(z_index-1).x(ray_index);
        yRay = scaleFactorX*rays(z_index-1).y(ray_index); 
        plot(xRay, yRay,'+','MarkerSize',10,'LineWidth',2,'color', color)   
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
