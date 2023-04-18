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
mu = 00;
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
radiusObs = sigmaLZ0/5;...4.5 sigmaLZ0/4.5;%%sigmaLZ0/4;%sigmaLo/4;        
%traslado
xt = 1.1*radiusObs;...2.5*radiusObs;... 0.0;...2.3*radiusObs;...1.1*radiusObs;%.15*sigmaLo;
yt = 0.01;   %xt=0;

[~,ro]      = cart2pol(X-xt,X'-yt);   
obstruction = double(ro<=radiusObs);
clear ro                         
g=g.*(1-obstruction);

total_rays = 9; 
hankel_type = 1;
raysH1 = get_init_rays_structure(xt, yt, radiusObs, total_rays, hankel_type);
hankel_type = 2;
raysH2 = get_init_rays_structure(xt, yt, radiusObs, total_rays, hankel_type);

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

difcart = [4*dx,4*dx,dz];

for z_index = 2:length(z) 

    % ray of previous step on z
    rayH1 = raysH1(z_index-1);
    rayH2 = raysH2(z_index-1);
    zi    = z(z_index-1);
    get_default_figure();
    fig = figure(6);
    fig.Position = [514 364 494 525];

%     AxesH = axes;
%     InSet = get(AxesH, 'TightInset');
%     set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])
    get_plot_intensity(xNormalized, xNormalized, g, pxy);
    hold on
    plot_rays(rayH1, scaleFactorX, 'r');
    plot_rays(rayH2, scaleFactorX, 'y');
    hold off
    title(['$ z/z_R = $',num2str(scaleFactorZ*zi)],Interpreter="latex")
    set(gca, 'FontSize', 14)

    path_file = ['images\Lm0_z_',num2str(scaleFactorZ*zi),'n.png'];
    print(path_file, '-dpng', '-r600')

    for ray_index = 1:total_rays
        
        % Ray traycing of H1
        %extract info of ray in previous step on z
        [ri, thi, xi, yi, hankel] = get_values_ray(rayH1, ray_index);
        pointcylindrical  = [ri,thi,zi];
        pointcart  = [xi,yi,zi];

%         hankel = 1;
        [H1r,H1th,H1z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                                 x, th, z, ...
                                                 ri, thi, zi, ...
                                                 units, ... 
                                                 hankel);

        [grad_f] = CylindricalGradient(unwrap(angle(H1r )), ...
                                       unwrap(angle(H1th)), ...
                                       unwrap(angle(H1z )), ...
                                       x, th, z, ...
                                       ri, thi, zi, ...
                                       dx, dth, dz);

        gradcart = cylindrical2cartesiangrad(grad_f, thi);
        newpoint = ray_tracing_cartesian(pointcart,gradcart,difcart);

       is_cross = is_cross_point_origin(pointcart, newpoint);
 
       if is_cross
            hankel = 2;
       else
            hankel =  hankel;
       end

        raysH1 = save_actual_values_on_ray(raysH1, z_index, ray_index, ...
                                           newpoint(1), newpoint(2), newpoint(3), hankel);

         % Ray traicing H2
        hankeltype2 = 2;
        [ri, thi, xi, yi, hankel] = get_values_ray(rayH2, ray_index);
        pointcylindrical  = [ri,thi,zi];
        pointcart  = [xi,yi,zi];     
        [H2r,H2th,H2z] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                                 x, th, z, ...
                                                 ri, thi, zi, ...
                                                 units, ... 
                                                 hankel);

        [grad_f] = CylindricalGradient(unwrap(angle(H2r )), ...
                                       unwrap(angle(H2th)), ...
                                       unwrap(angle(H2z )), ...
                                       x, th, z, ...
                                       ri, thi, zi, ...
                                       dx, dth, dz);
        gradcart = cylindrical2cartesiangrad(grad_f, thi);
        newpoint = ray_tracing_cartesian(pointcart,gradcart,difcart);
        
        raysH2 = save_actual_values_on_ray(raysH2, z_index, ray_index, ...
                                           newpoint(1), newpoint(2), newpoint(3), hankel);

    end

    G = fftshift(fft2(g));
    g =ifft2(fftshift(G.*prop));
    gx(:,z_index)   = g(Nx/2+1,:);
    gy(:,z_index)   = g(:,Nx/2+1);
    gg(:,z_index,:) = g';

    pause(.01)
end

close(v);


function is_cross = is_cross_point_origin(point, new_point)
    is_sign_changed_x = check_sign_change(point(1), new_point(1));
    is_sign_changed_y = check_sign_change(point(2), new_point(2));
    is_cross = is_sign_changed_x && is_sign_changed_y;
end

function is_sign_changed = check_sign_change(x, y)
% This function returns true if the signs of x and y are different, and
% false otherwise.

    % Compute the signs of x and y
    x_sign = sign(x);
    y_sign = sign(y);
    
    % Check if the signs are different
    if x_sign ~= y_sign
        is_sign_changed = true;
    else
        is_sign_changed = false;
    end
end

function ray = save_actual_values_on_ray(ray, z_index, ray_index, x, y, z, hankel_type)
    [th,r] =cart2pol(x,y);
    ray(z_index).r(ray_index)      = r;        
    ray(z_index).th(ray_index)     = th;
    ray(z_index).x(ray_index)      = x;
    ray(z_index).y(ray_index)      = y;
    ray(z_index).z(ray_index)      = z;
    ray(z_index).hankel(ray_index) = hankel_type;
end

function [r, th, x, y, hankel] = get_values_ray(ray,ray_index)
    r  = ray.r(ray_index);
    th = ray.th(ray_index);
    x  = ray.x(ray_index);
    y = ray.y(ray_index);
    hankel =ray.hankel(ray_index);
end


function ray = get_default_ray(x, y, z, hankel_type)
    [th,r] = cart2pol(x,y);
    ray    = struct();
    ray.x  = x;
    ray.y  = y;
    ray.r  = r;
    ray.th = th;
    ray.z  = z;
%     ray.cartpoint = [x,y,z];
    ray.hankel = hankel_type;
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

function rays_structure = get_init_rays_structure(xt,yt,radiusObs,number_rays, hankel_type)
    rays_structure = ([]);

    for ray_index = 1:number_rays
        xi = xt + radiusObs*cos(ray_index*(2*pi)/(number_rays)); 
        yi = yt + radiusObs*sin(ray_index*(2*pi)/(number_rays));
        zi = 0;
    
        ray = get_default_ray(xi, yi, zi, hankel_type);
        
        [rays_structure] = mergeStructure(rays_structure,ray);
    
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
    xlabel('$x/w_0$',Interpreter='latex')
    ylabel('$y/w_0$',Interpreter='latex')
    g(1,1) = pxyz;
    colormap(mapgreen)
    set(gca,'YDir','normal')
    axis square


end

function get_default_figure()

% Defaults for this blog post
width = 8;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 30;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% set(0,'FontSize',     fsz);
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1)-40, defpos(2)-40, width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

end

function resize_fcn
set(gcf,'units','pixels');
set(gca,'units','pixels');
w_pos = get(gcf, 'position');
set(gca, 'position', [0 0 w_pos(3) w_pos(4)]);
end
