%Start working with this file, last implementation on structures that are
%working.
%Save 3d matrix, and rays for work in other scripts 
% rays transversal x,z it doesnt work chek in otehr scripts or information
% of rays maybe is wrong its important for verify, its missing lines of plot xz and yz rays with optical field 
% verify another files and write a new folder with this last information of
% propagation of laguerre, this is a good base for scripts repo
% 2^9 x 2^9 x 2 ^8 matrix aprox 1 gb
% find blue cyan color in scripts for plots on laguerre.


folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));



%--------Propagacion paraxial (unidades fisicas) por espectro angular-----%



mapgreen = AdvancedColormap('kbcw',256,[0 30 100 255]/255);  %color del haz
%---------------------------indices Laguerre Gauss------------------------%
nu = 20;
mu = 10;
%% Physical parameters [microns]
initialWaist     = 100;
wavelength       = 0.6328;

k                = 2*pi/wavelength;
RayleighDistance = k*initialWaist^2/2;
% El zls(elegante) es zl/2

%% sampling of vectors 
%First, we estimate samplig in z-direction with propagation distance 
% z-direction
Dz = 0.5*RayleighDistance; % z-window (propagation distance)
Nz = 2^7;                  % number of points in z-direction
dz = Dz/Nz;                % Resolution in z
nz = 0:Nz-1;               % vector with N-points with resolution 1
z  = nz*dz;                % z-vector z of propagation 

%Second, we estimate sampling in x,y-direction in terms of waist of Laguerre
% Guassian Beam

% waist of Gaussian Beam to Dz distance on z
waistZ  = initialWaist*sqrt(1+Dz.^2/RayleighDistance^2);
% waist of Laguerre Gaussian Beam to distancec Dz
sigmaLZ = waistZ*sqrt((2*nu+mu+1));

% y,x-direction
Nx  =  2^9;                % Number of points in x,y axis
n   = -Nx/2+.05:Nx/2-1+.05; % vector with N-points with resolution 1
Dx  = (2*sigmaLZ)*1.37;      % Size of window 
dx  = Dx/Nx;                % Resolution
x   = n*dx;                 % Vector
y   = x;
[X] = meshgrid(x);

% Muestreo del vector de frecuencia
Du     = 1/dx;         % Tamaño de la ventana de Fourier
du     = 1/Dx;         % Resolución
u      = n*du;         % Vector
[U]    = meshgrid(u);
% Vectores kx,ky
kx     = 2*pi*u;
[Kx]   = meshgrid(kx); 
% Coordenadas polares para el Laguerre
[TH,R] = cart2pol(X,X');

%%
Nth =  2^10;                % Number of points in x,y axis
nth   = -Nth/2:Nth/2-1; % vector with N-points with resolution 1
Dth = 2*pi;      % Size of window 
dth = Dth/Nth;                % Resolution
th  = nth*dth;                 % Vector


%%
% Normalization of vectors

scaleFactorX = 1/(sqrt(2)*initialWaist);
scaleFactorZ = 2/RayleighDistance;
xNormalized = scaleFactorX*x;
zNormalized = scaleFactorZ*z;


%-----------------Laguerre Gauss sin ostrucción en z=0--------------------%
zDistance = 0;

% Las dos soluciones del Laguerre Gauss
Ln = exp(1i*mu*TH).*  laguerregz(nu,mu,initialWaist,RayleighDistance,R,zDistance); 
Xn = exp(1i*mu*TH).*pxlaguerregz(nu,mu,initialWaist,RayleighDistance,R,zDistance);

% funciones Hankel
H1 = Ln+1i*Xn;
H2 = Ln-1i*Xn; 

% Función a propagar (Recuerde normalizar respecto al máximo de la función
% en cuestión)
g = Ln;
% Graficando Laguerre
figure(1)
pcolor(xNormalized,xNormalized,abs(g).^2)
axis square
shading flat
colormap(mapgreen)
% xlabel('x') 
% ylabel('y') 
pxy=max(max(g));

%---------------Funcion a propagar con obstrucción en z=0-----------------%
sigmaLZ0  = initialWaist*sqrt((2*nu+mu+1)); % cintura de Laguerre en z=0
radiusObs = sigmaLZ0/4;%%sigmaLZ0/4;%sigmaLo/4;          % tamaño de la obstrucción (radio de la obstruccion)
%traslado
xt = 1.1*radiusObs;...1.1*radiusObs;%.15*sigmaLo;
yt = 0;   %xt=0;

[~,ro]      = cart2pol(X-xt,X'-yt);   %aplicando la traslación a coordendas xy
obstruction = double(ro<=radiusObs);  %creando la obstrucción   
clear ro                              %limpiando coordenadas trasladadas
% Función a propagar con obstrucción
g=g.*(1-obstruction);

figure(2)
pcolor(xNormalized,xNormalized,abs(g))
axis square
shading flat
colormap(mapgreen)
axis1=gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)

%---------------calculo de los rayos en el punto (rx,z=0)-----------------%
%--------construccion de rayos por ciclo for
pn=50;        %numero de puntos que queremos muestrear alrededor del circulo

%structures para rayos ahi guardaremos la informacion para cada rayo
rayo = ([]);
for ray_index = 1:pn
    %matrices de ceros para colocar las coordenadas x,y de las Hankel
    rayo(ray_index).xH1 = zeros(1,length(z));
    rayo(ray_index).xH2 = zeros(1,length(z));
    rayo(ray_index).yH1 = zeros(1,length(z));
    rayo(ray_index).yH2 = zeros(1,length(z));
end

%para el gradiente en z=0 las funciones gradiente z=cte no cambian para
%cualquier rayo por lo que las calculamos antes del ciclo for
zDistance = z(1);


fTH = exp(1i*mu*TH);


for ray_index = 1:pn
    %en que punto de la circunferencia estamos en x,y
    xi = xt + radiusObs*cos(ray_index*(2*pi)/(pn)); 
    yi = yt + radiusObs*sin(ray_index*(2*pi)/(pn));
    zi = 0;
    %guardando los puntos de la circunferencia asociado a cada rayo
    rayo(ray_index).xc = xi; 
    rayo(ray_index).yc = yi;
    %calculando el radio al origen de cada punto
    [thi,ri]=cart2pol(xi,yi);
    %calculo de las pendientes en el punto (r,z)  Hankel 1
    % para r=cte;
    %azimuth part
    fth  = exp(1i*mu*thi);

    hankeltype = 1;
    [H1r,H1th,H1z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,x,th,z,ri,thi,zi,hankeltype); 
    % pendiente para el frente de onda para cada punto y guardandola
    [mzrH1,mzthH1,mrthH1] = gradientrthz(unwrap(angle(H1r )), ...
                                         unwrap(angle(H1th)), ...
                                         unwrap(angle(H1z )),k,dx,dth,dz,ri,thi,zi);
    rayo(ray_index).mH1 = mzrH1;
    rayo(ray_index).mzthH1 = mzthH1;
    rayo(ray_index).mrthH1 = mrthH1;


    %calculo de las pendientes en el punto (r,z)  Hankel 2


    hankeltype = 2;
    [H2r,H2th,H2z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,x,th,z,ri,thi,zi,hankeltype); 
    
    [mzrH2,mzthH2,mrthH2] = gradientrthz(unwrap(angle(H2r)), ...
                                         unwrap(angle(H2th)), ...
                                         unwrap(angle(H2z)),k,dx,dth,dz,ri,thi,zi); 
    rayo(ray_index).mH2 = mzrH2; 
    rayo(ray_index).mzthH2 = mzthH2; 
    rayo(ray_index).mrthH2 = mrthH2; 


end

% Graficando los puntos que se propagaran
figure(3)
pcolor(xNormalized,xNormalized,abs(g))
axis square
shading flat
colormap(mapgreen)
axis1 = gca;
set(axis1,'FontSize',13);
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$y$','Interpreter','latex','FontSize',18)
hold on
for ray_index = 1:pn
    plot(scaleFactorX*rayo(ray_index).xc,scaleFactorX*rayo(ray_index).yc, '+','MarkerSize',10,'LineWidth',2,'color','y')
end
hold off

% Condiciones iniciales para el loop, las dos Hankel parten del mismo r
% rH1=zeros(1,pn);
for ray_index = 1:pn
    %guardando los datos de los rayos de sus componentes x,y
    rayo(ray_index).xH1(1) = rayo(ray_index).xc;
    rayo(ray_index).xH2(1) = rayo(ray_index).xc;
    rayo(ray_index).yH1(1) = rayo(ray_index).yc;
    rayo(ray_index).yH2(1) = rayo(ray_index).yc;
    %calculando el radio para cada punto
    rayo(ray_index).rH1 = sqrt((rayo(ray_index).xH1(1))^2+(rayo(ray_index).yH1(1))^2);
    rayo(ray_index).rH2 = rayo(ray_index).rH1;    %los radios debidos a Hankel 1 y 2 parten del mismo punto
    %para pasar de coordenadas polares a cartesianas requerimos el angulo
    rayo(ray_index).thH1 = atan2(rayo(ray_index).yc,rayo(ray_index).xc); %atan2 da el signo correcto 
    rayo(ray_index).thH2 = atan2(rayo(ray_index).yc,rayo(ray_index).xc); %atan2 da el signo correcto 

end
%-----------------------------fin de rayos--------------------------------%

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
    plotRays(rayo,pn,scaleFactorX,z_index)
    text(-6,6,[' z = ',num2str(scaleFactorZ*z(z_index-1))],'Color','y','FontSize',16)
    hold off
    g(1,1) = pxyz;
%     writeVideo(vidObj1, getframe(gca));
    %-------------------------Calculo de Rayos----------------------------%   
    % Dado z encontramos z+dz y la posición de r que le corresponde a este
    % punto z+dz con la pendiente del frente de onda que hay entre estos 
    % dos puntos

    for ray_index = 1:pn
        
        mH1 = rayo(ray_index).mH1;
        mzthH1 = rayo(ray_index).mzthH1;
        mrthH1 = rayo(ray_index).mrthH1;

        % increase r
        rH1 = rayo(ray_index).rH1;
        rH1 = ((1/mH1)*(dz)+rH1);
        % new r on ray
        rayo(ray_index).rH1 = rH1;        

        %change th
        thH1 = rayo(ray_index).thH1;
%         thH1 = ((1/mzthH1)*(dz)+thH1);
%         if thH1 > pi
%             thH1 = thH1 - pi;
%         elseif thH1 < pi
%             thH1 = thH1 + pi;
%         end
        
        rayo(ray_index).thH1 = thH1;

        % Reescribiendo los valores nuevos de x,y dado este radio y
        % guardandolos en los rayos
        rayo(ray_index).xH1(z_index) = rH1*cos(thH1);
        rayo(ray_index).yH1(z_index) = rH1*sin(thH1);
        %-----------Calculando las pendientes del siguiente paso----------%
        %calculo de la pendiente en (r,z(ii)) de Hankel 1
        % para r=cte;

        hankeltype = 1;
        [H1r,H1th,H1z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,x,th,z,rH1,thH1,zDistance,hankeltype); 

        [mzrH1,mzthH1,mrthH1] = gradientrthz(unwrap(angle(H1r )), ...
                                             unwrap(angle(H1th)), ...
                                             unwrap(angle(H1z )),k,dx,dth,dz,ri,thi,zDistance);
        mH1 = mzrH1;
        rayo(ray_index).mH1 = mH1;
        rayo(ray_index).mzthH1 = mzthH1;

        %calculo de la pendiente en (r,z(ii)) de Hankel 2
        %Aqui el rayo que entra debido a Hankel 1 despues se rige por Hankel 2
        %cuando pasa por el origen, para ello calculamos el angulo con
        %atan2 para saber en que cuadrante estamos y si así paso por el
        %origen
        
        mH2 = rayo(ray_index).mH2;
        mzthH2 = rayo(ray_index).mzthH2;

        % increase r
        rH2 = rayo(ray_index).rH2;
        rH2 = ((1/mH2)*(dz)+rH2);
        % new r on ray
        rayo(ray_index).rH2 = rH2;
        
        %change th
        thH2 = rayo(ray_index).thH2;
%         thH2 = ((1/mzthH2)*(dz)+thH2);
%         
%         if thH2 > pi
%             thH2 = thH2 - pi;
%         elseif thH1 < pi
%             thH2 = thH2 + pi;
%         end


        rayo(ray_index).thH2 = thH2;


       

        rayo(ray_index).xH2(z_index)=rH2*cos(thH2);
        rayo(ray_index).yH2(z_index)=rH2*sin(thH2);
        
        
        if  rH2<0
            hankeltype = 1;
            [H2r,H2th,H2z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,x,th,z,rH2,thH2,zDistance,hankeltype); 
        else
            hankeltype = 2;
            [H2r,H2th,H2z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,x,th,z,rH2,thH2,zDistance,hankeltype); 
        end

        [mzrH2,mzthH2,mrthH2] = gradientrthz(unwrap(angle(H2r)), ...
                                             unwrap(angle(H2th)), ...
                                             unwrap(angle(H2z)),k,dx,dth,dz,ri,thi,zDistance); 

        mH2 = mzrH2;
        rayo(ray_index).mH2 = mH2;
        rayo(ray_index).mzthH2 = mzthH2;

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

%graficando el campo en z
pxyz=g(1,1);
g(1,1)=pxy;
fig=figure(6);
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
    for ray_index=1:pn
        xRay = scaleFactorX*rayo(ray_index).xH1(z_index);
        yRay = scaleFactorX*rayo(ray_index).yH1(z_index); 
        plot(xRay, yRay,'+','MarkerSize',10,'LineWidth',2,'color','r')
        
        xRay = scaleFactorX*rayo(ray_index).xH2(z_index);
        yRay = scaleFactorX*rayo(ray_index).yH2(z_index); 
        plot(xRay, yRay,'+','MarkerSize',10,'LineWidth',2,'color','y')
    end
% %     text(-5,5,[' z = ',num2str(scaleFactorZ*Dz)],'Color','y','FontSize',16)
hold off
g(1,1)=pxyz;

% writeVideo(vidObj1, getframe(gca));
% close(vidObj1);
%%
figure(9)
% set(gca,'un','n','pos',[0,0,1,1])
pxy=gy(1,1); gy(1,1)=1; % Para H2 (Recuerde normalizar respecto al máximo de la función en cuestión)
% Para nu=2
pcolor(zNormalized,xNormalized,abs(gx))
% % Para nu=12
% pcolor(z/(zl/2),x/(sqrt(2)*wo),abs(gx))
shading interp
caxis([0,5])
caxis([0 2.75])
caxis([0 3.3])
caxis([0 2])
caxis([0 1.2])
% % Para Pacific Rim
% caxis([0 2.25])
% caxis([0 1.2])
% brighten(0.3)
% axis off
% fig.MenuBar='none';
%     fig.Position=[1,1,1354.0,733];
colormap(jet)
colormap(mapgreen)
hold on 
xlabel('$z/zR$','interpreter','latex','FontSize',20)
ylabel('$x/w_0$','interpreter','latex','FontSize',30)
ax = gca;
ax.FontSize = 18;

plot(scaleFactorZ*z,scaleFactorX*rayo(1).xH2,'linewidth',3,'color','y')
plot(scaleFactorZ*z,scaleFactorX*rayo(25).xH2,'linewidth',3,'color','y')

hold off
