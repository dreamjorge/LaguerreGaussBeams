close(figure(100))
figure (100)
% mapgreen = AdvancedColormap('kbcw',256,[0 245 250 255]/255);
mapgreen = AdvancedColormap('kbcw',256,[0 5 30 255]/255);  %color del haz
% mapgreen = AdvancedColormap('kbcw',256,[0 60 200 255]/255); 
bb = gg;
%iindex = find(abs(bb)<=.0599585);
iindex = find(abs(bb)<=.00599585);

bb(iindex) = NaN;


hold on 
plotRays3dLaguerre(rayo,z)

daspect([0.7,.1,.1]) 
daspect([1.1,.1,.1]) 

axis tight 
% view(36,20) 
% view(45,1) 
view(33,2) 

camzoom(1.4) 
% camproj perspective
axis off
% 
% slz = slice(z,x,y,abs(bb),[],[0],[]); 
% 
% rotate(slz,[-1,0,0],-45)
% xd = get(slz,'XData');
% yd = get(slz,'YData');
% zd = get(slz,'ZData');
% h = slice(z,x,y,abs(bb),xd,yd,zd);
% h.FaceColor = 'interp';
% h.EdgeColor = 'none';
% h.DiffuseStrength = 0.8;
%%

% zprop  = [0, (1/5)*RayleighDistance, (1/4)*RayleighDistance, (1/3)*RayleighDistance, 1.4*(1/3)*RayleighDistance];
zprop  = [0, 0.75*Dz/3 , 1.25*Dz/2, 0.99*Dz];
% zprop  = [0, 0.8*(1/5)*RayleighDistance, (1/3)*RayleighDistance, 1.48*(1/3)*RayleighDistance];

% hzprop = slice(z,x,y,abs(bb),[],[0.01],[]); 
% colormap(mapgreen)
% 
% % plotRays3dLaguerre(rayo,z)
% hzprop(1).FaceColor = 'interp'; 
% hzprop(1).EdgeColor = 'none'; 
%
% hold off

%%
hzprop = slice(z,x,y,abs(bb),zprop,[0],[]); 
hzprop(1).FaceColor = 'interp'; 
hzprop(2).FaceColor = 'interp';
hzprop(3).FaceColor = 'interp';
hzprop(4).FaceColor = 'interp';
hzprop(5).FaceColor = 'interp';
hzprop(5).FaceAlpha = '.6';

hzprop(1).EdgeColor = 'none'; 
hzprop(2).EdgeColor = 'none'; 
hzprop(3).EdgeColor = 'none'; 
hzprop(4).EdgeColor = 'none'; 
hzprop(5).EdgeColor = 'none'; 

colormap(mapgreen)

% plotRays3dLaguerre(rayo,z)

hold off