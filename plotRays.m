function plotRays(rayo,total_rays,scaleFactorX,z_index)

   for ray_index = 1:total_rays
        
        xRay = scaleFactorX*rayo(ray_index).xH1(z_index-1);
        yRay = scaleFactorX*rayo(ray_index).yH1(z_index-1); 
        plot(xRay, yRay,'+','MarkerSize',10,'LineWidth',2,'color','r')
        
        xRay = scaleFactorX*rayo(ray_index).xH2(z_index-1);
        yRay = scaleFactorX*rayo(ray_index).yH2(z_index-1); 
        plot(xRay, yRay,'+','MarkerSize',10,'LineWidth',2,'color','y')
   end


end