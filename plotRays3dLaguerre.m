function plotRays3dLaguerre(rays,z_distance)

    for ray_index = 1:length(rays)
%           p  = patchline(z_distance,rays(ray_index).xH1,rays(ray_index).yH1,'linestyle','-','edgecolor','r','linewidth',3,'edgealpha',0.35);
%         w = plot3(z_distance, rays(ray_index).yH1 ,rays(ray_index).xH1,'r','LineWidth',2,'FaceAlpha',.3,'EdgeAlpha',.3);
%         he_mh = w.MarkerHandle;
% %         he_mh.FaceColorType = 'truecoloralpha';
%         he_mh.FaceColorData = uint8(255*[1;0;0;0.3]); 
% %         he.CapSize = 0;
%         w = plot3(z_distance, rays(ray_index).yH2 ,rays(ray_index).xH2,'y','LineWidth',2,'FaceAlpha',.3,'EdgeAlpha',.3);
          p  = patchline(z_distance,rays(ray_index).xH2,rays(ray_index).yH2,'linestyle','-','edgecolor','y','linewidth',3,'edgealpha',0.20535);

        alpha(.5)

    end

    
end