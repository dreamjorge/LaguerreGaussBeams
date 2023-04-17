function newpoint = ray_tracing_cartesian(pointcart,gradcart,difcart) 


    newpoint = pointcart - gradcart.*difcart;
%     x_new = x - gradcart(1)*df;
%     y_new = y - gradcart(2)*dy;
%     z_new = z - gradcart(3)*dz;

end

