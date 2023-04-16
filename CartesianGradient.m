function [grad] = CartesianGradient(fx,fy,fz, ...
                                    x, y, z, ...
                                    x0, y0, z0, ...
                                    dx,dy,dz)

% Find the indices of the closest point in the meshgrid
    [~,x0_idx] = min(abs(x-x0));
    [~,y0_idx] = min(abs(y-y0));
    [~,z0_idx] = min(abs(z-z0));

    gx = gradient(fx)/dx;
    gy = gradient(fy)/dy;
    gz = gradient(fz)/dz;

    grad = [gx,gy,gz]; 

end