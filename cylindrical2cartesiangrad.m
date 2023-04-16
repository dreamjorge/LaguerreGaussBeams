function [grad_x, grad_y, grad_z] = cylindrical2cartesiangrad(grad_rho, grad_phi, grad_z, x, y, z)
% Convert gradient from cylindrical to Cartesian coordinates

% Compute rho and phi from x, y
rho = sqrt(x.^2 + y.^2);
phi = atan2(y, x);

% Compute conversion factors
cos_phi = cos(phi);
sin_phi = sin(phi);

% Convert gradient to Cartesian coordinates
grad_x = grad_rho .* cos_phi - grad_phi .* sin_phi;
grad_y = grad_rho .* sin_phi + grad_phi .* cos_phi;
grad_z = grad_z;
end