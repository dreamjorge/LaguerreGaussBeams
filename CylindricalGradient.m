function [grad_f] = CylindricalGradient(frho, fphi, fz, ...
                                        r, th, z, ...
                                        ri, thi, zi, ...
                                        dr, dth, dz)



df_drho = gradient(frho)/dr; % Compute derivative of f with respect to rho
df_dphi = (1/ri)*gradient(fphi)/dth; % Compute derivative of f with respect to phi
df_dz   = gradient(fz)/dz; % Compute derivative of f with respect to z

% Find the indices of the closest point in the meshgrid
[~,r_idx]  = min(abs(r-ri));
[~,th_idx] = min(abs(th-thi));
[~,z_idx]  = min(abs(z-zi));

df_drho = df_drho (r_idx);
df_dphi = df_dphi(th_idx);
df_dz  = df_dz (z_idx);

grad_f = [df_drho, df_dphi, df_dz];
%grad_magnitude = sqrt(df_drho.^2 + df_dphi.^2 + df_dz.^2);
% grad_f = grad_f./grad_magnitude;

end