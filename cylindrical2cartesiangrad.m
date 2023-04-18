function gradcart = cylindrical2cartesiangrad(gradcyl, thi)
% Convert gradient from cylindrical to Cartesian coordinates

    grad_x = gradcyl(1)*cos(thi) - gradcyl(2)*sin(thi);
    grad_y = gradcyl(1)*sin(thi) + gradcyl(2)*cos(thi);
    grad_z = gradcyl(3);
    gradcart = [grad_x, grad_y, grad_z];
    gradcart = gradcart/norm(gradcart);

end

