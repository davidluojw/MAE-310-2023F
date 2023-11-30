clear all; clc;

kappa = 1.0; % isotropic homogeneous heat conductivity

% manufactured solution
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) -2*x*(x-1)-2*y*(y-1);

% FEM mesh settings
n_en = 4;

n_el_x = 2;
n_el_y = 3;
n_el   = n_el_x * n_el_y;

n_np_x = n_el_x + 1;
n_np_y = n_el_y + 1;
n_np   = n_np_x * n_np_y;

x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = 1 / n_el_x;
hh_y = 1 / n_el_y;

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx;
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y;
  end
end

% setup of IEN array
IEN = zeros(n_en, n_el);
for ey = 1 : n_el_y
  for ex = 1 : n_el_x
    ee = (ey-1)*n_el_x + ex;
    IEN(1,ee) = (ey-1)*n_np_x + ex;
    IEN(2,ee) = (ey-1)*n_np_x + ex + 1;
    IEN(3,ee) = ey * n_np_x + ex + 1;
    IEN(4,ee) = ey * n_np_x + ex;
  end
end










% EOF