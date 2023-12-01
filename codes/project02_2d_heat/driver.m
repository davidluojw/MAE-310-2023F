clear all; clc;

kappa = 1.0; % isotropic homogeneous heat conductivity

% manufactured solution and source term
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) -2*x*(x-1)-2*y*(y-1);

% FEM mesh settings
n_en = 4; % 4-node quadrilateral element

n_el_x = 2;               % number of element in x-direction
n_el_y = 3;               % number of element in y-direction
n_el   = n_el_x * n_el_y; % total number of element in 2D domain

n_np_x = n_el_x + 1;      % number of node points in x-direction
n_np_y = n_el_y + 1;      % number of node points in y-direction
n_np   = n_np_x * n_np_y; % total number of node points in 2D domain

% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = 1 / n_el_x;        % mesh size in the x-direction
hh_y = 1 / n_el_y;        % mesh size in the y-direction

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index for the (nx, ny) node
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y;
  end
end

% setup the IEN array for element with local node numbering as
% a=4 ------- a=3
% |           |
% |           |
% |           |
% |           |
% a=1 ------- a=2
IEN = zeros(n_en, n_el);
for ey = 1 : n_el_y
  for ex = 1 : n_el_x
    ee = (ey-1)*n_el_x + ex;
    IEN(1,ee) = (ey-1)* n_np_x + ex;
    IEN(2,ee) = (ey-1)* n_np_x + ex + 1;
    IEN(3,ee) =  ey   * n_np_x + ex + 1;
    IEN(4,ee) =  ey   * n_np_x + ex;
  end
end










% EOF