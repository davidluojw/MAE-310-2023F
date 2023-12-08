clear all; clc;

kappa = 1.0; % isotropic homogeneous heat conductivity

% manufactured solution and source term
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) -2*x*(x-1)-2*y*(y-1);

%Dirichlet BC
g = @(x, y) 0.1 * sin((x + y) * 2 * pi);

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% FEM mesh settings
n_en = 4; % 4-node quadrilateral element

n_el_x = 200;               % number of element in x-direction
n_el_y = 200;               % number of element in y-direction
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

% ID array
ID = zeros(n_np, 1);
counter = 1;
for ny = 2 : n_np_y-1
  for nx = 2 : n_np_x-1
    ID( (ny-1)*n_np_x + nx ) = counter;
    counter = counter + 1;
  end
end

n_eq = n_np - n_np_x * 2 - n_np_y * 2 + 4;

LM = ID(IEN);

% Start the assembly procedure
K = spalloc(n_eq, n_eq, 9*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_en, n_en);
   f_ele = zeros(n_en, 1);

   x_ele = zeros(n_en, 1);
   y_ele = x_ele;
   for aa = 1 : n_en
     x_ele(aa) = x_coor( IEN(aa,ee) );
     y_ele(aa) = y_coor( IEN(aa,ee) );
   end

   % loop over quadrature points   
   for ll = 1 : n_int
     x_l = 0.0; y_l = 0.0;
     dx_dxi = 0.0; dx_deta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0;
     for aa = 1 : n_en
        x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
        y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
        [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
        dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
        dx_deta = dx_deta + x_ele(aa) * Na_eta;
        dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
        dy_deta = dy_deta + y_ele(aa) * Na_eta;
     end

     detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

     for aa = 1 : n_en
       [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
       Na_x = (Na_xi * dy_deta    - Na_eta * dy_dxi) / detJ;
       Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;

       f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Quad(aa, xi(ll), eta(ll));
       for bb = 1 : n_en
         [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
         Nb_x = (Nb_xi * dy_deta    - Nb_eta * dy_dxi) / detJ;
         Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;

         k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa *...
           ( Na_x * Nb_x + Na_y * Nb_y);
       end % end of bb-loop
     end % end of aa-loop
   end % end of quadrature loop

   % global assembly
   for aa = 1 : n_en
     PP = LM(aa, ee);
     if PP > 0
       F(PP) = F(PP) + f_ele(aa);
       for bb = 1 : n_en
         QQ = LM(bb, ee);
         if QQ > 0
           K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
         else
           % do something for non-zero g boundary condition
           F(PP) = F(PP) - k_ele(aa, bb) * g(x_ele(bb), y_ele(bb));
         end
       end
     end
   end
end % end of element loop

d_temp = K \ F;
disp = zeros(n_np, 1);

for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = d_temp(index);
  end
end

% plot the solution
[X, Y] = meshgrid( 0:hh_x:1, 0:hh_y:1 );
Z = reshape(disp, n_np_x, n_np_y);
surf(X, Y, Z');

% postprocess the solution by calculating the error measured in L2 norm
errorL2 = 0.0; bottomL2 = 0.0;
errorH1 = 0.0; bottomH1 = 0.0;
for ee = 1 : n_el
  x_ele = x_coor( IEN(1:n_en, ee) );
  y_ele = y_coor( IEN(1:n_en, ee) );
  u_ele = disp(   IEN(1:n_en, ee) );

  for ll = 1 : n_int
    x_l = 0.0; y_l = 0.0; u_l = 0.0;
    u_l_xi = 0.0; u_l_eta = 0.0;
    dx_dxi = 0.0; dy_dxi = 0.0; dx_deta = 0.0; dy_deta = 0.0;
    for aa = 1 : n_en
      x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
      y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
      u_l = u_l + u_ele(aa) * Quad(aa, xi(ll), eta(ll));
      [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
      u_l_xi  = u_l_xi  + u_ele(aa) * Na_xi;
      u_l_eta = u_l_eta + u_ele(aa) * Na_eta;
      dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
      dx_deta = dx_deta + x_ele(aa) * Na_eta;
      dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
      dy_deta = dy_deta + y_ele(aa) * Na_eta;
    end
    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    u_l_x = (u_l_xi * dy_deta - u_l_eta * dy_dxi) / detJ;
    u_l_y = (u_l_xi * (-dx_deta) + u_l_eta * dx_dxi) / detJ;

    errorL2 = errorL2 + weight(ll) * detJ * (u_l - exact(x_l, y_l))^2;
    errorH1 = errorH1 + weight(ll) * detJ *...
      (( u_l_x- exact_x(x_l,y_l))^2 + ( u_l_y - exact_y(x_l,y_l))^2);
    bottomL2 = bottomL2 + weight(ll) * detJ * exact(x_l, y_l)^2;
    bottomH1 = bottomH1 + weight(ll) * detJ * (exact_x(x_l,y_l)^2 + exact_y(x_l,y_l)^2);
  end
end

errorL2 = sqrt(errorL2) / sqrt(bottomL2);
errorH1 = sqrt(errorH1) / sqrt(bottomH1);

% EOF