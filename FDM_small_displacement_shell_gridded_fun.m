function [s1,s2,s12,z_m,R1,R2] = FDM_small_displacement_shell_gridded_fun(x,y,z,h,R1,R2,gamma,P,vv_o,a,b,plotting)

% This function is written to work on the gridded data extracted from the DIC
%
% Inputs:   x, y and z are 2D arrays defining the surface.
%           h can be 1D array or a single value defining thickness of the
%           surface. If h_inp is given as array, the size of h should
%           match x, y and z.
%           R1 and R2 are the 2D arrays containing principal radii of
%           curvature.
%           gamma is the 2D arrays containing anlge representing the direction of 
%           pricipal curvature in Cartesian coordinate system.
%           P is the uniform normal pressure.
%           vv_o is the Poisson's ratio of the matrial.
%           a and b are the major and minor radius of the ellipse respectively.
%           Multiple graphs can be obtained by inputting 'plotting' variable 1.
%
% Outputs:  s1, s2 and s12 are the stresses of middle surface in principal
%           curvature direction, output as 2D arrays corresponding to x, y,
%           z.
%           x, y, z are the 2D arrays definng surface as gridded data.
%           R1 and R2 in output are 2D arrays for radii of principal curvature after smoothing.
%           z_m is the z coordinate of the appex.

global  dx  dy

%% Input

E = 1e6;

if length(h) == 1               % thickness calculation
    h = ones(size(x)).*h;
end

[n,m] = size(z);
n=n-5;
m=m-5;

x1 = -2*dx:dx:(n+2)*dx;
y1 = -2*dy:dy:(m+2)*dy;

xc = 0:0.0001:a;
yc = sqrt((1-xc.^2/a^2)*b^2);


vv = vv_o;

z(:,1) = z(:,5);
z(:,2) = z(:,4);
 
z(1,:) = z(5,:);
z(2,:) = z(4,:);



%% Finding Principal directions

[dzdy,dzdx] = gradient(z,dy,dx);

for j=1:length(y1)-2
    for i=1:length(x1)
        ex(i,j) = x(i,j)-sqrt((1-y(i,j).^2/b^2)*a^2);
        ey(i,j) = y(i,j)-sqrt((1-x(i,j).^2/a^2)*b^2);
            if abs(ex(i,j)) < dx/2 || abs(ey(i,j)) < dy/2
                nj(j) = i;
                ni(i) = j;
            end
    end
end

nj = [nj 0 0];
ni = [ni 0 0];

c=1;
for j=1:length(y1)-2
    for i=nj(j):-1:1
        if j<length(y1)-2
            if  i+1 > nj(j+1) 
                nbc(j) = i;                               % -1 for only 0 boundary condition
                xbc1(c) = x(i,j);
                ybc1(c) = y(i,j);
                c=c+1;
            end
        elseif j==length(y1)-2
                nbc(j) = 1; 
                xbc1(c) = x(i,j);
                ybc1(c) = y(i,j);
                c=c+1;
        end
    end
end

for i=1:length(x1)
    d2zdy(i,:) = derivate2(z(i,:),y(i,:));
end

for j=1:length(y1)
    d2zdx(:,j) = derivate2(z(:,j),x(:,j));
end

[d2zdx_dy,red] = gradient(dzdx,dy,dx);

%% pricipal directions

% Geometry 

a11 = 1 + dzdx.^2;
a22 = 1 + dzdy.^2;
a12 = dzdx.*dzdy;

a3_1 = (-1).*dzdx.*(1+dzdy.^2+dzdx.^2).^(-1/2);
a3_2 = (-1).*dzdy.*(1+dzdy.^2+dzdx.^2).^(-1/2);
a3_3 = (1+dzdy.^2+dzdx.^2).^(-1/2);

th_1 = gamma;

% Dealing with singularity (optional)

% [K,H,Pmax,Pmin] = surfature(x,y,z);
%  
% term = sqrt(H.^2 - K);
% [min_term,c] = min(term(1:end-4,3));
% % clear c;
% [TF,Prom] = islocalmin(term(3:end-4,3),'MinProminence',1);  %2.5
% if max(Prom) > 0.5
%     [min_prom,c] = max(Prom);
% else
%     c = length(term(3:end-4,3));
% end
% %c = find(TF);
% c=c+2;
% 
% for j=1:m
%     for i=1:n
% %         if b/a < 0.5  && (x(i,j)==0 || y(i,j)==0)        % on for y(i,j) for low aspect ratio
% %             th_1(i,j) = 0;
% %         else
% %             
% %         end
%         if dbda(i,j)*cos(gamma(i,j)) < -1        % condition to go lower than -90 degree
%             th_1(i,j) = th_1(i,j) - pi;
%         end
%     end
% end
% 
%th_1(3,:) = 0;                  % at x=0
% if b/a < 0.5
%     th_1(:,3) = 0;              % at y=0
% else
%     th_1(1:c,3) = 0;            % at y=0
%     th_1(c+1:end,3) = -pi/2;
%     
% end

th_1(:,1) = -th_1(:,5);
th_1(:,2) = -th_1(:,4);

th_1(1,:) = -th_1(5,:);
th_1(2,:) = -th_1(4,:);

for j=1:size(z,2)
    for i=1:size(z,1)
        if th_1(i,j) > 45*pi/180         % condition to remain in possible domain
            th_2(i,j) = th_1(i,j) - pi/2;
        else
            th_2(i,j) = th_1(i,j) + pi/2;
        end
    end
end
 
for j=1:size(z,2)
    for i=1:size(z,1)
        Rot_1 = [cos(th_1(i,j)) + a3_1(i,j).^2*(1-cos(th_1(i,j)))                   ,  a3_1(i,j)*a3_2(i,j)*(1-cos(th_1(i,j)))-a3_3(i,j)*sin(th_1(i,j))  , a3_1(i,j)*a3_3(i,j)*(1-cos(th_1(i,j)))+a3_2(i,j)*sin(th_1(i,j)) ; ...
             a3_2(i,j)*a3_1(i,j)*(1-cos(th_1(i,j))) + a3_3(i,j)*sin(th_1(i,j))  ,  cos(th_1(i,j)) + a3_2(i,j).^2*(1-cos(th_1(i,j)))                 , a3_2(i,j)*a3_3(i,j)*(1-cos(th_1(i,j)))-a3_1(i,j)*sin(th_1(i,j)) ; ...
             a3_3(i,j)*a3_1(i,j)*(1-cos(th_1(i,j))) - a3_2(i,j)*sin(th_1(i,j))  ,  a3_3(i,j)*a3_2(i,j)*(1-cos(th_1(i,j)))+a3_1(i,j)*sin(th_1(i,j))  , cos(th_1(i,j)) + a3_3(i,j).^2*(1-cos(th_1(i,j))) ];
         
         Rot_2 = [cos(th_2(i,j)) + a3_1(i,j).^2*(1-cos(th_2(i,j)))                   ,  a3_1(i,j)*a3_2(i,j)*(1-cos(th_2(i,j)))-a3_3(i,j)*sin(th_2(i,j))  , a3_1(i,j)*a3_3(i,j)*(1-cos(th_2(i,j)))+a3_2(i,j)*sin(th_2(i,j)) ; ...
             a3_2(i,j)*a3_1(i,j)*(1-cos(th_2(i,j))) + a3_3(i,j)*sin(th_2(i,j))  ,  cos(th_2(i,j)) + a3_2(i,j).^2*(1-cos(th_2(i,j)))                 , a3_2(i,j)*a3_3(i,j)*(1-cos(th_2(i,j)))-a3_1(i,j)*sin(th_2(i,j)) ; ...
             a3_3(i,j)*a3_1(i,j)*(1-cos(th_2(i,j))) - a3_2(i,j)*sin(th_2(i,j))  ,  a3_3(i,j)*a3_2(i,j)*(1-cos(th_2(i,j)))+a3_1(i,j)*sin(th_2(i,j))  , cos(th_2(i,j)) + a3_3(i,j).^2*(1-cos(th_2(i,j))) ];
         
         p1 = Rot_1*[1;0;dzdx(i,j)]; 
         p1_1(i,j) = p1(1)/norm(p1);
         p1_2(i,j) = p1(2)/norm(p1);
         p1_3(i,j) = p1(3)/norm(p1);
         
         p2 = Rot_2*[1;0;dzdx(i,j)]; 
         p2_1(i,j) = p2(1)/norm(p2);
         p2_2(i,j) = p2(2)/norm(p2);
         p2_3(i,j) = p2(3)/norm(p2);
         
    end
end


a1_1 = p1_1;
a1_2 = p1_2;
a1_3 = p1_3;
a2_1 = p2_1;
a2_2 = p2_2;
a2_3 = p2_3;

%%

A1 = a1_1.*a1_1 + a1_2.*a1_2 + a1_3.*a1_3;
B1 = a2_1.*a2_1 + a2_2.*a2_2 + a2_3.*a2_3;

d3zdx = derivative_2d(x,y,z,3,0);
d3zd2x_dy = derivative_2d(x,y,z,2,1);
d3zdx_2dy = derivative_2d(x,y,z,1,2);
d3zdy = derivative_2d(x,y,z,0,3);

d4zdx = derivative_2d(x,y,z,4,0);
d4zd3x_dy = derivative_2d(x,y,z,3,1);
d4zd2x_2dy = derivative_2d(x,y,z,2,2);
d4zdx_3dy = derivative_2d(x,y,z,1,3);
d4zdy = derivative_2d(x,y,z,0,4);

dAdx = derivative_2d(x,y,A1,1,0);
dBdx = derivative_2d(x,y,B1,1,0);
dAdy = derivative_2d(x,y,A1,0,1);
dBdy = derivative_2d(x,y,B1,0,1);

dR1dx = derivative_2d(x,y,R1,1,0);
dR2dx = derivative_2d(x,y,R2,1,0);
dR1dy = derivative_2d(x,y,R1,0,1);
dR2dy = derivative_2d(x,y,R2,0,1);

d2Adx = derivative_2d(x,y,A1,2,0);
d2Bdx = derivative_2d(x,y,B1,2,0);
d2Ady = derivative_2d(x,y,A1,0,2);
d2Bdy = derivative_2d(x,y,B1,0,2);
d2Adx_dy = derivative_2d(x,y,A1,1,1);
d2Bdx_dy = derivative_2d(x,y,B1,1,1);

d2R1dx = derivative_2d(x,y,R1,2,0);
d2R2dx = derivative_2d(x,y,R2,2,0);
d2R1dy = derivative_2d(x,y,R1,0,2);
d2R2dy = derivative_2d(x,y,R2,0,2);
d2R1dx_dy = derivative_2d(x,y,R1,1,1);
d2R2dx_dy = derivative_2d(x,y,R2,1,1);

d3Adx = derivative_2d(x,y,A1,3,0);
d3Bdx = derivative_2d(x,y,B1,3,0);
d3Ady = derivative_2d(x,y,A1,0,3);
d3Bdy = derivative_2d(x,y,B1,0,3);
d3Ad2x_dy = derivative_2d(x,y,A1,2,1);
d3Bd2x_dy = derivative_2d(x,y,B1,2,1);
d3Adx_d2y = derivative_2d(x,y,A1,1,2);
d3Bdx_d2y = derivative_2d(x,y,B1,1,2);

d3R1dx = derivative_2d(x,y,R1,3,0);
d3R2dx = derivative_2d(x,y,R2,3,0);
d3R1dy = derivative_2d(x,y,R1,0,3);
d3R2dy = derivative_2d(x,y,R2,0,3);
d3R1d2x_dy = derivative_2d(x,y,R1,2,1);
d3R2d2x_dy = derivative_2d(x,y,R2,2,1);
d3R1dx_d2y = derivative_2d(x,y,R1,1,2);
d3R2dx_d2y = derivative_2d(x,y,R2,1,2);

da1_1dx = derivative_2d(x,y,a1_1,1,0);
da1_2dx = derivative_2d(x,y,a1_2,1,0);
da2_1dx = derivative_2d(x,y,a2_1,1,0);
da2_2dx = derivative_2d(x,y,a2_2,1,0);
da1_1dy = derivative_2d(x,y,a1_1,0,1);
da1_2dy = derivative_2d(x,y,a1_2,0,1);
da2_1dy = derivative_2d(x,y,a2_1,0,1);
da2_2dy = derivative_2d(x,y,a2_2,0,1);

d2a1_1dx = derivative_2d(x,y,a1_1,2,0);
d2a1_2dx = derivative_2d(x,y,a1_2,2,0);
d2a2_1dx = derivative_2d(x,y,a2_1,2,0);
d2a2_2dx = derivative_2d(x,y,a2_2,2,0);
d2a1_1dy = derivative_2d(x,y,a1_1,0,2);
d2a1_2dy = derivative_2d(x,y,a1_2,0,2);
d2a2_1dy = derivative_2d(x,y,a2_1,0,2);
d2a2_2dy = derivative_2d(x,y,a2_2,0,2);
d2a1_1dx_dy = derivative_2d(x,y,a1_1,1,1);
d2a1_2dx_dy = derivative_2d(x,y,a1_2,1,1);
d2a2_1dx_dy = derivative_2d(x,y,a2_1,1,1);
d2a2_2dx_dy = derivative_2d(x,y,a2_2,1,1);

d3a1_1dx = derivative_2d(x,y,a1_1,3,0);
d3a1_2dx = derivative_2d(x,y,a1_2,3,0);
d3a2_1dx = derivative_2d(x,y,a2_1,3,0);
d3a2_2dx = derivative_2d(x,y,a2_2,3,0);
d3a1_1dy = derivative_2d(x,y,a1_1,0,3);
d3a1_2dy = derivative_2d(x,y,a1_2,0,3);
d3a2_1dy = derivative_2d(x,y,a2_1,0,3);
d3a2_2dy = derivative_2d(x,y,a2_2,0,3);

d3a1_1d2x_dy = derivative_2d(x,y,a1_1,2,1);
d3a1_2d2x_dy = derivative_2d(x,y,a1_2,2,1);
d3a2_1d2x_dy = derivative_2d(x,y,a2_1,2,1);
d3a2_2d2x_dy = derivative_2d(x,y,a2_2,2,1);

d3a1_1dx_d2y = derivative_2d(x,y,a1_1,1,2);
d3a1_2dx_d2y = derivative_2d(x,y,a1_2,1,2);
d3a2_1dx_d2y = derivative_2d(x,y,a2_1,1,2);
d3a2_2dx_d2y = derivative_2d(x,y,a2_2,1,2);

dhdx = derivative_2d(x,y,h,1,0);
dhdy = derivative_2d(x,y,h,0,1);

d2hdx = derivative_2d(x,y,h,2,0);
d2hdy = derivative_2d(x,y,h,0,2);
d2hdx_dy = derivative_2d(x,y,h,1,1);

d3hdx = derivative_2d(x,y,h,3,0);
d3hdy = derivative_2d(x,y,h,0,3);
d3hd2x_dy = derivative_2d(x,y,h,2,1);
d3hdx_d2y = derivative_2d(x,y,h,1,2);

x = x(3:n+3,3:m+3);
y = y(3:n+3,3:m+3);
z = z(3:n+3,3:m+3);

dzdx = dzdx(3:n+3,3:m+3);
dzdy = dzdy(3:n+3,3:m+3);
d2zdx = d2zdx(3:n+3,3:m+3);
d2zdy = d2zdy(3:n+3,3:m+3);
d2zdx_dy = d2zdx_dy(3:n+3,3:m+3);


a1_1 = a1_1(3:n+3,3:m+3);
a1_2 = a1_2(3:n+3,3:m+3);
a1_3 = a1_3(3:n+3,3:m+3);
a2_1 = a2_1(3:n+3,3:m+3);
a2_2 = a2_2(3:n+3,3:m+3);
a2_3 = a2_3(3:n+3,3:m+3);
a3_1 = a3_1(3:n+3,3:m+3);
a3_2 = a3_2(3:n+3,3:m+3);
a3_3 = a3_3(3:n+3,3:m+3);

A1 = A1(3:n+3,3:m+3);
B1 = B1(3:n+3,3:m+3);
R1 = R1(3:n+3,3:m+3);
R2 = R2(3:n+3,3:m+3);

h = h(3:n+3,3:m+3);


%% Calculating cooefficients

dxdu1 = a1_1;
dydu1 = a1_2;
dxdu2 = a2_1;
dydu2 = a2_2;

d2xdu1 = da1_1dx.*dxdu1+da1_1dy.*dydu1;
d2xdu1_du2 = da1_1dy.*dydu2+dxdu2.*da1_1dx;
d2xdu2 = da2_1dy.*dydu2+dxdu2.*da2_1dx;

d2ydu1 = da1_2dx.*dxdu1+da1_2dy.*dydu1;
d2ydu1_du2 = da1_2dy.*dydu2+dxdu2.*da1_2dx;
d2ydu2 = da2_2dy.*dydu2+dxdu2.*da2_2dx;

d3xdu1 = d2a1_1dy.*dydu1.^2 + 2.*dxdu1.*dydu1.*d2a1_1dx_dy + dxdu1.^2.*d2a1_1dx + da1_1dx.*d2xdu1 + da1_1dy.*d2ydu1;
d3xd2u1_du2 = dydu2.*d2a1_1dy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2a1_1dx_dy+da1_1dx.*d2xdu1_du2+da1_1dy.*d2ydu1_du2+dxdu2.*dxdu1.*d2a1_1dx;
d3xdu1_2du2 = dydu2.^2.*d2a1_1dy+da1_1dy.*d2ydu2+d2xdu2.*da1_1dx+2.*dxdu2.*dydu2.*d2a1_1dx_dy+dxdu2.^2.*d2a1_1dx;
d3xdu2 = dydu2.^2.*d2a2_1dy+da2_1dy.*d2ydu2+d2xdu2.*da2_1dx+2.*dxdu2.*dydu2.*d2a2_1dx_dy+dxdu2.^2.*d2a2_1dx;

d3ydu1 = d2a1_2dy.*dydu1.^2+2.*dxdu1.*dydu1.*d2a1_2dx_dy+dxdu1.^2.*d2a1_2dx+da1_2dx.*d2xdu1+da1_2dy.*d2ydu1;
d3yd2u1_du2 = dydu2.*d2a1_2dy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2a1_2dx_dy+da1_2dx.*d2xdu1_du2+da1_2dy.*d2ydu1_du2+dxdu2.*dxdu1.*d2a1_2dx;
d3ydu1_2du2 = dydu2.^2.*d2a1_2dy+da1_2dy.*d2ydu2+d2xdu2.*da1_2dx+2.*dxdu2.*dydu2.*d2a1_2dx_dy+dxdu2.^2.*d2a1_2dx;
d3ydu2 = dydu2.^2.*d2a2_2dy+da2_2dy.*d2ydu2+d2xdu2.*da2_2dx+2.*dxdu2.*dydu2.*d2a2_2dx_dy+dxdu2.^2.*d2a2_2dx;

d4xdu1 = d3a1_1dy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3a1_1dx_d2y+3.*dxdu1.*d2a1_1dx.*d2xdu1+3.*d2a1_1dy.*dydu1.*d2ydu1+d2a1_1dx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3a1_1d2x_dy+dxdu1.^3.*d3a1_1dx+da1_1dx.*d3xdu1+da1_1dy.*d3ydu1;
d4xd2u1_2du2 = dydu2.^2.*d3a1_1dy.*dydu1+d2a1_1dx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2a1_1dy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3a1_1dx_d2y+da1_1dx.*d3xdu1_2du2+da1_1dy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2a1_1dx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3a1_1d2x_dy+dxdu2.^2.*dxdu1.*d3a1_1dx;
d4xdu2 = 3.*dydu2.*d2a2_1dy.*d2ydu2+dydu2.^3.*d3a2_1dy+da2_1dy.*d3ydu2+d3xdu2.*da2_1dx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2a2_1dx_dy+3.*dxdu2.*dydu2.^2.*d3a2_1dx_d2y+3.*dxdu2.*d2xdu2.*d2a2_1dx+3.*dxdu2.^2.*dydu2.*d3a2_1d2x_dy+dxdu2.^3.*d3a2_1dx;

d4ydu1 = d3a1_2dy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3a1_2dx_d2y+3.*dxdu1.*d2a1_2dx.*d2xdu1+3.*d2a1_2dy.*dydu1.*d2ydu1+d2a1_2dx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3a1_2d2x_dy+dxdu1.^3.*d3a1_2dx+da1_2dx.*d3xdu1+da1_2dy.*d3ydu1;
d4yd2u1_2du2 = dydu2.^2.*d3a1_2dy.*dydu1+d2a1_2dx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2a1_2dy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3a1_2dx_d2y+da1_2dx.*d3xdu1_2du2+da1_2dy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2a1_2dx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3a1_2d2x_dy+dxdu2.^2.*dxdu1.*d3a1_2dx;
d4ydu2 = 3.*dydu2.*d2a2_2dy.*d2ydu2+dydu2.^3.*d3a2_2dy+da2_2dy.*d3ydu2+d3xdu2.*da2_2dx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2a2_2dx_dy+3.*dxdu2.*dydu2.^2.*d3a2_2dx_d2y+3.*dxdu2.*d2xdu2.*d2a2_2dx+3.*dxdu2.^2.*dydu2.*d3a2_2d2x_dy+dxdu2.^3.*d3a2_2dx;

dAdu1 = dAdx.*dxdu1 + dAdy.*dydu1;
dAdu2 = dAdy.*dydu2+dxdu2.*dAdx;
d2Adu1 = d2Ady.*dydu1.^2+2.*dxdu1.*dydu1.*d2Adx_dy+dxdu1.^2.*d2Adx+dAdx.*d2xdu1+dAdy.*d2ydu1;
d2Adu1_du2 = dydu2.*d2Ady.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2Adx_dy+dAdx.*d2xdu1_du2+dAdy.*d2ydu1_du2+dxdu2.*dxdu1.*d2Adx;
d2Adu2 = dydu2.^2.*d2Ady+dAdy.*d2ydu2+d2xdu2.*dAdx+2.*dxdu2.*dydu2.*d2Adx_dy+dxdu2.^2.*d2Adx;
d3Adu1 = d3Ady.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3Adx_d2y+3.*dxdu1.*d2Adx.*d2xdu1+3.*d2Ady.*dydu1.*d2ydu1+d2Adx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3Ad2x_dy+dxdu1.^3.*d3Adx+dAdx.*d3xdu1+dAdy.*d3ydu1;
d3Adu1_2du2 = dydu2.^2.*d3Ady.*dydu1+d2Adx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2Ady.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3Adx_d2y+dAdx.*d3xdu1_2du2+dAdy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2Adx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3Ad2x_dy+dxdu2.^2.*dxdu1.*d3Adx;
d3Ad2u1_du2 = dydu2.*d3Ady.*dydu1.^2+(2.*dydu2.*dxdu1.*dydu1+dxdu2.*dydu1.^2).*d3Adx_d2y+d2Adx.*(2.*dxdu1.*d2xdu1_du2+dxdu2.*d2xdu1)+d2Adx_dy.*(2.*dydu1.*d2xdu1_du2+2.*dxdu1.*d2ydu1_du2+dydu2.*d2xdu1+dxdu2.*d2ydu1)+d2Ady.*(2.*dydu1.*d2ydu1_du2+dydu2.*d2ydu1)+(dydu2.*dxdu1.^2+2.*dxdu2.*dxdu1.*dydu1).*d3Ad2x_dy+dAdx.*d3xd2u1_du2+dAdy.*d3yd2u1_du2+dxdu2.*dxdu1.^2.*d3Adx;
d3Adu2 = 3.*dydu2.*d2Ady.*d2ydu2+dydu2.^3.*d3Ady+dAdy.*d3ydu2+d3xdu2.*dAdx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2Adx_dy+3.*dxdu2.*dydu2.^2.*d3Adx_d2y+3.*dxdu2.*d2xdu2.*d2Adx+3.*dxdu2.^2.*dydu2.*d3Ad2x_dy+dxdu2.^3.*d3Adx;

dBdu1 = dBdx.*dxdu1+dBdy.*dydu1;
dBdu2 = dBdy.*dydu2+dxdu2.*dBdx;
d2Bdu1 = d2Bdy.*dydu1.^2+2.*dxdu1.*dydu1.*d2Bdx_dy+dxdu1.^2.*d2Bdx+dBdx.*d2xdu1+dBdy.*d2ydu1;
d2Bdu1_du2 = dydu2.*d2Bdy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2Bdx_dy+dBdx.*d2xdu1_du2+dBdy.*d2ydu1_du2+dxdu2.*dxdu1.*d2Bdx;
d2Bdu2 = dydu2.^2.*d2Bdy+dBdy.*d2ydu2+d2xdu2.*dBdx+2.*dxdu2.*dydu2.*d2Bdx_dy+dxdu2.^2.*d2Bdx;
d3Bdu1 = d3Bdy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3Bdx_d2y+3.*dxdu1.*d2Bdx.*d2xdu1+3.*d2Bdy.*dydu1.*d2ydu1+d2Bdx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3Bd2x_dy+dxdu1.^3.*d3Bdx+dBdx.*d3xdu1+dBdy.*d3ydu1;
d3Bdu1_2du2 = dydu2.^2.*d3Bdy.*dydu1+d2Bdx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2Bdy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3Bdx_d2y+dBdx.*d3xdu1_2du2+dBdy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2Bdx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3Bd2x_dy+dxdu2.^2.*dxdu1.*d3Bdx;
d3Bd2u1_du2 = dydu2.*d3Bdy.*dydu1.^2+(2.*dydu2.*dxdu1.*dydu1+dxdu2.*dydu1.^2).*d3Bdx_d2y+d2Bdx.*(2.*dxdu1.*d2xdu1_du2+dxdu2.*d2xdu1)+d2Bdx_dy.*(2.*dydu1.*d2xdu1_du2+2.*dxdu1.*d2ydu1_du2+dydu2.*d2xdu1+dxdu2.*d2ydu1)+d2Bdy.*(2.*dydu1.*d2ydu1_du2+dydu2.*d2ydu1)+(dydu2.*dxdu1.^2+2.*dxdu2.*dxdu1.*dydu1).*d3Bd2x_dy+dBdx.*d3xd2u1_du2+dBdy.*d3yd2u1_du2+dxdu2.*dxdu1.^2.*d3Bdx;
d3Bdu2 = 3.*dydu2.*d2Bdy.*d2ydu2+dydu2.^3.*d3Bdy+dBdy.*d3ydu2+d3xdu2.*dBdx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2Bdx_dy+3.*dxdu2.*dydu2.^2.*d3Bdx_d2y+3.*dxdu2.*d2xdu2.*d2Bdx+3.*dxdu2.^2.*dydu2.*d3Bd2x_dy+dxdu2.^3.*d3Bdx;

dR1du1 = dR1dx.*dxdu1+dR1dy.*dydu1;
dR1du2 = dR1dy.*dydu2+dxdu2.*dR1dx;
d2R1du1 = d2R1dy.*dydu1.^2+2.*dxdu1.*dydu1.*d2R1dx_dy+dxdu1.^2.*d2R1dx+dR1dx.*d2xdu1+dR1dy.*d2ydu1;
d2R1du1_du2 = dydu2.*d2R1dy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2R1dx_dy+dR1dx.*d2xdu1_du2+dR1dy.*d2ydu1_du2+dxdu2.*dxdu1.*d2R1dx;
d2R1du2 = dydu2.^2.*d2R1dy+dR1dy.*d2ydu2+d2xdu2.*dR1dx+2.*dxdu2.*dydu2.*d2R1dx_dy+dxdu2.^2.*d2R1dx;
d3R1du1 = d3R1dy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3R1dx_d2y+3.*dxdu1.*d2R1dx.*d2xdu1+3.*d2R1dy.*dydu1.*d2ydu1+d2R1dx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3R1d2x_dy+dxdu1.^3.*d3R1dx+dR1dx.*d3xdu1+dR1dy.*d3ydu1;
d3R1du1_2du2 = dydu2.^2.*d3R1dy.*dydu1+d2R1dx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2R1dy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3R1dx_d2y+dR1dx.*d3xdu1_2du2+dR1dy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2R1dx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3R1d2x_dy+dxdu2.^2.*dxdu1.*d3R1dx;
d3R1d2u1_du2 = dydu2.*d3R1dy.*dydu1.^2+(2.*dydu2.*dxdu1.*dydu1+dxdu2.*dydu1.^2).*d3R1dx_d2y+d2R1dx.*(2.*dxdu1.*d2xdu1_du2+dxdu2.*d2xdu1)+d2R1dx_dy.*(2.*dydu1.*d2xdu1_du2+2.*dxdu1.*d2ydu1_du2+dydu2.*d2xdu1+dxdu2.*d2ydu1)+d2R1dy.*(2.*dydu1.*d2ydu1_du2+dydu2.*d2ydu1)+(dydu2.*dxdu1.^2+2.*dxdu2.*dxdu1.*dydu1).*d3R1d2x_dy+dR1dx.*d3xd2u1_du2+dR1dy.*d3yd2u1_du2+dxdu2.*dxdu1.^2.*d3R1dx;
d3R1du2 = 3.*dydu2.*d2R1dy.*d2ydu2+dydu2.^3.*d3R1dy+dR1dy.*d3ydu2+d3xdu2.*dR1dx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2R1dx_dy+3.*dxdu2.*dydu2.^2.*d3R1dx_d2y+3.*dxdu2.*d2xdu2.*d2R1dx+3.*dxdu2.^2.*dydu2.*d3R1d2x_dy+dxdu2.^3.*d3R1dx;

dR2du1 = dR2dx.*dxdu1+dR2dy.*dydu1;
dR2du2 = dR2dy.*dydu2+dxdu2.*dR2dx;
d2R2du1 = d2R2dy.*dydu1.^2+2.*dxdu1.*dydu1.*d2R2dx_dy+dxdu1.^2.*d2R2dx+dR2dx.*d2xdu1+dR2dy.*d2ydu1;
d2R2du1_du2 = dydu2.*d2R2dy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2R2dx_dy+dR2dx.*d2xdu1_du2+dR2dy.*d2ydu1_du2+dxdu2.*dxdu1.*d2R2dx;
d2R2du2 = dydu2.^2.*d2R2dy+dR2dy.*d2ydu2+d2xdu2.*dR2dx+2.*dxdu2.*dydu2.*d2R2dx_dy+dxdu2.^2.*d2R2dx;
d3R2du1 = d3R2dy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3R2dx_d2y+3.*dxdu1.*d2R2dx.*d2xdu1+3.*d2R2dy.*dydu1.*d2ydu1+d2R2dx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3R2d2x_dy+dxdu1.^3.*d3R2dx+dR2dx.*d3xdu1+dR2dy.*d3ydu1;
d3R2du1_2du2 = dydu2.^2.*d3R2dy.*dydu1+d2R2dx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2R2dy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3R2dx_d2y+dR2dx.*d3xdu1_2du2+dR2dy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2R2dx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3R2d2x_dy+dxdu2.^2.*dxdu1.*d3R2dx;
d3R2d2u1_du2 = dydu2.*d3R2dy.*dydu1.^2+(2.*dydu2.*dxdu1.*dydu1+dxdu2.*dydu1.^2).*d3R2dx_d2y+d2R2dx.*(2.*dxdu1.*d2xdu1_du2+dxdu2.*d2xdu1)+d2R2dx_dy.*(2.*dydu1.*d2xdu1_du2+2.*dxdu1.*d2ydu1_du2+dydu2.*d2xdu1+dxdu2.*d2ydu1)+d2R2dy.*(2.*dydu1.*d2ydu1_du2+dydu2.*d2ydu1)+(dydu2.*dxdu1.^2+2.*dxdu2.*dxdu1.*dydu1).*d3R2d2x_dy+dR2dx.*d3xd2u1_du2+dR2dy.*d3yd2u1_du2+dxdu2.*dxdu1.^2.*d3R2dx;
d3R2du2 = 3.*dydu2.*d2R2dy.*d2ydu2+dydu2.^3.*d3R2dy+dR2dy.*d3ydu2+d3xdu2.*dR2dx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2R2dx_dy+3.*dxdu2.*dydu2.^2.*d3R2dx_d2y+3.*dxdu2.*d2xdu2.*d2R2dx+3.*dxdu2.^2.*dydu2.*d3R2d2x_dy+dxdu2.^3.*d3R2dx;

dhdu1 = dhdx.*dxdu1 + dhdy.*dydu1;
dhdu2 = dhdy.*dydu2+dxdu2.*dhdx;
d2hdu1 = d2hdy.*dydu1.^2+2.*dxdu1.*dydu1.*d2hdx_dy+dxdu1.^2.*d2hdx+dhdx.*d2xdu1+dhdy.*d2ydu1;
d2hdu1_du2 = dydu2.*d2hdy.*dydu1+(dydu2.*dxdu1+dxdu2.*dydu1).*d2hdx_dy+dhdx.*d2xdu1_du2+dhdy.*d2ydu1_du2+dxdu2.*dxdu1.*d2hdx;
d2hdu2 = dydu2.^2.*d2hdy+dhdy.*d2ydu2+d2xdu2.*dhdx+2.*dxdu2.*dydu2.*d2hdx_dy+dxdu2.^2.*d2hdx;
d3hdu1 = d3hdy.*dydu1.^3+3.*dxdu1.*dydu1.^2.*d3hdx_d2y+3.*dxdu1.*d2hdx.*d2xdu1+3.*d2hdy.*dydu1.*d2ydu1+d2hdx_dy.*(3.*dydu1.*d2xdu1+3.*dxdu1.*d2ydu1)+3.*dxdu1.^2.*dydu1.*d3hd2x_dy+dxdu1.^3.*d3hdx+dhdx.*d3xdu1+dhdy.*d3ydu1;
d3hdu1_2du2 = dydu2.^2.*d3hdy.*dydu1+d2hdx_dy.*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*dydu2.*d2xdu1_du2+2.*dxdu2.*d2ydu1_du2)+d2hdy.*(d2ydu2.*dydu1+2.*dydu2.*d2ydu1_du2)+(dydu2.^2.*dxdu1+2.*dxdu2.*dydu2.*dydu1).*d3hdx_d2y+dhdx.*d3xdu1_2du2+dhdy.*d3ydu1_2du2+(d2xdu2.*dxdu1+2.*dxdu2.*d2xdu1_du2).*d2hdx+(2.*dxdu2.*dydu2.*dxdu1+dxdu2.^2.*dydu1).*d3hd2x_dy+dxdu2.^2.*dxdu1.*d3hdx;
d3hd2u1_du2 = dydu2.*d3hdy.*dydu1.^2+(2.*dydu2.*dxdu1.*dydu1+dxdu2.*dydu1.^2).*d3hdx_d2y+d2hdx.*(2.*dxdu1.*d2xdu1_du2+dxdu2.*d2xdu1)+d2hdx_dy.*(2.*dydu1.*d2xdu1_du2+2.*dxdu1.*d2ydu1_du2+dydu2.*d2xdu1+dxdu2.*d2ydu1)+d2hdy.*(2.*dydu1.*d2ydu1_du2+dydu2.*d2ydu1)+(dydu2.*dxdu1.^2+2.*dxdu2.*dxdu1.*dydu1).*d3hd2x_dy+dhdx.*d3xd2u1_du2+dhdy.*d3yd2u1_du2+dxdu2.*dxdu1.^2.*d3hdx;
d3hdu2 = 3.*dydu2.*d2hdy.*d2ydu2+dydu2.^3.*d3hdy+dhdy.*d3ydu2+d3xdu2.*dhdx+(3.*dydu2.*d2xdu2+3.*dxdu2.*d2ydu2).*d2hdx_dy+3.*dxdu2.*dydu2.^2.*d3hdx_d2y+3.*dxdu2.*d2xdu2.*d2hdx+3.*dxdu2.^2.*dydu2.*d3hd2x_dy+dxdu2.^3.*d3hdx;

C_01 	=	(1/12).*dxdu1.*R1.^(-1).*A1.^(-2).*B1.^(-1).*((-1).*dxdu2.^2.*((-2)+vv).*A1.^2+dxdu1.^2.*B1.^2).*h.^2;
C_02	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*((-1).*dxdu2.*(dxdu2.*dydu1+2.*dxdu1.*dydu2).*((-2)+vv).*A1.^2+3.*dxdu1.^2.*dydu1.*B1.^2).*h.^2;
C_03	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*((-1).*dydu2.*(2.*dxdu2.*dydu1+dxdu1.*dydu2).*((-2)+vv).*A1.^2+3.*dxdu1.*dydu1.^2.*B1.^2).*h.^2;
C_04	=	(1/12).*dydu1.*R1.^(-1).*A1.^(-2).*B1.^(-1).*((-1).*dydu2.^2.*((-2)+vv).*A1.^2+dydu1.^2.*B1.^2).*h.^2;
C_05	=	(1/12).*R1.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dxdu2.^2.*((-1)+vv).*A1.^2.*(3.*R1.^2+h.^2)+dxdu1.^2.*B1.^2.*(12.*R1.^2+h.^2));
C_06	=	(1/6).*R1.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dxdu2.*dydu2.*((-1)+vv).*A1.^2.*(3.*R1.^2+h.^2)+dxdu1.*dydu1.*B1.^2.*(12.*R1.^2+h.^2));
C_07	=	(1/12).*R1.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dydu2.^2.*((-1)+vv).*A1.^2.*(3.*R1.^2+h.^2)+dydu1.^2.*B1.^2.*(12.*R1.^2+h.^2));
C_08	=	(1/12).*dxdu1.*dxdu2.*R1.^(-1).*R2.^(-1).*(6.*R1.*R2.*(1+vv)+(-1).*((-2)+vv).*h.^2);
C_09	=	(-1/12).*(dxdu2.*dydu1+dxdu1.*dydu2).*R1.^(-1).*R2.^(-1).*((-6).*R1.*R2.*(1+vv)+((-2)+vv).*h.^2);
C_10	=	(1/12).*dydu1.*dydu2.*R1.^(-1).*R2.^(-1).*(6.*R1.*R2.*(1+vv)+(-1).*((-2)+vv).*h.^2);
C_11	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-3).*B1.^(-2).*h.*((-1).*dAdu2.*dxdu1.*dxdu2.*(R2.*(1+(-2).*vv)+2.*R1.*((-1)+vv)).*A1.^2.*B1.*h+(-3).*dAdu1.*dxdu1.^2.*R2.*B1.^3.*h+dxdu1.*R2.*A1.*B1.^2.*(3.*dhdu1.*dxdu1.*B1+(dBdu1.*dxdu1+3.*d2xdu1.*B1).*h)+R2.*A1.^3.*(3.*dxdu2.*((-2).*dhdu2.*dxdu1.*((-1)+vv)+dhdu1.*dxdu2.*vv).*B1+(dBdu1.*dxdu2.^2.*((-3)+vv)+(-1).*d2xdu2.*dxdu1.*((-2)+vv).*B1+dxdu2.*((-2)+vv).*(dBdu2.*dxdu1+(-2).*d2xdu1_du2.*B1)).*h));
C_12	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-3).*B1.^(-2).*h.*((-1).*dAdu2.*(dxdu2.*dydu1+dxdu1.*dydu2).*(R2.*(1+(-2).*vv)+2.*R1.*((-1)+vv)).*A1.^2.*B1.*h+(-6).*dAdu1.*dxdu1.*dydu1.*R2.*B1.^3.*h+R2.*A1.*B1.^2.*(6.*dhdu1.*dxdu1.*dydu1.*B1+(2.*dBdu1.*dxdu1.*dydu1+3.*(d2ydu1.*dxdu1+d2xdu1.*dydu1).*B1).*h)+R2.*A1.^3.*(6.*((-1).*dhdu2.*dxdu1.*dydu2.*((-1)+vv)+dxdu2.*((-1).*dhdu2.*dydu1.*((-1)+vv)+dhdu1.*dydu2.*vv)).*B1+(((-2)+vv).*(dBdu2.*dxdu1.*dydu2+(-1).*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*d2xdu1_du2.*dydu2).*B1)+dxdu2.*(2.*dBdu1.*dydu2.*((-3)+vv)+((-2)+vv).*(dBdu2.*dydu1+(-2).*d2ydu1_du2.*B1))).*h));
C_13	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-3).*B1.^(-2).*h.*((-1).*dAdu2.*dydu1.*dydu2.*(R2.*(1+(-2).*vv)+2.*R1.*((-1)+vv)).*A1.^2.*B1.*h+(-3).*dAdu1.*dydu1.^2.*R2.*B1.^3.*h+dydu1.*R2.*A1.*B1.^2.*(3.*dhdu1.*dydu1.*B1+(dBdu1.*dydu1+3.*d2ydu1.*B1).*h)+R2.*A1.^3.*(3.*dydu2.*((-2).*dhdu2.*dydu1.*((-1)+vv)+dhdu1.*dydu2.*vv).*B1+(dBdu1.*dydu2.^2.*((-3)+vv)+(-1).*d2ydu2.*dydu1.*((-2)+vv).*B1+dydu2.*((-2)+vv).*(dBdu2.*dydu1+(-2).*d2ydu1_du2.*B1)).*h));
C_14	=	(1/12).*R1.^(-3).*((-1).*dAdu1.*dxdu1.*R1.*A1.^(-2).*B1.*(12.*R1.^2+h.^2)+(-2).*dAdu2.*dxdu2.*R1.^2.*R2.^(-1).*((-1)+vv).*B1.^(-1).*(3.*R1.*R2+h.^2)+A1.^(-1).*h.^(-1).*(12.*dhdu1.*dxdu1.*R1.^3.*B1+12.*R1.^3.*(dBdu1.*dxdu1+d2xdu1.*B1).*h+3.*dhdu1.*dxdu1.*R1.*B1.*h.^2+((-2).*dR1du1.*dxdu1.*B1+R1.*(dBdu1.*dxdu1+d2xdu1.*B1)).*h.^3)+2.*((-1)+vv).*A1.*B1.^(-2).*h.^(-1).*((-3).*dhdu2.*dxdu2.*R1.^3.*B1+3.*R1.^3.*(dBdu2.*dxdu2+(-1).*d2xdu2.*B1).*h+(-3).*dhdu2.*dxdu2.*R1.*B1.*h.^2+(dR1du2.*dxdu2.*B1+R1.*(dBdu2.*dxdu2+(-1).*d2xdu2.*B1)).*h.^3));
C_15	=	(1/12).*R1.^(-3).*((-1).*dAdu1.*dydu1.*R1.*A1.^(-2).*B1.*(12.*R1.^2+h.^2)+(-2).*dAdu2.*dydu2.*R1.^2.*R2.^(-1).*((-1)+vv).*B1.^(-1).*(3.*R1.*R2+h.^2)+A1.^(-1).*h.^(-1).*(12.*dhdu1.*dydu1.*R1.^3.*B1+12.*R1.^3.*(dBdu1.*dydu1+d2ydu1.*B1).*h+3.*dhdu1.*dydu1.*R1.*B1.*h.^2+((-2).*dR1du1.*dydu1.*B1+R1.*(dBdu1.*dydu1+d2ydu1.*B1)).*h.^3)+2.*((-1)+vv).*A1.*B1.^(-2).*h.^(-1).*((-3).*dhdu2.*dydu2.*R1.^3.*B1+3.*R1.^3.*(dBdu2.*dydu2+(-1).*d2ydu2.*B1).*h+(-3).*dhdu2.*dydu2.*R1.*B1.*h.^2+(dR1du2.*dydu2.*B1+R1.*(dBdu2.*dydu2+(-1).*d2ydu2.*B1)).*h.^3));
C_16	=	(1/12).*R1.^(-1).*R2.^(-2).*A1.^(-1).*B1.^(-1).*h.^(-1).*(dAdu2.*dxdu1.*B1.*h.*((-6).*R1.*R2.^2.*((-3)+vv)+(R2+(-2).*R1.*((-1)+vv)).*h.^2)+A1.*(6.*R1.*R2.^2.*((-1).*dhdu2.*dxdu1.*((-1)+vv)+2.*dhdu1.*dxdu2.*vv).*B1+6.*R1.*R2.^2.*(dBdu1.*dxdu2.*((-3)+vv)+d2xdu1_du2.*(1+vv).*B1).*h+3.*R2.*((-2).*dhdu2.*dxdu1.*((-1)+vv)+dhdu1.*dxdu2.*vv).*B1.*h.^2+((dR2du2.*dxdu1.*((-2)+vv)+(-1).*dR2du1.*dxdu2.*vv).*B1+R2.*(dBdu1.*dxdu2.*((-3)+2.*vv)+(-1).*d2xdu1_du2.*((-2)+vv).*B1)).*h.^3));
C_17	=	(1/12).*R1.^(-1).*R2.^(-2).*A1.^(-1).*B1.^(-1).*h.^(-1).*(dAdu2.*dydu1.*B1.*h.*((-6).*R1.*R2.^2.*((-3)+vv)+(R2+(-2).*R1.*((-1)+vv)).*h.^2)+A1.*(6.*R1.*R2.^2.*((-1).*dhdu2.*dydu1.*((-1)+vv)+2.*dhdu1.*dydu2.*vv).*B1+6.*R1.*R2.^2.*(dBdu1.*dydu2.*((-3)+vv)+d2ydu1_du2.*(1+vv).*B1).*h+3.*R2.*((-2).*dhdu2.*dydu1.*((-1)+vv)+dhdu1.*dydu2.*vv).*B1.*h.^2+((dR2du2.*dydu1.*((-2)+vv)+(-1).*dR2du1.*dydu2.*vv).*B1+R2.*(dBdu1.*dydu2.*((-3)+2.*vv)+(-1).*d2ydu1_du2.*((-2)+vv).*B1)).*h.^3));
C_18	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-4).*B1.^(-3).*(3.*dAdu1.^2.*dxdu1.*R2.*B1.^4.*h.^2+(-1).*R2.*A1.*B1.^3.*h.*(3.*dAdu1.*dhdu1.*dxdu1.*B1+(d2Adu1.*dxdu1.*B1+dAdu1.*(dBdu1.*dxdu1.*(1+vv)+3.*d2xdu1.*B1)).*h)+A1.^2.*B1.^2.*h.*(3.*dhdu1.*R2.*B1.*(dBdu1.*dxdu1.*vv+d2xdu1.*B1)+(2.*dAdu2.^2.*dxdu1.*R1.*((-1)+vv)+R2.*((-1).*dBdu1.^2.*dxdu1+(-1).*dAdu1.*dAdu2.*dxdu2+(-2).*dAdu2.^2.*dxdu1.*((-1)+vv)+d2xdu1.*dBdu1.*B1+d2Bdu1.*dxdu1.*vv.*B1+d3xdu1.*B1.^2)).*h)+(-1).*A1.^3.*B1.*h.*((-3).*dAdu2.*R2.*(dhdu1.*dxdu2+2.*dhdu2.*dxdu1.*((-1)+vv)).*B1+((-2).*dAdu2.*R1.*((-1)+vv).*(dBdu1.*dxdu2+(-1).*d2xdu1_du2.*B1)+R2.*((-1).*(d2Adu1_du2.*dxdu2+2.*d2Adu2.*dxdu1.*((-1)+vv)).*B1+dAdu2.*(2.*dBdu2.*dxdu1.*((-1)+vv)+dBdu1.*dxdu2.*(1+vv)+d2xdu1_du2.*(1+(-2).*vv).*B1))).*h)+A1.^4.*((-12).*dxdu1.*(R2+R1.*vv).*B1.^4+3.*R2.*B1.*(2.*dhdu2.*((-1)+vv).*(dBdu1.*dxdu2+(-1).*d2xdu1_du2.*B1)+dhdu1.*vv.*((-1).*dBdu2.*dxdu2+d2xdu2.*B1)).*h+R2.*(dBdu2.*(dBdu1.*dxdu2.*(5+(-2).*vv)+d2xdu1_du2.*((-2)+vv).*B1)+B1.*(d2xdu2.*dBdu1.*((-3)+vv)+((-2)+vv).*(d2Bdu1_du2.*dxdu2+(-1).*d3xdu1_2du2.*B1))).*h.^2));
C_19	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-4).*B1.^(-3).*(3.*dAdu1.^2.*dydu1.*R2.*B1.^4.*h.^2+(-1).*R2.*A1.*B1.^3.*h.*(3.*dAdu1.*dhdu1.*dydu1.*B1+(d2Adu1.*dydu1.*B1+dAdu1.*(dBdu1.*dydu1.*(1+vv)+3.*d2ydu1.*B1)).*h)+A1.^2.*B1.^2.*h.*(3.*dhdu1.*R2.*B1.*(dBdu1.*dydu1.*vv+d2ydu1.*B1)+(2.*dAdu2.^2.*dydu1.*R1.*((-1)+vv)+R2.*((-1).*dBdu1.^2.*dydu1+(-1).*dAdu1.*dAdu2.*dydu2+(-2).*dAdu2.^2.*dydu1.*((-1)+vv)+d2ydu1.*dBdu1.*B1+d2Bdu1.*dydu1.*vv.*B1+d3ydu1.*B1.^2)).*h)+(-1).*A1.^3.*B1.*h.*((-3).*dAdu2.*R2.*(dhdu1.*dydu2+2.*dhdu2.*dydu1.*((-1)+vv)).*B1+((-2).*dAdu2.*R1.*((-1)+vv).*(dBdu1.*dydu2+(-1).*d2ydu1_du2.*B1)+R2.*((-1).*(d2Adu1_du2.*dydu2+2.*d2Adu2.*dydu1.*((-1)+vv)).*B1+dAdu2.*(2.*dBdu2.*dydu1.*((-1)+vv)+dBdu1.*dydu2.*(1+vv)+d2ydu1_du2.*(1+(-2).*vv).*B1))).*h)+A1.^4.*((-12).*dydu1.*(R2+R1.*vv).*B1.^4+3.*R2.*B1.*(2.*dhdu2.*((-1)+vv).*(dBdu1.*dydu2+(-1).*d2ydu1_du2.*B1)+dhdu1.*vv.*((-1).*dBdu2.*dydu2+d2ydu2.*B1)).*h+R2.*(dBdu2.*(dBdu1.*dydu2.*(5+(-2).*vv)+d2ydu1_du2.*((-2)+vv).*B1)+B1.*(d2ydu2.*dBdu1.*((-3)+vv)+((-2)+vv).*(d2Bdu1_du2.*dydu2+(-1).*d3ydu1_2du2.*B1))).*h.^2));
C_20	=	(1/12).*R1.^(-4).*R2.^(-1).*A1.^(-2).*B1.^(-2).*h.^(-1).*(dAdu1.*R1.*R2.*B1.^2.*h.*((-12).*dBdu1.*R1.^3.*vv+((-1).*dBdu1.*R1.*vv+dR1du1.*B1).*h.^2)+(-2).*R1.*R2.*((-1)+vv).*A1.^2.*((-3).*dAdu2.*dhdu2.*R1.^3.*B1+3.*R1.^3.*(dAdu2.*dBdu2+(-1).*d2Adu2.*B1).*h+(-3).*dAdu2.*dhdu2.*R1.*B1.*h.^2+(dAdu2.*dR1du2.*B1+R1.*(dAdu2.*dBdu2+(-1).*d2Adu2.*B1)).*h.^3)+(-1).*A1.*B1.*((-12).*dBdu1.*dhdu1.*R1.^4.*R2.*vv.*B1+(-6).*R1.^4.*R2.*((-2).*dBdu1.^2+dAdu2.^2.*((-1)+vv)+2.*d2Bdu1.*vv.*B1).*h+3.*dhdu1.*R1.*R2.*B1.*((-1).*dBdu1.*R1.*vv+dR1du1.*B1).*h.^2+((-2).*dAdu2.^2.*R1.^3.*((-1)+vv)+(-2).*dR1du1.^2.*R2.*B1.^2+R1.*R2.*B1.*(dBdu1.*dR1du1+d2R1du1.*B1)+R1.^2.*R2.*(dBdu1.^2+(-1).*d2Bdu1.*vv.*B1)).*h.^3));
C_21	=	(1/12).*R1.^(-1).*R2.^(-3).*A1.^(-2).*B1.^(-2).*h.^(-1).*((-2).*dBdu1.*dBdu2.*R2.^2.*((-1)+vv).*A1.^2.*h.*(3.*R1.*R2+h.^2)+R2.*A1.*B1.*(6.*dBdu1.*dhdu2.*R1.*R2.^2.*((-1)+vv).*A1+(-6).*R1.*R2.^2.*(dAdu2.*dBdu1.*(1+vv)+(-1).*d2Bdu1_du2.*((-1)+vv).*A1).*h+6.*dBdu1.*dhdu2.*R2.*((-1)+vv).*A1.*h.^2+(2.*dAdu2.*dBdu1.*R1.*((-1)+vv)+dBdu1.*dR2du2.*(3+(-2).*vv).*A1+(-1).*R2.*(dAdu2.*dBdu1.*vv+(-2).*d2Bdu1_du2.*((-1)+vv).*A1)).*h.^3)+(-1).*B1.^2.*((-12).*dAdu2.*dhdu1.*R1.*R2.^3.*A1+12.*R1.*R2.^3.*(dAdu1.*dAdu2+(-1).*d2Adu1_du2.*A1).*h+3.*dhdu1.*R2.*A1.*((-1).*dAdu2.*R2+dR2du2.*vv.*A1).*h.^2+((-2).*dR2du1.*dR2du2.*vv.*A1.^2+R2.^2.*(dAdu1.*dAdu2+(-1).*d2Adu1_du2.*A1)+R2.*A1.*(dAdu2.*dR2du1+d2R2du1_du2.*vv.*A1)).*h.^3));
C_22	=	R1.^(-2).*R2.^(-2).*h.^(-1).*((-1).*dhdu1.*R1.*R2.*(R2+R1.*vv).*B1+(dBdu1.*R1.*R2.^2.*((-1)+vv)+dR1du1.*R2.^2.*B1+R1.^2.*((-1).*dBdu1.*R2.*((-1)+vv)+dR2du1.*vv.*B1)).*h);

C_23	=	(1/12).*dxdu2.*R2.^(-1).*A1.^(-1).*B1.^(-2).*(dxdu2.^2.*A1.^2+(-1).*dxdu1.^2.*((-2)+vv).*B1.^2).*h.^2;
C_24	=	(1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(3.*dxdu2.^2.*dydu2.*A1.^2+(-1).*dxdu1.*(2.*dxdu2.*dydu1+dxdu1.*dydu2).*((-2)+vv).*B1.^2).*h.^2;
C_25	=	(1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(3.*dxdu2.*dydu2.^2.*A1.^2+(-1).*dydu1.*(dxdu2.*dydu1+2.*dxdu1.*dydu2).*((-2)+vv).*B1.^2).*h.^2;
C_26	=	(1/12).*dydu2.*R2.^(-1).*A1.^(-1).*B1.^(-2).*(dydu2.^2.*A1.^2+(-1).*dydu1.^2.*((-2)+vv).*B1.^2).*h.^2;
C_27	=	(1/12).*dxdu1.*dxdu2.*R1.^(-1).*R2.^(-1).*(6.*R1.*R2.*(1+vv)+(-1).*((-2)+vv).*h.^2);
C_28	=	(-1/12).*(dxdu2.*dydu1+dxdu1.*dydu2).*R1.^(-1).*R2.^(-1).*((-6).*R1.*R2.*(1+vv)+((-2)+vv).*h.^2);
C_29	=	(1/12).*dydu1.*dydu2.*R1.^(-1).*R2.^(-1).*(6.*R1.*R2.*(1+vv)+(-1).*((-2)+vv).*h.^2);
C_30	=	(1/12).*R2.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dxdu1.^2.*((-1)+vv).*B1.^2.*(3.*R2.^2+h.^2)+dxdu2.^2.*A1.^2.*(12.*R2.^2+h.^2));
C_31	=	(1/6).*R2.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dxdu1.*dydu1.*((-1)+vv).*B1.^2.*(3.*R2.^2+h.^2)+dxdu2.*dydu2.*A1.^2.*(12.*R2.^2+h.^2));
C_32	=	(1/12).*R2.^(-2).*A1.^(-1).*B1.^(-1).*((-2).*dydu1.^2.*((-1)+vv).*B1.^2.*(3.*R2.^2+h.^2)+dydu2.^2.*A1.^2.*(12.*R2.^2+h.^2));
C_33	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-2).*B1.^(-3).*h.*(dAdu2.*dxdu2.^2.*R1.*A1.^2.*B1.*h+dxdu1.*R1.*(dAdu2.*dxdu1.*((-3)+vv)+dAdu1.*dxdu2.*((-2)+vv)).*B1.^3.*h+3.*dxdu2.*R1.*A1.^3.*(dhdu2.*dxdu2.*B1+((-1).*dBdu2.*dxdu2+d2xdu2.*B1).*h)+A1.*B1.^2.*(3.*dxdu1.*R1.*((-2).*dhdu1.*dxdu2.*((-1)+vv)+dhdu2.*dxdu1.*vv).*B1+((-2).*dBdu1.*dxdu1.*dxdu2.*R2.*((-1)+vv)+R1.*((-2).*d2xdu1_du2.*dxdu1.*((-2)+vv).*B1+dxdu2.*(dBdu1.*dxdu1.*((-1)+2.*vv)+(-1).*d2xdu1.*((-2)+vv).*B1))).*h));
C_34	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-2).*B1.^(-3).*h.*(2.*dAdu2.*dxdu2.*dydu2.*R1.*A1.^2.*B1.*h+R1.*(dydu1.*(2.*dAdu2.*dxdu1.*((-3)+vv)+dAdu1.*dxdu2.*((-2)+vv))+dAdu1.*dxdu1.*dydu2.*((-2)+vv)).*B1.^3.*h+3.*R1.*A1.^3.*(2.*dhdu2.*dxdu2.*dydu2.*B1+((-2).*dBdu2.*dxdu2.*dydu2+(d2ydu2.*dxdu2+d2xdu2.*dydu2).*B1).*h)+(-1).*A1.*B1.^2.*((-6).*R1.*((-1).*dhdu1.*dxdu1.*dydu2.*((-1)+vv)+dydu1.*((-1).*dhdu1.*dxdu2.*((-1)+vv)+dhdu2.*dxdu1.*vv)).*B1+(2.*dBdu1.*(dxdu2.*dydu1+dxdu1.*dydu2).*R2.*((-1)+vv)+R1.*(2.*(d2ydu1_du2.*dxdu1+d2xdu1_du2.*dydu1).*((-2)+vv).*B1+dydu2.*(dBdu1.*dxdu1.*(1+(-2).*vv)+d2xdu1.*((-2)+vv).*B1)+dxdu2.*(dBdu1.*dydu1.*(1+(-2).*vv)+d2ydu1.*((-2)+vv).*B1))).*h));
C_35	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-2).*B1.^(-3).*h.*(dAdu2.*dydu2.^2.*R1.*A1.^2.*B1.*h+dydu1.*R1.*(dAdu2.*dydu1.*((-3)+vv)+dAdu1.*dydu2.*((-2)+vv)).*B1.^3.*h+3.*dydu2.*R1.*A1.^3.*(dhdu2.*dydu2.*B1+((-1).*dBdu2.*dydu2+d2ydu2.*B1).*h)+A1.*B1.^2.*(3.*dydu1.*R1.*((-2).*dhdu1.*dydu2.*((-1)+vv)+dhdu2.*dydu1.*vv).*B1+((-2).*dBdu1.*dydu1.*dydu2.*R2.*((-1)+vv)+R1.*((-2).*d2ydu1_du2.*dydu1.*((-2)+vv).*B1+dydu2.*(dBdu1.*dydu1.*((-1)+2.*vv)+(-1).*d2ydu1.*((-2)+vv).*B1))).*h));
C_36	=	(1/12).*R1.^(-2).*R2.^(-1).*A1.^(-1).*B1.^(-1).*h.^(-1).*(dAdu2.*dxdu1.*R1.*B1.*h.*(6.*R1.*R2.*((-3)+vv)+((-3)+2.*vv).*h.^2)+A1.*(6.*R1.^2.*R2.*((-1).*dhdu1.*dxdu2.*((-1)+vv)+2.*dhdu2.*dxdu1.*vv).*B1+6.*R1.^2.*R2.*((-1).*dBdu1.*dxdu2.*((-3)+vv)+d2xdu1_du2.*(1+vv).*B1).*h+3.*R1.*((-2).*dhdu1.*dxdu2.*((-1)+vv)+dhdu2.*dxdu1.*vv).*B1.*h.^2+((-2).*dBdu1.*dxdu2.*R2.*((-1)+vv)+(dR1du1.*dxdu2.*((-2)+vv)+(-1).*dR1du2.*dxdu1.*vv).*B1+R1.*(dBdu1.*dxdu2+(-1).*d2xdu1_du2.*((-2)+vv).*B1)).*h.^3));
C_37	=	(1/12).*R1.^(-2).*R2.^(-1).*A1.^(-1).*B1.^(-1).*h.^(-1).*(dAdu2.*dydu1.*R1.*B1.*h.*(6.*R1.*R2.*((-3)+vv)+((-3)+2.*vv).*h.^2)+A1.*(6.*R1.^2.*R2.*((-1).*dhdu1.*dydu2.*((-1)+vv)+2.*dhdu2.*dydu1.*vv).*B1+6.*R1.^2.*R2.*((-1).*dBdu1.*dydu2.*((-3)+vv)+d2ydu1_du2.*(1+vv).*B1).*h+3.*R1.*((-2).*dhdu1.*dydu2.*((-1)+vv)+dhdu2.*dydu1.*vv).*B1.*h.^2+((-2).*dBdu1.*dydu2.*R2.*((-1)+vv)+(dR1du1.*dydu2.*((-2)+vv)+(-1).*dR1du2.*dydu1.*vv).*B1+R1.*(dBdu1.*dydu2+(-1).*d2ydu1_du2.*((-2)+vv).*B1)).*h.^3));
C_38	=	(1/12).*R1.^(-1).*R2.^(-3).*A1.^(-2).*B1.^(-2).*h.^(-1).*(2.*dAdu1.*dxdu1.*R1.*R2.*((-1)+vv).*B1.^3.*h.*(3.*R2.^2+h.^2)+dAdu2.*dxdu2.*R1.*R2.*A1.^2.*B1.*h.*(12.*R2.^2+h.^2)+(-2).*((-1)+vv).*A1.*B1.^2.*(3.*dhdu1.*dxdu1.*R1.*R2.^3.*B1+3.*R1.*R2.^3.*(dBdu1.*dxdu1+d2xdu1.*B1).*h+3.*dhdu1.*dxdu1.*R1.*R2.*B1.*h.^2+(dBdu1.*dxdu1.*R2.^2+(-1).*dR2du1.*dxdu1.*R1.*B1+d2xdu1.*R1.*R2.*B1).*h.^3)+R1.*A1.^3.*(12.*dhdu2.*dxdu2.*R2.^3.*B1+(-12).*R2.^3.*(dBdu2.*dxdu2+(-1).*d2xdu2.*B1).*h+3.*dhdu2.*dxdu2.*R2.*B1.*h.^2+((-2).*dR2du2.*dxdu2.*B1+R2.*((-1).*dBdu2.*dxdu2+d2xdu2.*B1)).*h.^3));
C_39	=	(1/12).*R1.^(-1).*R2.^(-3).*A1.^(-2).*B1.^(-2).*h.^(-1).*(2.*dAdu1.*dydu1.*R1.*R2.*((-1)+vv).*B1.^3.*h.*(3.*R2.^2+h.^2)+dAdu2.*dydu2.*R1.*R2.*A1.^2.*B1.*h.*(12.*R2.^2+h.^2)+(-2).*((-1)+vv).*A1.*B1.^2.*(3.*dhdu1.*dydu1.*R1.*R2.^3.*B1+3.*R1.*R2.^3.*(dBdu1.*dydu1+d2ydu1.*B1).*h+3.*dhdu1.*dydu1.*R1.*R2.*B1.*h.^2+(dBdu1.*dydu1.*R2.^2+(-1).*dR2du1.*dydu1.*R1.*B1+d2ydu1.*R1.*R2.*B1).*h.^3)+R1.*A1.^3.*(12.*dhdu2.*dydu2.*R2.^3.*B1+(-12).*R2.^3.*(dBdu2.*dydu2+(-1).*d2ydu2.*B1).*h+3.*dhdu2.*dydu2.*R2.*B1.*h.^2+((-2).*dR2du2.*dydu2.*B1+R2.*((-1).*dBdu2.*dydu2+d2ydu2.*B1)).*h.^3));
C_40	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-3).*B1.^(-4).*(dAdu1.*dAdu2.*dxdu1.*R1.*(5+(-2).*vv).*B1.^4.*h.^2+R1.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dxdu2.*vv.*B1+(d2Adu2.*dxdu2.*vv.*B1+(-1).*dAdu2.*(dBdu2.*dxdu2.*(1+vv)+(-1).*d2xdu2.*B1)).*h)+(-1).*A1.^2.*B1.^2.*h.*((-3).*R1.*B1.*(2.*dBdu1.*dhdu1.*dxdu2.*((-1)+vv)+(-2).*d2xdu1_du2.*dhdu1.*((-1)+vv).*B1+dhdu2.*(dBdu1.*dxdu1+d2xdu1.*vv.*B1))+((-2).*dBdu1.*R2.*((-1)+vv).*(dBdu1.*dxdu2+(-1).*d2xdu1_du2.*B1)+R1.*(dBdu1.*dBdu2.*dxdu1+dAdu2.^2.*dxdu2+d2xdu1_du2.*dBdu1.*B1+(-1).*d2Bdu1_du2.*dxdu1.*B1+(-2).*d2xdu1_du2.*dBdu1.*vv.*B1+(-2).*d3xd2u1_du2.*B1.^2+d3xd2u1_du2.*vv.*B1.^2+2.*dxdu2.*((-1)+vv).*(dBdu1.^2+(-1).*d2Bdu1.*B1))).*h)+A1.*B1.^3.*h.*((-3).*dxdu1.*R1.*((-2).*dAdu2.*dhdu1.*((-1)+vv)+dAdu1.*dhdu2.*vv).*B1+(2.*dAdu2.*dBdu1.*dxdu1.*R2.*((-1)+vv)+R1.*((-2).*dAdu1.*dBdu1.*dxdu2.*((-1)+vv)+(d2xdu1_du2.*dAdu1+d2Adu1_du2.*dxdu1).*((-2)+vv).*B1+(-1).*dAdu2.*(dBdu1.*dxdu1.*(1+vv)+(-1).*d2xdu1.*((-3)+vv).*B1))).*h)+(-1).*A1.^4.*(12.*dxdu2.*(R1+R2.*vv).*B1.^4+(-3).*dBdu2.^2.*dxdu2.*R1.*h.^2+(-1).*R1.*B1.^2.*h.*(3.*d2xdu2.*dhdu2+d3xdu2.*h)+R1.*B1.*h.*(d2Bdu2.*dxdu2.*h+3.*dBdu2.*(dhdu2.*dxdu2+d2xdu2.*h))));
C_41	=	(1/12).*R1.^(-1).*R2.^(-1).*A1.^(-3).*B1.^(-4).*(dAdu1.*dAdu2.*dydu1.*R1.*(5+(-2).*vv).*B1.^4.*h.^2+R1.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dydu2.*vv.*B1+(d2Adu2.*dydu2.*vv.*B1+(-1).*dAdu2.*(dBdu2.*dydu2.*(1+vv)+(-1).*d2ydu2.*B1)).*h)+(-1).*A1.^2.*B1.^2.*h.*((-3).*R1.*B1.*(2.*dBdu1.*dhdu1.*dydu2.*((-1)+vv)+(-2).*d2ydu1_du2.*dhdu1.*((-1)+vv).*B1+dhdu2.*(dBdu1.*dydu1+d2ydu1.*vv.*B1))+((-2).*dBdu1.*R2.*((-1)+vv).*(dBdu1.*dydu2+(-1).*d2ydu1_du2.*B1)+R1.*(dBdu1.*dBdu2.*dydu1+dAdu2.^2.*dydu2+d2ydu1_du2.*dBdu1.*B1+(-1).*d2Bdu1_du2.*dydu1.*B1+(-2).*d2ydu1_du2.*dBdu1.*vv.*B1+(-2).*d3yd2u1_du2.*B1.^2+d3yd2u1_du2.*vv.*B1.^2+2.*dydu2.*((-1)+vv).*(dBdu1.^2+(-1).*d2Bdu1.*B1))).*h)+A1.*B1.^3.*h.*((-3).*dydu1.*R1.*((-2).*dAdu2.*dhdu1.*((-1)+vv)+dAdu1.*dhdu2.*vv).*B1+(2.*dAdu2.*dBdu1.*dydu1.*R2.*((-1)+vv)+R1.*((-2).*dAdu1.*dBdu1.*dydu2.*((-1)+vv)+(d2ydu1_du2.*dAdu1+d2Adu1_du2.*dydu1).*((-2)+vv).*B1+(-1).*dAdu2.*(dBdu1.*dydu1.*(1+vv)+(-1).*d2ydu1.*((-3)+vv).*B1))).*h)+(-1).*A1.^4.*(12.*dydu2.*(R1+R2.*vv).*B1.^4+(-3).*dBdu2.^2.*dydu2.*R1.*h.^2+(-1).*R1.*B1.^2.*h.*(3.*d2ydu2.*dhdu2+d3ydu2.*h)+R1.*B1.*h.*(d2Bdu2.*dydu2.*h+3.*dBdu2.*(dhdu2.*dydu2+d2ydu2.*h))));
C_42	=	(-1/12).*R1.^(-3).*R2.^(-1).*A1.^(-2).*B1.^(-2).*h.^(-1).*(dBdu1.*dBdu2.*R1.^2.*A1.^2.*h.*(12.*R1.*R2+h.^2)+R1.*A1.*B1.*((-12).*dBdu1.*dhdu2.*R1.^2.*R2.*A1+6.*R1.^2.*R2.*(dAdu2.*dBdu1.*(1+vv)+(-2).*d2Bdu1_du2.*A1).*h+(-3).*dBdu1.*dhdu2.*R1.*A1.*h.^2+(R1.*(dAdu2.*dBdu1.*vv+(-1).*d2Bdu1_du2.*A1)+dBdu1.*((-2).*dAdu2.*R2.*((-1)+vv)+dR1du2.*A1)).*h.^3)+B1.^2.*((-6).*dAdu2.*dhdu1.*R1.^3.*R2.*((-1)+vv).*A1+6.*R1.^3.*R2.*((-1)+vv).*(dAdu1.*dAdu2+(-1).*d2Adu1_du2.*A1).*h+3.*R1.*A1.*((-2).*dAdu2.*dhdu1.*R1.*((-1)+vv)+dhdu2.*dR1du1.*vv.*A1).*h.^2+((-2).*dR1du1.*dR1du2.*vv.*A1.^2+2.*R1.^2.*((-1)+vv).*(dAdu1.*dAdu2+(-1).*d2Adu1_du2.*A1)+R1.*A1.*(dAdu2.*dR1du1.*((-3)+2.*vv)+d2R1du1_du2.*vv.*A1)).*h.^3));
C_43	=	(1/12).*R1.^(-1).*R2.^(-4).*A1.^(-2).*B1.^(-2).*h.^(-1).*(R1.*A1.^3.*h.^2.*((-3).*dhdu2.*dR2du2.*R2.*B1+(2.*dR2du2.^2.*B1+R2.*(dBdu2.*dR2du2+(-1).*d2R2du2.*B1)).*h)+(-2).*dAdu1.*dBdu1.*R1.*R2.^2.*((-1)+vv).*B1.^2.*h.*(3.*R2.^2+h.^2)+R1.*R2.*A1.^2.*(12.*dAdu2.*dhdu2.*R2.^3.*vv.*B1+(-12).*R2.^3.*vv.*(dAdu2.*dBdu2+(-1).*d2Adu2.*B1).*h+3.*dAdu2.*dhdu2.*R2.*vv.*B1.*h.^2+(-1).*(dAdu2.*dR2du2.*B1+R2.*vv.*(dAdu2.*dBdu2+(-1).*d2Adu2.*B1)).*h.^3)+R2.*A1.*B1.*(6.*dBdu1.*dhdu1.*R1.*R2.^3.*((-1)+vv).*B1+(-6).*R1.*R2.^3.*(2.*dAdu2.^2+(-1).*((-1)+vv).*(dBdu1.^2+d2Bdu1.*B1)).*h+6.*dBdu1.*dhdu1.*R1.*R2.*((-1)+vv).*B1.*h.^2+(-1).*((-2).*dBdu1.^2.*R2.^2.*((-1)+vv)+R1.*(2.*dBdu1.*dR2du1.*((-1)+vv).*B1+R2.*(dAdu2.^2+(-2).*d2Bdu1.*((-1)+vv).*B1))).*h.^3));
C_44	=	R1.^(-2).*R2.^(-2).*h.^(-1).*((-1).*dhdu2.*R1.*R2.*(R1+R2.*vv).*A1+((-1).*dAdu2.*R1.*R2.^2.*((-1)+vv)+dR1du2.*R2.^2.*vv.*A1+R1.^2.*(dAdu2.*R2.*((-1)+vv)+dR2du2.*A1)).*h);

C_45	=	(-1/12).*A1.^(-3).*B1.^(-3).*(dxdu2.^2.*A1.^2+dxdu1.^2.*B1.^2).^2.*h.^2;
C_46	=	(-1/3).*A1.^(-3).*B1.^(-3).*(dxdu2.^2.*A1.^2+dxdu1.^2.*B1.^2).*(dxdu2.*dydu2.*A1.^2+dxdu1.*dydu1.*B1.^2).*h.^2;
C_47	=	(-1/6).*A1.^(-3).*B1.^(-3).*(3.*dxdu2.^2.*dydu2.^2.*A1.^4+(dxdu2.^2.*dydu1.^2+4.*dxdu1.*dxdu2.*dydu1.*dydu2+dxdu1.^2.*dydu2.^2).*A1.^2.*B1.^2+3.*dxdu1.^2.*dydu1.^2.*B1.^4).*h.^2;
C_48	=	(-1/3).*A1.^(-3).*B1.^(-3).*(dxdu2.*dydu2.*A1.^2+dxdu1.*dydu1.*B1.^2).*(dydu2.^2.*A1.^2+dydu1.^2.*B1.^2).*h.^2;
C_49	=	(-1/12).*A1.^(-3).*B1.^(-3).*(dydu2.^2.*A1.^2+dydu1.^2.*B1.^2).^2.*h.^2;
C_50	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*(dxdu1.*dxdu2.^2.*((-2)+vv).*A1.^2+(-1).*dxdu1.^3.*B1.^2).*h.^2;
C_51	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*(dxdu2.*(dxdu2.*dydu1+2.*dxdu1.*dydu2).*((-2)+vv).*A1.^2+(-3).*dxdu1.^2.*dydu1.*B1.^2).*h.^2;
C_52	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*(dydu2.*(2.*dxdu2.*dydu1+dxdu1.*dydu2).*((-2)+vv).*A1.^2+(-3).*dxdu1.*dydu1.^2.*B1.^2).*h.^2;
C_53	=	(1/12).*R1.^(-1).*A1.^(-2).*B1.^(-1).*(dydu1.*dydu2.^2.*((-2)+vv).*A1.^2+(-1).*dydu1.^3.*B1.^2).*h.^2;
C_54	=	(-1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(dxdu2.^3.*A1.^2+(-1).*dxdu1.^2.*dxdu2.*((-2)+vv).*B1.^2).*h.^2;
C_55	=	(-1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(3.*dxdu2.^2.*dydu2.*A1.^2+(-1).*dxdu1.*(2.*dxdu2.*dydu1+dxdu1.*dydu2).*((-2)+vv).*B1.^2).*h.^2;
C_56	=	(-1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(3.*dxdu2.*dydu2.^2.*A1.^2+(-1).*dydu1.*(dxdu2.*dydu1+2.*dxdu1.*dydu2).*((-2)+vv).*B1.^2).*h.^2;
C_57	=	(-1/12).*R2.^(-1).*A1.^(-1).*B1.^(-2).*(dydu2.^3.*A1.^2+(-1).*dydu1.^2.*dydu2.*((-2)+vv).*B1.^2).*h.^2;
C_58	=	(-1/6).*A1.^(-4).*B1.^(-4).*h.*(dAdu2.*dxdu2.^3.*A1.^4.*B1.*h+(-1).*dxdu1.*dxdu2.*(dAdu2.*dxdu1+dAdu1.*dxdu2).*A1.^2.*B1.^3.*h+(-3).*dAdu1.*dxdu1.^3.*B1.^5.*h+dxdu1.^2.*A1.*B1.^4.*(3.*dhdu1.*dxdu1.*B1+(dBdu1.*dxdu1+3.*d2xdu1.*B1).*h)+3.*dxdu2.^2.*A1.^5.*(dhdu2.*dxdu2.*B1+((-1).*dBdu2.*dxdu2+d2xdu2.*B1).*h)+A1.^3.*B1.^2.*(3.*dxdu1.*dxdu2.*(dhdu2.*dxdu1+dhdu1.*dxdu2).*B1+(d2xdu2.*dxdu1.^2.*B1+dxdu2.^2.*((-1).*dBdu1.*dxdu1+d2xdu1.*B1)+dxdu1.*dxdu2.*((-1).*dBdu2.*dxdu1+4.*d2xdu1_du2.*B1)).*h));
C_59	=	(-1/6).*A1.^(-4).*B1.^(-4).*h.*(3.*dAdu2.*dxdu2.^2.*dydu2.*A1.^4.*B1.*h+(-1).*(dAdu1.*dxdu2.^2.*dydu1+dAdu2.*dxdu1.^2.*dydu2+2.*dxdu1.*dxdu2.*(dAdu2.*dydu1+dAdu1.*dydu2)).*A1.^2.*B1.^3.*h+(-9).*dAdu1.*dxdu1.^2.*dydu1.*B1.^5.*h+3.*dxdu1.*A1.*B1.^4.*(3.*dhdu1.*dxdu1.*dydu1.*B1+(dBdu1.*dxdu1.*dydu1+(d2ydu1.*dxdu1+2.*d2xdu1.*dydu1).*B1).*h)+3.*dxdu2.*A1.^5.*(3.*dhdu2.*dxdu2.*dydu2.*B1+((-3).*dBdu2.*dxdu2.*dydu2+(d2ydu2.*dxdu2+2.*d2xdu2.*dydu2).*B1).*h)+A1.^3.*B1.^2.*(3.*(dhdu1.*dxdu2.^2.*dydu1+dhdu2.*dxdu1.^2.*dydu2+2.*dxdu1.*dxdu2.*(dhdu2.*dydu1+dhdu1.*dydu2)).*B1+(dxdu2.^2.*((-1).*dBdu1.*dydu1+d2ydu1.*B1)+dxdu1.*((-1).*dBdu2.*dxdu1.*dydu2+(d2ydu2.*dxdu1+2.*d2xdu2.*dydu1+4.*d2xdu1_du2.*dydu2).*B1)+dxdu2.*((-2).*dBdu2.*dxdu1.*dydu1+4.*(d2ydu1_du2.*dxdu1+d2xdu1_du2.*dydu1).*B1+dydu2.*((-2).*dBdu1.*dxdu1+2.*d2xdu1.*B1))).*h));
C_60	=	(-1/6).*A1.^(-4).*B1.^(-4).*h.*(3.*dAdu2.*dxdu2.*dydu2.^2.*A1.^4.*B1.*h+(-1).*(dAdu2.*dxdu2.*dydu1.^2+2.*(dAdu2.*dxdu1+dAdu1.*dxdu2).*dydu1.*dydu2+dAdu1.*dxdu1.*dydu2.^2).*A1.^2.*B1.^3.*h+(-9).*dAdu1.*dxdu1.*dydu1.^2.*B1.^5.*h+3.*dydu1.*A1.*B1.^4.*(3.*dhdu1.*dxdu1.*dydu1.*B1+(dBdu1.*dxdu1.*dydu1+(2.*d2ydu1.*dxdu1+d2xdu1.*dydu1).*B1).*h)+3.*dydu2.*A1.^5.*(3.*dhdu2.*dxdu2.*dydu2.*B1+((-3).*dBdu2.*dxdu2.*dydu2+(2.*d2ydu2.*dxdu2+d2xdu2.*dydu2).*B1).*h)+A1.^3.*B1.^2.*(3.*(dhdu2.*dxdu2.*dydu1.^2+2.*(dhdu2.*dxdu1+dhdu1.*dxdu2).*dydu1.*dydu2+dhdu1.*dxdu1.*dydu2.^2).*B1+(dydu2.^2.*((-1).*dBdu1.*dxdu1+d2xdu1.*B1)+dydu1.*((-1).*dBdu2.*dxdu2.*dydu1+(2.*d2ydu2.*dxdu1+4.*d2ydu1_du2.*dxdu2+d2xdu2.*dydu1).*B1)+dydu2.*((-2).*dBdu2.*dxdu1.*dydu1+4.*(d2ydu1_du2.*dxdu1+d2xdu1_du2.*dydu1).*B1+dxdu2.*((-2).*dBdu1.*dydu1+2.*d2ydu1.*B1))).*h));
C_61	=	(-1/6).*A1.^(-4).*B1.^(-4).*h.*(dAdu2.*dydu2.^3.*A1.^4.*B1.*h+(-1).*dydu1.*dydu2.*(dAdu2.*dydu1+dAdu1.*dydu2).*A1.^2.*B1.^3.*h+(-3).*dAdu1.*dydu1.^3.*B1.^5.*h+dydu1.^2.*A1.*B1.^4.*(3.*dhdu1.*dydu1.*B1+(dBdu1.*dydu1+3.*d2ydu1.*B1).*h)+3.*dydu2.^2.*A1.^5.*(dhdu2.*dydu2.*B1+((-1).*dBdu2.*dydu2+d2ydu2.*B1).*h)+A1.^3.*B1.^2.*(3.*dydu1.*dydu2.*(dhdu2.*dydu1+dhdu1.*dydu2).*B1+(d2ydu2.*dydu1.^2.*B1+dydu2.^2.*((-1).*dBdu1.*dydu1+d2ydu1.*B1)+dydu1.*dydu2.*((-1).*dBdu2.*dydu1+4.*d2ydu1_du2.*B1)).*h));
C_62	=	(-1/12).*R1.^(-2).*A1.^(-3).*B1.^(-2).*h.*((-1).*dAdu2.*dxdu1.*dxdu2.*R1.*A1.^2.*B1.*h+(-3).*dAdu1.*dxdu1.^2.*R1.*B1.^3.*h+dxdu1.*A1.*B1.^2.*(6.*dhdu1.*dxdu1.*R1.*B1+((-3).*dR1du1.*dxdu1.*B1+R1.*(2.*dBdu1.*dxdu1+3.*d2xdu1.*B1)).*h)+A1.^3.*(6.*dxdu2.*R1.*(dhdu2.*dxdu1+(-1).*dhdu1.*dxdu2.*((-1)+vv)).*B1+(dxdu2.*((-2).*dR1du2.*dxdu1+dR1du1.*dxdu2.*((-2)+vv)).*B1+R1.*(dBdu1.*dxdu2.^2+(-1).*d2xdu2.*dxdu1.*((-2)+vv).*B1+dxdu2.*((-2)+vv).*(dBdu2.*dxdu1+(-2).*d2xdu1_du2.*B1))).*h));
C_63	=	(1/12).*R1.^(-2).*A1.^(-3).*B1.^(-2).*h.*(dAdu2.*(dxdu2.*dydu1+dxdu1.*dydu2).*R1.*A1.^2.*B1.*h+6.*dAdu1.*dxdu1.*dydu1.*R1.*B1.^3.*h+(-1).*A1.*B1.^2.*(12.*dhdu1.*dxdu1.*dydu1.*R1.*B1+((-6).*dR1du1.*dxdu1.*dydu1.*B1+R1.*(4.*dBdu1.*dxdu1.*dydu1+3.*(d2ydu1.*dxdu1+d2xdu1.*dydu1).*B1)).*h)+A1.^3.*(6.*R1.*((-1).*dhdu2.*dxdu1.*dydu2+dxdu2.*((-1).*dhdu2.*dydu1+2.*dhdu1.*dydu2.*((-1)+vv))).*B1+(2.*(dR1du2.*dxdu1.*dydu2+dxdu2.*(dR1du2.*dydu1+(-1).*dR1du1.*dydu2.*((-2)+vv))).*B1+R1.*((-1).*((-2)+vv).*(dBdu2.*dxdu1.*dydu2+(-1).*(d2ydu2.*dxdu1+d2xdu2.*dydu1+2.*d2xdu1_du2.*dydu2).*B1)+(-1).*dxdu2.*(2.*dBdu1.*dydu2+((-2)+vv).*(dBdu2.*dydu1+(-2).*d2ydu1_du2.*B1)))).*h));
C_64	=	(-1/12).*R1.^(-2).*A1.^(-3).*B1.^(-2).*h.*((-1).*dAdu2.*dydu1.*dydu2.*R1.*A1.^2.*B1.*h+(-3).*dAdu1.*dydu1.^2.*R1.*B1.^3.*h+dydu1.*A1.*B1.^2.*(6.*dhdu1.*dydu1.*R1.*B1+((-3).*dR1du1.*dydu1.*B1+R1.*(2.*dBdu1.*dydu1+3.*d2ydu1.*B1)).*h)+A1.^3.*(6.*dydu2.*R1.*(dhdu2.*dydu1+(-1).*dhdu1.*dydu2.*((-1)+vv)).*B1+(dydu2.*((-2).*dR1du2.*dydu1+dR1du1.*dydu2.*((-2)+vv)).*B1+R1.*(dBdu1.*dydu2.^2+(-1).*d2ydu2.*dydu1.*((-2)+vv).*B1+dydu2.*((-2)+vv).*(dBdu2.*dydu1+(-2).*d2ydu1_du2.*B1))).*h));
C_65	=	(-1/12).*R2.^(-2).*A1.^(-2).*B1.^(-3).*h.*(2.*dAdu2.*dxdu2.^2.*R2.*A1.^2.*B1.*h+dxdu1.*R2.*(dAdu2.*dxdu1+dAdu1.*dxdu2.*((-2)+vv)).*B1.^3.*h+(-3).*dxdu2.*A1.^3.*((-2).*dhdu2.*dxdu2.*R2.*B1+(dR2du2.*dxdu2.*B1+R2.*(dBdu2.*dxdu2+(-1).*d2xdu2.*B1)).*h)+(-1).*A1.*B1.^2.*(6.*dxdu1.*R2.*((-1).*dhdu1.*dxdu2+dhdu2.*dxdu1.*((-1)+vv)).*B1+(dxdu1.*(2.*dR2du1.*dxdu2+(-1).*dR2du2.*dxdu1.*((-2)+vv)).*B1+R2.*(2.*d2xdu1_du2.*dxdu1.*((-2)+vv).*B1+dxdu2.*(dBdu1.*dxdu1+d2xdu1.*((-2)+vv).*B1))).*h));
C_66	=	(-1/12).*R2.^(-2).*A1.^(-2).*B1.^(-3).*h.*(4.*dAdu2.*dxdu2.*dydu2.*R2.*A1.^2.*B1.*h+R2.*(dydu1.*(2.*dAdu2.*dxdu1+dAdu1.*dxdu2.*((-2)+vv))+dAdu1.*dxdu1.*dydu2.*((-2)+vv)).*B1.^3.*h+3.*A1.^3.*(4.*dhdu2.*dxdu2.*dydu2.*R2.*B1+((-2).*dR2du2.*dxdu2.*dydu2.*B1+R2.*((-2).*dBdu2.*dxdu2.*dydu2+(d2ydu2.*dxdu2+d2xdu2.*dydu2).*B1)).*h)+(-1).*A1.*B1.^2.*((-6).*R2.*(dhdu1.*dxdu1.*dydu2+dydu1.*(dhdu1.*dxdu2+(-2).*dhdu2.*dxdu1.*((-1)+vv))).*B1+(2.*(dR2du1.*dxdu1.*dydu2+dydu1.*(dR2du1.*dxdu2+(-1).*dR2du2.*dxdu1.*((-2)+vv))).*B1+R2.*(2.*(d2ydu1_du2.*dxdu1+d2xdu1_du2.*dydu1).*((-2)+vv).*B1+dydu2.*(dBdu1.*dxdu1+d2xdu1.*((-2)+vv).*B1)+dxdu2.*(dBdu1.*dydu1+d2ydu1.*((-2)+vv).*B1))).*h));
C_67	=	(-1/12).*R2.^(-2).*A1.^(-2).*B1.^(-3).*h.*(2.*dAdu2.*dydu2.^2.*R2.*A1.^2.*B1.*h+dydu1.*R2.*(dAdu2.*dydu1+dAdu1.*dydu2.*((-2)+vv)).*B1.^3.*h+(-3).*dydu2.*A1.^3.*((-2).*dhdu2.*dydu2.*R2.*B1+(dR2du2.*dydu2.*B1+R2.*(dBdu2.*dydu2+(-1).*d2ydu2.*B1)).*h)+(-1).*A1.*B1.^2.*(6.*dydu1.*R2.*((-1).*dhdu1.*dydu2+dhdu2.*dydu1.*((-1)+vv)).*B1+(dydu1.*(2.*dR2du1.*dydu2+(-1).*dR2du2.*dydu1.*((-2)+vv)).*B1+R2.*(2.*d2ydu1_du2.*dydu1.*((-2)+vv).*B1+dydu2.*(dBdu1.*dydu1+d2ydu1.*((-2)+vv).*B1))).*h));
C_68	=	(-1/12).*A1.^(-5).*B1.^(-5).*(15.*dAdu1.^2.*dxdu1.^2.*B1.^6.*h.^2+(-1).*dxdu1.*A1.*B1.^5.*h.*(21.*dAdu1.*dhdu1.*dxdu1.*B1+(4.*d2Adu1.*dxdu1.*B1+dAdu1.*(dBdu1.*dxdu1.*(7+vv)+18.*d2xdu1.*B1)).*h)+dxdu2.*A1.^5.*B1.*h.*(3.*dAdu2.*dhdu2.*dxdu2.*(2+vv).*B1+(d2Adu2.*dxdu2.*(1+vv).*B1+(-1).*dAdu2.*(dBdu2.*dxdu2.*(7+vv)+(-6).*d2xdu2.*B1)).*h)+(-1).*A1.^3.*B1.^3.*h.*(3.*(3.*dAdu2.*dhdu2.*dxdu1.^2+dAdu1.*dhdu1.*dxdu2.^2.*vv+2.*dxdu1.*dxdu2.*(dAdu1.*dhdu2+(-1).*dAdu2.*dhdu1.*vv)).*B1+(dAdu1.*dBdu1.*dxdu2.^2.*((-3)+vv)+dxdu1.*(dAdu2.*dBdu2.*dxdu1.*((-3)+vv)+(2.*d2xdu2.*dAdu1+4.*d2xdu1_du2.*dAdu2+(-1).*d2Adu2.*dxdu1.*((-3)+vv)).*B1)+dxdu2.*((-2).*dAdu1.*dBdu2.*dxdu1+4.*d2xdu1_du2.*dAdu1.*B1+2.*dAdu2.*(3.*dBdu1.*dxdu1+d2xdu1.*B1))).*h)+A1.^2.*B1.^4.*(6.*dhdu1.^2.*dxdu1.^2.*B1.^2+3.*dxdu1.*B1.*(dBdu1.*dhdu1.*dxdu1.*(2+vv)+(6.*d2xdu1.*dhdu1+d2hdu1.*dxdu1).*B1).*h+(4.*dAdu2.^2.*dxdu1.^2+(-1).*dBdu1.^2.*dxdu1.^2+2.*dAdu1.*dAdu2.*dxdu1.*dxdu2+6.*d2xdu1.*dBdu1.*dxdu1.*B1+B1.*(d2Bdu1.*dxdu1.^2.*(1+vv)+3.*d2xdu1.^2.*B1+4.*d3xdu1.*dxdu1.*B1)).*h.^2)+A1.^6.*(6.*dhdu2.^2.*dxdu2.^2.*B1.^2+3.*dxdu2.*B1.*((-7).*dBdu2.*dhdu2.*dxdu2+(6.*d2xdu2.*dhdu2+d2hdu2.*dxdu2).*B1).*h+(15.*dBdu2.^2.*dxdu2.^2+(-18).*d2xdu2.*dBdu2.*dxdu2.*B1+B1.*((-4).*d2Bdu2.*dxdu2.^2+3.*d2xdu2.^2.*B1+4.*d3xdu2.*dxdu2.*B1)).*h.^2)+A1.^4.*B1.^2.*(6.*((-2).*dhdu1.*dhdu2.*dxdu1.*dxdu2.*((-1)+vv)+dhdu2.^2.*dxdu1.^2.*vv+dhdu1.^2.*dxdu2.^2.*vv).*B1.^2+3.*B1.*(dxdu2.^2.*((-3).*dBdu1.*dhdu1+d2hdu1.*vv.*B1)+dxdu1.*((-1).*dBdu2.*dhdu2.*dxdu1.*vv+(2.*d2xdu2.*dhdu1+4.*d2xdu1_du2.*dhdu2+d2hdu2.*dxdu1.*vv).*B1)+2.*dxdu2.*((-1).*dBdu2.*dhdu1.*dxdu1+(2.*d2xdu1_du2.*dhdu1+(-1).*d2hdu1_du2.*dxdu1.*((-1)+vv)).*B1+dhdu2.*(dBdu1.*dxdu1.*vv+d2xdu1.*B1))).*h+((-1).*dAdu2.^2.*dxdu2.^2+dxdu2.^2.*(4.*dBdu1.^2+d2Bdu1.*((-3)+vv).*B1)+2.*B1.*((-2).*d2xdu1_du2.*dBdu2.*dxdu1+2.*(d2xdu1_du2.^2+d3xdu1_2du2.*dxdu1).*B1+d2xdu2.*((-1).*dBdu1.*dxdu1+d2xdu1.*B1))+2.*dxdu2.*(dBdu2.*(dBdu1.*dxdu1+(-1).*d2xdu1.*B1)+2.*B1.*((-1).*d2xdu1_du2.*dBdu1+d3xd2u1_du2.*B1))).*h.^2));
C_69	=	(-1/6).*A1.^(-5).*B1.^(-5).*(15.*dAdu1.^2.*dxdu1.*dydu1.*B1.^6.*h.^2+(-1).*A1.^3.*B1.^3.*h.*(3.*(dxdu1.*(dhdu2.*(3.*dAdu2.*dydu1+dAdu1.*dydu2)+(-1).*dAdu2.*dhdu1.*dydu2.*vv)+dxdu2.*(dAdu1.*dhdu1.*dydu2.*vv+dydu1.*(dAdu1.*dhdu2+(-1).*dAdu2.*dhdu1.*vv))).*B1+(3.*dAdu2.*dBdu1.*dxdu1.*dydu2+dBdu2.*dxdu1.*((-1).*dAdu1.*dydu2+dAdu2.*dydu1.*((-3)+vv))+d2ydu2.*dAdu1.*dxdu1.*B1+2.*d2ydu1_du2.*dAdu2.*dxdu1.*B1+d2xdu2.*dAdu1.*dydu1.*B1+2.*d2xdu1_du2.*dAdu2.*dydu1.*B1+3.*d2Adu2.*dxdu1.*dydu1.*B1+2.*d2xdu1_du2.*dAdu1.*dydu2.*B1+d2xdu1.*dAdu2.*dydu2.*B1+(-1).*d2Adu2.*dxdu1.*dydu1.*vv.*B1+dxdu2.*(3.*dAdu2.*dBdu1.*dydu1+(-1).*dAdu1.*dBdu2.*dydu1+dAdu1.*dBdu1.*dydu2.*((-3)+vv)+2.*d2ydu1_du2.*dAdu1.*B1+d2ydu1.*dAdu2.*B1)).*h)+(-1).*A1.*B1.^5.*h.*(21.*dAdu1.*dhdu1.*dxdu1.*dydu1.*B1+(4.*d2Adu1.*dxdu1.*dydu1.*B1+dAdu1.*(dBdu1.*dxdu1.*dydu1.*(7+vv)+9.*(d2ydu1.*dxdu1+d2xdu1.*dydu1).*B1)).*h)+A1.^5.*B1.*h.*(3.*dAdu2.*dhdu2.*dxdu2.*dydu2.*(2+vv).*B1+(d2Adu2.*dxdu2.*dydu2.*(1+vv).*B1+dAdu2.*((-1).*dBdu2.*dxdu2.*dydu2.*(7+vv)+3.*(d2ydu2.*dxdu2+d2xdu2.*dydu2).*B1)).*h)+A1.^2.*B1.^4.*(6.*dhdu1.^2.*dxdu1.*dydu1.*B1.^2+3.*B1.*(dBdu1.*dhdu1.*dxdu1.*dydu1.*(2+vv)+(3.*d2xdu1.*dhdu1.*dydu1+dxdu1.*(3.*d2ydu1.*dhdu1+d2hdu1.*dydu1)).*B1).*h+(4.*dAdu2.^2.*dxdu1.*dydu1+(-1).*dBdu1.^2.*dxdu1.*dydu1+dAdu1.*dAdu2.*(dxdu2.*dydu1+dxdu1.*dydu2)+3.*dBdu1.*(d2ydu1.*dxdu1+d2xdu1.*dydu1).*B1+B1.*((3.*d2xdu1.*d2ydu1+2.*d3xdu1.*dydu1).*B1+dxdu1.*(d2Bdu1.*dydu1.*(1+vv)+2.*d3ydu1.*B1))).*h.^2)+A1.^6.*(6.*dhdu2.^2.*dxdu2.*dydu2.*B1.^2+3.*B1.*((-7).*dBdu2.*dhdu2.*dxdu2.*dydu2+(3.*d2xdu2.*dhdu2.*dydu2+dxdu2.*(3.*d2ydu2.*dhdu2+d2hdu2.*dydu2)).*B1).*h+(15.*dBdu2.^2.*dxdu2.*dydu2+(-9).*dBdu2.*(d2ydu2.*dxdu2+d2xdu2.*dydu2).*B1+B1.*((3.*d2xdu2.*d2ydu2+2.*d3xdu2.*dydu2).*B1+dxdu2.*((-4).*d2Bdu2.*dydu2+2.*d3ydu2.*B1))).*h.^2)+A1.^4.*B1.^2.*(6.*(dhdu2.*dxdu1.*((-1).*dhdu1.*dydu2.*((-1)+vv)+dhdu2.*dydu1.*vv)+dhdu1.*dxdu2.*((-1).*dhdu2.*dydu1.*((-1)+vv)+dhdu1.*dydu2.*vv)).*B1.^2+3.*B1.*((-1).*dBdu2.*dhdu1.*dxdu1.*dydu2+d2ydu2.*dhdu1.*dxdu1.*B1+d2xdu2.*dhdu1.*dydu1.*B1+2.*d2xdu1_du2.*dhdu1.*dydu2.*B1+d2hdu1_du2.*dxdu1.*dydu2.*B1+d2hdu2.*dxdu1.*dydu1.*vv.*B1+(-1).*d2hdu1_du2.*dxdu1.*dydu2.*vv.*B1+dhdu2.*((-1).*dBdu2.*dxdu1.*dydu1.*vv+2.*(d2ydu1_du2.*dxdu1+d2xdu1_du2.*dydu1).*B1+dydu2.*(dBdu1.*dxdu1.*vv+d2xdu1.*B1))+dxdu2.*((-1).*dBdu2.*dhdu1.*dydu1+2.*d2ydu1_du2.*dhdu1.*B1+d2hdu1_du2.*dydu1.*B1+(-1).*d2hdu1_du2.*dydu1.*vv.*B1+dhdu2.*(dBdu1.*dydu1.*vv+d2ydu1.*B1)+dydu2.*((-3).*dBdu1.*dhdu1+d2hdu1.*vv.*B1))).*h+(dBdu1.*dBdu2.*dxdu1.*dydu2+(-1).*dAdu2.^2.*dxdu2.*dydu2+(-1).*d2ydu2.*dBdu1.*dxdu1.*B1+(-2).*d2ydu1_du2.*dBdu2.*dxdu1.*B1+(-1).*d2xdu2.*dBdu1.*dydu1.*B1+(-2).*d2xdu1_du2.*dBdu2.*dydu1.*B1+(-2).*d2xdu1_du2.*dBdu1.*dydu2.*B1+(-1).*d2xdu1.*dBdu2.*dydu2.*B1+d2xdu2.*d2ydu1.*B1.^2+4.*d2xdu1_du2.*d2ydu1_du2.*B1.^2+d2xdu1.*d2ydu2.*B1.^2+2.*d3ydu1_2du2.*dxdu1.*B1.^2+2.*d3xdu1_2du2.*dydu1.*B1.^2+2.*d3xd2u1_du2.*dydu2.*B1.^2+dxdu2.*(dBdu2.*(dBdu1.*dydu1+(-1).*d2ydu1.*B1)+2.*B1.*((-1).*d2ydu1_du2.*dBdu1+d3yd2u1_du2.*B1)+dydu2.*(4.*dBdu1.^2+d2Bdu1.*((-3)+vv).*B1))).*h.^2));
C_70	=	(-1/12).*A1.^(-5).*B1.^(-5).*(15.*dAdu1.^2.*dydu1.^2.*B1.^6.*h.^2+(-1).*dydu1.*A1.*B1.^5.*h.*(21.*dAdu1.*dhdu1.*dydu1.*B1+(4.*d2Adu1.*dydu1.*B1+dAdu1.*(dBdu1.*dydu1.*(7+vv)+18.*d2ydu1.*B1)).*h)+dydu2.*A1.^5.*B1.*h.*(3.*dAdu2.*dhdu2.*dydu2.*(2+vv).*B1+(d2Adu2.*dydu2.*(1+vv).*B1+(-1).*dAdu2.*(dBdu2.*dydu2.*(7+vv)+(-6).*d2ydu2.*B1)).*h)+(-1).*A1.^3.*B1.^3.*h.*(3.*(3.*dAdu2.*dhdu2.*dydu1.^2+dAdu1.*dhdu1.*dydu2.^2.*vv+2.*dydu1.*dydu2.*(dAdu1.*dhdu2+(-1).*dAdu2.*dhdu1.*vv)).*B1+(dAdu1.*dBdu1.*dydu2.^2.*((-3)+vv)+dydu1.*(dAdu2.*dBdu2.*dydu1.*((-3)+vv)+(2.*d2ydu2.*dAdu1+4.*d2ydu1_du2.*dAdu2+(-1).*d2Adu2.*dydu1.*((-3)+vv)).*B1)+dydu2.*((-2).*dAdu1.*dBdu2.*dydu1+4.*d2ydu1_du2.*dAdu1.*B1+2.*dAdu2.*(3.*dBdu1.*dydu1+d2ydu1.*B1))).*h)+A1.^2.*B1.^4.*(6.*dhdu1.^2.*dydu1.^2.*B1.^2+3.*dydu1.*B1.*(dBdu1.*dhdu1.*dydu1.*(2+vv)+(6.*d2ydu1.*dhdu1+d2hdu1.*dydu1).*B1).*h+(4.*dAdu2.^2.*dydu1.^2+(-1).*dBdu1.^2.*dydu1.^2+2.*dAdu1.*dAdu2.*dydu1.*dydu2+6.*d2ydu1.*dBdu1.*dydu1.*B1+B1.*(d2Bdu1.*dydu1.^2.*(1+vv)+3.*d2ydu1.^2.*B1+4.*d3ydu1.*dydu1.*B1)).*h.^2)+A1.^6.*(6.*dhdu2.^2.*dydu2.^2.*B1.^2+3.*dydu2.*B1.*((-7).*dBdu2.*dhdu2.*dydu2+(6.*d2ydu2.*dhdu2+d2hdu2.*dydu2).*B1).*h+(15.*dBdu2.^2.*dydu2.^2+(-18).*d2ydu2.*dBdu2.*dydu2.*B1+B1.*((-4).*d2Bdu2.*dydu2.^2+3.*d2ydu2.^2.*B1+4.*d3ydu2.*dydu2.*B1)).*h.^2)+A1.^4.*B1.^2.*(6.*((-2).*dhdu1.*dhdu2.*dydu1.*dydu2.*((-1)+vv)+dhdu2.^2.*dydu1.^2.*vv+dhdu1.^2.*dydu2.^2.*vv).*B1.^2+3.*B1.*(dydu2.^2.*((-3).*dBdu1.*dhdu1+d2hdu1.*vv.*B1)+dydu1.*((-1).*dBdu2.*dhdu2.*dydu1.*vv+(2.*d2ydu2.*dhdu1+4.*d2ydu1_du2.*dhdu2+d2hdu2.*dydu1.*vv).*B1)+2.*dydu2.*((-1).*dBdu2.*dhdu1.*dydu1+(2.*d2ydu1_du2.*dhdu1+(-1).*d2hdu1_du2.*dydu1.*((-1)+vv)).*B1+dhdu2.*(dBdu1.*dydu1.*vv+d2ydu1.*B1))).*h+((-1).*dAdu2.^2.*dydu2.^2+dydu2.^2.*(4.*dBdu1.^2+d2Bdu1.*((-3)+vv).*B1)+2.*B1.*((-2).*d2ydu1_du2.*dBdu2.*dydu1+2.*(d2ydu1_du2.^2+d3ydu1_2du2.*dydu1).*B1+d2ydu2.*((-1).*dBdu1.*dydu1+d2ydu1.*B1))+2.*dydu2.*(dBdu2.*(dBdu1.*dydu1+(-1).*d2ydu1.*B1)+2.*B1.*((-1).*d2ydu1_du2.*dBdu1+d3yd2u1_du2.*B1))).*h.^2));
C_71	=	(1/12).*R1.^(-3).*R2.^(-1).*A1.^(-4).*B1.^(-3).*((-3).*dAdu1.^2.*dxdu1.*R1.^2.*R2.*B1.^4.*h.^2+R1.*R2.*A1.*B1.^3.*h.*(9.*dAdu1.*dhdu1.*dxdu1.*R1.*B1+((-6).*dAdu1.*dR1du1.*dxdu1.*B1+R1.*(d2Adu1.*dxdu1.*B1+dAdu1.*(dBdu1.*dxdu1.*(3+vv)+3.*d2xdu1.*B1))).*h)+R1.*R2.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dxdu1.*R1.*(3+(-2).*vv).*B1+((-1).*dAdu2.*(dR1du1.*dxdu2+dR1du2.*dxdu1.*(3+(-2).*vv)).*B1+R1.*(d2Adu2.*dxdu1.*(3+(-2).*vv).*B1+dAdu2.*((-1).*dBdu1.*dxdu2.*((-2)+vv)+dBdu2.*dxdu1.*((-3)+2.*vv)+d2xdu1_du2.*B1))).*h)+(-1).*R2.*A1.^2.*B1.^2.*(6.*dhdu1.^2.*dxdu1.*R1.^2.*B1.^2+3.*R1.*B1.*((-4).*dhdu1.*dR1du1.*dxdu1.*B1+R1.*(dBdu1.*dhdu1.*dxdu1.*(2+vv)+(2.*d2xdu1.*dhdu1+d2hdu1.*dxdu1).*B1)).*h+(6.*dR1du1.^2.*dxdu1.*B1.^2+(-1).*R1.*B1.*(4.*dBdu1.*dR1du1.*dxdu1+3.*(d2xdu1.*dR1du1+d2R1du1.*dxdu1).*B1)+R1.^2.*(dAdu2.^2.*dxdu1+(-1).*dBdu1.^2.*dxdu1+2.*d2xdu1.*dBdu1.*B1+B1.*(d2Bdu1.*dxdu1.*(1+vv)+d3xdu1.*B1))).*h.^2)+A1.^4.*(6.*R1.^2.*B1.^2.*(2.*dxdu1.*R1.*vv.*B1.^2+R2.*(2.*dhdu1.*dhdu2.*dxdu2.*((-1)+vv)+(-1).*dhdu2.^2.*dxdu1.*vv+2.*dxdu1.*B1.^2))+3.*R1.*R2.*B1.*(2.*(dhdu2.*dR1du1.*dxdu2+dR1du2.*((-1).*dhdu1.*dxdu2.*((-1)+vv)+dhdu2.*dxdu1.*vv)).*B1+(-1).*R1.*(2.*dBdu2.*dhdu1.*dxdu2.*((-1)+vv)+((-2).*d2xdu2.*dhdu1.*((-1)+vv)+(-2).*d2hdu1_du2.*dxdu2.*((-1)+vv)+d2hdu2.*dxdu1.*vv).*B1+dhdu2.*(2.*dBdu1.*dxdu2+(-1).*dBdu2.*dxdu1.*vv+2.*d2xdu1_du2.*B1))).*h+R2.*((-2).*dR1du2.*(2.*dR1du1.*dxdu2+dR1du2.*dxdu1.*vv).*B1.^2+R1.*B1.*(dBdu2.*dR1du1.*dxdu2.*((-2)+vv)+(2.*d2R1du1_du2.*dxdu2+(-1).*d2xdu2.*dR1du1.*((-2)+vv)+d2R1du2.*dxdu1.*vv).*B1+dR1du2.*(2.*dBdu1.*dxdu2+(-1).*dBdu2.*dxdu1.*vv+2.*d2xdu1_du2.*B1))+R1.^2.*(dBdu2.*(3.*dBdu1.*dxdu2+(-1).*d2xdu1_du2.*((-2)+vv).*B1)+B1.*((-1).*d2xdu2.*dBdu1+(-2).*d2Bdu1_du2.*dxdu2+d3xdu1_2du2.*((-2)+vv).*B1))).*h.^2));
C_72	=	(1/12).*R1.^(-3).*R2.^(-1).*A1.^(-4).*B1.^(-3).*((-3).*dAdu1.^2.*dydu1.*R1.^2.*R2.*B1.^4.*h.^2+R1.*R2.*A1.*B1.^3.*h.*(9.*dAdu1.*dhdu1.*dydu1.*R1.*B1+((-6).*dAdu1.*dR1du1.*dydu1.*B1+R1.*(d2Adu1.*dydu1.*B1+dAdu1.*(dBdu1.*dydu1.*(3+vv)+3.*d2ydu1.*B1))).*h)+R1.*R2.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dydu1.*R1.*(3+(-2).*vv).*B1+((-1).*dAdu2.*(dR1du1.*dydu2+dR1du2.*dydu1.*(3+(-2).*vv)).*B1+R1.*(d2Adu2.*dydu1.*(3+(-2).*vv).*B1+dAdu2.*((-1).*dBdu1.*dydu2.*((-2)+vv)+dBdu2.*dydu1.*((-3)+2.*vv)+d2ydu1_du2.*B1))).*h)+(-1).*R2.*A1.^2.*B1.^2.*(6.*dhdu1.^2.*dydu1.*R1.^2.*B1.^2+3.*R1.*B1.*((-4).*dhdu1.*dR1du1.*dydu1.*B1+R1.*(dBdu1.*dhdu1.*dydu1.*(2+vv)+(2.*d2ydu1.*dhdu1+d2hdu1.*dydu1).*B1)).*h+(6.*dR1du1.^2.*dydu1.*B1.^2+(-1).*R1.*B1.*(4.*dBdu1.*dR1du1.*dydu1+3.*(d2ydu1.*dR1du1+d2R1du1.*dydu1).*B1)+R1.^2.*(dAdu2.^2.*dydu1+(-1).*dBdu1.^2.*dydu1+2.*d2ydu1.*dBdu1.*B1+B1.*(d2Bdu1.*dydu1.*(1+vv)+d3ydu1.*B1))).*h.^2)+A1.^4.*(6.*R1.^2.*B1.^2.*(2.*dydu1.*R1.*vv.*B1.^2+R2.*(2.*dhdu1.*dhdu2.*dydu2.*((-1)+vv)+(-1).*dhdu2.^2.*dydu1.*vv+2.*dydu1.*B1.^2))+3.*R1.*R2.*B1.*(2.*(dhdu2.*dR1du1.*dydu2+dR1du2.*((-1).*dhdu1.*dydu2.*((-1)+vv)+dhdu2.*dydu1.*vv)).*B1+(-1).*R1.*(2.*dBdu2.*dhdu1.*dydu2.*((-1)+vv)+((-2).*d2ydu2.*dhdu1.*((-1)+vv)+(-2).*d2hdu1_du2.*dydu2.*((-1)+vv)+d2hdu2.*dydu1.*vv).*B1+dhdu2.*(2.*dBdu1.*dydu2+(-1).*dBdu2.*dydu1.*vv+2.*d2ydu1_du2.*B1))).*h+R2.*((-2).*dR1du2.*(2.*dR1du1.*dydu2+dR1du2.*dydu1.*vv).*B1.^2+R1.*B1.*(dBdu2.*dR1du1.*dydu2.*((-2)+vv)+(2.*d2R1du1_du2.*dydu2+(-1).*d2ydu2.*dR1du1.*((-2)+vv)+d2R1du2.*dydu1.*vv).*B1+dR1du2.*(2.*dBdu1.*dydu2+(-1).*dBdu2.*dydu1.*vv+2.*d2ydu1_du2.*B1))+R1.^2.*(dBdu2.*(3.*dBdu1.*dydu2+(-1).*d2ydu1_du2.*((-2)+vv).*B1)+B1.*((-1).*d2ydu2.*dBdu1+(-2).*d2Bdu1_du2.*dydu2+d3ydu1_2du2.*((-2)+vv).*B1))).*h.^2));
C_73	=	(1/12).*R1.^(-1).*R2.^(-3).*A1.^(-3).*B1.^(-4).*(3.*dAdu1.*dAdu2.*dxdu1.*R1.*R2.^2.*B1.^4.*h.^2+(-1).*R1.*R2.*A1.*B1.^3.*h.*((-3).*R2.*((-2).*dxdu1.*(dAdu2.*dhdu1+dAdu1.*dhdu2.*((-1)+vv))+dAdu1.*dhdu1.*dxdu2.*vv).*B1+(((-1).*dxdu1.*(2.*dAdu2.*dR2du1+dAdu1.*dR2du2.*((-2)+vv))+dAdu1.*dR2du1.*dxdu2.*vv).*B1+R2.*(dAdu1.*dBdu1.*dxdu2.*(3+(-2).*vv)+(2.*d2Adu1_du2.*dxdu1+d2xdu1_du2.*dAdu1.*((-2)+vv)).*B1+dAdu2.*(dBdu1.*dxdu1.*((-2)+vv)+d2xdu1.*B1))).*h)+(-1).*R1.*R2.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dxdu2.*R2.*(2+vv).*B1+((-4).*dAdu2.*dR2du2.*dxdu2.*B1+R2.*(d2Adu2.*dxdu2.*(1+vv).*B1+(-1).*dAdu2.*(dBdu2.*dxdu2.*(3+vv)+(-2).*d2xdu2.*B1))).*h)+R1.*A1.^2.*B1.^2.*((-6).*dhdu1.*R2.^2.*((-2).*dhdu2.*dxdu1.*((-1)+vv)+dhdu1.*dxdu2.*vv).*B1.^2+(-3).*R2.*B1.*((-2).*(dxdu1.*(dhdu1.*dR2du2+(-1).*dhdu2.*dR2du1.*((-1)+vv))+dhdu1.*dR2du1.*dxdu2.*vv).*B1+R2.*(2.*(d2xdu1_du2.*dhdu1+(-1).*d2xdu1.*dhdu2.*((-1)+vv)+(-1).*d2hdu1_du2.*dxdu1.*((-1)+vv)).*B1+dxdu2.*(dBdu1.*dhdu1.*((-3)+2.*vv)+d2hdu1.*vv.*B1))).*h+((-2).*dR2du1.*(2.*dR2du2.*dxdu1+dR2du1.*dxdu2.*vv).*B1.^2+R2.*B1.*(2.*(d2xdu1_du2.*dR2du1+d2R2du1_du2.*dxdu1).*B1+(-1).*dR2du2.*(dBdu1.*dxdu1+d2xdu1.*((-2)+vv).*B1)+dxdu2.*(dBdu1.*dR2du1.*((-3)+2.*vv)+d2R2du1.*vv.*B1))+R2.^2.*(dAdu2.^2.*dxdu2+B1.*(d2xdu1_du2.*dBdu1+d3xd2u1_du2.*((-2)+vv).*B1)+(-1).*dxdu2.*(dBdu1.^2+d2Bdu1.*((-3)+2.*vv).*B1))).*h.^2)+A1.^4.*(12.*dxdu2.*R2.^2.*(R1+R2.*vv).*B1.^4+(-3).*dBdu2.^2.*dxdu2.*R1.*R2.^2.*h.^2+R1.*R2.*B1.*h.*((-6).*dBdu2.*dR2du2.*dxdu2.*h+R2.*(d2Bdu2.*dxdu2.*h+3.*dBdu2.*(3.*dhdu2.*dxdu2+d2xdu2.*h)))+(-1).*R1.*B1.^2.*(6.*dR2du2.^2.*dxdu2.*h.^2+(-3).*R2.*h.*(4.*dhdu2.*dR2du2.*dxdu2+(d2xdu2.*dR2du2+d2R2du2.*dxdu2).*h)+R2.^2.*(6.*dhdu2.^2.*dxdu2+6.*d2xdu2.*dhdu2.*h+h.*(3.*d2hdu2.*dxdu2+d3xdu2.*h)))));
C_74	=	(1/12).*R1.^(-1).*R2.^(-3).*A1.^(-3).*B1.^(-4).*(3.*dAdu1.*dAdu2.*dydu1.*R1.*R2.^2.*B1.^4.*h.^2+(-1).*R1.*R2.*A1.*B1.^3.*h.*((-3).*R2.*((-2).*dydu1.*(dAdu2.*dhdu1+dAdu1.*dhdu2.*((-1)+vv))+dAdu1.*dhdu1.*dydu2.*vv).*B1+(((-1).*dydu1.*(2.*dAdu2.*dR2du1+dAdu1.*dR2du2.*((-2)+vv))+dAdu1.*dR2du1.*dydu2.*vv).*B1+R2.*(dAdu1.*dBdu1.*dydu2.*(3+(-2).*vv)+(2.*d2Adu1_du2.*dydu1+d2ydu1_du2.*dAdu1.*((-2)+vv)).*B1+dAdu2.*(dBdu1.*dydu1.*((-2)+vv)+d2ydu1.*B1))).*h)+(-1).*R1.*R2.*A1.^3.*B1.*h.*(3.*dAdu2.*dhdu2.*dydu2.*R2.*(2+vv).*B1+((-4).*dAdu2.*dR2du2.*dydu2.*B1+R2.*(d2Adu2.*dydu2.*(1+vv).*B1+(-1).*dAdu2.*(dBdu2.*dydu2.*(3+vv)+(-2).*d2ydu2.*B1))).*h)+R1.*A1.^2.*B1.^2.*((-6).*dhdu1.*R2.^2.*((-2).*dhdu2.*dydu1.*((-1)+vv)+dhdu1.*dydu2.*vv).*B1.^2+(-3).*R2.*B1.*((-2).*(dydu1.*(dhdu1.*dR2du2+(-1).*dhdu2.*dR2du1.*((-1)+vv))+dhdu1.*dR2du1.*dydu2.*vv).*B1+R2.*(2.*(d2ydu1_du2.*dhdu1+(-1).*d2ydu1.*dhdu2.*((-1)+vv)+(-1).*d2hdu1_du2.*dydu1.*((-1)+vv)).*B1+dydu2.*(dBdu1.*dhdu1.*((-3)+2.*vv)+d2hdu1.*vv.*B1))).*h+((-2).*dR2du1.*(2.*dR2du2.*dydu1+dR2du1.*dydu2.*vv).*B1.^2+R2.*B1.*(2.*(d2ydu1_du2.*dR2du1+d2R2du1_du2.*dydu1).*B1+(-1).*dR2du2.*(dBdu1.*dydu1+d2ydu1.*((-2)+vv).*B1)+dydu2.*(dBdu1.*dR2du1.*((-3)+2.*vv)+d2R2du1.*vv.*B1))+R2.^2.*(dAdu2.^2.*dydu2+B1.*(d2ydu1_du2.*dBdu1+d3yd2u1_du2.*((-2)+vv).*B1)+(-1).*dydu2.*(dBdu1.^2+d2Bdu1.*((-3)+2.*vv).*B1))).*h.^2)+A1.^4.*(12.*dydu2.*R2.^2.*(R1+R2.*vv).*B1.^4+(-3).*dBdu2.^2.*dydu2.*R1.*R2.^2.*h.^2+R1.*R2.*B1.*h.*((-6).*dBdu2.*dR2du2.*dydu2.*h+R2.*(d2Bdu2.*dydu2.*h+3.*dBdu2.*(3.*dhdu2.*dydu2+d2ydu2.*h)))+(-1).*R1.*B1.^2.*(6.*dR2du2.^2.*dydu2.*h.^2+(-3).*R2.*h.*(4.*dhdu2.*dR2du2.*dydu2+(d2ydu2.*dR2du2+d2R2du2.*dydu2).*h)+R2.^2.*(6.*dhdu2.^2.*dydu2+6.*d2ydu2.*dhdu2.*h+h.*(3.*d2hdu2.*dydu2+d3ydu2.*h)))));
C_75	=	(1/12).*A1.^(-6).*B1.^(-6).*(15.*dAdu1.^3.*dxdu1.*B1.^7.*h.^2+(-1).*dAdu1.*A1.*B1.^6.*h.*(21.*dAdu1.*dhdu1.*dxdu1.*B1+(10.*d2Adu1.*dxdu1.*B1+dAdu1.*(dBdu1.*dxdu1.*(7+4.*vv)+15.*d2xdu1.*B1)).*h)+A1.^2.*B1.^5.*(6.*dAdu1.*dhdu1.^2.*dxdu1.*B1.^2+3.*B1.*(2.*d2Adu1.*dhdu1.*dxdu1.*B1+dAdu1.*(2.*dBdu1.*dhdu1.*dxdu1.*(1+2.*vv)+(7.*d2xdu1.*dhdu1+d2hdu1.*dxdu1).*B1)).*h+(9.*dAdu1.*dAdu2.^2.*dxdu1+(-3).*dAdu1.^2.*dAdu2.*dxdu2+B1.*(d2Adu1.*dBdu1.*dxdu1.*(2+vv)+(4.*d2Adu1.*d2xdu1+d3Adu1.*dxdu1).*B1)+dAdu1.*((-3).*dBdu1.^2.*dxdu1+d2xdu1.*dBdu1.*(7+vv).*B1+B1.*(d2Bdu1.*dxdu1.*(1+4.*vv)+6.*d3xdu1.*B1))).*h.^2)+(-1).*A1.^6.*B1.*(6.*dAdu2.*dhdu2.^2.*dxdu2.*vv.*B1.^2+3.*B1.*(2.*d2Adu2.*dhdu2.*dxdu2.*vv.*B1+dAdu2.*((-2).*dBdu2.*dhdu2.*dxdu2.*(1+2.*vv)+(d2hdu2.*dxdu2.*vv+d2xdu2.*dhdu2.*(2+vv)).*B1)).*h+(B1.*((-1).*d2Adu2.*dBdu2.*dxdu2.*(1+4.*vv)+(d3Adu2.*dxdu2.*vv+d2Adu2.*d2xdu2.*(1+vv)).*B1)+dAdu2.*(dBdu2.^2.*dxdu2.*(7+4.*vv)+(-1).*d2xdu2.*dBdu2.*(7+vv).*B1+B1.*((-1).*d2Bdu2.*dxdu2.*(2+vv)+2.*d3xdu2.*B1))).*h.^2)+A1.^4.*B1.^3.*(6.*((-1).*dAdu2.*dhdu1.*(dhdu1.*dxdu2+2.*dhdu2.*dxdu1.*((-1)+vv))+dAdu1.*dhdu2.^2.*dxdu1.*vv).*B1.^2+3.*B1.*(2.*dAdu2.*dBdu1.*dhdu1.*dxdu2+dAdu2.*dBdu1.*dhdu1.*dxdu2.*vv+(-1).*dBdu2.*dhdu1.*((-2).*dAdu2.*dxdu1.*((-1)+vv)+dAdu1.*dxdu2.*vv)+2.*d2hdu1_du2.*dAdu2.*dxdu1.*B1+2.*d2Adu2.*dhdu1.*dxdu1.*B1+(-1).*d2hdu1.*dAdu2.*dxdu2.*B1+(-2).*d2Adu1_du2.*dhdu1.*dxdu2.*B1+d2xdu2.*dAdu1.*dhdu1.*vv.*B1+(-2).*d2xdu1_du2.*dAdu2.*dhdu1.*vv.*B1+d2hdu2.*dAdu1.*dxdu1.*vv.*B1+(-2).*d2hdu1_du2.*dAdu2.*dxdu1.*vv.*B1+(-2).*d2Adu2.*dhdu1.*dxdu1.*vv.*B1+dhdu2.*(2.*dAdu2.*dBdu1.*dxdu1+2.*dAdu1.*dBdu1.*dxdu2.*((-1)+vv)+dAdu2.*dBdu1.*dxdu1.*vv+(-1).*dAdu1.*dBdu2.*dxdu1.*vv+2.*d2xdu1_du2.*dAdu1.*B1+3.*d2xdu1.*dAdu2.*B1+2.*d2Adu1_du2.*dxdu1.*B1)).*h+((-1).*dAdu2.^3.*dxdu2+dBdu2.*(dAdu1.*dBdu1.*dxdu2.*(5+(-2).*vv)+((-2).*d2xdu1_du2.*dAdu1+d2Adu1_du2.*dxdu1.*((-2)+vv)).*B1)+B1.*((-2).*d2Bdu1_du2.*dAdu1.*dxdu2+2.*d2Adu1_du2.*dBdu1.*dxdu2+d2xdu2.*dAdu1.*dBdu1.*((-3)+vv)+d2Bdu1_du2.*dAdu1.*dxdu2.*vv+d2Adu1_du2.*dBdu1.*dxdu2.*vv+2.*d3xdu1_2du2.*dAdu1.*B1+2.*d3Adu1_2du2.*dxdu1.*B1+(-1).*d3Ad2u1_du2.*dxdu2.*B1+(-1).*d3Adu1_2du2.*dxdu1.*vv.*B1+d2Adu2.*(dBdu1.*dxdu1.*(1+vv)+(-1).*d2xdu1.*((-3)+vv).*B1))+dAdu2.*(B1.*(6.*d2xdu1_du2.*dBdu1+d2Bdu1_du2.*dxdu1.*(2+vv)+2.*d3xd2u1_du2.*B1)+dBdu2.*((-1).*dBdu1.*dxdu1.*(3+2.*vv)+d2xdu1.*((-3)+vv).*B1)+dxdu2.*((-4).*dBdu1.^2+d2Bdu1.*(1+vv).*B1))).*h.^2)+A1.^7.*(6.*dhdu2.^2.*B1.^2.*(dBdu2.*dxdu2+(-1).*d2xdu2.*B1)+(-3).*B1.*(7.*dBdu2.^2.*dhdu2.*dxdu2+(-1).*dBdu2.*(7.*d2xdu2.*dhdu2+d2hdu2.*dxdu2).*B1+B1.*(d2hdu2.*d2xdu2.*B1+dhdu2.*((-2).*d2Bdu2.*dxdu2+2.*d3xdu2.*B1))).*h+h.^2.*(15.*dBdu2.^3.*dxdu2+(-15).*d2xdu2.*dBdu2.^2.*B1+2.*dBdu2.*B1.*((-5).*d2Bdu2.*dxdu2+3.*d3xdu2.*B1)+B1.^2.*(4.*d2Bdu2.*d2xdu2+d3Bdu2.*dxdu2+(-1).*B1.*d4xdu2)))+(-1).*A1.^5.*B1.^2.*(6.*B1.^2.*(2.*dhdu1.*dhdu2.*((-1)+vv).*(dBdu1.*dxdu2+(-1).*d2xdu1_du2.*B1)+dhdu1.^2.*vv.*((-1).*dBdu2.*dxdu2+d2xdu2.*B1)+dhdu2.^2.*(dBdu1.*dxdu1+d2xdu1.*vv.*B1))+3.*B1.*(5.*dBdu1.*dBdu2.*dhdu1.*dxdu2+(-1).*dAdu2.^2.*dhdu2.*dxdu2+(-3).*d2xdu2.*dBdu1.*dhdu1.*B1+(-2).*d2xdu1_du2.*dBdu2.*dhdu1.*B1+d2hdu2.*dBdu1.*dxdu1.*B1+(-2).*d2hdu1_du2.*dBdu1.*dxdu2.*B1+(-2).*d2Bdu1_du2.*dhdu1.*dxdu2.*B1+2.*d2hdu1_du2.*dBdu1.*dxdu2.*vv.*B1+(-1).*d2hdu1.*dBdu2.*dxdu2.*vv.*B1+2.*d2hdu1_du2.*d2xdu1_du2.*B1.^2+2.*d3xdu1_2du2.*dhdu1.*B1.^2+d2hdu2.*d2xdu1.*vv.*B1.^2+(-2).*d2hdu1_du2.*d2xdu1_du2.*vv.*B1.^2+d2hdu1.*d2xdu2.*vv.*B1.^2+dhdu2.*((-2).*dxdu2.*((-1)+vv).*(dBdu1.^2+(-1).*d2Bdu1.*B1)+2.*B1.*(d2Bdu1_du2.*dxdu1+d2xdu1_du2.*dBdu1.*vv+d3xd2u1_du2.*B1)+(-1).*dBdu2.*(3.*dBdu1.*dxdu1+d2xdu1.*vv.*B1))).*h+h.^2.*(3.*dBdu1.*dBdu2.^2.*dxdu1+(-2).*d2Adu2.*dAdu2.*dxdu2.*B1+dAdu2.^2.*(3.*dBdu2.*dxdu2+(-1).*d2xdu2.*B1)+(-1).*dBdu2.*(B1.*((-2).*d2xdu1_du2.*dBdu1+3.*d2Bdu1_du2.*dxdu1+2.*d3xd2u1_du2.*B1)+dxdu2.*(9.*dBdu1.^2+d2Bdu1.*((-5)+2.*vv).*B1))+B1.*((-1).*d2Bdu2.*dBdu1.*dxdu1+5.*d2Bdu1_du2.*dBdu1.*dxdu2+(-2).*d3xdu1_2du2.*dBdu1.*B1+d3Bdu1_2du2.*dxdu1.*B1+(-2).*d3Bd2u1_du2.*dxdu2.*B1+d3Bd2u1_du2.*dxdu2.*vv.*B1+d2xdu2.*(4.*dBdu1.^2+d2Bdu1.*((-3)+vv).*B1)+2.*B1.^2.*d4xd2u1_2du2)))+(-1).*A1.^3.*B1.^4.*(6.*dhdu1.^2.*B1.^2.*(dBdu1.*dxdu1.*vv+d2xdu1.*B1)+3.*B1.*((-1).*dBdu1.^2.*dhdu1.*dxdu1+dAdu1.*dAdu2.*(5.*dhdu2.*dxdu1+(-3).*dhdu1.*dxdu2)+(-2).*dAdu2.^2.*dhdu1.*dxdu1.*((-1)+vv)+dBdu1.*(d2hdu1.*dxdu1.*vv+d2xdu1.*dhdu1.*(2+vv)).*B1+B1.*(d2hdu1.*d2xdu1.*B1+2.*dhdu1.*(d2Bdu1.*dxdu1.*vv+d3xdu1.*B1))).*h+h.^2.*(dBdu1.^3.*dxdu1+(-1).*(dBdu1.*(d2xdu1.*dBdu1+2.*d2Bdu1.*dxdu1)+3.*d2Adu1_du2.*dAdu1.*dxdu2+d2Adu2.*dAdu1.*dxdu1.*((-5)+2.*vv)).*B1+(2.*d3xdu1.*dBdu1+d3Bdu1.*dxdu1.*vv+d2Bdu1.*d2xdu1.*(1+vv)).*B1.^2+4.*dAdu2.^2.*(dBdu1.*dxdu1+d2xdu1.*B1)+dAdu2.*(dAdu1.*dBdu2.*dxdu1.*((-5)+2.*vv)+(2.*d2xdu1_du2.*dAdu1+5.*d2Adu1_du2.*dxdu1).*B1+dxdu2.*(dAdu1.*dBdu1.*(3+2.*vv)+(-1).*d2Adu1.*B1))+B1.^3.*d4xdu1)));
C_76	=	(1/12).*A1.^(-6).*B1.^(-6).*(15.*dAdu1.^3.*dydu1.*B1.^7.*h.^2+(-1).*dAdu1.*A1.*B1.^6.*h.*(21.*dAdu1.*dhdu1.*dydu1.*B1+(10.*d2Adu1.*dydu1.*B1+dAdu1.*(dBdu1.*dydu1.*(7+4.*vv)+15.*d2ydu1.*B1)).*h)+A1.^2.*B1.^5.*(6.*dAdu1.*dhdu1.^2.*dydu1.*B1.^2+3.*B1.*(2.*d2Adu1.*dhdu1.*dydu1.*B1+dAdu1.*(2.*dBdu1.*dhdu1.*dydu1.*(1+2.*vv)+(7.*d2ydu1.*dhdu1+d2hdu1.*dydu1).*B1)).*h+(9.*dAdu1.*dAdu2.^2.*dydu1+(-3).*dAdu1.^2.*dAdu2.*dydu2+B1.*(d2Adu1.*dBdu1.*dydu1.*(2+vv)+(4.*d2Adu1.*d2ydu1+d3Adu1.*dydu1).*B1)+dAdu1.*((-3).*dBdu1.^2.*dydu1+d2ydu1.*dBdu1.*(7+vv).*B1+B1.*(d2Bdu1.*dydu1.*(1+4.*vv)+6.*d3ydu1.*B1))).*h.^2)+(-1).*A1.^6.*B1.*(6.*dAdu2.*dhdu2.^2.*dydu2.*vv.*B1.^2+3.*B1.*(2.*d2Adu2.*dhdu2.*dydu2.*vv.*B1+dAdu2.*((-2).*dBdu2.*dhdu2.*dydu2.*(1+2.*vv)+(d2hdu2.*dydu2.*vv+d2ydu2.*dhdu2.*(2+vv)).*B1)).*h+(B1.*((-1).*d2Adu2.*dBdu2.*dydu2.*(1+4.*vv)+(d3Adu2.*dydu2.*vv+d2Adu2.*d2ydu2.*(1+vv)).*B1)+dAdu2.*(dBdu2.^2.*dydu2.*(7+4.*vv)+(-1).*d2ydu2.*dBdu2.*(7+vv).*B1+B1.*((-1).*d2Bdu2.*dydu2.*(2+vv)+2.*d3ydu2.*B1))).*h.^2)+A1.^4.*B1.^3.*(6.*((-1).*dAdu2.*dhdu1.*(dhdu1.*dydu2+2.*dhdu2.*dydu1.*((-1)+vv))+dAdu1.*dhdu2.^2.*dydu1.*vv).*B1.^2+3.*B1.*(2.*dAdu2.*dBdu1.*dhdu1.*dydu2+dAdu2.*dBdu1.*dhdu1.*dydu2.*vv+(-1).*dBdu2.*dhdu1.*((-2).*dAdu2.*dydu1.*((-1)+vv)+dAdu1.*dydu2.*vv)+2.*d2hdu1_du2.*dAdu2.*dydu1.*B1+2.*d2Adu2.*dhdu1.*dydu1.*B1+(-1).*d2hdu1.*dAdu2.*dydu2.*B1+(-2).*d2Adu1_du2.*dhdu1.*dydu2.*B1+d2ydu2.*dAdu1.*dhdu1.*vv.*B1+(-2).*d2ydu1_du2.*dAdu2.*dhdu1.*vv.*B1+d2hdu2.*dAdu1.*dydu1.*vv.*B1+(-2).*d2hdu1_du2.*dAdu2.*dydu1.*vv.*B1+(-2).*d2Adu2.*dhdu1.*dydu1.*vv.*B1+dhdu2.*(2.*dAdu2.*dBdu1.*dydu1+2.*dAdu1.*dBdu1.*dydu2.*((-1)+vv)+dAdu2.*dBdu1.*dydu1.*vv+(-1).*dAdu1.*dBdu2.*dydu1.*vv+2.*d2ydu1_du2.*dAdu1.*B1+3.*d2ydu1.*dAdu2.*B1+2.*d2Adu1_du2.*dydu1.*B1)).*h+((-1).*dAdu2.^3.*dydu2+dBdu2.*(dAdu1.*dBdu1.*dydu2.*(5+(-2).*vv)+((-2).*d2ydu1_du2.*dAdu1+d2Adu1_du2.*dydu1.*((-2)+vv)).*B1)+B1.*((-2).*d2Bdu1_du2.*dAdu1.*dydu2+2.*d2Adu1_du2.*dBdu1.*dydu2+d2ydu2.*dAdu1.*dBdu1.*((-3)+vv)+d2Bdu1_du2.*dAdu1.*dydu2.*vv+d2Adu1_du2.*dBdu1.*dydu2.*vv+2.*d3ydu1_2du2.*dAdu1.*B1+2.*d3Adu1_2du2.*dydu1.*B1+(-1).*d3Ad2u1_du2.*dydu2.*B1+(-1).*d3Adu1_2du2.*dydu1.*vv.*B1+d2Adu2.*(dBdu1.*dydu1.*(1+vv)+(-1).*d2ydu1.*((-3)+vv).*B1))+dAdu2.*(B1.*(6.*d2ydu1_du2.*dBdu1+d2Bdu1_du2.*dydu1.*(2+vv)+2.*d3yd2u1_du2.*B1)+dBdu2.*((-1).*dBdu1.*dydu1.*(3+2.*vv)+d2ydu1.*((-3)+vv).*B1)+dydu2.*((-4).*dBdu1.^2+d2Bdu1.*(1+vv).*B1))).*h.^2)+A1.^7.*(6.*dhdu2.^2.*B1.^2.*(dBdu2.*dydu2+(-1).*d2ydu2.*B1)+(-3).*B1.*(7.*dBdu2.^2.*dhdu2.*dydu2+(-1).*dBdu2.*(7.*d2ydu2.*dhdu2+d2hdu2.*dydu2).*B1+B1.*(d2hdu2.*d2ydu2.*B1+dhdu2.*((-2).*d2Bdu2.*dydu2+2.*d3ydu2.*B1))).*h+h.^2.*(15.*dBdu2.^3.*dydu2+(-15).*d2ydu2.*dBdu2.^2.*B1+2.*dBdu2.*B1.*((-5).*d2Bdu2.*dydu2+3.*d3ydu2.*B1)+B1.^2.*(4.*d2Bdu2.*d2ydu2+d3Bdu2.*dydu2+(-1).*B1.*d4ydu2)))+(-1).*A1.^5.*B1.^2.*(6.*B1.^2.*(2.*dhdu1.*dhdu2.*((-1)+vv).*(dBdu1.*dydu2+(-1).*d2ydu1_du2.*B1)+dhdu1.^2.*vv.*((-1).*dBdu2.*dydu2+d2ydu2.*B1)+dhdu2.^2.*(dBdu1.*dydu1+d2ydu1.*vv.*B1))+3.*B1.*(5.*dBdu1.*dBdu2.*dhdu1.*dydu2+(-1).*dAdu2.^2.*dhdu2.*dydu2+(-3).*d2ydu2.*dBdu1.*dhdu1.*B1+(-2).*d2ydu1_du2.*dBdu2.*dhdu1.*B1+d2hdu2.*dBdu1.*dydu1.*B1+(-2).*d2hdu1_du2.*dBdu1.*dydu2.*B1+(-2).*d2Bdu1_du2.*dhdu1.*dydu2.*B1+2.*d2hdu1_du2.*dBdu1.*dydu2.*vv.*B1+(-1).*d2hdu1.*dBdu2.*dydu2.*vv.*B1+2.*d2hdu1_du2.*d2ydu1_du2.*B1.^2+2.*d3ydu1_2du2.*dhdu1.*B1.^2+d2hdu2.*d2ydu1.*vv.*B1.^2+(-2).*d2hdu1_du2.*d2ydu1_du2.*vv.*B1.^2+d2hdu1.*d2ydu2.*vv.*B1.^2+dhdu2.*((-2).*dydu2.*((-1)+vv).*(dBdu1.^2+(-1).*d2Bdu1.*B1)+2.*B1.*(d2Bdu1_du2.*dydu1+d2ydu1_du2.*dBdu1.*vv+d3yd2u1_du2.*B1)+(-1).*dBdu2.*(3.*dBdu1.*dydu1+d2ydu1.*vv.*B1))).*h+h.^2.*(3.*dBdu1.*dBdu2.^2.*dydu1+(-2).*d2Adu2.*dAdu2.*dydu2.*B1+dAdu2.^2.*(3.*dBdu2.*dydu2+(-1).*d2ydu2.*B1)+(-1).*dBdu2.*(B1.*((-2).*d2ydu1_du2.*dBdu1+3.*d2Bdu1_du2.*dydu1+2.*d3yd2u1_du2.*B1)+dydu2.*(9.*dBdu1.^2+d2Bdu1.*((-5)+2.*vv).*B1))+B1.*((-1).*d2Bdu2.*dBdu1.*dydu1+5.*d2Bdu1_du2.*dBdu1.*dydu2+(-2).*d3ydu1_2du2.*dBdu1.*B1+d3Bdu1_2du2.*dydu1.*B1+(-2).*d3Bd2u1_du2.*dydu2.*B1+d3Bd2u1_du2.*dydu2.*vv.*B1+d2ydu2.*(4.*dBdu1.^2+d2Bdu1.*((-3)+vv).*B1)+2.*B1.^2.*d4yd2u1_2du2)))+(-1).*A1.^3.*B1.^4.*(6.*dhdu1.^2.*B1.^2.*(dBdu1.*dydu1.*vv+d2ydu1.*B1)+3.*B1.*((-1).*dBdu1.^2.*dhdu1.*dydu1+dAdu1.*dAdu2.*(5.*dhdu2.*dydu1+(-3).*dhdu1.*dydu2)+(-2).*dAdu2.^2.*dhdu1.*dydu1.*((-1)+vv)+dBdu1.*(d2hdu1.*dydu1.*vv+d2ydu1.*dhdu1.*(2+vv)).*B1+B1.*(d2hdu1.*d2ydu1.*B1+2.*dhdu1.*(d2Bdu1.*dydu1.*vv+d3ydu1.*B1))).*h+h.^2.*(dBdu1.^3.*dydu1+(-1).*(dBdu1.*(d2ydu1.*dBdu1+2.*d2Bdu1.*dydu1)+3.*d2Adu1_du2.*dAdu1.*dydu2+d2Adu2.*dAdu1.*dydu1.*((-5)+2.*vv)).*B1+(2.*d3ydu1.*dBdu1+d3Bdu1.*dydu1.*vv+d2Bdu1.*d2ydu1.*(1+vv)).*B1.^2+4.*dAdu2.^2.*(dBdu1.*dydu1+d2ydu1.*B1)+dAdu2.*(dAdu1.*dBdu2.*dydu1.*((-5)+2.*vv)+(2.*d2ydu1_du2.*dAdu1+5.*d2Adu1_du2.*dydu1).*B1+dydu2.*(dAdu1.*dBdu1.*(3+2.*vv)+(-1).*d2Adu1.*B1))+B1.^3.*d4ydu1)));
C_77	=	(1/12).*R1.^(-4).*R2.^(-1).*A1.^(-4).*B1.^(-4).*(3.*dAdu1.^2.*R1.^2.*R2.*B1.^4.*((-1).*dBdu1.*R1.*vv+dR1du1.*B1).*h.^2+(-1).*R1.*R2.*A1.*B1.^3.*h.*(9.*dAdu1.*dhdu1.*R1.*B1.*((-1).*dBdu1.*R1.*vv+dR1du1.*B1)+((-6).*dAdu1.*dR1du1.^2.*B1.^2+R1.*B1.*(d2Adu1.*dR1du1.*B1+dAdu1.*(dBdu1.*dR1du1.*(3+vv)+3.*d2R1du1.*B1))+R1.^2.*((-1).*d2Adu1.*dBdu1.*vv.*B1+dAdu1.*(2.*dBdu1.^2+(-3).*d2Bdu1.*vv.*B1))).*h)+R2.*A1.^2.*B1.^2.*(6.*dhdu1.^2.*R1.^2.*B1.^2.*((-1).*dBdu1.*R1.*vv+dR1du1.*B1)+3.*R1.*B1.*((-4).*dhdu1.*dR1du1.^2.*B1.^2+R1.*B1.*(dBdu1.*dhdu1.*dR1du1.*(2+vv)+(2.*d2R1du1.*dhdu1+d2hdu1.*dR1du1).*B1)+R1.^2.*(dBdu1.^2.*dhdu1+2.*dAdu1.*dAdu2.*dhdu2.*((-1)+vv)+(-1).*d2hdu1.*dBdu1.*vv.*B1+(-2).*d2Bdu1.*dhdu1.*vv.*B1)).*h+(6.*dR1du1.^3.*B1.^3+(-2).*dR1du1.*R1.*B1.^2.*(2.*dBdu1.*dR1du1+3.*d2R1du1.*B1)+R1.^3.*((-1).*dBdu1.^3+dAdu2.^2.*dBdu1.*((-2)+vv)+(-2).*dAdu1.*dAdu2.*dBdu2.*((-1)+vv)+2.*(d2Bdu1.*dBdu1+d2Adu2.*dAdu1.*((-1)+vv)).*B1+(-1).*d3Bdu1.*vv.*B1.^2)+R1.^2.*B1.*(dAdu2.^2.*dR1du1+(-1).*dBdu1.^2.*dR1du1+(-2).*dAdu1.*dAdu2.*dR1du2.*((-1)+vv)+2.*d2R1du1.*dBdu1.*B1+B1.*(d2Bdu1.*dR1du1.*(1+vv)+d3R1du1.*B1))).*h.^2)+R1.*R2.*A1.^3.*B1.*((-12).*dAdu2.*dhdu1.*dhdu2.*R1.^2.*((-1)+vv).*B1.^2+3.*R1.*B1.*(dAdu2.*(2.*dhdu1.*dR1du2.*((-1)+vv)+dhdu2.*dR1du1.*((-3)+2.*vv)).*B1+R1.*((-2).*(d2Adu2.*dhdu1+d2Adu1_du2.*dhdu2).*((-1)+vv).*B1+dAdu2.*(dBdu1.*dhdu2.*vv+2.*((-1)+vv).*(dBdu2.*dhdu1+(-1).*d2hdu1_du2.*B1)))).*h+(2.*dAdu2.*dR1du1.*dR1du2.*(3+(-2).*vv).*B1.^2+R1.*B1.*((2.*d2Adu1_du2.*dR1du2.*((-1)+vv)+d2Adu2.*dR1du1.*((-3)+2.*vv)).*B1+(-1).*dAdu2.*(dBdu1.*dR1du2.*vv+((-3)+2.*vv).*(dBdu2.*dR1du1+(-1).*d2R1du1_du2.*B1)))+R1.^2.*(dAdu2.*((-2).*dBdu1.*dBdu2.*vv+d2Bdu1_du2.*vv.*B1)+B1.*(d2Adu2.*dBdu1.*vv+2.*((-1)+vv).*(d2Adu1_du2.*dBdu2+(-1).*d3Adu1_2du2.*B1)))).*h.^2)+A1.^4.*(12.*dBdu1.*R1.^3.*(R1+R2.*vv).*B1.^4+(-3).*dBdu1.*dBdu2.^2.*R1.^3.*R2.*h.^2+R1.^2.*R2.*B1.*h.*((-3).*dBdu1.*dBdu2.*dR1du2.*h+R1.*(d2Bdu2.*dBdu1.*h+3.*dBdu2.*(3.*dBdu1.*dhdu2+d2Bdu1_du2.*h)))+(-1).*R1.*R2.*B1.^2.*(2.*dR1du2.*(dBdu1.*dR1du2+(-1).*dBdu2.*dR1du1.*vv).*h.^2+R1.*h.*(dhdu2.*((-6).*dBdu1.*dR1du2+3.*dBdu2.*dR1du1.*vv)+((-1).*d2R1du2.*dBdu1+(-2).*d2Bdu1_du2.*dR1du2+d2R1du1_du2.*dBdu2.*vv).*h)+R1.^2.*(6.*dBdu1.*dhdu2.^2+6.*d2Bdu1_du2.*dhdu2.*h+h.*(3.*d2hdu2.*dBdu1+d3Bdu1_2du2.*h)))+R2.*vv.*B1.^3.*(6.*dR1du1.*dR1du2.^2.*h.^2+(-2).*R1.*h.*(6.*dhdu2.*dR1du1.*dR1du2+(d2R1du2.*dR1du1+2.*d2R1du1_du2.*dR1du2).*h)+R1.^2.*(6.*dhdu2.^2.*dR1du1+6.*d2R1du1_du2.*dhdu2.*h+h.*(3.*d2hdu2.*dR1du1+d3R1du1_2du2.*h)))));
C_78	=	(1/12).*R1.^(-1).*R2.^(-4).*A1.^(-4).*B1.^(-4).*((-3).*dAdu1.^2.*dAdu2.*R1.*R2.^3.*B1.^4.*h.^2+R1.*R2.^2.*A1.*B1.^3.*h.*(9.*dAdu1.*dAdu2.*dhdu1.*R2.*B1+((-3).*dAdu1.*dAdu2.*dR2du1.*B1+R2.*(3.*d2Adu1_du2.*dAdu1.*B1+dAdu2.*((-2).*dAdu1.*dBdu1.*vv+d2Adu1.*B1))).*h)+R1.*A1.^5.*(6.*dhdu2.^2.*dR2du2.*R2.^2.*B1.^2+3.*R2.*B1.*((-4).*dhdu2.*dR2du2.^2.*B1+R2.*((-3).*dBdu2.*dhdu2.*dR2du2+(2.*d2R2du2.*dhdu2+d2hdu2.*dR2du2).*B1)).*h+(6.*dR2du2.^3.*B1.^2+6.*dR2du2.*R2.*B1.*(dBdu2.*dR2du2+(-1).*d2R2du2.*B1)+R2.^2.*(3.*dBdu2.^2.*dR2du2+(-3).*d2R2du2.*dBdu2.*B1+(-1).*d2Bdu2.*dR2du2.*B1+d3R2du2.*B1.^2)).*h.^2)+(-1).*R1.*R2.*A1.^2.*B1.^2.*(6.*dAdu2.*dhdu1.^2.*R2.^2.*B1.^2+3.*R2.*B1.*(dhdu1.*((-2).*dAdu2.*dR2du1+dAdu1.*dR2du2.*vv).*B1+R2.*((-2).*dAdu1.*dBdu1.*dhdu2.*((-1)+vv)+2.*d2Adu1_du2.*dhdu1.*B1+dAdu2.*((-1).*dBdu1.*dhdu1.*vv+d2hdu1.*B1))).*h+(2.*dR2du1.*(dAdu2.*dR2du1+(-1).*dAdu1.*dR2du2.*vv).*B1.^2+R2.*B1.*(dAdu1.*dBdu1.*dR2du2.*((-3)+2.*vv)+((-2).*d2Adu1_du2.*dR2du1+d2R2du1_du2.*dAdu1.*vv).*B1+dAdu2.*(dBdu1.*dR2du1.*vv+(-1).*d2R2du1.*B1))+R2.^2.*(dAdu2.^3+2.*dAdu1.*dBdu1.*dBdu2.*((-1)+vv)+B1.*((-2).*d2Bdu1_du2.*dAdu1.*((-1)+vv)+(-1).*d2Adu1_du2.*dBdu1.*vv+d3Ad2u1_du2.*B1)+(-1).*dAdu2.*(dBdu1.^2.*((-2)+vv)+d2Bdu1.*vv.*B1))).*h.^2)+R1.*A1.^3.*B1.*(6.*dhdu1.*R2.^2.*B1.^2.*((-2).*dBdu1.*dhdu2.*R2.*((-1)+vv)+dhdu1.*dR2du2.*vv.*B1)+3.*R2.*B1.*((-4).*dhdu1.*dR2du1.*dR2du2.*vv.*B1.^2+R2.^2.*(dAdu2.^2.*dhdu2+2.*((-1)+vv).*(dBdu1.*dBdu2.*dhdu1+(-1).*(d2hdu1_du2.*dBdu1+d2Bdu1_du2.*dhdu1+d2Bdu1.*dhdu2).*B1))+R2.*B1.*(dR2du2.*(dBdu1.*dhdu1.*((-3)+2.*vv)+d2hdu1.*vv.*B1)+2.*(dBdu1.*dhdu2.*dR2du1.*((-1)+vv)+d2R2du1_du2.*dhdu1.*vv.*B1))).*h+(6.*dR2du1.^2.*dR2du2.*vv.*B1.^3+(-2).*R2.^3.*(dAdu2.^2.*dBdu2+(-1).*d2Adu2.*dAdu2.*B1+((-1)+vv).*B1.*((-1).*d2Bdu1.*dBdu2+d3Bd2u1_du2.*B1))+(-2).*R2.*B1.^2.*(2.*d2R2du1_du2.*dR2du1.*vv.*B1+dR2du2.*(dBdu1.*dR2du1.*((-3)+2.*vv)+d2R2du1.*vv.*B1))+R2.^2.*B1.*(2.*dBdu1.*dBdu2.*dR2du1+(-1).*dAdu2.^2.*dR2du2+(-2).*dBdu1.*dBdu2.*dR2du1.*vv+(-3).*d2R2du1_du2.*dBdu1.*B1+(-2).*d2Bdu1_du2.*dR2du1.*B1+2.*d2R2du1_du2.*dBdu1.*vv.*B1+2.*d2Bdu1_du2.*dR2du1.*vv.*B1+d3R2d2u1_du2.*vv.*B1.^2+dR2du2.*(dBdu1.^2+d2Bdu1.*((-3)+2.*vv).*B1))).*h.^2)+R2.*A1.^4.*(12.*dAdu2.*R2.^2.*(R2+R1.*vv).*B1.^4+(-3).*dAdu2.*dBdu2.^2.*R1.*R2.^2.*vv.*h.^2+R1.*R2.*B1.*h.*((-1).*dAdu2.*dBdu2.*dR2du2.*(3+vv).*h+R2.*vv.*(3.*d2Adu2.*dBdu2.*h+dAdu2.*(9.*dBdu2.*dhdu2+d2Bdu2.*h)))+R1.*B1.^2.*((-4).*dAdu2.*dR2du2.^2.*h.^2+R2.*h.*(d2Adu2.*dR2du2.*(1+vv).*h+dAdu2.*(3.*dhdu2.*dR2du2.*(2+vv)+2.*d2R2du2.*h))+(-1).*R2.^2.*vv.*(3.*dAdu2.*(2.*dhdu2.^2+d2hdu2.*h)+h.*(6.*d2Adu2.*dhdu2+d3Adu2.*h)))));
C_79	=	(-1).*R1.^(-2).*R2.^(-2).*(R1.^2+R2.^2+2.*R1.*R2.*vv).*A1.*B1;
C_80	=	(-1).*E.^(-1).*h.^(-1).*P.*((-1)+vv.^2).*A1.*B1;


%% Finding boundary
clear nj nbc;

for j=1:m+1
    for i=1:n+1
        ex(i,j) = x(i,j)-sqrt((1-y(i,j).^2/b^2)*a^2);
        ey(i,j) = y(i,j)-sqrt((1-x(i,j).^2/a^2)*b^2);
            if abs(ex(i,j)) < dx/2 || abs(ey(i,j)) < dy/2
                nj(j) = i;
            end
    end
end

%nj(1:m+1) = n+1;                % On for rectangular boundary

c=1;
for j=1:m+1
    for i=nj(j):-1:1
        if j<m+1
            if  i+1 > nj(j+1) 
                nbc(j) = i;                               % -1 for only 0 boundary condition
                xbc1(c) = x(i,j);
                ybc1(c) = y(i,j);
                c=c+1;
            end
        elseif j==m+1
                nbc(j) = 1; 
                xbc1(c) = x(i,j);
                ybc1(c) = y(i,j);
                c=c+1;
        end
    end
end

nbcc = cumsum(nbc-1);
nbcc = [0 nbcc];
nbcc(end) = nbcc(end)+5;

A=zeros(nbcc(end)*3);

%% Creating matrixes
for j=1:m
    for i=1:nbc(j)-1
        
                                                 %% corner
        if j == 1 && i==1                                               
            
                % 1st equation

                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+1 ) = 1;
                
                B( (i+nbcc(j)-1)*3+1 ) = 0;
                
                
                % 2nd equation

                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+2 ) = 1;
                
                B( (i+nbcc(j)-1)*3+2 ) = 0;
                
                % 3rd equation

                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+1 ) = 	(-1).*dx.^(-3).*C_50(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(1/2).*dx.^(-1).*C_71(i,j)...
                                                                     -( dx.^(-3).*C_50(i,j)+dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(-1/2).*dx.^(-1).*C_71(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+1 ) = 	(1/2).*dx.^(-3).*C_50(i,j)...
                                                                    -( (-1/2).*dx.^(-3).*C_50(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+1 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1).*dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(1/2).*dy.^(-1).*C_72(i,j)...
                                                                     + dx.^(-2).*dy.^(-1).*C_51(i,j)+dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(-1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    -( (-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ) ...
                                                                    + (-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    - ( (1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ) ;
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+1 ) = 	(1/2).*dy.^(-3).*C_53(i,j)...
                                                                    +(-1/2).*dy.^(-3).*C_53(i,j);
                
                
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+2 ) = 	(-1).*dx.^(-3).*C_54(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(1/2).*dx.^(-1).*C_73(i,j) ...
                                                                    + dx.^(-3).*C_54(i,j)+dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(-1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+2 ) = 	(1/2).*dx.^(-3).*C_54(i,j) ...
                                                                    + (-1/2).*dx.^(-3).*C_54(i,j) ;
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+2 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1).*dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(1/2).*dy.^(-1).*C_74(i,j) ...
                                                                    -( dx.^(-2).*dy.^(-1).*C_55(i,j)+dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(-1/2).*dy.^(-1).*C_74(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ...
                                                                    -( (-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ) ...
                                                                    -( (-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ) ...
                                                                    + (1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+2 ) = 	(1/2).*dy.^(-3).*C_57(i,j) ...
                                                                     - ( (-1/2).*dy.^(-3).*C_57(i,j) );

                
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1).*dx.^(-3).*C_58(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(1/2).*dx.^(-1).*C_75(i,j) ...
                                                                    + (-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+dx.^(-3).*C_58(i,j)+dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(-1/2).*dx.^(-1).*C_75(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(1/2).*dx.^(-3).*C_58(i,j) ...
                                                                    + dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+(-1).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1).*dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(1/2).*dy.^(-1).*C_76(i,j) ...
                                                                    + (-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+dx.^(-2).*dy.^(-1).*C_59(i,j)+dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(-1/2).*dy.^(-1).*C_76(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ;
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3 ) = 	(1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(1/2).*dy.^(-3).*C_61(i,j) ...
                                                                    + dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3 ) = 	(1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                
                B( (i+nbcc(j)-1)*3+3 ) = - C_80(i,j);
            
        
        
        elseif j==1 &&  i > 1 % && i+1 < nbc(j) && i < nbc(j+1)                                           % Lower boundary
            
                % 1st equation

                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-2).*C_05(i,j)+(-1/2).*  dx.^(-1).*C_14(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_05(i,j)+(-2).*  dy.^(-2).*C_07(i,j)+C_20(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-2).*C_05(i,j)+  (1/2).*dx.^(-1).*C_14(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+1 ) = 	(-1/4).*dx.^(-1).*dy.^(-1).*  C_06(i,j) ...
                                                                    + (1/4).*dx.^(-1).*dy.^(-1).*C_06(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+1 ) = 	dy.^(-2).*C_07(i,j)+(1/2).*dy.^(-1).*C_15(i,j) ...
                                                                    + dy.^(-2).*C_07(i,j)+(-1/2).*dy.^(-1).*C_15(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_06(i,j) ...
                                                                    + (-1/4).*dx.^(-1)  .*dy.^(-1).*C_06(i,j);

                
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-2).*C_08(i,j)+(-1/2).* dx.^(-1).*C_16(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_08(i,j)+(-2).*  dy.^(-2).*C_10(i,j)+C_21(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-2).*C_08(i,j)  +(1/2).*dx.^(-1).*C_16(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+2 ) = 	(-1/4).*dx.^(-1).*dy.^(-1).*  C_09(i,j) ...
                                                                    - ( (1/4).*dx.^(-1).*dy.^(-1).*C_09(i,j) );
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+2 ) = 	dy.^(-2).*C_10(i,j)+(1/2).*dy.^(-1).*C_17(i,j) ...
                                                                    - ( dy.^(-2).*C_10(i,j)+(-1/2).*dy.^(-1).*C_17(i,j) );
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_09(i,j) ...
                                                                    - ( (-1/4).*dx.^(-1)  .*dy.^(-1).*C_09(i,j) );

                
                if i==2
                    A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j) ...
                                                                        + (-1/2).*dx.^(-3).*C_01(i,j);
                else
                    A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j);
                    A( (i+nbcc(j)-1)*3+1 , (i-2 + nbcc(j  )-1)*3+3 ) = 	(-1/2).*dx.^(-3).*C_01(i,j);
                end
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+3 ) = 	dx.^(-3).*C_01(i,j)+dx.^(-1).*dy.^(-2).*C_03(i,j)+  dx.^(-2).*C_11(i,j)+(-1/2).*dx.^(-1).*C_18(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-1).*dx.^(-3).*C_01(i,j)+(-1).*dx.^(-1).*  dy.^(-2).*C_03(i,j)+dx.^(-2).*C_11(i,j)+(1/2).*dx.^(-1).*C_18(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+2 + nbcc(j  )-1)*3+3 ) = 	(1/2).*dx.^(-3).*C_01(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*dx.^(  -2).*dy.^(-1).*C_02(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_12(i,j) ...
                                                                    + (-1/2).*dx.^(-2).*dy.^(-1)  .*C_02(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(1/4)  .*dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-1).*  dx.^(-2).*dy.^(-1).*C_02(i,j)+(-1).*dy.^(-3).*C_04(i,j)+  dy.^(-2).*C_13(i,j)+(1/2).*dy.^(-1).*C_19(i,j) ...
                                                                    + dx.^(-2).*dy.^(-1).*C_02(i,j)+dy.^(-3).*C_04(i,j)+dy.^(-2).*C_13(i,j)+  (-1/2).*dy.^(-1).*C_19(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*  dx.^(-2).*dy.^(-1).*C_02(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_12(i,j) ...
                                                                    + (-1/2).*dx.^(-2).*dy.^(-1).*C_02(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(-1/4).*  dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+2)-1)*3+3 ) = 	(1/2).*dy.^(-3).*C_04(i,j) ...
                                                                    + (-1/2).*dy.^(-3).*C_04(i,j);

                B( (i+nbcc(j)-1)*3+1 ) = 0;
                
                
                % 2nd equation
              
                
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+2 ) = 		1;
                                
                B( (i+nbcc(j)-1)*3+2 ) = 0;
                
                
                % 3rd equation

                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j) ...
                                                                         - ( (-1/2).*dx.^(-3).*C_50(i,j) );
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j) ;
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+1 ) = 	(-1/2).*dx.^(-3).*C_50(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-3).*C_50(i,j)+dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(-1/2).*dx.^(-1).*C_71(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+1 ) = 	(-1).*dx.^(-3).*C_50(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(1/2).*dx.^(-1).*C_71(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+1 ) = 	(1/2).*dx.^(-3).*C_50(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    + 	(-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+1 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1).*dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(1/2).*dy.^(-1).*C_72(i,j) ...
                                                                    + 	dx.^(-2).*dy.^(-1).*C_51(i,j)+dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(-1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    + (-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+1 ) = 	(1/2).*dy.^(-3).*C_53(i,j) ...
                                                                    + 	(-1/2).*dy.^(-3).*C_53(i,j);
                
                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j) ...
                                                                        + (-1/2).*dx.^(-3).*C_54(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+2 ) = 	(-1/2).*dx.^(-3).*C_54(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-3).*C_54(i,j)+dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(-1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+2 ) = 	(-1).*dx.^(-3).*C_54(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+2 ) = 	(1/2).*dx.^(-3).*C_54(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ...
                                                                    - ( (-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j));
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+2 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1).*dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(1/2).*dy.^(-1).*C_74(i,j) ...
                                                                    - ( dx.^(-2).*dy.^(-1).*C_55(i,j)+dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(-1/2).*dy.^(-1).*C_74(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ...
                                                                    - ( (-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+2 ) = 	(1/2).*dy.^(-3).*C_57(i,j) ...
                                                                    - ( (-1/2).*dy.^(-3).*C_57(i,j) );

                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j) ...
                                                                        + dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+(-1).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1).*dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(1/2).*dy.^(-1).*C_76(i,j) ...
                                                                        + (-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+dx.^(-2).*dy.^(-1).*C_59(i,j)+dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(-1/2).*dy.^(-1).*C_76(i,j) ...
                                                                        + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                        + (1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+(-1).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1).*dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(1/2).*dy.^(-1).*C_76(i,j) ...
                                                                        + (-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+dx.^(-2).*dy.^(-1).*C_59(i,j)+dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(-1/2).*dy.^(-1).*C_76(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j+1)-1)*3+3 ) = 	(-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                        + (1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+dx.^(-3).*C_58(i,j)+dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(-1/2).*dx.^(-1).*C_75(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1).*dx.^(-3).*C_58(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(1/2).*dx.^(-1).*C_75(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(1/2).*dx.^(-3).*C_58(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3 ) = 	(1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3 ) = 	(-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(1/2).*dy.^(-3).*C_61(i,j) ...
                                                                    + dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3 ) = 	(1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);

                
                B( (i+nbcc(j)-1)*3+3 ) = -C_80(i,j);
                

        elseif j > 1  && i == 1 %&& i+1 < nbc(j) && i < nbc(j+1)                                   % Left boundary
                    
                % 1st equation
                
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+1 ) = 1;
                
                B( (i+nbcc(j)-1)*3+1 ) = 0;
                
                % 2nd equation
                                
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+1 ) = 		dy.^(-2).*C_29(i,j)+(-1/2).*dy.^(-1).*C_37(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+1 ) = 		(-1/4).*dx.^(-1)  .*dy.^(-1).*C_28(i,j) ...
                                                                        - ( (1/4).*dx.^(-1).*dy.^(-1).*C_28(i,j) );
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+1 ) = 		(-2).*dx.^(-2).*C_27(i,j)+(-2).*  dy.^(-2).*C_29(i,j)+C_42(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+1 ) = 		dx.^(-2).*C_27(i,j)+(1/2).*dx.^(-1).*C_36(i,j) ...
                                                                        - ( dx.^(-2).*C_27(i,j)+(-1/2).*  dx.^(-1).*C_36(i,j) );
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+1 ) = 		dy.^(-2).*C_29(i,j)+(1/2).*dy.^(-1).*  C_37(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+1 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_28(i,j) ...
                                                                        - ( (-1/4).*dx.^(-1).*dy.^(  -1).*C_28(i,j) );
                
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+2 ) = 		dy.^(-2).* C_32(i,j)+(-1/2).*dy.^(-1).*C_39(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+2 ) = 		(-1/4).*dx.^(-1)  .*dy.^(-1).*C_31(i,j) ...
                                                                        + (1/4).*dx.^(-1).*dy.^(-1).*C_31(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+2 ) = 		(-2).*dx.^(-2).*C_30(i,j)+(-2).*  dy.^(-2).*C_32(i,j)+C_43(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+2 ) = 		dx.^(-2).*C_30(i,j)+(1/2).*dx.^(-1).*C_38(i,j) ...
                                                                        + dx.^(-2).*C_30(i,j)+(-1/2).*  dx.^(-1).*C_38(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+2 ) = 		dy.^(-2).*C_32(i,j)+(1/2).*dy.^(-1).*  C_39(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+2 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_31(i,j) ...
                                                                        + (-1/4).*dx.^(-1).*dy.^(  -1).*C_31(i,j);
                
                if j==2
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) = 		(-2).*dx.^(-2).*C_33(i,j)+(-2).*dy.^(-2).*  C_35(i,j)+C_44(i,j) ...
                                                                            + (-1/2).*dy.^(-3).*C_26(i,j);
                else
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) = 		(-2).*dx.^(-2).*C_33(i,j)+(-2).*dy.^(-2).*  C_35(i,j)+C_44(i,j);
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-2)-1)*3+3 ) = 		(-1/2).*dy.^(-3).*C_26(i,j);
                end
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+3 ) = 		dx.^(-2).*dy.^(-1).*C_24(i,j)+dy.^(-3).*C_26(i,j)+dy.^(-2).*C_35(i,j)+(-1/2).*dy.^(-1).*C_41(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+3 ) = 		(-1/2).*dx.^(-2).*dy.^(-1).*  C_24(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_25(i,j)+(-1/4)  .*dx.^(-1).*dy.^(-1).*C_34(i,j) ...
                                                                        + (-1/2).*dx.^(-2).*dy.^(  -1).*C_24(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_25(i,j)+(  1/4).*dx.^(-1).*dy.^(-1).*C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+3 ) = 		(-1).*dx.^(-3).*C_23(i,j)+(  -1).*dx.^(-1).*dy.^(-2).*C_25(i,j)+dx.^(-2).*C_33(i,j)  +(1/2).*dx.^(-1).*C_40(i,j) ...
                                                                        + dx.^(-3).*C_23(i,j)+dx.^(-1).*dy.^(-2).*C_25(i,j)+dx.^(-2).*C_33(i,j)+(-1/2).*dx.^(-1).*C_40(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+2 + nbcc(j  )-1)*3+3 ) = 		(1/2).*dx.^(-3).*C_23(i,j) ...
                                                                        + (-1/2).*dx.^(-3).*C_23(i,j);   
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+3 ) = 		(-1).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(-1).*dy.^(-3).*C_26(i,j)+dy.^(-2).*C_35(i,j)+(1/2).*dy.^(-1).*C_41(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+3 ) = 		(1/2).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(1/2).*  dx.^(-1).*dy.^(-2).*C_25(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*  C_34(i,j) ...
                                                                        + (1/2).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(-1/2).*dx.^(-1).*  dy.^(-2).*C_25(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+2)-1)*3+3 ) = 		(1/2).*dy.^(-3).*C_26(i,j);

                B( (i+nbcc(j)-1)*3+2 ) = 0;
                
                
                % 3rd equation

                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j)...
                                                                        + (-1/2).*dy.^(-3).*C_53(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+1 ) = 	(-1/2).*dy.^(-3).*C_53(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+1 ) = 	dx.^(-2).*dy.^(-1).*C_51(i,j)+dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(-1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+1 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    - ( (-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+1 ) = 	(-1).*dx.^(-3).*C_50(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(1/2).*dx.^(-1).*C_71(i,j) ...
                                                                    - ( dx.^(-3).*C_50(i,j)+dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(-1/2).*dx.^(-1).*C_71(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+1 ) = 	(1/2).*dx.^(-3).*C_50(i,j) ...
                                                                    - ( (-1/2).*dx.^(-3).*C_50(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+1 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1).*dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) ...
                                                                    - ( (1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j) );
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+1 ) = 	(1/2).*dy.^(-3).*C_53(i,j);
                
                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j) ...
                                                                        - ( (-1/2).*dy.^(-3).*C_57(i,j) );
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+2 ) = 	(-1/2).*dy.^(-3).*C_57(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+2 ) = 	dx.^(-2).*dy.^(-1).*C_55(i,j)+dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(-1/2).*dy.^(-1).*C_74(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+2 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ...
                                                                    + (-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+2 ) = 	(-1).*dx.^(-3).*C_54(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(1/2).*dx.^(-1).*C_73(i,j) ...
                                                                    + dx.^(-3).*C_54(i,j)+dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(-1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+2 ) = 	(1/2).*dx.^(-3).*C_54(i,j) ... 
                                                                    + (-1/2).*dx.^(-3).*C_54(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+2 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1).*dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(1/2).*dy.^(-1).*C_74(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j) ...
                                                                    + (1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+2 ) = 	(1/2).*dy.^(-3).*C_57(i,j);
                
                
                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j) ...
                                                                        + dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1).*dx.^(-3).*C_58(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(1/2).*dx.^(-1).*C_75(i,j) ...
                                                                        + (-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+dx.^(-3).*C_58(i,j)+dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(-1/2).*dx.^(-1).*C_75(i,j) ...
                                                                        + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                        + (1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1).*dx.^(-3).*C_58(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(1/2).*dx.^(-1).*C_75(i,j) ...
                                                                        + (-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+dx.^(-3).*C_58(i,j)+dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(-1/2).*dx.^(-1).*C_75(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-2)-1)*3+3 ) = 	(-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                        + (1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+dx.^(-2).*dy.^(-1).*C_59(i,j)+dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(-1/2).*dy.^(-1).*C_76(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+3 ) = 	(1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3 ) = 	(-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(1/2).*dx.^(-3).*C_58(i,j) ...
                                                                    + dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+(-1).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1).*dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(1/2).*dy.^(-1).*C_76(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j) ...
                                                                    + (1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3 ) = 	(1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j) ...
                                                                    + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(1/2).*dy.^(-3).*C_61(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3 ) = 	(1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j) ...
                                                                    + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                                                                
                B( (i+nbcc(j)-1)*3+3 ) = -C_80(i,j);
                

        else                        %% middle region                                                 
                % 1st equation
                
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j-1)-1)*3+1 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_06(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j-1)-1)*3+1 ) = 	dy.^(-2).*C_07(i,j)+(-1/2).*dy.^(-1).*C_15(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+1 ) = 	(-1/4).*dx.^(-1)  .*dy.^(-1).*C_06(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-2).*C_05(i,j)+(-1/2).*  dx.^(-1).*C_14(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_05(i,j)+(-2).*  dy.^(-2).*C_07(i,j)+C_20(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-2).*C_05(i,j)+  (1/2).*dx.^(-1).*C_14(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+1 ) = 	(-1/4).*dx.^(-1).*dy.^(-1).*  C_06(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+1 ) = 	dy.^(-2).*C_07(i,j)+(1/2).*dy.^(-1).*C_15(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_06(i,j);
                
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j-1)-1)*3+2 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_09(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j-1)-1)*3+2 ) = 	dy.^(-2).*C_10(i,j)+(-1/2).*dy.^(-1).*C_17(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+2 ) = 	(-1/4).*dx.^(-1)  .*dy.^(-1).*C_09(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-2).*C_08(i,j)+(-1/2).* dx.^(-1).*C_16(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_08(i,j)+(-2).*  dy.^(-2).*C_10(i,j)+C_21(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-2).*C_08(i,j)  +(1/2).*dx.^(-1).*C_16(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+2 ) = 	(-1/4).*dx.^(-1).*dy.^(-1).*  C_09(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+2 ) = 	dy.^(-2).*C_10(i,j)+(1/2).*dy.^(-1).*C_17(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/4).*dx.^(-1).*dy.^(-1).*C_09(i,j);
                
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j);
                if i==2
                    A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 )...
                                                                            + (-1/2).*dx.^(-3).*C_01(i,j);
                else
                    %A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j);
                    A( (i+nbcc(j)-1)*3+1 , (i-2 + nbcc(j  )-1)*3+3 ) = 	(-1/2).*dx.^(-3).*C_01(i,j);
                end
                if j==2
                    A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 )...
                                                                        + (-1/2).*dy.^(-3).*C_04(i,j);
                else
                    %A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j);
                    A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j-2)-1)*3+3 ) = 	(-1/2).*dy.^(-3).*C_04(i,j);
                end
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j-1)-1)*3+3 ) = 	(-1/2).*dx.^(-2).*dy.^(-1)  .*C_02(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(1/4)  .*dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j-1)-1)*3+3 ) = 	dx.^(-2).*dy.^(-1).*C_02(i,j)+dy.^(-3).*C_04(i,j)+dy.^(-2).*C_13(i,j)+  (-1/2).*dy.^(-1).*C_19(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+3 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_02(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(-1/4).*  dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j  )-1)*3+3 ) = 	dx.^(-3).*C_01(i,j)+dx.^(-1).*dy.^(-2).*C_03(i,j)+  dx.^(-2).*C_11(i,j)+(-1/2).*dx.^(-1).*C_18(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j  )-1)*3+3 ) = 	(-2).*  dx.^(-2).*C_11(i,j)+(-2).*dy.^(-2).*C_13(i,j)+C_22(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-1).*dx.^(-3).*C_01(i,j)+(-1).*dx.^(-1).*  dy.^(-2).*C_03(i,j)+dx.^(-2).*C_11(i,j)+(1/2).*dx.^(-1).*C_18(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+2 + nbcc(j  )-1)*3+3 ) = 	(1/2).*dx.^(-3).*C_01(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*dx.^(  -2).*dy.^(-1).*C_02(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-1).*  dx.^(-2).*dy.^(-1).*C_02(i,j)+(-1).*dy.^(-3).*C_04(i,j)+  dy.^(-2).*C_13(i,j)+(1/2).*dy.^(-1).*C_19(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*  dx.^(-2).*dy.^(-1).*C_02(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_03(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_12(i,j);
                A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+2)-1)*3+3 ) = 	(1/2).*dy.^(-3).*C_04(i,j);

                B( (i+nbcc(j)-1)*3+1 ) = 0;
                
                %% 2nd equation
                
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j-1)-1)*3+1 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_28(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+1 ) = 		dy.^(-2).*C_29(i,j)+(-1/2).*dy.^(-1).*C_37(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+1 ) = 		(-1/4).*dx.^(-1)  .*dy.^(-1).*C_28(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j  )-1)*3+1 ) = 		dx.^(-2).*C_27(i,j)+(-1/2).*  dx.^(-1).*C_36(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+1 ) = 		(-2).*dx.^(-2).*C_27(i,j)+(-2).*  dy.^(-2).*C_29(i,j)+C_42(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+1 ) = 		dx.^(-2).*C_27(i,j)+(1/2).*dx.^(-1).*C_36(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+1 ) = 		(-1/4).*dx.^(-1).*dy.^(  -1).*C_28(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+1 ) = 		dy.^(-2).*C_29(i,j)+(1/2).*dy.^(-1).*  C_37(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+1 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_28(i,j);
                
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j-1)-1)*3+2 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_31(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+2 ) = 		dy.^(-2).* C_32(i,j)+(-1/2).*dy.^(-1).*C_39(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+2 ) = 		(-1/4).*dx.^(-1)  .*dy.^(-1).*C_31(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j  )-1)*3+2 ) = 		dx.^(-2).*C_30(i,j)+(-1/2).*  dx.^(-1).*C_38(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+2 ) = 		(-2).*dx.^(-2).*C_30(i,j)+(-2).*  dy.^(-2).*C_32(i,j)+C_43(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+2 ) = 		dx.^(-2).*C_30(i,j)+(1/2).*dx.^(-1).*C_38(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+2 ) = 		(-1/4).*dx.^(-1).*dy.^(  -1).*C_31(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+2 ) = 		dy.^(-2).*C_32(i,j)+(1/2).*dy.^(-1).*  C_39(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+2 ) = 		(1/4).*dx.^(-1).*dy.^(-1).*C_31(i,j);
                
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) = 		(-2).*dx.^(-2).*C_33(i,j)+(-2).*dy.^(-2).*  C_35(i,j)+C_44(i,j);
                if i==2
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) = 		A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) ...
                                                                            + (-1/2).*dx.^(-3).*C_23(i,j);
                else
                    A( (i+nbcc(j)-1)*3+2 , (i-2 + nbcc(j  )-1)*3+3 ) = 		(-1/2).*dx.^(-3).*C_23(i,j);
                end
                
                if j==2
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) =     A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j  )-1)*3+3 ) ...
                                                                            + (-1/2).*dy.^(-3).*C_26(i,j);
                else
                    A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-2)-1)*3+3 ) = 		(-1/2).*dy.^(-3).*C_26(i,j);
                end
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j-1)-1)*3+3 ) = 		(-1/2).*dx.^(-2).*dy.^(  -1).*C_24(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_25(i,j)+(  1/4).*dx.^(-1).*dy.^(-1).*C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j-1)-1)*3+3 ) = 		dx.^(-2).*dy.^(-1).*C_24(i,j)+dy.^(-3).*C_26(i,j)+dy.^(-2).*C_35(i,j)+(-1/2).*dy.^(-1).*C_41(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+3 ) = 		(-1/2).*dx.^(-2).*dy.^(-1).*  C_24(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_25(i,j)+(-1/4)  .*dx.^(-1).*dy.^(-1).*C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j  )-1)*3+3 ) = 		dx.^(-3).*C_23(i,j)+dx.^(-1).*dy.^(-2).*C_25(i,j)+dx.^(-2).*C_33(i,j)+(-1/2).*dx.^(-1).*C_40(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+3 ) = 		(-1).*dx.^(-3).*C_23(i,j)+(  -1).*dx.^(-1).*dy.^(-2).*C_25(i,j)+dx.^(-2).*C_33(i,j)  +(1/2).*dx.^(-1).*C_40(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+2 + nbcc(j  )-1)*3+3 ) = 		(1/2).*dx.^(-3).*C_23(i,j);  
                A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+3 ) = 		(1/2).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(-1/2).*dx.^(-1).*  dy.^(-2).*C_25(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+3 ) = 		(-1).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(-1).*dy.^(-3).*C_26(i,j)+dy.^(-2).*C_35(i,j)+(1/2).*dy.^(-1).*C_41(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+3 ) = 		(1/2).*dx.^(-2).*dy.^(-1).*C_24(i,j)+(1/2).*  dx.^(-1).*dy.^(-2).*C_25(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*  C_34(i,j);
                A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+2)-1)*3+3 ) = 		(1/2).*dy.^(-3).*C_26(i,j);

                B( (i+nbcc(j)-1)*3+2 ) = 0;
                
                %% 3rd equation
                
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	(-2).*dx.^(-2).*C_62(i,j)+(-2).*dy.^(-2).*C_64(i,j)+C_77(i,j);
                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) ...
                                                                        - ( (-1/2).*dx.^(-3).*C_50(i,j) );
                else
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+1 ) = 	(-1/2).*dx.^(-3).*C_50(i,j);
                end
                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+1 ) ...
                                                                        + (-1/2).*dy.^(-3).*C_53(i,j) ;
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+1 ) = 	(-1/2).*dy.^(-3).*C_53(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j-1)-1)*3+1 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+1 ) = 	dx.^(-2).*dy.^(-1).*C_51(i,j)+dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(-1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+1 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+1 ) = 	dx.^(-3).*C_50(i,j)+dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(-1/2).*dx.^(-1).*C_71(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+1 ) = 	(-1).*dx.^(-3).*C_50(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_52(i,j)+dx.^(-2).*C_62(i,j)+(1/2).*dx.^(-1).*C_71(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+1 ) = 	(1/2).*dx.^(-3).*C_50(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+1 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(-1).*dy.^(-3).*C_53(i,j)+dy.^(-2).*C_64(i,j)+(1/2).*dy.^(-1).*C_72(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+1 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_51(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_52(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_63(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+1 ) = 	(1/2).*dy.^(-3).*C_53(i,j);
                
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	(-2).*dx.^(-2).*C_65(i,j)+(-2).*dy.^(-2).*C_67(i,j)+C_78(i,j);
                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 )...
                                                                        + (-1/2).*dx.^(-3).*C_54(i,j) ;
                else
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+2 ) = 	(-1/2).*dx.^(-3).*C_54(i,j);
                end
                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+2 )...
                                                                        - ( (-1/2).*dy.^(-3).*C_57(i,j) ) ;                    
                else
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+2 ) = 	(-1/2).*dy.^(-3).*C_57(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j-1)-1)*3+2 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+2 ) = 	dx.^(-2).*dy.^(-1).*C_55(i,j)+dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(-1/2).*dy.^(-1).*C_74(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+2 ) = 	(-1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+2 ) = 	dx.^(-3).*C_54(i,j)+dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(-1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+2 ) = 	(-1).*dx.^(-3).*C_54(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_56(i,j)+dx.^(-2).*C_65(i,j)+(1/2).*dx.^(-1).*C_73(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+2 ) = 	(1/2).*dx.^(-3).*C_54(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+2 ) = 	(-1).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(-1).*dy.^(-3).*C_57(i,j)+dy.^(-2).*C_67(i,j)+(1/2).*dy.^(-1).*C_74(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+2 ) = 	(1/2).*dx.^(-2).*dy.^(-1).*C_55(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_56(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_66(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+2 ) = 	(1/2).*dy.^(-3).*C_57(i,j);
                
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+dx.^(-3).*C_58(i,j)+dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(-1/2).*dx.^(-1).*C_75(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	6.*dx.^(-4).*C_45(i,j)+4.*dx.^(-2).*dy.^(-2).*C_47(i,j)+6.*dy.^(-4).*C_49(i,j)+(-2).*dx.^(-2).*C_68(i,j)+(-2).*dy.^(-2).*C_70(i,j)+C_79(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	(-4).*dx.^(-4).*C_45(i,j)+(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1).*dx.^(-3).*C_58(i,j)+(-1).*dx.^(-1).*dy.^(-2).*C_60(i,j)+dx.^(-2).*C_68(i,j)+(1/2).*dx.^(-1).*C_75(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+dx.^(-2).*dy.^(-1).*C_59(i,j)+dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(-1/2).*dy.^(-1).*C_76(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	(-2).*dx.^(-2).*dy.^(-2).*C_47(i,j)+(-4).*dy.^(-4).*C_49(i,j)+(-1).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1).*dy.^(-3).*C_61(i,j)+dy.^(-2).*C_70(i,j)+(1/2).*dy.^(-1).*C_76(i,j);
                if i==2
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+3 )...
                                                                        + (1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) ...
                                                                        + dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3 ) ...
                                                                        + (-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j-1)-1)*3+3 ) = 	(1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(-1/2).*dx.^(-3).*C_58(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j+1)-1)*3+3 ) = 	(-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                end
                if j==2
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j  )-1)*3+3 ) ...
                                                                        + (1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j  )-1)*3+3 ) ...
                                                                        + dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 	A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 )...
                                                                        + (-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                else
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j-2)-1)*3+3 ) = 	(1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(-1/2).*dy.^(-3).*C_61(i,j);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-2)-1)*3+3 ) = 	(-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                end
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j-1)-1)*3+3 ) = 	(-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+3 ) = 	(1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(-1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3 ) = 	(-1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j  )-1)*3+3 ) = 	dx.^(-4).*C_45(i,j)+(1/2).*dx.^(-3).*C_58(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+3 ) = 	(1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(-1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(-1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+3 ) = 	(-1/2).*dx.^(-3).*dy.^(-1).*C_46(i,j)+dx.^(-2).*dy.^(-2).*C_47(i,j)+(-1/2).*dx.^(-1).*dy.^(-3).*C_48(i,j)+(1/2).*dx.^(-2).*dy.^(-1).*C_59(i,j)+(1/2).*dx.^(-1).*dy.^(-2).*C_60(i,j)+(1/4).*dx.^(-1).*dy.^(-1).*C_69(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3 ) = 	(1/4).*dx.^(-3).*dy.^(-1).*C_46(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3 ) = 	(-1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+2)-1)*3+3 ) = 	dy.^(-4).*C_49(i,j)+(1/2).*dy.^(-3).*C_61(i,j);
                A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3 ) = 	(1/4).*dx.^(-1).*dy.^(-3).*C_48(i,j);

                B( (i+nbcc(j)-1)*3+3 ) = -C_80(i,j);
                
        end
    end
end

for j=1:m
    for i=1:nbc(j)-1
        
        if j<m
            if i > 1
                if i-1 >= nbc(j+2) &&  i-1 <= nj(j+2)   % for on       (1)
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3) = 0;
                elseif i-1 > nj(j+2)                    % for outside (above)
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j)-1)*3+3) =  A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3);
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3) = 0; 
                end
            end
            if i >= nbc(j+2) &&  i <= nj(j+2)       % for on   (2)
               A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+1) = 0;
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+2) = 0;
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+3) = 0;
            elseif i > nj(j+2)                      % for outside (above)    
               A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) - A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+1);       % u
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+1) = 0;       % u
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) - A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+2);       % v
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+2) = 0;       % v
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+3) = 0;
            end
            if i+1 >= nbc(j+2) &&  i+1 <= nj(j+2)           % for on    (3)
               A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
            elseif  i+1 > nj(j+2)                                   % for outside
                if i-1 < nbc(j+2)                    % additional conditions
                    A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
                elseif i+1 < nbc(j)                % additional conditions
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
                else                              % additional conditions
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
                end
            end            
        elseif j == m
            if i > 1                        % for outside  (1b)
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j)-1)*3+3) =  A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3);
                A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+2)-1)*3+3) = 0;
            end
                                             % for outside  (2b)
               A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j+2)-1)*3+3) = 0;

               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) - A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+1);       % u
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+1) = 0;       % u
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) - A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+2);       % v
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+2) = 0;       % v
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+3);
               A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+2)-1)*3+3) = 0;
               
               if i+1 < nbc(j)                      % for outside  (3b) additional conditions
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3);
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
               else
                    A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+2)-1)*3+3) = 0;
               end
        end
        
        if i > 2
            if i-2 >= nbc(j+1) &&  i-2 <= nj(j+1)       % for on only    (4)
                A( (i+nbcc(j)-1)*3+3 , (i-2 + nbcc(j+1)-1)*3+3) = 0;
            end        
        end
        if i > 1
            if i-1 >= nbc(j+1) && i-1 <= nj(j+1)         % for on only   (5)
               A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+1 ) = 0;
               A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+2 ) = 0;
               A( (i+nbcc(j)-1)*3+1 , (i-1 + nbcc(j+1)-1)*3+3 ) = 0;

               A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+1 ) = 0;
               A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+2 ) = 0;
               A( (i+nbcc(j)-1)*3+2 , (i-1 + nbcc(j+1)-1)*3+3 ) = 0;

               A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+1 ) = 0;
               A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+2 ) = 0;
               A( (i+nbcc(j)-1)*3+3 , (i-1 + nbcc(j+1)-1)*3+3 ) = 0;
            end
        end
        if i >= nbc(j+1) && i <= nj(j+1)         % for on only          (6)
           A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i   + nbcc(j+1)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i   + nbcc(j+1)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j+1)-1)*3+3) = 0;
        end
        if i+1 >= nbc(j+1) && i+1 <= nj(j+1)         % for on only      (7)
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j+1)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j+1)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j+1)-1)*3+3) = 0;
        end
        if i+2 >= nbc(j+1) &&  i+2 <= nj(j+1)   % for on             (8)
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3) = 0;
        elseif i+2 > nj(j+1)                    % for outside
            if      i < nbc(j+1)                % additional conditions
                A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+1)-1)*3+3) =  A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j+1)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3) = 0; 
            elseif  j>1 && i+2 < nbc(j-1)              % additional conditions
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3) = A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3);
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3) = 0; 
            else                                % additional conditions
                A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j+1)-1)*3+3) = 0;
            end
        end
        
        if i+1 == nbc(j)                                % for on      (9)
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+1 ) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+2 ) = 0;
           A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j  )-1)*3+3 ) = 0;

           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+1 ) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+2 ) = 0;
           A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j  )-1)*3+3 ) = 0;

           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+1 ) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+2 ) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j  )-1)*3+3 ) = 0;
        end
        if i+2 >= nbc(j) &&  i+2 <= nj(j)           % for on           (10)
           A( (i+nbcc(j)-1)*3+1 , (i+2 + nbcc(j)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+2 , (i+2 + nbcc(j)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+1) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+2) = 0;
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+3) = 0;
        elseif i+2 > nj(j)                                      % for outside (right)
           A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) =  A( (i+nbcc(j)-1)*3+1 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+1 , (i+2 + nbcc(j)-1)*3+3);
           A( (i+nbcc(j)-1)*3+1 , (i+2 + nbcc(j)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) =  A( (i+nbcc(j)-1)*3+2 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+2 , (i+2 + nbcc(j)-1)*3+3);
           A( (i+nbcc(j)-1)*3+2 , (i+2 + nbcc(j)-1)*3+3) = 0;

           A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) =  A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+1) - A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+1);   % u
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+1) = 0;   % u
           A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) =  A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+2) - A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+2);   % v
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+2) = 0;   % v
           A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) =  A( (i+nbcc(j)-1)*3+3 , (i + nbcc(j)-1)*3+3) + A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+3);
           A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j)-1)*3+3) = 0; 
        end
        
        if j > 1 
            if i+1 == nbc(j-1) %|| i+1 > nbc(j-1)           % on       (11)
               A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+1 ) = 0;            % u
               A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+2 ) = 0;            % v
               A( (i+nbcc(j)-1)*3+1 , (i+1 + nbcc(j-1)-1)*3+3 ) = 0;            % w
               
               A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+1 ) = 0;            % u
               A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+2 ) = 0;            % v
               A( (i+nbcc(j)-1)*3+2 , (i+1 + nbcc(j-1)-1)*3+3 ) = 0;            % w
               
               A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+1 ) = 0;            % u
               A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+2 ) = 0;            % v
               A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-1)-1)*3+3 ) = 0;            % w
            end
            
            if i+2 >= nbc(j-1) &&  i+2 <= nj(j-1)       % on           (12)
               A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3) = 0;
            elseif i+2 > nj(j-1)                        % outside
               A( (i+nbcc(j)-1)*3+3 , (i   + nbcc(j-1)-1)*3+3) = A((i+nbcc(j)-1)*3+3 ,(i   + nbcc(j-1)-1)*3+3) + A((i+nbcc(j)-1)*3+3 ,(i+2 + nbcc(j-1)-1)*3+3); 
               A( (i+nbcc(j)-1)*3+3 , (i+2 + nbcc(j-1)-1)*3+3) = 0;
            end  
        end
        
        if j>2
            if i+1 == nbc(j-2)                   % on               (13)
               A( (i+nbcc(j)-1)*3+3 , (i+1 + nbcc(j-2)-1)*3+3) = 0;
            end
        end
       
    end
end

% A = A(:,1:nbcc(end)*3);
A = A(1:nbcc(end-1)*3,1:nbcc(end-1)*3);

disp = A\B';

%disp = gmres(A,B',[],[],1000);
% As = sparse(A);
% Bs = sparse(B);
% whos As
% whos Bs
% clear A;
% 
% disp = As\Bs';

%res = A*disp-B';

c=1;
for j=1:m+1
    for i=1:nj(j)
        if j==m+1 || i == nbc(j) || i > nbc(j)
            u(i,j) = 0;
            v(i,j) = 0;
            w(i,j) = 0;
        else
            u(i,j) = disp((i+nbcc(j)-1)*3+1);
            v(i,j) = disp((i+nbcc(j)-1)*3+2);
            w(i,j) = disp((i+nbcc(j)-1)*3+3);
            c=c+1;
        end
        
        a1(i,j) = norm([a1_1(i,j) a1_2(i,j) a1_3(i,j)]);
        a2(i,j) = norm([a2_1(i,j) a2_2(i,j) a2_3(i,j)]);
        a3(i,j) = norm([a3_1(i,j) a3_2(i,j) a3_3(i,j)]);
        
    end
end

[dudy,dudx] = gradient(u,dy,dx);
[dvdy,dvdx] = gradient(v,dy,dx);
[dwdy,dwdx] = gradient(w,dy,dx);

dudx(1,:) = u(2,:)/dx;
dvdx(1,:) = 0;
dwdx(1,:) = 0;

dudy(:,1) = 0;
dvdy(:,1) = v(:,2)/dy;
dwdy(:,1) = 0;

e1 = 1./A1.*(dudx.*dxdu1 + dudy.*dydu1) + 1./(A1.*B1).*(dAdx.*dxdu2 + dAdy.*dydu2).*v - w./R1;
e2 = 1./B1.*(dvdx.*dxdu2 + dvdy.*dydu2) + 1./(A1.*B1).*(dBdx.*dxdu1 + dBdy.*dydu1).*u - w./R2;
e12 = 1./(A1.*B1).*( B1.*(dvdx.*dxdu1 + dvdy.*dydu1) - v.*dBdu1 + A1.*(dudx.*dxdu2 + dudy.*dydu2) - u.*dAdu2 );

s1 = E./(1-vv.^2).*(e1 + vv.*e2);
s2 = E./(1-vv.^2).*(e2 + vv.*e1);
s12 = E./2./(1+vv).*e12;

s1_o = s1(1,1);
s2_o = s2(1,1);

R1_o = R1(1,1);
R2_o = R2(1,1);

s1_x = s1(:,1);
s2_x = s2(:,1);
s1_y = s1(1,:);
s2_y = s2(1,:);

z_m = max(max(z));

u_new = u.*a1_1./a1 + v.*a2_1./a2 + w.*a3_1./a3;
v_new = u.*a1_2./a1 + v.*a2_2./a2 + w.*a3_2./a3;
w_new = u.*a1_3./a1 + v.*a2_3./a2 + w.*a3_3./a3;


if plotting ==1

%% 3D plots

    figure;
    title('Undeformed shell showing 1st and 2nd principal curvature directions')
    surf(x,y,z)
    hold on
    quiver3(x,y,z, a1_1, a1_2, a1_3)
    quiver3(x,y,z, a2_1, a2_2, a2_3)
    hold off
    legend('Undeformed surface','\alpha direction','\beta direction')
    legend boxoff
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    xlim([0 inf])
    ylim([0 inf])
    zlim([0 inf])

    figure
    plot(xc,yc,xbc1,ybc1,'-',x,y,'.')
    xlabel('x (m)')
    ylabel('y (m)')
    title('Finite difference grid')
    xlim([0 inf])
    ylim([0 inf])

    figure;
    surf(x,y,R1)
    title('Principal radius of curvature R_1')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('R_1 (m)')
    xlim([0 inf])
    ylim([0 inf])

    figure;
    surf(x,y,R2)
    title('Principal radius of curvature R_2')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('R_2 (m)')
    xlim([0 inf])
    ylim([0 inf])
    

end


end
