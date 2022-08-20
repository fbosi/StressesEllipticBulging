function [p1_1,p1_2,p1_3,p2_1,p2_2,p2_3,a3_1,a3_2,a3_3,R1,R2,c,Prom] = Principal_curvature_fun (x,y,z,dzdx,dzdy,d2zdx,d2zdx_dy,d2zdy,a,b)

dx = x(2,1)-x(1,1);
dy = y(1,2)-y(1,1);
[n,m]=size(z);

%% Geometry 

a11 = 1 + dzdx.^2;
a22 = 1 + dzdy.^2;
a12 = dzdx.*dzdy;

a3_1 = (-1).*dzdx.*(1+dzdy.^2+dzdx.^2).^(-1/2);

a3_2 = (-1).*dzdy.*(1+dzdy.^2+dzdx.^2).^(-1/2);

a3_3 = (1+dzdy.^2+dzdx.^2).^(-1/2);

b11 = (1+dzdy.^2+dzdx.^2).^(-1/2).*d2zdx;
b22 = d2zdy.*(1+dzdy.^2+dzdx.^2).^(-1/2);
b12 = (1+dzdy.^2+dzdx.^2).^(-1/2).*d2zdx_dy;

lam_1 = (-(b11.*a22 - b22.*a11) + sqrt((b11.*a22 - a11.*b22).^2 - 4*(b12.*a22 - b22.*a12).*(b11.*a12 - b12.*a11)) ) ./ ( 2*(b12.*a22 - b22.*a12) );
dbda = lam_1;

sin_a = dzdx./(sqrt(1 + dzdx.^2));
sin_b = dzdy./(sqrt(1 + dzdy.^2));
cos_a = 1./(sqrt(1 + dzdx.^2));
cos_b = 1./(sqrt(1 + dzdy.^2));
sin_g = sqrt(1 - sin_a.^2.*sin_b.^2);
cos_g = sin_a.*sin_b;

gamma = acos(cos_g);

th_1 = atan(dbda.*sin(pi - gamma)./(1 - dbda.*cos(pi - gamma)));   %zeros(size(a11));


% for j=1:m
%     for i=1:n
%         if th_1(i,j)<0
%             th_1(i,j)=th_1(i,j)+pi/2;
%         end
%     end
% end


%% Dealing with singularity

[K,H,Pmax,Pmin] = surfature(x,y,z);
 
term = sqrt(H.^2 - K);
[min_term,c] = min(term(1:end-4,3));
% clear c;
[TF,Prom] = islocalmin(term(3:end-4,3),'MinProminence',1);  %2.5
if max(Prom) > 0.5
    [min_prom,c] = max(Prom);
else
    c = length(term(3:end-4,3));
end

c=c+2;

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

if a > b
    th_1(3,:) = 0;                  % at x=0
    if b/a < 0.5
        th_1(:,3) = 0;              % at y=0
    else
        th_1(1:c,3) = 0;            % at y=0
        th_1(c+1:end,3) = -pi/2;
    end
else
    th_1(3,:) = pi/2;   % at x=0
    th_1(:,3) = pi/2;  
    for j=1:m
        for i=1:n
            if th_1(i,j) > 135*pi/180 ||  th_1(i,j) < -45*pi/180        % condition to remain in possible domain
                th_1(i,j) = th_1(i,j) + pi;
            end
        end
    end
end
%%

th_1(:,1) = -th_1(:,5);
th_1(:,2) = -th_1(:,4);

th_1(1,:) = -th_1(5,:);
th_1(2,:) = -th_1(4,:);

for j=1:m
    for i=1:n
        if th_1(i,j) > 45*pi/180         % condition to remain in possible domain
            th_2(i,j) = th_1(i,j) - pi/2;
        else
            th_2(i,j) = th_1(i,j) + pi/2;
        end
        
    end
end


for j=1:m
    for i=1:n
        Rot_1 = [cos(th_1(i,j)) + a3_1(i,j).^2*(1-cos(th_1(i,j)))                   ,  a3_1(i,j)*a3_2(i,j)*(1-cos(th_1(i,j)))-a3_3(i,j)*sin(th_1(i,j))  , a3_1(i,j)*a3_3(i,j)*(1-cos(th_1(i,j)))+a3_2(i,j)*sin(th_1(i,j)) ; ...
             a3_2(i,j)*a3_1(i,j)*(1-cos(th_1(i,j))) + a3_3(i,j)*sin(th_1(i,j))  ,  cos(th_1(i,j)) + a3_2(i,j).^2*(1-cos(th_1(i,j)))                 , a3_2(i,j)*a3_3(i,j)*(1-cos(th_1(i,j)))-a3_1(i,j)*sin(th_1(i,j)) ; ...
             a3_3(i,j)*a3_1(i,j)*(1-cos(th_1(i,j))) - a3_2(i,j)*sin(th_1(i,j))  ,  a3_3(i,j)*a3_2(i,j)*(1-cos(th_1(i,j)))+a3_1(i,j)*sin(th_1(i,j))  , cos(th_1(i,j)) + a3_3(i,j).^2*(1-cos(th_1(i,j))) ];
         
         Rot_2 = [cos(th_2(i,j)) + a3_1(i,j).^2*(1-cos(th_2(i,j)))                   ,  a3_1(i,j)*a3_2(i,j)*(1-cos(th_2(i,j)))-a3_3(i,j)*sin(th_2(i,j))  , a3_1(i,j)*a3_3(i,j)*(1-cos(th_2(i,j)))+a3_2(i,j)*sin(th_2(i,j)) ; ...
             a3_2(i,j)*a3_1(i,j)*(1-cos(th_2(i,j))) + a3_3(i,j)*sin(th_2(i,j))  ,  cos(th_2(i,j)) + a3_2(i,j).^2*(1-cos(th_2(i,j)))                 , a3_2(i,j)*a3_3(i,j)*(1-cos(th_2(i,j)))-a3_1(i,j)*sin(th_2(i,j)) ; ...
             a3_3(i,j)*a3_1(i,j)*(1-cos(th_2(i,j))) - a3_2(i,j)*sin(th_2(i,j))  ,  a3_3(i,j)*a3_2(i,j)*(1-cos(th_2(i,j)))+a3_1(i,j)*sin(th_2(i,j))  , cos(th_2(i,j)) + a3_3(i,j).^2*(1-cos(th_2(i,j))) ];
         
         p1 = Rot_1*[1;0;dzdx(i,j)]; 
         p1_1(i,j) = p1(1);
         p1_2(i,j) = p1(2);
         p1_3(i,j) = p1(3);
         
         p2 = Rot_2*[1;0;dzdx(i,j)]; 
         p2_1(i,j) = p2(1);
         p2_2(i,j) = p2(2);
         p2_3(i,j) = p2(3);
         
    end
end

for j=1:m
    for i=1:n
        if i==1
              dp1_1dx(i,j) = ( p1_1(i+1,j) - p1_1(i,j) ) / dx;      %forward difference
              dp1_2dx(i,j) = ( p1_2(i+1,j) - p1_2(i,j) ) / dx;      
              dp1_3dx(i,j) = ( p1_3(i+1,j) - p1_3(i,j) ) / dx;      
              
              dp2_1dx(i,j) = ( p2_1(i+1,j) - p2_1(i,j) ) / dx;      
              dp2_2dx(i,j) = ( p2_2(i+1,j) - p2_2(i,j) ) / dx;      
              dp2_3dx(i,j) = ( p2_3(i+1,j) - p2_3(i,j) ) / dx;      
        elseif i==n
              dp1_1dx(i,j) = ( p1_1(i,j) - p1_1(i-1,j) ) / dx;      %backward difference
              dp1_2dx(i,j) = ( p1_2(i,j) - p1_2(i-1,j) ) / dx;
              dp1_3dx(i,j) = ( p1_3(i,j) - p1_3(i-1,j) ) / dx;
              
              dp2_1dx(i,j) = ( p2_1(i,j) - p2_1(i-1,j) ) / dx;
              dp2_2dx(i,j) = ( p2_2(i,j) - p2_2(i-1,j) ) / dx;
              dp2_3dx(i,j) = ( p2_3(i,j) - p2_3(i-1,j) ) / dx;
        else
              dp1_1dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p1_1((-1)+i,j) + p1_1(1+i,j) );
              dp1_2dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p1_2((-1)+i,j) + p1_2(1+i,j) );
              dp1_3dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p1_3((-1)+i,j) + p1_3(1+i,j) );
         
              dp2_1dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p2_1((-1)+i,j) + p2_1(1+i,j) );
              dp2_2dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p2_2((-1)+i,j) + p2_2(1+i,j) );
              dp2_3dx(i,j) = (1/2).*dx.^(-1).*( (-1).*p2_3((-1)+i,j) + p2_3(1+i,j) );
        end
        
        if j==1
              dp1_1dy(i,j) = ( p1_1(i,j+1) - p1_1(i,j) ) / dy;      %forward difference
              dp1_2dy(i,j) = ( p1_2(i,j+1) - p1_2(i,j) ) / dy;      
              dp1_3dy(i,j) = ( p1_3(i,j+1) - p1_3(i,j) ) / dy;      
              
              dp2_1dy(i,j) = ( p2_1(i,j+1) - p2_1(i,j) ) / dy;      
              dp2_2dy(i,j) = ( p2_2(i,j+1) - p2_2(i,j) ) / dy;      
              dp2_3dy(i,j) = ( p2_3(i,j+1) - p2_3(i,j) ) / dy;      
        elseif j==m
              dp1_1dy(i,j) = ( p1_1(i,j) - p1_1(i,j-1) ) / dy;      %backward difference
              dp1_2dy(i,j) = ( p1_2(i,j) - p1_2(i,j-1) ) / dy;
              dp1_3dy(i,j) = ( p1_3(i,j) - p1_3(i,j-1) ) / dy;
              
              dp2_1dy(i,j) = ( p2_1(i,j) - p2_1(i,j-1) ) / dy;
              dp2_2dy(i,j) = ( p2_2(i,j) - p2_2(i,j-1) ) / dy;
              dp2_3dy(i,j) = ( p2_3(i,j) - p2_3(i,j-1) ) / dy;
        else
              dp1_1dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p1_1(i,(-1)+j) + p1_1(i,1+j) );
              dp1_2dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p1_2(i,(-1)+j) + p1_2(i,1+j) );
              dp1_3dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p1_3(i,(-1)+j) + p1_3(i,1+j) );

              dp2_1dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p2_1(i,(-1)+j) + p2_1(i,1+j) );
              dp2_2dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p2_2(i,(-1)+j) + p2_2(i,1+j) );
              dp2_3dy(i,j) = (1/2).*dy.^(-1).*( (-1).*p2_3(i,(-1)+j) + p2_3(i,1+j) );
        end
      
    end
end

b11_p = (a3_1.*dp1_1dx + a3_2.*dp1_2dx + a3_3.*dp1_3dx).*p1_1 + (a3_1.*dp1_1dy + a3_2.*dp1_2dy + a3_3.*dp1_3dy).*p1_2 ; 
b22_p = (a3_1.*dp2_1dx + a3_2.*dp2_2dx + a3_3.*dp2_3dx).*p2_1 + (a3_1.*dp2_1dy + a3_2.*dp2_2dy + a3_3.*dp2_3dy).*p2_2 ; 

p11 = p1_1.*p1_1 + p1_2.*p1_2 + p1_3.*p1_3;
p22 = p2_1.*p2_1 + p2_2.*p2_2 + p2_3.*p2_3;

p1_1 = p1_1./sqrt(p11);
p1_2 = p1_2./sqrt(p11);
p1_3 = p1_3./sqrt(p11);

p2_1 = p2_1./sqrt(p22);
p2_2 = p2_2./sqrt(p22);
p2_3 = p2_3./sqrt(p22);

R1 = (p11./b11_p);
R2 = (p22./b22_p);

end