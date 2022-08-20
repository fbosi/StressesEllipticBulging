clear;

global dx  dy

ho = 100e-6;           % thickness in meters
v = 0.5; 

%% Inflation test of elliptical membrane - Experiments - with syncing

path = pwd;
filename = 'Grid_data_Inflation-0150_0.csv';

% Points for extracting data from DIC in portrait orientation
point_x = 83;     % for ellipse  b/a 0.5
point_y = 143;
a = 0.03;           % in meters
b = a*2;            % in meters


%% Data extract DIC


data_grid_dic = xlsread( [path '\' filename] );  
P=0.9175*6894.76;

[row,col]=size(data_grid_dic);
data_grid_dic(row+1:point_x*point_y,:) = nan;
x_dic=reshape(data_grid_dic(:,1),[point_x,point_y]);
y_dic=reshape(data_grid_dic(:,2),[point_x,point_y]);
z_dic=reshape(data_grid_dic(:,3),[point_x,point_y]);
e1_dic=reshape(data_grid_dic(:,10),[point_x,point_y]);
e2_dic=reshape(data_grid_dic(:,11),[point_x,point_y]);
k1_dic=reshape(data_grid_dic(:,end-6),[point_x,point_y]);  
k2_dic=reshape(data_grid_dic(:,end-5),[point_x,point_y]);   
gamma_dic = reshape(data_grid_dic(:,end-4),[point_x,point_y]);  
phi_dic = reshape(data_grid_dic(:,end-3),[point_x,point_y]);   

e3_dic = -e1_dic-e2_dic;
h_dic = exp(e3_dic)*ho;          % in meters

% unit conversion
x_dic=x_dic*1e-3;
y_dic=y_dic*1e-3;
z_dic=z_dic*1e-3;
R1_dic = 1./k1_dic/2/1000;        % mm
R2_dic = 1./k2_dic/2/1000;        % mm

% dx = nanmean(nanmean(diff(x)));   % optional
% dy = nanmean(nanmean(diff(y')));

dx = 0.0005;
dy = 0.0005;

% grid creation for FDM
x1 = -0.01:dx:a+2*dx;
y1 = -0.01:dy:b+2*dy;
[x,y] = meshgrid(x1,y1);
x = x';
y = y';

if abs(sum(sum(x-x_dic,'omitnan'))) > 1e-5
    disp('Something wrong')
    pause;
end

% Extrapolation
[q,r] = size(x);
d=0.7;
counter =1;
for j=1:r
    for i=1:q
        if x(i,j)^2/(d*a)^2 + y(i,j)^2/(d*b)^2 > 1 && ~isnan(z_dic(i,j))
            x_out(counter) = x(i,j);
            y_out(counter) = y(i,j);
            z_out(counter) = z_dic(i,j);
            R1_out(counter) = R1_dic(i,j);
            R2_out(counter) = R2_dic(i,j);
            h_out(counter) = h_dic(i,j);
            gamma_out(counter) = gamma_dic(i,j);
            phi_out(counter) = phi_dic(i,j);
            counter = counter + 1;
        end
    end
end

F_R1 = scatteredInterpolant(x_out',y_out',R1_out','natural','linear');
F_R2 = scatteredInterpolant(x_out',y_out',R2_out','natural','linear');
F_h = scatteredInterpolant(x_out',y_out',h_out','natural','linear');
F_gamma = scatteredInterpolant(x_out',y_out',gamma_out','natural','linear');
F_phi = scatteredInterpolant(x_out',y_out',phi_out','natural','linear');

xc = min(x1):0.0001:a;
yc = sqrt((1-xc.^2/a^2)*b^2);

PolySurf = fit([x_out',y_out'],z_out','poly44');
Zcoord = PolySurf(xc(1),yc(1));

x_out = [x_out xc];
y_out = [y_out yc];
z_out = [z_out-Zcoord zeros(size(xc))];

F_z = scatteredInterpolant(x_out',y_out',z_out','natural','linear');

for j=1:r
    for i=1:q
        if isnan(z_dic(i,j))
            z(i,j) = F_z(x(i,j),y(i,j));
            R1(i,j) = F_R1(x(i,j),y(i,j));
            R2(i,j) = F_R2(x(i,j),y(i,j));
            h(i,j) = F_h(x(i,j),y(i,j));
            gamma(i,j) = F_gamma(x(i,j),y(i,j));
            phi(i,j) = F_phi(x(i,j),y(i,j));
        else
            z(i,j) = z_dic(i,j)-Zcoord;
            R1(i,j) = R1_dic(i,j);
            R2(i,j) = R2_dic(i,j);
            h(i,j) = h_dic(i,j);
            gamma(i,j) = gamma_dic(i,j);
            phi(i,j) = phi_dic(i,j);
        end
    end
end

x_array = x_dic(:);
y_array = y_dic(:);
h_array = h_dic(:);
z_array = z_dic(:);

x_array(isnan(x_array)) = [];   % to remove nan
y_array(isnan(y_array)) = [];   % to remove nan
z_array(isnan(z_array)) = [];   % to remove nan

% Removing outliers
for j=1:size(x,2)
    for i=1:size(x,1)
        if gamma(i,j) > 135*pi/180 ||  gamma(i,j) < -45*pi/180        % condition to remain in possible domain
            gamma(i,j) = gamma(i,j) + pi;
        end
    end
end

x_array_outlier = x(:);
y_array_outlier = y(:);
gamma_array_outlier = gamma(:);
TF = isoutlier(gamma_array_outlier,'median','ThresholdFactor',3);
gamma_array_outlier = gamma_array_outlier(~TF);
x_array_outlier = x_array_outlier(~TF);
y_array_outlier = y_array_outlier(~TF);
%gamma_array(isoutlier(gamma_array)) = [];

F_gamma_outlier = scatteredInterpolant(x_array_outlier,y_array_outlier,gamma_array_outlier,'natural','linear');
TF = isoutlier(gamma,'median','ThresholdFactor',3);
c=0;
for j=1:r
    for i=1:q
        if TF(i,j)
            c=c+1;
            gamma(i,j) = F_gamma_outlier(x(i,j),y(i,j));
        end
    end
end

x_array = [x_array ; xc'];
y_array = [y_array ; yc'];
z_array = [z_array ; zeros(size(xc))'];

% gid reduction
x = x(19:end,19:end);
y = y(19:end,19:end);
z = z(19:end,19:end);
h = h(19:end,19:end);
R1 = R1(19:end,19:end);
R2 = R2(19:end,19:end);
gamma = gamma(19:end,19:end);
phi = phi(19:end,19:end);

%%%  Swapping x and y - required for potrait orientation only
temp = x;
x=y;
y=temp;
temp = a;
a=b;
b=temp;
temp = R1;
R1=R2;
R2=temp;
x=x';
y=y';
z=z';
h=h';
R1=R1';
R2=R2';
gamma=gamma';
gamma = gamma-pi/2;

%%% smoothing (optional)
[sp,R1,rho]=spaps({x(:,1),y(1,:)},R1,mean(mean(R1))/1000,1,3);
[sp,R2,rho]=spaps({x(:,1),y(1,:)},R2,mean(mean(R2))/1000,1,3);
[sp,gamma,rho]=spaps({x(:,1),y(1,:)},gamma,mean(mean(gamma))/1000,1,3);


%% Solution

[s1,s2,s12,z_m,R1,R2] = FDM_small_displacement_shell_gridded_fun(x,y,z,h,R1,R2,gamma,P,v,a,b,1);

S1 = (s1 + s2)/2 + sqrt( ((s1-s2)/2).^2 + s12.^2 );
S2 = (s1 + s2)/2 - sqrt( ((s1-s2)/2).^2 + s12.^2 );
S_eqv = 1/sqrt(2)*sqrt( (S1-S2).^2 + S1.^2 + S2.^2 );

x = x(3:end-2,3:end-2);
y = y(3:end-2,3:end-2);
z = z(3:end-2,3:end-2);

%% Contour plots

z_dic = z_dic(21:end-2,21:end-2)';
        
figure;
hold on
fig1=surfc(x,y,z-max(min(z_dic)),S1/1e6,'EdgeColor','none','EdgeColor','none','FaceColor','interp');
title('1st principal stress \sigma_I (MPa)')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
cmap=[255 0 0
255 92 0
255 185 0
231 255 0
139 255 0
46 255 0
0  255 46
0  255 139
0  255 231
0  185 255
0  92  255
0  0   255]/255;
colormap(flip(cmap));
colorbar
grid off
zlim([-max(min(z_dic)) inf])
caxis([S1(end,1) S1(1,1)]/1e6)

figure;
hold on
fig2=surfc(x,y,z - max(min(z_dic)),S2/1e6,'EdgeColor','none','EdgeColor','none','FaceColor','interp');
title('2nd principal stress \sigma_{II} (MPa)')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
cmap=[255 0 0
255 92 0
255 185 0
231 255 0
139 255 0
46 255 0
0  255 46
0  255 139
0  255 231
0  185 255
0  92  255
0  0   255]/255;
colormap(flip(cmap));
colorbar
grid off
zlim([-max(min(z_dic)) inf])
caxis([S2(end,1) S2(1,1)]/1e6)

figure;
hold on
fig3=surfc(x,y,z - max(min(z_dic)),S_eqv/1e6,'EdgeColor','none','FaceColor','interp');
title('Equivalent stress $\bar{\sigma}$ (MPa)','Interpreter','Latex')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
cmap=[255 0 0
255 92 0
255 185 0
231 255 0
139 255 0
46 255 0
0  255 46
0  255 139
0  255 231
0  185 255
0  92  255
0  0   255]/255;
colormap(flip(cmap));
colorbar
grid off
zlim([-max(min(z_dic)) inf])
caxis([S_eqv(end,1) S_eqv(1,1)]/1e6)


