clear; clc; close all;
%% read in the simulated strain, compute stress, and plot on Parkfield Section
% Haoran Meng, Jan 22, 2020
% Simulation code: SPECFEM3D
% 3D velocity model: CVMS4.26
% Finite fault model: Yue et al. 2020, in preparation 
%   FFM inverted from hr-GPS/static-GPS/Strong-ground-motion/InSAR/


%% Mainshock
load pk_str_RC71.mat; flag = 1;
% RidgeCrest event 2019/ 7/ 6 Mag:7.10
% Initial time at 3:19:53; Hypocenter at  35.77 N, -117.60 E, 8 km

%% Forshock
% load pk_str_RC64.mat; flag = 2;
% RidgeCrest event 2019/ 7/ 4 Mag:6.40
% Initial time at 17:33:49; Hypocenter at  35.69 N, -117.50 E, 11 km

%% time samples, fs = 5 Hz, start from 8 sec before the origin time
t = 1:2000; fs = 5; t = -8 + (t-1)/fs;

%% compute stress using lame parameters
% x: east; y: north; z: up
% Tensile positive, compression negative
[NX,NZ] = size(lambda);
for i = 1:NX
    for j = 1:NZ
         temp = lambda(i,j)*(rcv_Exx(i,j,:) + rcv_Eyy(i,j,:) + rcv_Ezz(i,j,:));
         rcv_Sxx(i,j,:) = temp + 2*mu(i,j)*rcv_Exx(i,j,:);
         rcv_Syy(i,j,:) = temp + 2*mu(i,j)*rcv_Eyy(i,j,:);
         rcv_Szz(i,j,:) = temp + 2*mu(i,j)*rcv_Ezz(i,j,:);
         rcv_Sxy(i,j,:) = 2*mu(i,j)*rcv_Exy(i,j,:);
         rcv_Sxz(i,j,:) = 2*mu(i,j)*rcv_Exz(i,j,:);
         rcv_Syz(i,j,:) = 2*mu(i,j)*rcv_Eyz(i,j,:);
    end
end

%% compute direction vectors and corresponding stress components
% pkst_lat & pkst_lon are the lat and lon of grids on Parkfield Section,
% 200 km long (NX = 41), 30 km deep(NZ = 7)
[arclen,az] = distance(pkst_lat(1),pkst_lon(1),pkst_lat(end),pkst_lon(end));
az = az + 90 - 360;
% normal to the vertical fault interface, towards NE
ex = sin(az/180*pi); ey = cos(az/180*pi); ez = 0;

% compute the traction
rcv_Tx = rcv_Sxx*ex + rcv_Sxy*ey + rcv_Sxz*ez;
rcv_Ty = rcv_Sxy*ex + rcv_Syy*ey + rcv_Syz*ez;
rcv_Tz = rcv_Sxz*ex + rcv_Syz*ey + rcv_Szz*ez;

% normal stress, tensile positive
rcv_Nrm = rcv_Tx*ex + rcv_Ty*ey + rcv_Tz*ez;

% along strike shear stress
az = az + 90; ex = sin(az/180*pi); ey = cos(az/180*pi); ez = 0;
rcv_Stk = rcv_Tx*ex + rcv_Ty*ey + rcv_Tz*ez;

% along dip shear stress
ex = 0; ey = 0; ez = 1;
rcv_Dp = rcv_Tx*ex + rcv_Ty*ey + rcv_Tz*ez;

%% lowpass filter results without phase shift, corner at 5 sec
fc = 1/5; [b,a] = butter(6,fc/(fs/2));
for i = 1:NX
    for j = 1:NZ
        temp = filtfilt(b,a,squeeze(rcv_Nrm(i,j,:))); rcv_Nrm(i,j,:) = temp;
        temp = filtfilt(b,a,squeeze(rcv_Stk(i,j,:))); rcv_Stk(i,j,:) = temp;
        temp = filtfilt(b,a,squeeze(rcv_Dp(i,j,:))); rcv_Dp(i,j,:) = temp;
    end
end

%% plot frames
for k = 1:200
    disp(k/200);
    n = k*5;
    h = figure('visible','off');
    
    % normal stress
    subplot(3,1,1)
    pcolor(hori,dep,transpose(squeeze(rcv_Nrm(:,:,n)/1e3)))
    axis equal; shading interp; set(gca, 'Xdir', 'reverse')
    text(4,5,'SE'); text(204,5,'NW'); text(-6,5,'(kPa)');
    xlim([0 200]); ylim([-30 0]);
    colormap(redblue);colorbar; hold on; 
    if flag == 1; caxis([-50 50]); elseif flag == 2; caxis([-10 10]); end
    xlabel('Along-fault distance (km)');ylabel('Depth (km)');
    title(['Normal stress at t = ' num2str(t(n)) ' s' ]); set(gca,'fontsize',12)
    
    % along strike shear stress
    subplot(3,1,2)
    pcolor(hori,dep,transpose(squeeze(rcv_Stk(:,:,n)/1e3)))
    axis equal; shading interp; set(gca, 'Xdir', 'reverse')
    text(4,5,'SE'); text(204,5,'NW'); text(-6,5,'(kPa)');
    xlim([0 200]); ylim([-30 0]);
    colormap(redblue);colorbar; hold on;
    if flag == 1; caxis([-50 50]); elseif flag == 2; caxis([-10 10]); end
    xlabel('Along-fault distance (km)');ylabel('Depth (km)');
    title('Along-strike shear stress'); set(gca,'fontsize',12)
    
    % along dip shear stress
    subplot(3,1,3)
    pcolor(hori,dep,transpose(squeeze(rcv_Dp(:,:,n)/1e3)))
    axis equal; shading interp; set(gca, 'Xdir', 'reverse')
    text(4,5,'SE'); text(204,5,'NW'); text(-6,5,'(kPa)');
    xlim([0 200]); ylim([-30 0]);
    colormap(redblue);colorbar; hold on;
    if flag == 1; caxis([-50 50]); elseif flag == 2; caxis([-10 10]); end
    xlabel('Along-fault distance (km)');ylabel('Depth (km)');
    title('Along-dip shear stress'); set(gca,'fontsize',12)
    
    % save as figures
    if flag == 1
        if k < 10
            saveas(h,['./figures/RC71_PK_STRESS_00' num2str(k) ],'png')
        elseif k<100 && k>=10
            saveas(h,['./figures/RC71_PK_STRESS_0' num2str(k) ],'png')
        else
            saveas(h,['./figures/RC71_PK_STRESS_' num2str(k) ],'png')
        end
    end
    if flag == 2
        if k < 10
            saveas(h,['./figures/RC64_PK_STRESS_00' num2str(k) ],'png')
        elseif k<100 && k>=10
            saveas(h,['./figures/RC64_PK_STRESS_0' num2str(k) ],'png')
        else
            saveas(h,['./figures/RC64_PK_STRESS_' num2str(k) ],'png')
        end
    end

end


function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG,
%   COLORMAP, RGBPLOT.
%   Adam Auton, 9th October 2009
if nargin < 1, m = size(get(gcf,'colormap'),1); end
if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end
c = [r g b];
end