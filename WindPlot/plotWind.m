% plot wind data stored in a .bts file based on the script of NREL
% readfile_BTS.m
% Written by Yang Yang in Ningbo University
% Date: April 27, 2021.
clear all;
close all;
% read wind data
Windfile = 'G:\StudyJeyla\Research\Training\OpenFAST_Training1\Lecture3\MLife\10MWRWT\WindData\OOStar_Demo_Kaimal_6.bts';  % change this for your own case
[velocity, twrVelocity, ygrid, zgrid, zTwr, nz, ny, dz, dy, dt, zHub, z1,mffws] = readfile_BTS(Windfile);
% ygrid: lateral gird; zgrid: vertical grid
NumPts = length(velocity(:,1,1,1)); % length of wind speed time series
TempVel = sqrt(velocity(:,1,:,:).^2+velocity(:,2,:,:).^2+velocity(:,3,:,:).^2);
ResVel = zeros(NumPts,ny,nz); % resulant wind speed
for i=1:NumPts
    for j=1:ny
         for k=1:nz
             ResVel(i,j,k) = TempVel(i,1,j,k);
         end
    end
end
% grid definition
TimeDur = dt*(NumPts-1); % Time duration 
Nhub = (nz-1)./2 + 1;    % Hub node
tgrid = 0:dt:TimeDur;    % Time vector
% hub height wind speed
Hubwind = zeros(length(tgrid),length(ygrid));
Hubwind(:,:)=ResVel(:,:,Nhub);
%% plot hub height wind speed
WindSpd_min = 0.25*mffws;   % for plot
WindSpd_max = 1.8*mffws;    % for plot
sizeFront = 20;             % for plot
[Ygrid,Tgrid] = meshgrid(ygrid,tgrid); % surface grid for hub-height wind
figure
surf(Ygrid,Tgrid,Hubwind);
xlabel('Horizontal position/m','FontWeight','bold', 'FontSize',sizeFront, 'FontName','times new roman');
ylabel('Time/s','FontWeight','bold', 'FontSize',sizeFront, 'FontName','times new roman');
zlabel('Wind speed/(m/s)','FontWeight','bold', 'FontSize',sizeFront, 'FontName','times new roman')
xlim([ygrid(1) ygrid(end)]);ylim([tgrid(1) tgrid(end)]);zlim([WindSpd_min WindSpd_max]);
view(102,16);
shading interp;
set(gca,'FontWeight','bold','FontSize',sizeFront,'FontName','times new roman');
colormap jet;
caxis([WindSpd_min,WindSpd_max]);
colorbar('FontWeight','bold','FontSize',sizeFront,'FontName','times new roman');
%% plot full field wind speed
dtSlice = 10; % s, delta time for slice
tSlice = 0:dtSlice:TimeDur;
%plotFFWind(ResVel,tgrid,ygrid,zgrid,tSlice)
[YG,TG,ZG]=meshgrid(ygrid,tgrid,zgrid);
figure
%slice(YG,TG,ZG,ResVel,ygrid,tSlice,zgrid);
h = slice(YG,TG,ZG,ResVel,ygrid(2:end-1),tSlice,zgrid(2:end-1));
set(h,'EdgeColor','none','FaceColor','interp','AlphaData',0.5); % trying to add transparency for each slice
colormap jet;
shading interp;
view(75,20);
caxis([WindSpd_min,WindSpd_max]);
colorbar('FontWeight','bold','FontSize',sizeFront,'FontName','times new roman');
ylabel({'Time/s'},'FontWeight','bold','FontSize',sizeFront,'FontName','Times New Roman',...
    'HorizontalAlignment','left');
ylim([tgrid(1) tgrid(end)])
xlabel({'Lateral position/m'},'FontWeight','bold','FontSize',sizeFront,'FontName','Times New Roman');
xlim([ygrid(1) ygrid(end)]);
zlabel({'Vertical/m'},'FontWeight','bold','FontSize',sizeFront,'FontName','Times New Roman');
zlim([zgrid(1) zgrid(end)]);
set(gca,'FontWeight','bold','FontSize',sizeFront,'FontName','times new roman');