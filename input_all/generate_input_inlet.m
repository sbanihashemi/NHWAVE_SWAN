% construct ESTUARY grid.
clear;clc

% Grid parameters
Mglob = 64; % x-direction
Nglob = 60; % y-direction 

% Depth_Type = CELL_GRID

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct dx_r and dy_r.  Integrate to get x_r and y_r.
i=[1:Mglob];
j=[1:Nglob];

mp=(Mglob+1)/2; % ????
i_length = length(i)-2;
% Determine river channel center point 
mid=ceil(Nglob/2);
pt_shift = 0;    % shift channel southward
channel_center = mid - pt_shift

% NHWAVE uses uniform grid
dx1 = 200;% minimum dx 
dx2 = 200;% maximum dx

i_center = 31;% the point offshore of river mouth in i direction
j_center = channel_center; %position of center point
dx_r = i*0 + dx1;

dy1 = 200;% minimum dy
dy2 = 200;% maximum dy
dy_r=j*0+dy1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need NO CHANGE

[pm,pn]=meshgrid(1./dx_r,1./dy_r);
% Calculate x_r and y_r
dx=[dx_r(1)./2 0.5.*(dx_r(1:end-1)+dx_r(2:end))];
dy=[dy_r(1)./2 0.5.*(dy_r(1:end-1)+dy_r(2:end))];
[x_r,y_r]=meshgrid(cumsum(dx),cumsum(dy));

% Shift grid so x_u(:,1)=0 and y_v(1,:)=0.
y_r = y_r - y_r(1,1) - (y_r(2,1)-y_r(1,1))/2;
x_r = x_r - x_r(1,1) - (x_r(1,2)-x_r(1,1))/2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct grid topography.
B_total = 10;      % breath in points = B_total+1
                   % low resolution, use B = 10
h_mouth = 4.05;    % river mouth (shoaling area) depth
slope_land = 0;    % landward slope
slope_sea  = 1/200;% seaward slope
%RiverLength_highres = 3e3;    % highly resolved river length
i_land = Mglob % land is the end of cross-shore coordinate
i_mouth = i_center + 1 % beginning point of river mouth in crossshore
RiverLength_total = x_r(1,end) - x_r(1,i_mouth)
dh(1:Mglob) = 0;   % channel-shoaling area difference for constant depth channel
 % change of river channel depth per grid distance along j direction
%j_mid=ceil(Mp/2);
j_1 = channel_center - B_total/2;
j_2 = channel_center + B_total/2;
Width = y_r(j_2,1) - y_r(j_1,1) % Width of river channel 
%bb = Width/2;       % tuning parameter (see Li and Valle-Levinson 99)
W = y_r(channel_center,1);


% ------ make convergent channel (bb = width/50 to width/100 over 20-pts) -----
% ------ centered at i = i_land + 10                         ------------------
total_pt = i_land - i_mouth +1;% grid points for river
width_para(1:Mglob) = 2;            % constant width, wider channel
bb = Width./width_para;

% ------ make depth -------------------------
for i=1:Mglob
for j=1:Nglob
    h_min(j,i) = h_mouth + slope_land * ( x_r(1,i_mouth) - x_r(1,i) )+ dh(i);%h_min= 4.05 everywhere
end
end

for j=1:Nglob
for i=1:Mglob

if ( (j>j_1) & (j<=j_2) )% river channel grid point in j direction

  %h(j,i) = h_mouth+dh(i)* exp(-( (y_r(j,i)-W )/bb(i))^2 ); %Gaussian shape
  
  h(j,i) = h_mouth + dh(i); % Rectangular shape

else
  h(j,i)=h_mouth;
  %h(j,i) = h_min(j,i) + dh(i)* exp(-( (y_r(j_1,i)-W )/bb(i) )^2 ); 

end
end
end
%return
% -------- this creates offshore slope --------------------------------------
h_mouth_shoal = h(j_1-1,i_mouth)      % shoal depth at the mouth 0.05
depth_max = 4.05;

for j=1:Nglob;
for i=1:i_mouth-1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %h_temp = h_mouth_shoal + slope_sea * (x_r(j,i_mouth)-x_r(j,i));
  %h(j,i) = max( h(j,i), h_temp); % constant slope beach
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %h_ebp = 0.08 * (x_r(j,i_mouth)-x_r(j,i))^(2/3);% EBP 
  %h(j,i) = max( h(j,i), h_ebp );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   h_tan0 = 25;%tanh-shaped beach profile
%   h_alpha = 0.05/h_tan0;
%   h_tan = h_tan0 * tanh(h_alpha*(x_r(j,i_mouth)-x_r(j,i)));
%   h(j,i) = max(h(j,i), h_tan); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h(j,i)=4.05;
end
end
%return
% i_max = max( find( h(1,:) > depth_max ) );
% h(1,i_max)

% for j=1:Mp
% for i=1:i_mouth
%   h(j,i) = min( h(j,i), h(j,i_max) );
% end
% end

disp(['max depth is ' num2str(max(max(h)))])
% ----------------------------------------------------------------------

dx=ones(Nglob,1)*dx_r;
dy=dy_r'*ones(1,Mglob);
hbar=sum(sum(h.*dx.*dy))./sum(sum(dx.*dy));
disp(['Hbar = ',num2str(hbar)]);

disp(['Western cross-sectional area = ',...
	num2str(sum(h(2:end-1,1).*dy(2:end-1,1)))]);
disp(['Eastern cross-sectional area = ',...
	num2str(sum(h(2:end-1,end).*dy(2:end-1,end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create masking grids used for SWAN input file
rmask=ones(size(x_r));

for i=i_mouth:i_mouth
    for j=1:(j_1)
        rmask(j,i)=0;
    end
    for j=(j_2+1):Nglob
        rmask(j,i)=0;
    end
end

for i=i_mouth:Mglob
    rmask(1,i)=0;
end

for i=i_mouth:Mglob
    rmask(Nglob,i)=0;
end

for j=1:Nglob
    rmask(j,Mglob)=0;
end

%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate SWAN computational grid and bathymetry

%Generate SWAN bathymetry
swan_bath=h;
swan_bath(rmask==0)=-10;
save swan_bathy.bot swan_bath -ASCII

%Generate SWAN computational grid and bathymetry
gridx = [0:Mglob-1]*dx1;
gridy = [0:Nglob-1]*dy1;
for j = 1:Nglob;
swan_gridx(((j-1)*Mglob+1):Mglob*j)=gridx;
end
MNglob = Mglob*Nglob;
for j = 1:Nglob;
    swan_gridy(((j-1)*Mglob+1):Mglob*j) = gridy(j)*ones(1,Mglob);
end
swan_grid2(1:MNglob) = swan_gridx;
swan_grid2(MNglob+1:2*MNglob) = swan_gridy;
swan_grid2 = swan_grid2';
save swan_grid_coord.grd swan_grid2 -ASCII

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate NHWAVE bathymetry
h_inlet = h;
h_inlet(rmask==0)=-9999;
save depth.txt h_inlet -ascii


% Plotting
plot1=1;

if plot1==1
figure;clf
rmask(find(rmask==0)) = nan;
h_m = h.*rmask;
pcolor(x_r,y_r,h_m)

%set(gca,'dataaspectratio',[1 1 1])
hc=colorbar('horiz');
hcl=get(hc,'ylabel');
set(hcl,'string','Depth (m)');
xlabel ('xi distance (m)');
ylabel ('eta distance (m)');
title ('Topography');


%return

% figure;clf
% plot (y_r(2:end-1,:),h(2:end-1,:),'-o');
end
return
figure;clf
subplot(2,1,1);
plot (x_r./1000,dx_r./1000);hold on;
plot (x_r./1000,dx_r./1000,'o');
xlabel('Along channel distance (km)');
ylabel('Along channel grid spacing (km)');
subplot(2,1,2);
plot (y_r./1000,dy_r./1000);hold on;
plot (y_r./1000,dy_r./1000,'o');
xlabel('Cross channel distance (km)');
ylabel('Cross channel grid spacing (km)');


