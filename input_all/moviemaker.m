clc;clear

dep = load('depth');
prof = depth(:,1);

size_x = 64;
size_y = 8;
size_z = 20;

len_x = 14300/64;
len_y = 10;
%len_z = ;

dx= len_x/size_x;
dy= len_y/size_y;
dz= prof/size_z;

x = [dx:dx:len_x];
y = [dy:dy:len_y];
%z = -[len_z:-dz:dz];

time=load('time');
time = round(time);
t_str = int2str(time);

nFrames=length(time);

% load data

for i=1:nFrames;
%for i = 45;

    fnum = sprintf('%.4d',i);
    uu=load(['u_' fnum]);
    vv=load(['v_' fnum]);
    ww=load(['w_' fnum]);
    eta=load(['eta_' fnum]);
        
% x-direction uniform do not need x-dir average
% 1 represents take only one column in x-direction

    u=reshape(uu(:,1),64,60)';
    v=reshape(vv(:,1),64,60)';
    w=reshape(ww(:,1),64,60)';

    umean(:,i) = mean(u,2);% average over y-direction
    vmean(:,i) = mean(v,2);
    wmean(:,i) = mean(w,2);
    
    umean_mid(i) = mean(umean(:,i));

    % plot mean profile u-velocity
    figure(1)
    clf
    plot(umean(:,i)/umean_mid(i),z,'LineWidth',2)
    axis([0.6 1.3 -len_z -dz])
    ylabel('Z (m)','FontSize',10);
    title(['U (Time = ',t_str(i,:),' s)'])
    pause(0.5)
    
    uflux(i)=sum(umean(:,i))*dz*dy;
    vflux(i)=sum(vmean(:,i))*dz*dx;
    wflux(i)=sum(wmean(:,i))*dx*dy;
    
   % u',v',w' are defined as u - <u>, where <.> is average over x,y direction 
    udiff(:,:,i) = u - umean(:,i)*ones(1,64);
    vdiff(:,:,i) = v - vmean(:,i)*ones(1,64);
    wdiff(:,:,i) = w - wmean(:,i)*ones(1,64);
    
end

udiff_max = max(max(max(udiff)));
udiff_min = min(min(min(udiff)));
    
vdiff_max = max(max(max(vdiff)));
vdiff_min = min(min(min(vdiff)));
    
wdiff_max = max(max(max(wdiff)));
wdiff_min = min(min(min(wdiff)));

umn_max = max(max(umean));
umn_min = min(min(umean));

vmn_max = max(max(vmean));
vmn_min = min(min(vmean));

wmn_max = max(max(wmean));
wmn_min = min(min(wmean));

%%%%%%%%%%%%%%%% start to plot figures %%%%%%%%%%%%%%%%%%%%%%%%
mov1(1:nFrames) = struct('cdata', [],'colormap', []);
mov2(1:nFrames) = struct('cdata', [],'colormap', []);
mov3(1:nFrames) = struct('cdata', [],'colormap', []);
 
for i=1:nFrames;   
    %% velocity difference contour
    figure(1);
    clf;
    subplot(3,1,1)
    contourf(y,z,udiff(:,:,i)/umean_mid(i),10);
   % caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.15 0.15])
    colorbar
    ylabel('Z (m)','FontSize',10);
    title(['U (Time = ',t_str(i,:),' s)'])
    
    subplot(3,1,2)
    contourf(y,z,vdiff(:,:,i)/umean_mid(i),10);
    %caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.15 0.15])
    colorbar
    ylabel('Z (m)','FontSize',10);
    title(['V (Time = ',t_str(i,:),' s)'])

    subplot(3,1,3)
    contourf(y,z,wdiff(:,:,i)/umean_mid(i),10);
    %caxis([udiff_min udiff_max]/umean_mid(i))
    caxis([-0.15 0.15])
    colorbar
    xlabel('Y (m)','FontSize',10)
    ylabel('Z (m)','FontSize',10);
    title(['W (Time = ',t_str(i,:),' s)'])

    set(gca,'nextplot','replacechildren');
    mov(i) = getframe(gcf);
    
   % mean velocity plot

   
    figure(2)
    clf
    subplot(3,1,1)
    plot(umean(:,i)/umean_mid(i),z,'LineWidth',2)
    axis([0.6 1.3 -len_z -dz])
    ylabel('Z (m)','FontSize',10);
    title(['U (Time = ',t_str(i,:),' s)'])
    
    subplot(3,1,2)
    plot(vmean(:,i)/umean_mid(i),z,'LineWidth',2)
   % axis([vmn_min vmn_max -len_z -dz])
    axis([-1e-1 1e-1 -len_z -dz])
    ylabel('Z (m)','FontSize',10);
    title(['V (Time = ',t_str(i,:),' s)']) 
    
    subplot(3,1,3)
    plot(wmean(:,i)/umean_mid(i),z,'LineWidth',2)
   % axis([wmn_min wmn_max -len_z -dz])
    axis([-1e-3 1e-3 -len_z -dz])
    xlabel('Y (m)','FontSize',10)
    ylabel('Z (m)','FontSize',10);
    title(['W (Time = ',t_str(i,:),' s)'])
    
    set(gca,'nextplot','replacechildren');
    mov2(i) = getframe(gcf);

%     % surface elevation plot 
%     
%     figure(3)
%     clf
%     plot(y,eta(:,1),'LineWidth',2)
%     hold on
%     plot(y,y*0,'k--','LineWidth',2)
%     axis([dy len_y -0.05 0.05])
%     xlabel('Y (m)','FontSize',10);
%     ylabel('Eta (m)','FontSize',10);
% 
%     set(gca,'nextplot','replacechildren');
%     mov3(i) = getframe(gcf);

end
%return
movie2avi(mov, 'uvw1hr.avi', 'FPS',2,'compression', 'Cinepak','quality',50,'keyframe',1);
movie2avi(mov2, 'uvwmean1hr.avi', 'FPS',2,'compression','Cinepak','quality',50,'keyframe',1);

%movie2avi(mov3, 'eta1hr.avi', 'FPS',2,'compression', 'Cinepak','quality',50,'keyframe',1);
return
% plot mean U profile comparison over time 

% figure(4);
% 
% u0=load('uprofstdy.mat'); % should be 50 hour run
% u0=u0.u;
% plot(u0,z,'LineWidth',2)
% hold on
% plot(umean(:,20),z,'r--','LineWidth',2)
% plot(umean(:,40),z,'k--','LineWidth',2)
% plot(umean(:,60),z,'g--','LineWidth',2)
% xlabel('U (m/s)','FontSize',14)
% ylabel('Z (m)','FontSize',14)
% h_legend = legend('No LC','LC 20 min','LC 40 min','LC 60 min');
% set(h_legend,'FontSize',14)


% plot flux in U, V, W direction

figure(5)
subplot(3,1,1)
plot(time,uflux./umean_mid/dy/len_z,'LineWidth',2)
ylabel('x-Flux','FontSize',10)

subplot(3,1,2)
plot(time,vflux./umean_mid/dx/len_z,'LineWidth',2)
axis([min(time) 4000 -2e-5 2e-5])
ylabel('y-Flux ','FontSize',10)

subplot(3,1,3)
plot(time,wflux./umean_mid/dy/dx/size_z,'LineWidth',2)
axis([min(time) 4000 -2e-5 2e-5])
xlabel('Time (s)','FontSize',10)
ylabel('z-Flux ','FontSize',10)
