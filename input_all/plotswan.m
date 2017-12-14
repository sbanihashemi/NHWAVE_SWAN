%%% this file is used to plot inlet test
%%% wave properties

clear;clc

size_x = 64;
size_y = 16;
size_z = 20;

dx= 14300/64;
dy= 1000;


grav = 9.8;
rho = 1000;

x_grid=[1:size_x]*dx;
y_grid=[1:size_y]*dy;

mask = zeros(size_y,size_x);

x_skip = [1:2:size_x];
y_skip = [1:2:size_y];

depth = load('depth');
time=load('time');
time = round(time);
t_str = int2str(time);
nFrames=length(time);
mov(1:nFrames) = struct('cdata', [],'colormap', []);

for i=nFrames:nFrames;
    f_number = i;
    fnum = sprintf('%.4d',f_number);
    
    hs=load(['hs_' fnum]);
    per=load(['per_' fnum]);
    wdir=load(['wdir_' fnum]);
    wdisbrk=load(['wdisbrk_' fnum]);
    wdisp=load(['wdisp_' fnum]);
    wdisfrc = load(['wdisfrc_' fnum]);
    wdiswcp = load(['wdiswcp_' fnum]);
    
    hs_prnt = ['hs',fnum];
    per_prnt = ['per',fnum];
    wdir_prnt = ['wdir',fnum];
    wdisp_prnt = ['wdisp',fnum];
    wdisbrk_prnt = ['wdisbrk',fnum];
    wdisfrc_prnt = ['wdisfrc',fnum];
    wdiswcp_prnt = ['wdiswcp',fnum];
    
    mask(depth<-1000)=1;
    hs(mask==1)=NaN;
    per(mask==1)=NaN;
    wdir(mask==1)=NaN;
    wdisbrk(mask==1)=NaN;
    wdisp(mask==1)=NaN;
    
    wdisfrc(mask==1)=NaN;
    wdiswcp(mask==1)=NaN;
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,hs,15);
    caxis([1.8 3]);colorbar
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    title('Significant Wave Height (m)','FontSize',12)
    print(hh,'-djpeg',hs_prnt)
    

end

    return 
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,per,15);
    caxis([0 10]);colorbar
    title('Peak Period (s)','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',per_prnt)
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,wdir,15);
    caxis([-40 40]);colorbar
    title('Peak Direction (\circ)','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',wdir_prnt)
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,wdisbrk,15);
    colorbar
    title('Wave Breaking Dissipation','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',wdisbrk_prnt)
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,wdisp,15);
    caxis([0 2.5e-4]);colorbar
    title('Total Dissipation','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',wdisp_prnt)
    
    
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,wdisfrc,15);
    caxis([0 2.5e-4]);colorbar
    title('Friction Dissipation','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',wdisfrc_prnt)
    
    
    hh=figure(1);
    clf;
    contourf(x_grid/1e3,y_grid/1e3,wdiswcp,15);
    caxis([0 5e-8]);colorbar
    title('Whitecapping Dissipation','FontSize',12)
    xlabel('Cross-shore (km)','FontSize',12)
    ylabel('Along-shore (km)','FontSize',12)
    print(hh,'-djpeg',wdiswcp_prnt)
    
%end