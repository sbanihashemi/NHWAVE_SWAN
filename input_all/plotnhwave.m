%This file is used to plot 3D output data from NHWAVE
clc;clear

% model setup 
size_x = 64;
size_y = 16;
size_z = 20;

dx= 14300/64;
dy= 1000;

x = [1:size_x]*dx;
y = [1:size_y]*dy;

% load depth 
h=load(['depth']); 

f_number =5;
fnum = sprintf('%.4d',f_number);

% eta, u, v, w
eta=load(['eta_' fnum]);
u=load(['u_' fnum]);
v=load(['v_' fnum]);
w=load(['w_' fnum]);

i=load(['i_' fnum]);
rho=load(['rho_' fnum]);

etaw=load(['etaw_' fnum]);


% %Stokes drift velocity
% 
ustk=load(['ustk_' fnum]);
vstk=load(['vstk_' fnum]);
wstk=load(['wstk_' fnum]);
% 
% %Lagrangian velocity
% %  U=u+ustk;
% %  V=v+vstk;
% %  W=w+wstk;
% 
% %Bernoulli head and gradient
% 
bhx=load(['bhx_' fnum]);
bhy=load(['bhy_' fnum]);
bhz=load(['bhz_' fnum]);

% wave pressure gradient
wpx=load(['wpx_' fnum]);
wpy=load(['wpy_' fnum]);
wpz=load(['wpz_' fnum]);

% 
% % %Vortex force 
vfx=load(['vortforcX_' fnum]);
vfy=load(['vortforcY_' fnum]);
vfz=load(['vortforcZ_' fnum]);
% 
% % Dissipation force
dissx = load(['dissx_' fnum]);
dissy = load(['dissy_' fnum]);
% 
% % Eddy viscosity 
% eddy= load(['c_' fnum]);

% total water depth
D=eta+h;
% construct z-grid, the relation is linear, we take uniform dz
% sigma = z+h/D
dz=D/size_z;

% For z-x view
% Specify location in y-direction for view
y_loc = 9;
dz_loc = dz(y_loc,:);

%NHWAVE output follows the rules:
%1. Each line is a vector in x-direction
%2. Second line follows in y-direction
%3. Z-direction is the slowest, from bottom (Kbeg) to surface (Kend)
%return

for i_z=1:size_z;
    for i_x=1:size_x;
        %have to jump size_y steps to get another z-level value
        u1(i_z,i_x)=u((i_z-1)*size_y+y_loc,i_x);
        v1(i_z,i_x)=v((i_z-1)*size_y+y_loc,i_x);
        w1(i_z,i_x)=w((i_z-1)*size_y+y_loc,i_x);
        
        i1(i_z,i_x)=i((i_z-1)*size_y+y_loc,i_x);
        rho1(i_z,i_x)=rho((i_z-1)*size_y+y_loc,i_x);
        
        ustk1(i_z,i_x)=ustk((i_z-1)*size_y+y_loc,i_x);
        vstk1(i_z,i_x)=vstk((i_z-1)*size_y+y_loc,i_x);
        wstk1(i_z,i_x)=wstk((i_z-1)*size_y+y_loc,i_x);
        
        bhx1(i_z,i_x)=bhx((i_z-1)*size_y+y_loc,i_x);
        bhy1(i_z,i_x)=bhy((i_z-1)*size_y+y_loc,i_x);
        bhz1(i_z,i_x)=bhz((i_z-1)*size_y+y_loc,i_x);
        
        wpx1(i_z,i_x)=wpx((i_z-1)*size_y+y_loc,i_x);
        wpy1(i_z,i_x)=wpy((i_z-1)*size_y+y_loc,i_x);
        wpz1(i_z,i_x)=wpz((i_z-1)*size_y+y_loc,i_x);
         
        vfx1(i_z,i_x)=vfx((i_z-1)*size_y+y_loc,i_x);
        vfy1(i_z,i_x)=vfy((i_z-1)*size_y+y_loc,i_x);
        vfz1(i_z,i_x)=vfz((i_z-1)*size_y+y_loc,i_x);
        
%         U1(i_z,i_x)=U((i_z-1)*size_y+y_loc,i_x);
%         V1(i_z,i_x)=V((i_z-1)*size_y+y_loc,i_x);
%         W1(i_z,i_x)=W((i_z-1)*size_y+y_loc,i_x);
        
        dissx1(i_z,i_x)=dissx((i_z-1)*size_y+y_loc,i_x);
        dissy1(i_z,i_x)=dissy((i_z-1)*size_y+y_loc,i_x);
%         
%         eddy1(i_z,i_x)=eddy((i_z-1)*size_y+y_loc,i_x);
          
    end;
end

for i_x=1:size_x;
    % wave dissipation
    dissx_flux(i_x) = sum(dissx1(:,i_x))*dz_loc(i_x);
    dissy_flux(i_x) = sum(dissy1(:,i_x))*dz_loc(i_x);
    
    % vortex force
    vfx_flux(i_x) = sum(vfx1(:,i_x))*dz_loc(i_x);
    vfy_flux(i_x) = sum(vfy1(:,i_x))*dz_loc(i_x);
    
    % current flux
    utrspt(i_x)=sum(u1(:,i_x))*dz_loc(i_x);
    vtrspt(i_x)=sum(v1(:,i_x))*dz_loc(i_x);
    wtrspt(i_x)=sum(w1(:,i_x))*dz_loc(i_x);
    
    u_bar(i_x)=utrspt(i_x)/D(y_loc,i_x);
    v_bar(i_x)=vtrspt(i_x)/D(y_loc,i_x);
    
%     % stokes flux
    ustk_trspt(i_x)=sum(ustk1(:,i_x))*dz_loc(i_x);
    vstk_trspt(i_x)=sum(vstk1(:,i_x))*dz_loc(i_x);
    wstk_trspt(i_x)=sum(wstk1(:,i_x))*dz_loc(i_x);
    ustk_bar(i_x)=ustk_trspt(i_x)/D(y_loc,i_x);
    
end

%return
for i=1:size_x;
    zz(:,i) = [-h(y_loc,i)+dz_loc(i):dz_loc(i):eta(y_loc,i)];
%     zz(:,i)=-dz_loc(i)*[size_z:-1:1];
end
zzint=[min(min(zz)):0.1:max(max(zz))];
len_int=length(zzint);
for i=1:size_x;
    u1int(:,i)=interp1(zz(:,i),u1(:,i),zzint,'spline');
    v1int(:,i)=interp1(zz(:,i),v1(:,i),zzint,'spline');
    w1int(:,i)=interp1(zz(:,i),w1(:,i),zzint,'spline');
    
    i1int(:,i)=interp1(zz(:,i),i1(:,i),zzint,'spline');
    rho1int(:,i)=interp1(zz(:,i),rho1(:,i),zzint,'spline');
    
    bhx1int(:,i)=interp1(zz(:,i),bhx1(:,i),zzint,'spline');
    bhy1int(:,i)=interp1(zz(:,i),bhy1(:,i),zzint,'spline');
    bhz1int(:,i)=interp1(zz(:,i),bhz1(:,i),zzint,'spline');
    
    wpx1int(:,i)=interp1(zz(:,i),wpx1(:,i),zzint,'spline');
    wpy1int(:,i)=interp1(zz(:,i),wpy1(:,i),zzint,'spline');
    wpz1int(:,i)=interp1(zz(:,i),wpz1(:,i),zzint,'spline');    
    
    vfx1int(:,i)=interp1(zz(:,i),vfx1(:,i),zzint,'spline');
    vfy1int(:,i)=interp1(zz(:,i),vfy1(:,i),zzint,'spline');
    vfz1int(:,i)=interp1(zz(:,i),vfz1(:,i),zzint,'spline');
    
    ustk1int(:,i)=interp1(zz(:,i),ustk1(:,i),zzint,'spline');
    vstk1int(:,i)=interp1(zz(:,i),vstk1(:,i),zzint,'spline');
    wstk1int(:,i)=interp1(zz(:,i),wstk1(:,i),zzint,'spline');
    
%     U1int(:,i)=interp1(zz(:,i),U1(:,i),zzint,'spline');
%     V1int(:,i)=interp1(zz(:,i),V1(:,i),zzint,'spline');
%     W1int(:,i)=interp1(zz(:,i),W1(:,i),zzint,'spline');
     
    dissx1int(:,i)=interp1(zz(:,i),dissx1(:,i),zzint,'spline');
    dissy1int(:,i)=interp1(zz(:,i),dissy1(:,i),zzint,'spline');
%     
%     eddy1int(:,i)=interp1(zz(:,i),eddy1(:,i),zzint,'v5cubic');
end
%return
for i=1:size_x;
    if(h(y_loc,i)>0);% water region
    for j=1:len_int;
        if zzint(j)>-h(y_loc,i);
            bot_loc(i)= j;
            u1int(1:j-1,i)=NaN;
            v1int(1:j-1,i)=NaN;
            w1int(1:j-1,i)=NaN;
            
            i1int(1:j-1,i)=NaN;
            rho1int(1:j-1,i)=NaN;
                  
%              U1int(1:j-1,i)=NaN;
%              V1int(1:j-1,i)=NaN;
%              W1int(1:j-1,i)=NaN;
            
            bhx1int(1:j-1,i)=NaN;
            bhy1int(1:j-1,i)=NaN;
            bhz1int(1:j-1,i)=NaN;
          
            wpx1int(1:j-1,i)=NaN;
            wpy1int(1:j-1,i)=NaN;
            wpz1int(1:j-1,i)=NaN;            
             
            vfx1int(1:j-1,i)=NaN;
            vfy1int(1:j-1,i)=NaN;
            vfz1int(1:j-1,i)=NaN;
            
            ustk1int(1:j-1,i)=NaN;
            vstk1int(1:j-1,i)=NaN;
            wstk1int(1:j-1,i)=NaN;

            dissx1int(1:j-1,i)=NaN;
            dissy1int(1:j-1,i)=NaN;
%             
%             eddy1int(1:j-1,i)=NaN;
            break;
        end;
    end;
    else
            bot_loc(i)= 120;
            u1int(:,i)=NaN;
            v1int(:,i)=NaN;
            w1int(:,i)=NaN;
            
            i1int(:,i)=NaN;
            rho1int(:,i)=NaN;
                                
%              U1int(:,i)=NaN;
%              V1int(:,i)=NaN;
%              W1int(:,i)=NaN;
            
            bhx1int(:,i)=NaN;
            bhy1int(:,i)=NaN;
            bhz1int(:,i)=NaN;
            
            wpx1int(:,i)=NaN;
            wpy1int(:,i)=NaN;
            wpz1int(:,i)=NaN;
            
            vfx1int(:,i)=NaN;
            vfy1int(:,i)=NaN;
            vfz1int(:,i)=NaN;
            
            ustk1int(:,i)=NaN;
            vstk1int(:,i)=NaN;
            wstk1int(:,i)=NaN;

            dissx1int(:,i)=NaN;
            dissy1int(:,i)=NaN;
%             
%             eddy1int(:,i)=NaN;
    end
    
end
return
% cross-shore profile

figure;
plot(x/1000,-h(y_loc,:),'LineWidth',2)
xlabel('Cross-shore (km)','FontS',12)
ylabel('Depth (m)','FontS',12)
title('Cross-shore profile','FontS',12)
hold on
plot(x/1000,eta(y_loc,:),'LineWidth',2)
axis([0 14.8 -100 5]);
for i=1:12
plot(x(i*5)/1000+u1int(end:-1:bot_loc(5*i),5*i)*0.5,zzint(end:-1:bot_loc(5*i)),'LineWidth',2)
plot(x(i*5)/1000+u1int(end:-1:bot_loc(5*i),5*i)*0,zzint(end:-1:bot_loc(5*i)),'k--','LineWidth',2)
end
umax=num2str(round(max(max(abs(u1)))*100)/100);
leg_name = sprintf('%s%s','u_{max}=',umax)
hh = legend(leg_name)
set(hh,'FontS',12)
%print -djpeg -r350 cross_shore_profileU_nwsw

% plot eta
figure
%contourf(x/1000,y/1000,eta)
pcolor(x/1000,y/1000,eta)
xlabel('Cross-shore (km)','FontS',12)
ylabel('Along-shore (km)','FontS',12)
title('Mean water level','FontS',12)
colorbar
%print -djpeg -r350 mean_water_level

%salinity
figure;
pcolor(x/1000,zzint,i1int);caxis([0 30]);colorbar;shading interp
title('Salinity');
xlabel('x (km)','FontS',12)
ylabel('z (m)','FontS',12)
print -djpeg -r350 salinity

figure;
pcolor(x/1000,zzint,rho1int);caxis([995 1035]);colorbar;shading interp
title('Density');
xlabel('x (km)','FontS',12)
ylabel('z (m)','FontS',12)
print -djpeg -r350 density

figure;
subplot(3,1,1)
%contourf(x,zzint,u1int);colorbar
pcolor(x/1000,zzint,u1int);caxis([-3.0 0]);colorbar;shading interp
title('u');
ylabel('z (m)','FontS',12)
subplot(3,1,2)
%contourf(x,zzint,v1int);colorbar
pcolor(x/1000,zzint,v1int);caxis([-1e-4 7e-4]);colorbar;shading interp
title('v');
ylabel('z (m)','FontS',12)
subplot(3,1,3)
%contourf(x,zzint,w1int);colorbar
pcolor(x/1000,zzint,w1int);caxis([-7e-3 4e-3]);colorbar;shading interp
title('w');
xlabel('x (km)','FontS',12)
ylabel('z (m)','FontS',12)
print -djpeg -r350 uvw_contour_nwsw

return 

% plot waveset down
figure
contourf(x/1000,y/1000,etaw)
xlabel('Cross-shore (km)','FontS',12)
ylabel('Along-shore (km)','FontS',12)
title('Wave setdown (m)','FontS',12)
colorbar
print -djpeg -r350 mean_water_level

return
% mean current
figure;
subplot(3,1,1)
plot(x,utrspt,'LineWidth',2)
hold on
plot(x,ustk_trspt,'r--','LineWidth',2)
legend('Current','Stokes drift')
% longshore current
subplot(3,1,2)
%v_bar(50)=0;
plot(x/1000,v_bar,'LineWidth',2)
xlabel('cross-shore (km)')
ylabel('Longshore current (m/s)')
title('Longshore current')
% wave setup
subplot(3,1,3)
plot(x(1:49),eta(1,1:49)-eta(1,1)+wave_setup(1,1),'--','LineWidth',2)
hold on
plot(x(1:49),wave_setup(1,1:49),'r--','LineWidth',2)
legend('model output','analytical')
title('Mean surface level')
return
print -djpeg -r350 u_v_eta

%return

% eddy viscosity 
figure;
contourf(x,zzint,eddy1int);colorbar;
xlabel('x (m)')
ylabel('z (m)')
print -djpeg -r350 eddy_viscosity
% 2D Eulerian current velocity



% Stokes drift velocity

figure;subplot(3,1,1)
%contourf(x,zzint,ustk1int);caxis([0.0015 0.025]);colorbar
pcolor(x/1000,zzint,ustk1int);caxis([0.0015 0.025]);colorbar;shading interp
title('u');
ylabel('z (m)','FontS',12)
subplot(3,1,2)
%contourf(x,zzint,vstk1int);caxis([0.0e-4 2.5e-4]);colorbar
pcolor(x/1000,zzint,vstk1int);caxis([0.0e-4 2.5e-4]);colorbar;shading interp
title('v');
ylabel('z (m)','FontS',12)
subplot(3,1,3)
%contourf(x,zzint,wstk1int);caxis([-2.5e-4 2e-4]);colorbar
pcolor(x/1000,zzint,wstk1int);caxis([-2.5e-4 2e-4]);colorbar;shading interp
title('w');
ylabel('z (m)','FontS',12)
xlabel('x (km)','FontS',12)

print -djpeg -r350 uvw_stokes_nwsw
% Lagrangian velocity

% figure;subplot(3,1,1)
% contourf(x,zzint,U1int);colorbar
% subplot(3,1,2)
% contourf(x,zzint,V1int);colorbar
% subplot(3,1,3)
% contourf(x,zzint,W1int);colorbar

% dissipation

figure;subplot(2,1,1)
contourf(x,zzint,dissx1int);colorbar
ylabel('z (m)','FontS',12)
subplot(2,1,2)
contourf(x,zzint,dissy1int);colorbar
xlabel('x (m)','FontS',12)
ylabel('z (m)','FontS',12)
print -djpeg -r350 diss_xy
% vortex force

figure;subplot(3,1,1)
%contourf(x,zzint,vfx1int);colorbar
pcolor(x/1000,zzint,vfx1int);caxis([-1e-4 1e-4]);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('Vortex Force X')
subplot(3,1,2)
%contourf(x,zzint,vfy1int);colorbar
pcolor(x/1000,zzint,vfy1int);caxis([-1e-4 1e-4]);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('Vortex Force Y')
subplot(3,1,3)
%contourf(x,zzint,vfz1int);colorbar;
pcolor(x/1000,zzint,vfz1int);caxis([-1.5e-3 3.5e-3]);colorbar;shading interp
xlabel('x (km)','FontS',12)
ylabel('z (m)','FontS',12)
title('Vortex Force Z')
print -djpeg -r350 vort_xyz_nwsw
% bernoulli head gradient

figure;subplot(3,1,1)
%contourf(x,zzint,bhx1int);caxis([-4e-5 3e-5]);colorbar
pcolor(x/1000,zzint,bhx1int);caxis([-4e-5 3e-5]);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('B-H Gradient X')
subplot(3,1,2)
%contourf(x,zzint,bhy1int);caxis([-1e-4 1e-4]);colorbar
pcolor(x/1000,zzint,bhy1int);caxis([-1e-4 1e-4]);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('B-H Gradient Y')
subplot(3,1,3)
%contourf(x,zzint,bhz1int);caxis([-1e-4 1e-4]);colorbar
pcolor(x/1000,zzint,bhz1int);caxis([-1e-4 1e-4]);colorbar;shading interp
xlabel('x (km)','FontS',12)
ylabel('z (m)','FontS',12)
title('B-H Gradient Z')
print -djpeg -r350 bh_xyz_nwsw

figure;subplot(3,1,1)
%contourf(x,zzint,bhx1int);caxis([-4e-5 3e-5]);colorbar
pcolor(x,zzint,wpx1int);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('Wave Pressure Gradient X')
subplot(3,1,2)
%contourf(x,zzint,bhy1int);caxis([-1e-4 1e-4]);colorbar
pcolor(x,zzint,wpy1int);colorbar;shading interp
ylabel('z (m)','FontS',12)
title('Wave Pressure Gradient Y')
subplot(3,1,3)
%contourf(x,zzint,bhz1int);caxis([-1e-4 1e-4]);colorbar
pcolor(x,zzint,wpz1int);colorbar;shading interp
xlabel('x (m)','FontS',12)
ylabel('z (m)','FontS',12)
title('Wave Pressure Gradient Z')
print -djpeg -r350 wp_xyz_nwsw
