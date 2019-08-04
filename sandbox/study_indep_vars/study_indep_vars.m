clear;clc;close all




r_0_vec = (25:25:300)*1e-6; % radii in meters
DT_vec  = 1:5;       % initial temperature difference between the air and sea
% in degrees Celsius, air temp is always 18
U_vec   = 10:10:100;  % wind speed in m/s


%% ---------- Begin ------------------------------

[r_0,DT] = meshgrid(r_0_vec,DT_vec);

V = 4/3*pi*r_0.^3;
rho_s = 1000; % kg m^-3
g     = 9.81; % m/s
c_ps  = 4189.5; % J kg^-1 K ^-1
L_v   = 2257000; % J kg^-1
tau_T = 0.3;
tau_r = 300;
r_eq = 2/3*r_0;
S0 = 34;

i = 1;
t_flight_vec = [.1];
for t_flight_ind = 1:length(t_flight_vec)
    for U = [10:10:100]
         t_flight = 2*U/g*t_flight_vec(t_flight_ind);
        

        
        Tf  = (-DT)*exp(-t_flight/tau_T)+DT;
        
        
        figure(3)
        rf  = (r_0-r_eq)*exp(-t_flight/tau_r)+r_eq;
        plot(r_0(1,:),rf(1,:))
        hold on
        %% Mechanical Energy ------------------------------
        KE    = 0.5.*rho_s.*V.*U.^2;
        PE = rho_s.*V.*g.*(U.^2./g);
        
        %% Thermal Energy ---------------------------------
        
        Q_s = rho_s.*V.*c_ps.*Tf;
        Q_L = rho_s.*V.*L_v.*(1 - (rf./r_0).^3);
        
        %%
        
        net_E(:,:,i) = max(Q_s-Q_L,0*zeros(size(Q_s))) - (KE+PE);
        
        figure(1)
        subplot(2,2,1)
        surf(r_0,DT,PE)
        title('PE')
        hold on
        
        subplot(2,2,2)
        surf(r_0,DT,KE)
        title('KE')
        hold on
        
        figure(4)
        plot(r_0(1,:),net_E(4,:,i))
        hold on
        
        figure(10)
        subplot(1,length(t_flight_vec),t_flight_ind)
        surf(r_0,DT,net_E(:,:,i).*(net_E(:,:,i)>0),U*ones(size(net_E(:,:,i).*(net_E(:,:,i)>0))))
        colormap(colorcube)
%         title(sprintf('net for t_f = %f ',t_flight))
        xlabel('r_0')
        ylabel('\Delta T')
        title('net E [J]')
           shading interp
        set(gca,'fontsize',18,'view',[40.900000000000006 26])
        colorbar
        caxis([10 100])
        set(gcf,'color','w','position',[597   232   832   573])
        i = i + 1;
        hold on
    end
    figure(1)
    
    subplot(2,2,3)
    surf(r_0,DT,Q_s)
    title('Q_s')

    subplot(2,2,4)
    surf(r_0,DT,Q_L)
    title('Q_L')
    for i = 1:4
        subplot(2,2,i)
        set(gca,'fontsize',18)
    end
    set(gcf,'color','w')
end
    
    
    
    %{
    addpath ~/Documents/MATLAB/util/export_fig/
 export_fig(gcf,'individ_E_comp.pdf')
 copyfile('individ_E_comp.pdf','~/MIT/Research/EmanuelGroup/GAMEStudentSeminar_15112017/')

    
 export_fig(gcf,'individ_E_comp.pdf')
 copyfile('individ_E_comp.pdf','~/MIT/Research/EmanuelGroup/GAMEStudentSeminar_15112017/')


 export_fig(gcf,sprintf('net_E_U%03d.pdf',U))
 copyfile(sprintf('net_E_U%03d.pdf',U),'~/MIT/Research/EmanuelGroup/GAMEStudentSeminar_15112017/')

%}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
