function [T_np1,Conv_top,Conv_bottom,Ra,delta_thick_top,delta_thick_bottom,is_convect,k_bar,Nu]=HP_Ice_Evolve_v8(T_n,k_i,rho_i,c_i,...
    dt,dz,T_surf,Base_flux,Tm,TTol,P_n,phasenew,Ra_c,eta_0,alpha,g_s,A,Mx,Rx)

% Improved function that takes the original vertical profiles of ice/water
% temperature (T_n), thermal conductivity (k_i), density
% (rho_i), heat capacity (c_i), and melting temperature (Tm) alongside time
% step (dt), vertical resolution (dz - in meters), surface temperature (in
% Kelvin), basal heat flux (in W/m^2), and error tolerances for the finite
% difference iterator (TTol) and calculates the thermal evolution of the 
% ice sheet.

% This function does thermal diffusion and phase change (via SeaFreeze) and
% parameterizes 1D convection using Ra/Nu # scaling

% This version also contains simulated convection in regions which exceed a
% predefined critical Rayleigh # (Ra_c).

matrix_dimension=length(T_n);
counter=0;

T_np1_km1=T_n;
%phi_np1_km2=Phi_n;
T_evolve=[];
%phi_evolve=[];
TErr=999;
%PhiErr=999;

%% OUTSIDE OF ITERATIVE LOOP - parameterized convective mixing of any regions capable of convecting
phase_top=phasenew(1);
top_hold=[1];               % top of ice phase boundaries
bottom_hold=[];             % bottom of ice phase boundaries
Ra=[];                      % Rayleigh #s
Nu=[];                      % Nusselt #s
eta=[];                     % viscosities
Conv_top=[];                % top of convection zones
Conv_bottom=[];             % bottom of convection zones
delta_thick_top=[];         % upper convective boundary layer thickness
delta_thick_bottom=[];      % lower convective boundary layer thickness
is_convect=[];              % binary matrix to tell if layer is convecting

% finding the phase change interfaces
for i=2:matrix_dimension
    phase=phasenew(i);
    if phase==phase_top
        depth=i;
    else
        top_hold=[top_hold i];
        bottom_hold=[bottom_hold i-1];
        phase_top=phase;
    end
end

bottom_hold=[bottom_hold matrix_dimension];

% calculating the viscosity, Rayleigh #, and Nusselt # profile
for i=1:length(top_hold)
    for j=top_hold(i):bottom_hold(i)
        eta(j)=eta_0(j)*exp(A*((Tm(j)/T_np1_km1(j))-1));
        Ra(j)=(alpha(j)*(rho_i(j)^2)*c_i(j)*g_s*(T_np1_km1(bottom_hold(i))-T_np1_km1(j))*((bottom_hold(i)-j)*dz)^3)/(k_i(j)*eta(j));
        if Ra(j)<0
            Ra(j)=0;
        else
        end
        Nu(j)=(Ra(j)/Ra_c)^(1/3);
    end
end

% checking to see if the different ice/water layers are convecting
for i=1:length(top_hold)
    for j=top_hold(i):bottom_hold(i)
        if Ra(j)>Ra_c
%             [T_ad,dT_dz_ad,P_ad,rho_ad,alpha_ad,Cp_ad,K_ad,phase_ad,g_s_ad,z_ad] = ...
%                 adiabat_profile(P_n(j),T_np1_km1(j),dz*(bottom_hold(i)-j+1),dz/1000,Mx,Rx);
            Conv_top=[Conv_top j];
            Conv_bottom=[Conv_bottom bottom_hold(i)];
%             Tprime=T_np1_km1';
%             En_hold=c_i(j:bottom_hold(i)).*Tprime(j:bottom_hold(i));%+L(j:bottom_hold(i)).*phi_np1_km1(j:bottom_hold(i));
%             En_ad=c_i(j:bottom_hold(i)).*T_ad;%+L(j:bottom_hold(i)).*phi_np1_km1(j:bottom_hold(i));
%             T_np1_km1(j:bottom_hold(i))=T_ad;
%             En_diff=sum(En_hold)-sum(En_ad);
%             En_cell=En_diff/(bottom_hold(i)-j+1);
%             T_cell=En_cell./c_i(j:bottom_hold(i));
%             T_np1_km1(j:bottom_hold(i))=T_np1_km1(j:bottom_hold(i))+T_cell';
            is_convect(i)=1;
            break
        else
            is_convect(i)=0;
        end
    end
end

% calculating convection boundary layer thicknesses (conductive regions)
for i=1:length(top_hold)
    if is_convect(i)>0 && Ra(i)>0
        delta_thick_top(i)=round((bottom_hold(i)-top_hold(i))*(Ra_c/Ra(top_hold(i)))^(1/3));
        delta_thick_bottom(i)=round((bottom_hold(i)-top_hold(i))*(0.28*Ra(top_hold(i))^-0.79)^(1/3)); % from Deschamps and Sotin (2001)
    else
        delta_thick_top(i)=0;
        delta_thick_bottom(i)=0;
    end
end

%% Thermal conduction iterations
while  TErr>TTol %|| PhiErr>PhiTol
    counter=counter+1;
    if counter>500
        counter
    else
    end

% %% Enthalpy Method
% Hs=c_i.*Tm;     % enthalpy of solid ice to compare in enthalpy method
% 
% % Value of Phi(n+1,k-1,j)
% phi_np1_km1=phi_np1_km2;
% % Latent Heat interpolation using SeaFreeze data stored in
% % Latent_Heat_Lookup.mat
% L=-interp1(Latent_Heat_Lookup(:,1),Latent_Heat_Lookup(:,2),P_n);
%     for i=1:matrix_dimension
% % Test to see if melting ice or freezing brine
%     En=c_i(i)*T_np1_km1(i)+L(i)*phi_np1_km2(i);
%         if En<Hs(i)
%             phi_np1_km1(i)=0;
%         elseif En>Hs(i)+L(i)
%             phi_np1_km1(i)=1;
%         else
%             phi_np1_km1(i)=((c_i(i)*T_np1_km1(i)+phi_np1_km2(i)*L(i))-Hs(i))/...
%             L(i);
%         end
%     end
% 
% % Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
% % iteration
% phi_np1_km2=phi_np1_km1;


%% Volume averaging conductivity, density*specific heat
% modifying thermal conductivities to account for convection and conduction
% boundary layers
for i=1:length(top_hold)
    for j=top_hold(i):bottom_hold(i)
        % if layer isn't convecting
        if is_convect(i)==0
            k_bar(j)=k_i(j);
        % if in conducting boundary layers
        elseif j-top_hold(i)+1<=delta_thick_top(i) || bottom_hold(i)-j+1<=delta_thick_bottom(i)
            k_bar(j)=k_i(j);
        % if the boundary layers take up the whole convecting region
        elseif bottom_hold(i)-top_hold(i)<=delta_thick_top(i)+delta_thick_bottom(i)
            k_bar(j)=k_i(j);
        % if layer is water
        elseif phasenew(j)==0
            k_bar(j)=k_i(j);
        else
            %k_bar(j)=max(Nu(top_hold(i)+delta_thick_top(i):bottom_hold(i)-delta_thick_bottom(i)))*k_i(j);
            k_bar(j)=max(Nu(top_hold(i):bottom_hold(i)))*k_i(j);
        end
    end
end

rho_c_bar=rho_i.*c_i;

% Creating shifted thermal conductivity vectors and up/down conductivity
% averages
k_jp1=circshift(k_bar,1);
k_jp1(1)=k_bar(1);
k_jm1=circshift(k_bar,-1);
k_jm1(end)=k_bar(end);

k_up_avg=(k_jp1+k_bar)./2;
k_down_avg=(k_jm1+k_bar)./2;

% Constant coefficient for solution matrix
C_0=dt./((dz^2)*rho_c_bar);

% Build diagonal vectors of solution matrix to put in Thomas Tridiagonal
% algorithm
a=zeros(1,length(T_np1_km1));
b=zeros(1,length(T_np1_km1));
c=zeros(1,length(T_np1_km1));

for i=1:length(T_np1_km1)
    a(i)=1+C_0(i)*(k_up_avg(i)+k_down_avg(i));
    b(i)=-C_0(i)*k_down_avg(i);
    c(i)=-C_0(i)*k_up_avg(i);
end

% To account for basal heat flux discretization
%c(end)=-2*C_0(end)*k_up_avg(end);
a(end)=1+C_0(i)*k_up_avg(i);

b(end)=[];
c(1)=[];

% Temp profile from previous time step (used in Thom_Trid algorithm) and
% inclusion of temperature boundary conditions
y=zeros(1,length(T_np1_km1));
y(1:end)=T_n;

% surface BC (surface temp)
y(1)=y(1)+C_0(1)*k_up_avg(1)*T_surf;
% basal BC (basal heat flux)
%y(end)=y(end)+2*C_0(end)*(Base_flux*dz);
y(end)=y(end)+C_0(end)*(Base_flux*dz);
%y(end)=y(end)+C_0(end)*k_down_avg(end)*(T_np1_km1(end)+Base_flux*dz/k_down_avg(end));

% Solve for new Temperature profile
T_np1_km1=Thomas_Trid(a,b,c,y)';

    %% Appending value to matrix to check for convergence
    T_evolve=[T_evolve T_np1_km1'];
    %phi_evolve=[phi_evolve phi_np1_km1'];
    
    if counter==1
        TErr=999;
        %PhiErr=999;
    else
        TErr=max(abs(T_evolve(:,counter)-T_evolve(:,counter-1)));
        %PhiErr=max(abs(phi_evolve(:,counter)-phi_evolve(:,counter-1)));
    end
end

%% Adiabatic temperature profile in any water layers
for i=1:length(top_hold)
        if phasenew(top_hold(i))==0 && Ra(top_hold(i))>Ra_c
            j=top_hold(i);
            [T_ad,dT_dz_ad,P_ad,rho_ad,alpha_ad,Cp_ad,K_ad,phase_ad,g_s_ad,z_ad] = ...
                adiabat_profile(P_n(j),T_np1_km1(j),dz*(bottom_hold(i)-j+1),dz/1000,Mx,Rx);
            Conv_top=[Conv_top j];
            Conv_bottom=[Conv_bottom bottom_hold(i)];
            Tprime=T_np1_km1';
            En_hold=c_i(j:bottom_hold(i)).*Tprime(j:bottom_hold(i));%+L(j:bottom_hold(i)).*phi_np1_km1(j:bottom_hold(i));
            En_ad=c_i(j:bottom_hold(i)).*T_ad;%+L(j:bottom_hold(i)).*phi_np1_km1(j:bottom_hold(i));
            T_np1_km1(j:bottom_hold(i))=T_ad;
            En_diff=sum(En_hold)-sum(En_ad);
            En_cell=En_diff/(bottom_hold(i)-j+1);
            T_cell=En_cell./c_i(j:bottom_hold(i));
            T_np1_km1(j:bottom_hold(i))=T_np1_km1(j:bottom_hold(i))+T_cell;
        else
        end
end

%% Final temperature profile
T_np1=T_np1_km1';


