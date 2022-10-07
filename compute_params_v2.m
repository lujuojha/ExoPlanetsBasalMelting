function [k_i,rho_i,c_i,alpha_i] = compute_params_v2(PT,phase_s)
parfor i=1:length(phase_s)
    Phase = phase_s(i);  
    if Phase==0;
        hold=SeaFreeze(PT(i,:),'water1');
        k_i(i)=.60;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    elseif Phase==1;
        hold=SeaFreeze(PT(i,:),'Ih');
        k_i(i)= 2.4;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    elseif Phase==2;
        hold=SeaFreeze(PT(i,:),'II');
        k_i(i)= 1.8;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    elseif Phase==3;
        hold=SeaFreeze(PT(i,:),'III');
        k_i(i)= 1.1;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    elseif Phase==5;
        hold=SeaFreeze(PT(i,:),'V');
        k_i(i)=1.5;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    elseif Phase==6;
        hold=SeaFreeze(PT(i,:),'VI');
        k_i(i)= 1.9;
        rho_i(i)=hold.rho;
        c_i(i)=hold.Cp;
        alpha_i(i)=hold.alpha;
    else
    end
end
end