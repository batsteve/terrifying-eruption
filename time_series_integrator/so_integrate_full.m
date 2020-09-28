function [NN, KK, ic, ind3] = so_integrate_full(param, DS, f_env, kappa)
% param -- struct containing certain scalar parameters
% DS -- struct containing the load signal
% f_env -- functional description of the coherent envelope
% kappa -- function related to ascending leg intersections

Zmax=0;
NN = -1;

t_max = 1*10^5;

% fatigue-crack parameters
dc=param.dc; sc=param.sc; da=param.da;

KK=Inf;
i=1;ic=0;
ind1=-1; % initial index to make sure we begin with positive load
ind2=0; %termination index indicating fracture
ind3=-1; % termination index for special failure

options = optimset('Display','none');

% problems with initialization?
KKu = inf;


dd = zeros(t_max, 1);
ddm = zeros(t_max, 1);
KK = zeros(t_max, 1);

while (i<t_max) && (ind2==0)
         
    if (i == 1)
        Dsp = 0;
    else
        Dsp = DS.sspp(i) - max(DS.ssnp(i-1), 0);
    end
    Zmax = DS.sspp(i);
    
    if (i == 1)
        dd(i)=lsqnonlin(@(d00)f_env(d00)-Zmax,exp(1)*sc/dc,0,dc,options);
        KK(i) = kappa(Zmax);
    else    % this is a normal loading cycle
        KK(i)=KKu-Dsp/da;
        dd(i)=ddm(i-1)-da*log(1-Dsp/da/KKu);
    end

    if (Zmax > f_env(dd(i)))
        if (dd(i)<dc) && (Zmax<sc)
            %dd(i)=lsqnonlin(@(d00)coh_env(d00,dc,sc)-Zmax,exp(1)*sc/dc,0,dc,options);
            dd(i)=lsqnonlin(@(d00)f_env(d00)-Zmax,exp(1)*sc/dc,0,dc,options);
            %KK(i)=Zmax/dd(i);
            KK(i) = kappa(Zmax);
            ic=ic+1; % number of catastrophic events for fracture / crossing of curve in the ascending part
        elseif (dd(i)>dc) && (Zmax<sc)
            ind2=1; %termination index
            NN=i; %number of cycles
            ind3=0; % catastrophic failure without exceedence of sc
        else
            ind2=1; %termination index
            NN=i; %number of cycles
            ind3=1; % exceedenve of sc
        end
        
    end
    
    % min of cycle
    

    Dsm = max(DS.ssnp(i), 0) -DS.sspp(i);
    
    KKm=DS.sspp(i)/dd(i); % Zmax still refers to the local max of the load
    ddm(i)=max(dd(i)+Dsm/KKm,0); % Zmax to make sure we don't get negative ddm
    KKu=KKm-exp(Dsm/da/KKm)*(KKm-KK(i));    
    
    %fprintf('KK = %0.2f.\n', KK(i));
    i = i+1;
end


