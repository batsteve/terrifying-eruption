function [NN, KK, ic, ind3] = so_decomp_synth(param, DS, f_env, kappa)
%% 
% a lot of initialization stuff here is leftover from the full SO
%  integration function, which this evolved out of.
%
% param -- struct containing certain scalar parameters
% DS -- struct containing the load signal
% f_env -- functional description of the coherent envelope
% kappa -- function related to ascending leg intersections

dc=param.dc; sc=param.sc; da=param.da;

%dd=0;
%KK=Inf;
%N=t_span_length;
%i=1;
ic=0; % j=0;
%ind1=-1; % initial index to make sure we begin with positive load
ind2=0; %termination index indicating fracture
%dcv=[0:dc/100:8*dc];
options = optimset('Display','none');

% problems with initialization?
KKu = inf;

K_slope_mu = param.K_slope_mu;


extreme_height_threshold = param.extreme_height_threshold;
[PP, BB] = findpeaks(DS.sspp, 'MinPeakHeight', extreme_height_threshold);
if isempty(BB)
    BB = [1];
    PP = [DS.sspp(1)];
elseif (BB(1) ~= 1)
    BB = [1, BB];
    PP = [DS.sspp(1), PP];
end

BB = [BB, length(DS.sspp)];
PP = [PP, DS.sspp(end)];

p_max = length(PP);
fprintf('Handling %d peaks specially.\n', p_max);
cp = 1;


dd = zeros(p_max, 1);
ddm = zeros(p_max, 1);
KK = zeros(p_max, 1);

while (ind2==0)
    
    i = BB(cp);
    
    % max of cycle
     
    if (i == 1)
        Dsp = 0;
    else
        Dsp = DS.sspp(i) - max(DS.ssnp(i-1), 0);
    end
    Zmax = DS.sspp(i);
    
    %if KK==Inf % this is the case during the first cycle
    if (i == 1)
        %dd(cp)=lsqnonlin(@(d00)coh_env(d00,dc,sc)-Zmax,exp(1)*sc/dc,0,dc,options);
        dd(cp)=lsqnonlin(@(d00)f_env(d00)-Zmax,exp(1)*sc/dc,0,dc,options);
        %KK(cp)=Zmax/dd(cp);
        KK(cp) = kappa(Zmax);
    else    % this is a normal loading cycle
        ddmm = 0;
        KK(cp)=KKu-Dsp/da;
        dd(cp)=ddmm-da*log(1-Dsp/da/KKu);
    end

    if (Zmax > f_env(dd(cp)))
        if (dd(cp)<dc) && (Zmax<sc)
            dd(cp)=lsqnonlin(@(d00)f_env(d00)-Zmax,exp(1)*sc/dc,0,dc,options);
            KK(cp) = min(kappa(Zmax),  KK(cp)); % change to handle weird failure cases 
            ic=ic+1; % number of catastrophic events for fracture / crossing of curve in the ascending part
        elseif (dd(cp)>dc) && (Zmax<sc)
            ind2=1; %termination index
            NN=i; %number of cycles
            ind3=0; % catastrophic failure without exceedence of sc
        else
            ind2=1; %termination index
            NN=i; %number of cycles
            ind3=1; % exceedenve of sc
        end
        
    end
    
    if (ind2 ~= 1)
   
        % min of cycle

        Dsm = max(DS.ssnp(i), 0) -DS.sspp(i);

        KKm=DS.sspp(i)/dd(cp); % Zmax still refers to the local max of the load
        ddm(cp)=max(dd(cp)+Dsm/KKm,0); % Zmax to make sure we don't get negative ddm
        KKu=KKm-exp(Dsm/da/KKm)*(KKm-KK(cp));



        n_gap = BB(cp+1) - (BB(cp) + 1);
        KKu_new = KKu + n_gap*K_slope_mu;
        if (KKu_new < 0)
            F = ceil(((KKu - 0)/(KKu - KKu_new))*n_gap);
            ind2=1; %termination index
            NN=BB(cp) + F; %number of cycles
            ind3=1; % exceedenve of sc
        else
            KKu = KKu_new;
        end


        if((KK(cp) < 0) && (ind2 ~= 1))
            fprintf('Somehow missed failure!')
            ind2=1; %termination index
            NN=BB(cp); %number of cycles
            ind3=1; % exceedenve of sc
        end


    %     fprintf('KK = %0.2f.\n', KK(cp));
        cp = cp + 1;
        
    end
end


