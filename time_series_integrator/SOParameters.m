classdef SOParameters < handle
    %SOPARAMETERS Serebrinsky-Ortiz model parameters
    
    properties
        
        dc;
        sc;
        da;
        
        NNc;
        QQe;
        n_SN;
        
        input_signal_length;

        K_slope_mu;
        extreme_height_threshold;
    end
    
    methods
        function param = SOParameters()
            %SOPARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            

            param.dc=1; param.sc=250; param.da=300;
            
            param.n_SN = 256;
            
            %[param.NNc,param.QQe] = run_counts(param); % note that param.QQe refers to amplitude - not the range as it is the case in the raincount alg
            
            param.input_signal_length = 8e4;
            
            param.K_slope_mu = -9.7e-3;
            param.extreme_height_threshold = 9.6;
        end
        
        
        function [NNc,QQe] = run_counts(param, e_forward)

            dc=param.dc; sc=param.sc; da=param.da;
            options = optimset('Display','none');

            %QQ=[0.25:0.25:param.sc];
            QQ = linspace(0, param.sc, param.n_SN + 1);
            QQ = QQ(2:end);
            
            ddm = zeros(length(QQ), 1);
            ddp = zeros(length(QQ), 1);
            NNc = zeros(length(QQ), 1);
            
            for i=1:length(QQ)
                ddm(i)=lsqnonlin(@(d00)e_forward(d00)-QQ(i),exp(1)*sc/dc,0,dc,options);
                ddp(i)=lsqnonlin(@(d00)e_forward(d00)-QQ(i),exp(1)*sc/dc,dc,10*dc,options);

                NNc(i)=da*(1/ddm(i)-1/ddp(i)); 
            end
            QQe=QQ;
        end
        
        
    end
    
end

