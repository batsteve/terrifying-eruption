classdef Plotter
    %PLOTTER Holder class for the various intermediate and debugging plots
    %   Figures 1 through 64 are reserved.
    
    properties(Constant)
        
        hold_all_figs = false;
        
    end
    
    methods(Static)
        function obj = Plotter()

        end
        
        function [ outcode ] = draw_spectrum_pdf(fS, S)
            Z = fS(S);
            zmin = Z(end);
            zmax = max(Z);

            figure(1);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(S, Z, 'LineWidth', 3);
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$f_\Sigma(\sigma)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([zmin, zmax]);
            title('probability density function of load', 'Interpreter', 'Latex')
            
            outcode = 1;
            
        end
        
        function [ outcode ] = draw_spectrum_over_time(fS, S, T)
            [SS, TT] = meshgrid(S, T);
            FFS = fS(SS);
            
            figure(2);
            clf;
            mesh(SS, TT, FFS);
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$f_\Sigma(\sigma; t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('joint probability density function of load across time', 'Interpreter', 'Latex')
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_phi(phi, S, T)
            [SS, TT] = meshgrid(S, T);
            PhiPhi = real(phi(SS, TT));
            
            figure(3);
            clf;
            mesh(SS, TT, PhiPhi);
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$\phi(\sigma, t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('surface plot of $\phi(\sigma, t)$', 'Interpreter', 'Latex')

            outcode = 1;
        end
        
        function [ outcode ] = draw_phi_truncated(phi, S, T, M)
            [SS, TT] = meshgrid(S, T);
            PhiPhi = real(phi(SS, TT));
            PhiPhi_trunc = max(PhiPhi, 0).*M(SS);
            
            figure(4);
            clf;
            mesh(TT', SS', PhiPhi_trunc');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            xlabel('$t$', 'Interpreter', 'Latex');
            zlabel('$\phi(\sigma, t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('truncated surface plot of $\phi(\sigma, t)$', 'Interpreter', 'Latex')

            outcode = 1;
        end
        
        function [ outcode ] = draw_coherent_envelope(e_inv_asc, e_inv_des, s_c, d_c)
            S = linspace(0, s_c, 129);
            Dasc = e_inv_asc(S);
            Ddes = e_inv_des(S);
            
            figure(5);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(Dasc, S, 'LineWidth', 3);
            plot(Ddes, S,  'LineWidth', 3);
            xlabel('$\delta$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([0, 4*d_c]);
            legend({'ascending', 'descending'}, 'Interpreter', 'Latex');
            title('coherent envelope', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_coherent_ascending_limb(e_inv_asc, S, d_c)
            Dasc = e_inv_asc(S);
            
            figure(27);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(Dasc, S, 'LineWidth', 3);
            xlabel('$\delta$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([0, d_c]);
            title('coherent envelope', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_coherent_descending_limb(e_inv_des, S, d_c)
            Ddes = e_inv_des(S);
            
            figure(28);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(Ddes, S, 'LineWidth', 3);
            xlabel('$\delta$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([d_c, 4*d_c]);
            title('coherent envelope', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_pdf_mdas_st(fMDAS, S, T)
            %Splot = [0, logspace(-3, log(max(S))/log(10), length(S)-1)];
            %Tplot = [0, logspace(0, log(max(T))/log(10), length(T)-1)];
            %[SS, TT] = meshgrid(Splot, Tplot);
            [SS, TT] = meshgrid(S, T);
%             FFM = zeros(size(SS));
%             for k = 1:size(SS, 1)
%                 for j = 1:size(SS, 2)
%                     FFM(k, j) = fMDAS(SS(k, j), TT(k, j));
%                 end
%             end
            
            FFM = fMDAS(SS, TT);
            
            figure(6);
            clf;
            mesh(SS, TT, real(FFM));
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$f_{\Sigma t}(\sigma, t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'XScale', 'log');
            %set(gca, 'YScale', 'log');
            set(gca, 'ZScale', 'log');
            zlim([1e-8, max(FFM(:))]);
            title('joint pdf of MDAS', 'Interpreter', 'Latex')

            outcode = 1;
        end
        
        function [ outcode ] = draw_spectrum_cdf(fS, S)
            FS = fS(S);
            GS = cumsum(FS);
            GS = GS./max(GS);
            
            figure(7);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(S, GS, 'LineWidth', 3)
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$F_\Sigma(\sigma)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('cumulative distribution function', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_q1_probability(Q1, P_Qplot, T)
            [PPq, TTq] = meshgrid(P_Qplot, T);
%             QQ_plot = zeros(size(PPq));
%             for k = 1:size(PPq, 1)
%                 for j = 1:size(PPq, 2)
%                     QQ_plot(k, j) = Q1(PPq(k, j), TTq(k, j));
%                 end
%             end

            QQ_plot = Q1(PPq, TTq);
            
            figure(8);
            clf;
            mesh(PPq, TTq, QQ_plot);
            xlabel('$\phi_0$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$Q_1 = P(\phi(\sigma, t) < \phi_0)$', 'Interpreter', 'Latex');
            %set(gca, 'XScale', 'log');
            %set(gca, 'YScale', 'log');
            set(gca, 'FontSize', 14);
            title('$Q_1$ -- probability that $\phi_0$ is not exceeded', 'Interpreter', 'Latex')
            
            figure(37);
            clf;
            mesh(PPq, TTq, 1-QQ_plot);
            xlabel('$\phi_0$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$1 - Q_1 = P(\phi(\sigma, t) > \phi_0)$', 'Interpreter', 'Latex');
            %set(gca, 'XScale', 'log');
            %set(gca, 'YScale', 'log');
            set(gca, 'ZScale', 'log');
            set(gca, 'FontSize', 14);
            title('$1 - Q_1$ -- probability that $\phi_0$ is exceeded', 'Interpreter', 'Latex')

            outcode = 1;
        end
        
        function [ outcode ] = draw_loqq1(logq, phi_max)
            P = logspace(0, log(phi_max)/log(10), 128);
            LL = zeros(length(P), 1);
            for k = 1:length(P)
                LL(k) = logq(P(k));
            end
            
            figure(32);
            clf;
            if ~Plotter.hold_all_figs, clf; end
            plot(P, LL, 'linewidth', 3);
            ylabel('$\int \log Q_1(\phi) dt$', 'Interpreter', 'Latex');
            xlabel('$\phi$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'XScale', 'log');
            title('logq vs $\phi$', 'Interpreter', 'Latex');
            
            figure(33);
            clf;
            if ~Plotter.hold_all_figs, clf; end
            plot(P, exp(LL), 'linewidth', 3);
            ylabel('$\exp(\int \log Q_1(\phi) dt)$', 'Interpreter', 'Latex');
            xlabel('$\phi$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'XScale', 'log');
            title('explogq vs $\phi$', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_kappa_family(kappa, kappa_inv, ...
                kappa_prime, kappa_prime_inv, S, K)
            
            K1 = kappa(S);
            Z2 = zeros(size(K));
            Z3 = kappa_prime(S);
            for k = 1:length(K)
                Z2(k) = kappa_inv(K(k));
            end
            Z4 = abs(kappa_prime_inv(K));
            
            figure(38);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            %plot(K1, S, 'LineWidth', 3);
            %plot(K, Z2, 'LineWidth', 3, 'LineStyle', ':'); %  'Marker', '*'
            plot(S, K1, 'LineWidth', 3);
            plot(Z2, K, 'LineWidth', 3, 'LineStyle', ':');
            plot(S, Z3, 'LineWidth', 3, 'LineStyle', '--')
            plot(Z4, K, 'LineWidth', 3, 'LineStyle', '-.')
            ylabel('$K_a^+$', 'Interpreter', 'Latex');
            xlabel('$\sigma_a$', 'Interpreter', 'Latex');
            K_max = max(abs(K));
            ylim([-K_max, K_max]);
            xlim([0, max(S)]);
            set(gca, 'FontSize', 14);
            legend({'forward', 'inverse', 'derivative', 'div-inv'}, 'Interpreter', 'Latex');
            title('kappa family', 'Interpreter', 'Latex');

            outcode = 1;
        end


function [ outcode ] = draw_eta_family(eta, eta_inv, ...
                eta_prime, eta_prime_inv, S)
            
            E = linspace(0, 1, 129);

            K1 = eta(S);
            Z2 = zeros(size(E));
            Z3 = eta_prime(S);
            for k = 1:length(E)
                Z2(k) = eta_inv(E(k));
            end
            Z4 = abs(eta_prime_inv(E));
            
            figure(51);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            %plot(K1, S, 'LineWidth', 3);
            %plot(K, Z2, 'LineWidth', 3, 'LineStyle', ':'); %  'Marker', '*'
            plot(S, K1, 'LineWidth', 3);
            plot(Z2, E, 'LineWidth', 3, 'LineStyle', ':');
            plot(S, Z3, 'LineWidth', 3, 'LineStyle', '--')
            plot(Z4, E, 'LineWidth', 3, 'LineStyle', '-.')
            ylabel('$K_a^+$', 'Interpreter', 'Latex');
            xlabel('$\sigma_a$', 'Interpreter', 'Latex');
            E_max = max(abs(E));
            ylim([-E_max, E_max]);
            xlim([0, max(S)]);
            set(gca, 'FontSize', 14);
            legend({'forward', 'inverse', 'derivative', 'div-inv'}, 'Interpreter', 'Latex');
            title('eta family', 'Interpreter', 'Latex');

            outcode = 1;
        end

        
        function [ outcode ] = draw_kappa(kappa,S)
            K1 = kappa(S);
            
            figure(9);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(S, K1, 'LineWidth', 3);
            ylabel('$K_a^+$', 'Interpreter', 'Latex');
            xlabel('$\sigma_a$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('kappa function', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_pdf_mdas_kt(fKt, K, T)
            [KKmdas, TTmdas] = meshgrid(K, T);
%             FFKt = zeros(size(KKmdas));
%             for k = 1:size(KKmdas, 1)
%                 for j = 1:size(KKmdas, 2)
%                     FFKt(k, j) = fKt(KKmdas(k, j), TTmdas(k, j));
%                 end
%             end
            
            FFKt = fKt(KKmdas, TTmdas);
            
            figure(10);
            clf;
            mesh(KKmdas, TTmdas, FFKt);
            xlabel('$K^+_t$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$f_{K t}(k, t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('joint pdf of MDAS', 'Interpreter', 'Latex')
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_psi_constant_K(psi, S, T, K_min, K_max, M)
            [SS, TT] = meshgrid(S, T);
            Kpsi = linspace(K_min, K_max, 4);
            
            figure(11);
            clf;
            figure(12);
            clf;
            for k = 1:4
                subplot(2, 2, k);
                ZZ = zeros(size(SS));
                for j = 1:size(SS, 1)
                    for l = 1:size(SS, 2)
                        ZZ(j, l) = psi(SS(j, l), TT(j, l), Kpsi(k));
                    end
                end
                figure(11);
                subplot(2, 2, k);
                mesh(SS, TT, ZZ);
                xlabel('$\sigma$', 'Interpreter', 'Latex');
                ylabel('$\Delta t$', 'Interpreter', 'Latex');
                zlabel('$\psi$', 'Interpreter', 'Latex');
                set(gca, 'FontSize', 14);
                title(sprintf('K = %0.2f', Kpsi(k)), 'Interpreter', 'Latex');

                figure(12);
                subplot(2, 2, k);
                mesh(SS, TT, ZZ.*(ZZ > 0).*M(SS));
                xlabel('$\sigma$', 'Interpreter', 'Latex');
                ylabel('$\Delta t$', 'Interpreter', 'Latex');
                zlabel('$\psi$', 'Interpreter', 'Latex');
                set(gca, 'FontSize', 14);
                title(sprintf('K = %0.2f', Kpsi(k)), 'Interpreter', 'Latex');
            end

            outcode = 1;
        end
        
        function [ outcode ] = draw_positive_anticipation(nDS, S, K, M)
            [SS_k, KK_k] = meshgrid(S, K);
            ZZ = zeros(size(SS_k));
            for k = 1:size(SS_k, 1)
                for j = 1:size(SS_k, 2)
                    ZZ(k, j) = nDS(SS_k(k, j), KK_k(k, j));
                end
            end
            
            figure(13);
            clf;
            mesh(SS_k, KK_k, ZZ.*(ZZ>0).*M(SS_k));
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$K^+_+$', 'Interpreter', 'Latex');
            zlabel('$\Delta t$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('Locus of positive anticipation', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_q2_probability(Q2, T, K)
            [TT_q2, KK_q2] = meshgrid(T, K);
            Q2_plot = zeros(size(KK_q2));
            for k = 1:size(KK_q2, 1)
                for j = 1:size(KK_q2, 2)
                    Q2_plot(k, j) = Q2(KK_q2(k, j), TT_q2(k, j));
                end
            end
            
            figure(14);
            clf;
            mesh(KK_q2, TT_q2, Q2_plot);
            xlabel('$K^+_a$', 'Interpreter', 'Latex');
            ylabel('$\Delta t$', 'Interpreter', 'Latex');
            zlabel('$Q_2 = P(\psi(K_a^+, \Delta t) < 0)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('$Q_2$ -- probability that $\psi$ is not positive', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_gap_pdf(fDS, T ,K)
            [TT_q2, KK_q2] = meshgrid(T, K);
            FDS = zeros(size(KK_q2));
            for k = 1:size(KK_q2, 1)
                for j = 1:size(KK_q2, 2)
                    FDS(k, j) = fDS(TT_q2(k, j), KK_q2(k, j));
                end
            end
            
            figure(15);
            clf;
            mesh(KK_q2, TT_q2, FDS);
            xlabel('$K^+_a$', 'Interpreter', 'Latex');
            ylabel('$\Delta t$', 'Interpreter', 'Latex');
            zlabel('$f_{\Delta t}(\Delta t; K^+_a)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('pdf for gap size before DS', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        % figure(16);
        % clf;
        % Fqui = zeros(size(T));
        % for k = 1:length(T)
        %     Fqui(k) = f_qui(T(k));
        % end
        % plot(T, FFail, 'LineWidth', 3);
        % xlabel('$t$', 'Interpreter', 'Latex');
        % ylabel('$f_T^{\mbox{QUI}}(t)$', 'Interpreter', 'Latex');
        % set(gca, 'FontSize', 14);
        % title('pdf for quiescent mode failure time', 'Interpreter', 'Latex');
        
        function [ outcode ] = draw_quiescent_failure_pdf(f_qui, T)
            P = f_qui(T);
            
            figure(17);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, P, 'LineWidth', 3);
            xlabel('$t$', 'Interpreter', 'Latex');
            ylabel('$f_T^{\mbox{QUI}}(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([min(P), max(P)]);
            title('pdf for quiescent mode failure time -- conv method', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_extreme_failure_pdf(f_ex, T)
            Zex = f_ex(T);
            
            figure(18);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, Zex, 'LineWidth', 3)
            xlabel('$t$', 'Interpreter', 'Latex');
            ylabel('$f_T^{\mbox{EX}}(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            title('pdf for extreme mode failure time', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_total_failure_pdf(f_full, T)
            P = f_full(T);
            
            figure(19);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, P, 'LineWidth', 3);
            xlabel('$t$', 'Interpreter', 'Latex');
            ylabel('$f_T(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            %ylim([min(F_ex), max(F_full)]);
            ylim([1e-8, 1e-2]);
            title('total pdf for failure time', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_failure_pdf_components(F_ex, G_ex, F_qui, G_qui, T)
            figure(20);
            clf;
            hold on
            plot(T, F_ex, 'LineWidth', 3);
            plot(T, G_ex, 'LineWidth', 3);
            plot(T, F_qui, 'LineWidth', 3);
            plot(T, G_qui, 'LineWidth', 3);
            xlabel('$t$', 'Interpreter', 'Latex');
            ylabel('$f_T(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            legend({'extreme pdf', 'extreme cdf', 'quiescent pdf', 'quiesecent cdf'})
            title('total failure time pdf components', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_mdas_marginal_pdf(fMDAS, S, T, ds)
            [SS, TT] = meshgrid(S, T);
            F_MDAS_TS = zeros(size(SS));
            for k_s = 1:size(SS, 1)
                for k_t = 1:size(SS, 2)
                    F_MDAS_TS(k_s, k_t) = fMDAS(SS(k_s, k_t), TT(k_s, k_t));
                end
            end
            F_MDAS_T = sum(F_MDAS_TS, 2)*ds;
            
            figure(21);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, F_MDAS_T, 'LineWidth', 3);
            xlabel('$t$', 'Interpreter', 'Latex');
            ylabel('$f_T^{MDAS}(t)$', 'Interpreter', 'Latex');
            ylim([max(1e-10, min(F_MDAS_T)), min(max(F_MDAS_T), 1)])
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            title('pdf for time of MDAS', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_mdas_damage_pdf_old(fMDAS, phi_inv, phiGrid, T, dt)
            %
            % Missing a transfer Jacobian!
            %
            [PP, TT] = meshgrid(phiGrid, T);
            FPt = zeros(size(PP));
            for k_p = 1:size(PP, 1)
                for k_t = 1:size(PP, 2)
                    s = phi_inv(PP(k_p, k_t), TT(k_p, k_t));
                    dp = 0.05;
                    J = (phi_inv(PP(k_p, k_t) + dp, TT(k_p, k_t)) - ...
                         phi_inv(PP(k_p, k_t) - dp, TT(k_p, k_t)))/(2*dp);
                    FPt(k_p, k_t) = fMDAS(s, TT(k_p, k_t))*abs(J);
                end
            end
            FPt = max(FPt, 0);
            FP = sum(FPt, 2)*dt;
            FP = FP/sum(FP);
            
            figure(26);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(phiGrid, FP, 'LineWidth', 3);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$f_{\Phi}(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            title('pdf for damage of MDAS--Missing Jacobian!', 'Interpreter', 'Latex');

            outcode = 1;
        end

        function [ outcode ] = draw_mdas_damage_pdf(f_phi, phi_max)
            P = linspace(0, phi_max, 256);
            Z = f_phi(P);
            
            figure(26);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(P, Z, 'LineWidth', 3);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$f_{\Phi}(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-8, 1e0]);
            title('pdf for damage of MDAS', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_delta_marginal_pdf(fKt, fDS, T, K, dk, dt)
            [TT_q2, KK_q2] = meshgrid(T, K);
            F_KT = zeros(size(TT_q2));
            for k_t = 1:size(TT_q2, 1)
                for k_k = 1:size(TT_q2, 2)
                    F_KT(k_t, k_k) = fKt(KK_q2(k_t, k_k), TT_q2(k_t, k_k));
                end
            end
            F_K = sum(F_KT, 2)*dt;
            F_DS_TK = zeros(size(TT_q2));
            for k_k = 1:size(TT_q2, 1)
                for k_t = 1:size(TT_q2, 2)
                    F_DS_TK(k_k, k_t) = fDS(TT_q2(k_k, k_t), KK_q2(k_k, k_t)).*F_K(k_k);
                end
            end
            F_DS_T = sum(F_DS_TK, 1)*dk;
            
            figure(22);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, F_DS_T, 'LineWidth', 3);
            xlabel('$dt$', 'Interpreter', 'Latex');
            ylabel('$f_T^{DS}(t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-15, 1]);
            title('pdf for time delta until DS', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_shifted_Q2(Q2, T_max, delta_K)
            n_plot = 256;
            T_squig = logspace(-1, log(T_max/2)/log(10), n_plot/2);
            T_squig = [-fliplr(T_squig), 0, T_squig];
            
            k_0 = 1.5;
            ZZ = Q2(k_0, T_squig + k_0/delta_K);
            
            figure(23);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T_squig, ZZ, 'LineWidth', 3);
            xlabel('$\tilde{T}$', 'Interpreter', 'Latex');
            ylabel('$Q_2$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            title('shifted representation of $Q_2$', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_shifted_gap_pdf(fDS, T_max, delta_K)
            n_plot = 256;
            T_squig = logspace(-1, log(T_max/2)/log(10), n_plot/2);
            T_squig = [-fliplr(T_squig), 0, T_squig];
            
            k_0 = 1.5;
            ZZ = fDS(T_squig + k_0/delta_K, k_0);
            
            figure(24);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T_squig, ZZ, 'LineWidth', 3);
            xlabel('$\tilde{T}$', 'Interpreter', 'Latex');
            ylabel('$f_{DS}$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            title('shifted representation of $f_{DS}$', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_anticipation_st_curve(psi, T, k_eps, s_c)
            
            k_0 = k_eps;
            Z = zeros(length(T), 1);
            for k = 1:length(T)
                Z(k) = fminbnd(@(s) psi(s, T(k), k_0).^2, 0, s_c);
            end
            
            figure(25);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, Z, 'LineWidth', 3)
            xlabel('$T$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('Locus of positive anticipation', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        
        function [ outcode ] = draw_phi_inverse(phi_inv, Phi_grid, T)
            [PP, TT] = meshgrid(Phi_grid, T);
            ZZ = phi_inv(PP, TT);
            
            figure(29);
            clf;
            mesh(PP, TT, ZZ);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$s = \phi^{-1}(\phi_0, t)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('inverse phi mapping', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
       
        
        function [ outcode ] = draw_mdas_exp_integral(Q1, fS_qui, phi, phi_grid, T )
            
            %fMDAS_exact = @(s, t) fS_qui(s) .* ...
            %    exp( integral( @(tau) log(Q1(phi(s, t), tau)'), 0, T_max));
            
            logq = @(P, tau) log(Q1(P, tau)');
            [PP, TT] = meshgrid(phi_grid, T);
            ZZ = zeros(size(PP));
            for k = 1:size(PP, 1)
                for j = 1:size(PP, 2)
                    ZZ(k, j) = logq(PP(k, j), TT(k, j));
                end
            end
            
            figure(31);
            clf;
            mesh(PP, TT, ZZ);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$\log Q$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('log Q', 'Interpreter', 'Latex');
            
            Z = zeros(length(phi_grid), 1);
            for k = 1:length(phi_grid)
                Z(k) = integral(@(t) logq(phi_grid(k), t), 1, max(T));
            end
            
            figure(32);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(phi_grid, Z, 'LineWidth', 3);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$\log Q$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            title('summed log Q', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_phi_jacobian(J, Phi_grid, T)
            
            [PP, TT] = meshgrid(Phi_grid, T);
            JJ = zeros(size(PP));
            for k = 1:size(PP, 1)
                for j = 1:size(PP, 2)
                    JJ(k, j) = J(PP(k, j), TT(k, j));
                end
            end
            
            figure(34);
            clf;
            mesh(PP, TT, JJ);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$J$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('$Kt--\phi$ Jacobian', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_phi_prime(phi_prime, phi_inv, Phi_grid, S, T, s_c)
            [SS, TTs] = meshgrid(S, T);
            ZZ1 = zeros(size(SS));
            for k = 1:size(SS, 1)
                for j = 1:size(SS, 2)
                    ZZ1(k, j) = phi_prime(SS(k, j), TTs(k, j));
                end
            end
            
            [PP, TTp] = meshgrid(Phi_grid, T);
            ZZ2 = zeros(size(PP));
            for k = 1:size(PP, 1)
                for j = 1:size(PP, 2)
                    ZZ2(k, j) = phi_prime(phi_inv(PP(k, j),  TTp(k, j)), TTp(k, j));
                end
            end
            
            ds = S(2)  - S(1);
            
            figure(35);
            clf;
            mesh(SS, TTs, ZZ1);
            xlabel('$\sigma$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$\phi''$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([2*ds, s_c  - 2*ds]);
            zlim([0, 2e5])
            title('$\phi''(\sigma, t)$', 'Interpreter', 'Latex');
            
            figure(36);
            clf;
            mesh(PP, TTp, ZZ2);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$t$', 'Interpreter', 'Latex');
            zlabel('$\phi''J$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('$\phi''(\phi, t)$', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_q1_theta(Q1_theta, th_max)
            TT = logspace(0, log(th_max)/log(10), 128);
            QQ = Q1_theta(TT);
            
            figure(39);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(TT, QQ, 'LineWidth', 3);
            xlabel('$\theta$', 'Interpreter', 'Latex');
            ylabel('$Q_1$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            set(gca, 'XScale', 'log');
            title('$Q_1(\theta)$ ', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_inverted_envelope(e_asc, e_des)
            S = linspace(0, 1, 128);
            Ea = e_asc(S);
            Ed = e_des(S);
            
            figure(40);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(Ea, S, 'LineWidth', 3);
            plot(Ed, S, 'LineWidth', 3);
            xlabel('$\delta$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            %set(gca, 'XScale', 'log');
            title('inverted coherent envelope', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_envelope_inversion_comparison(e_fit_asc, e_fit_des, e_inv_asc,e_inv_des)
            S = linspace(0, 1, 128);
            Efa = e_fit_asc(S);
            Efd = e_fit_des(S);
            
            Eta = e_inv_asc(S);
            Etd = e_inv_des(S);
            
            figure(41);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(Efa, S, 'LineWidth', 3, 'Color', 'blue', 'LineStyle', '-');
            plot(Efd, S, 'LineWidth', 3, 'Color', 'blue', 'LineStyle', '-');
            plot(Eta, S, 'LineWidth', 3, 'Color', 'red', 'LineStyle', ':');
            plot(Etd, S, 'LineWidth', 3, 'Color', 'red', 'LineStyle', ':');
            xlabel('$\delta$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'YScale', 'log');
            %set(gca, 'XScale', 'log');
            legend({'fit ascending', 'fit descending', 'true ascending', 'true descending'})
            title('coherent envelope (recovered vs true)', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = debug_draw_psi_obs(S_parts, T_parts, psi_parts)
            M1 = (psi_parts > 0);
            M2 = (psi_parts <= 0);
            
            figure(42);
            clf;
            hold on
            scatter(T_parts(M1), S_parts(M1), 'o');
            scatter(T_parts(M2), S_parts(M2), '+');
            xlabel('$T$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            
            title('coherent envelope (recovered vs true)', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = debug_draw_psi_obs_with_boundary(S_parts, T_parts, psi_parts, f)
            M1 = (psi_parts > 0);
            M2 = (psi_parts <= 0);
            
            S_max = max(1, max(S_parts));
            S = linspace(0, S_max, 128);
            T = f(S);
            
            figure(43);
            clf;
            hold on
            scatter(T_parts(M1), S_parts(M1), 'o');
            scatter(T_parts(M2), S_parts(M2), '+');
            plot(T, S, 'LineWidth', 3, 'LineStyle', ':');
            xlabel('$T$', 'Interpreter', 'Latex');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([0, max(T_parts)]);
            
            title('coherent envelope (recovered vs true)', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        function [ outcode ] = draw_nu_experiment(nu, s_c)
            SS = linspace(0, s_c, 129);
            SS = SS(2:end);
            NN = nu(SS);
            
            figure(44);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            %plot(SS, NN, 'LineWidth', 3);
            %xlabel('$\sigma$', 'Interpreter', 'Latex');
            %ylabel('$N$', 'Interpreter', 'Latex');
            plot(NN, SS, 'LineWidth', 3);
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            xlabel('$N$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            title('number of harmonic cycles until failure', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_reduced_gap_pdf(fDS, T_max)
            
            TT = linspace(-T_max, T_max, 257);
            
            figure(45);
            clf;
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(TT, fDS(TT, 0), 'LineWidth', 3);
            set(gca, 'YScale', 'log');
            set(gca, 'FontSize', 14);
            xlabel('$\tilde{T}$', 'Interpreter', 'Latex');
            ylabel('$f_{\tilde{T}}(d\tilde{t})$', 'Interpreter', 'Latex');
            title('reduced pdf of DS gap, $K_a = 0$', 'Interpreter', 'Latex');
            
            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_kappa_pdf(fKappa, FKappa, K)
            ffK = fKappa(K);
            FFK = FKappa(K);
            
            figure(45);
            clf
            hold on
            plot(K, ffK, 'LineWidth', 3);
            plot(K, FFK, 'LineWidth', 3);
            xlabel('$K(\sigma)$', 'Interpreter', 'Latex');
            ylabel('$f_K(k)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([min(K), max(K)]);
            set(gca, 'YScale', 'log');
            legend({'pdf', 'cdf'}, 'Interpreter', 'Latex');
            title('pdf of $\kappa(\sigma)$', 'Interpreter', 'Latex');

            outcode = 1;
        end



        function [ outcode ] = draw_eta_pdf(fEta, FEta)
            E = linspace(0, 1, 129);
            ffE = fEta(E);
            FFE = FEta(E);
            
            figure(52);
            clf
            subplot(2, 1, 1);
            hold on
            plot(E, ffE, 'LineWidth', 3);
            xlabel('$H(\sigma)$', 'Interpreter', 'Latex');
            ylabel('$f_H(h)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([min(E), max(E)]);
            set(gca, 'YScale', 'log');
            title('pdf of $\eta(\sigma)$', 'Interpreter', 'Latex');
            
            %figure(61);
            %clf
            subplot(2, 1, 2);
            hold on
            plot(E, FFE, 'LineWidth', 3);
            plot(E, 1-FFE, 'LineWidth', 3);
            xlabel('$H(\sigma)$', 'Interpreter', 'Latex');
            ylabel('$f_H(h)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([min(E), max(E)]);
            set(gca, 'YScale', 'log');
            legend({'cdf', '$1-$cdf'}, 'Interpreter', 'Latex');
            title('cdf of $\eta(\sigma)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_fphin_conditional(fphin, T)
            [TT1, TT2] = meshgrid(T, T);
            FPN = fphin(TT1, TT2);
            
            figure(47);
            clf;
            mesh(TT1, TT2, FPN);
            xlabel('$n$', 'Interpreter', 'Latex');
            ylabel('$\phi$', 'Interpreter', 'Latex');
            zlabel('$f_{\phi_n}(\phi_n; n)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            title('pdf of $\phi_n$ for fixed $n$', 'Interpreter', 'Latex')
            
            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_q_prod(q_prod, phi_max)
            K = linspace(0, phi_max, 129);
            Q = q_prod(K);
            
            figure(48);
            if ~Plotter.hold_all_figs, clf; end
            subplot(2, 1, 1)
            hold on
            plot(K, Q, 'LineWidth', 3);
            xlabel('$y$', 'Interpreter', 'Latex');
            ylabel('$F_K(k)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([min(K), max(K)]);
            set(gca, 'YScale', 'log');
            title('$Q_\kappa(y)$', 'Interpreter', 'Latex');
            
            subplot(2, 1, 2)
            hold on
            plot(K, 1 - Q, 'LineWidth', 3);
            xlabel('$y$', 'Interpreter', 'Latex');
            ylabel('$1 - F_K(k)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([min(K), max(K)]);
            set(gca, 'YScale', 'log');
            %title('$Q_\kappa(y)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_W(W, K_max)
            K = linspace(0, K_max, 129);
            WW = W(K);
            
            figure(49);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, WW, 'LineWidth', 3);
            xlabel('$K$', 'Interpreter', 'Latex');
            ylabel('$W(k)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %xlim([min(K_fix), max(K_fix)]);
            set(gca, 'YScale', 'log');
            title('$W_\kappa(k)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_fphi_comparison(fphi_1, fphi_2, T)
            P1 = fphi_1(T);
            P2 = fphi_2(T);
            
            figure(50);
            clf
            hold on
            plot(T, P1, 'LineWidth', 3);
            plot(T, P2, 'LineWidth', 3);
            xlabel('$\phi$', 'Interpreter', 'Latex');
            ylabel('$f_{\phi}(n)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %xlim([min(K_fix), max(K_fix)]);
            legend({'from phi marginalization', 'from Q-W integral'})
            set(gca, 'YScale', 'log');
            ylim([1e-15, 1e-1])
            title('Comparison of $f_{\phi(n)}$ computations', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_Veta(Veta, K_max)
            K = linspace(K_max - 1, K_max, 129);
            [K1, K2] = meshgrid(K, K);
            VV = zeros(size(K1));
            for k = 1:length(K)
                for j = 1:length(K)
                    VV(k, j) = Veta(K1(k, j), K2(k, j));
                end
            end
            
            figure(60);
            clf;
            mesh(K1, K2, real(VV));
            xlabel('$x_l$', 'Interpreter', 'Latex');
            ylabel('$x_u$', 'Interpreter', 'Latex');
            zlabel('$V(x_l, x_u)$', 'Interpreter', 'Latex');
            %zlim([-1, 1]);
            set(gca, 'FontSize', 14);
            title('helper integral $V(x_l, x_u)$', 'Interpreter', 'Latex')
            
            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_Vexp(Vexp, K_max)
            K = linspace(K_max - 1, K_max, 129);
            [K1, K2] = meshgrid(K, K);
            VV = zeros(size(K1));
            for k = 1:length(K)
                for j = 1:length(K)
                    VV(k, j) = Vexp(K1(k, j), K2(k, j));
                end
            end
            
            figure(64);
            clf;
            mesh(K1, K2, real(VV));
            xlabel('$x_l$', 'Interpreter', 'Latex');
            ylabel('$x_u$', 'Interpreter', 'Latex');
            zlabel('$exp(V(x_l, x_u))$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            %set(gca, 'ZScale', 'log');
            title('helper integral $\exp(V(x_l, x_u))$', 'Interpreter', 'Latex')
            
            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_Szeta(Szeta, K_max)
            K = linspace(0, K_max, 129);
            S = Szeta(K);
            
            figure(54);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, S, 'LineWidth', 3);
            xlabel('$\zeta$', 'Interpreter', 'Latex');
            ylabel('$S_{\kappa}(\zeta)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-10, 1e0])
            title('helper function $S_{\kappa}(\zeta)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_q_prod_prime(q_prod_prime, phi_max)
            K = linspace(0, phi_max, 129);
            Qp = q_prod_prime(K);
            
            figure(55);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, Qp, 'LineWidth', 3);
            xlabel('$\zeta$', 'Interpreter', 'Latex');
            ylabel('$Q''_{\kappa}(\zeta)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-4, 1e2]);
            title('pdf of $Q''_{\kappa}(\zeta)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_hvqs_quiescent_failure(f_n, T)
            P = f_n(T);
            
            figure(56);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, P, 'LineWidth', 3);
            xlabel('$n$', 'Interpreter', 'Latex');
            ylabel('$f_N(n)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-8, 1e-2]);
            title('pdf of quiescent failure time (hvqs method)', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        function [ outcode ] = draw_envelope_kappa_eta_transformation(kappa, eta, s_max)
            S = linspace(0, s_max, 129);
            K = kappa(S);
            H = eta(S);
            
            figure(57);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, S, 'LineWidth', 3);
            plot(H, S, 'LineWidth', 3);
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            xlabel('$\mathcal{N}_{\kappa,\eta}[ \frac{e(\sigma)}{\sigma} ]$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            legend({'$\kappa$', '$\eta$'}, 'Interpreter', 'Latex')
            title('$(\kappa, \eta)$ transformation of coherent envelope', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        function [ outcode ] = draw_envelope_kappa_eta_comparison(e_inv_asc,e_inv_des, kappa, eta, s_max)
            S = linspace(0, s_max, 129);
            K = kappa(S);
            H = eta(S);
            
            A = e_inv_asc(S);
            D = e_inv_des(S);
            
            figure(58);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            [p1] = plot(A, S, 'LineWidth', 3, 'Color', 'Blue');
            plot(D, S, 'LineWidth', 3, 'Color', 'Blue');
            [p3] = plot(K, S, 'LineWidth', 3, 'Color', 'Red');
            plot(H, S, 'LineWidth', 3, 'Color', 'Red');
            ylabel('$\sigma$', 'Interpreter', 'Latex');
            xlabel('$\delta$ OR $K$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            xlim([0, 6])
            legend([p1, p3], {'envelope', '$(\kappa, \eta)$'}, 'Interpreter', 'Latex')
            title('comparison of coherent envelope and $(\kappa, \eta)$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_Veta_integrand(integrand, K_max)
            K = linspace(K_max - 1, K_max, 257);
            P = integrand(K);
            
            figure(59);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, P, 'LineWidth', 3);
            xlabel('$s$', 'Interpreter', 'Latex');
            ylabel('$\log(1 - F_{\eta}(K_0^+ - s))$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([-1, -1e-10]);
            title('integrand for $V_{\eta}$ helper function', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_q_integrand(Q_integrand, K_max)
            K = linspace(0, K_max, 129);
            P = abs(Q_integrand(K));
            
            figure(62);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(K, P, 'LineWidth', 3);
            xlabel('$s$', 'Interpreter', 'Latex');
            ylabel('$\log(1 - F_{\kappa}(K_0^+ - s))$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            title('integrand for $Q_{\kappa}$ helper function', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        
        
        function [ outcode ] = draw_fna(fna, T)
            P = fna(T);
            
            figure(63);
            if ~Plotter.hold_all_figs, clf; end
            hold on
            plot(T, P, 'LineWidth', 3);
            xlabel('$a$', 'Interpreter', 'Latex');
            ylabel('$f_{n_a}(a)$', 'Interpreter', 'Latex');
            set(gca, 'FontSize', 14);
            set(gca, 'YScale', 'log');
            ylim([1e-15, 1]);
            title('distribution of last ascending upcrossing $n_a$', 'Interpreter', 'Latex');

            outcode = 1;
        end
        
        % next figure(65)
        
        function [ outcode ] = draw()
            outcode = 1;
        end
        
    end
end

