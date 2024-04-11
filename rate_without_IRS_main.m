%==========================================================================

%==========================================================================
% clear all;
% clc;
disp('_________________IRS and Relay Combination________________________');
% taking input for the number of iterations
NS=input('input number of samples for estimation (1000s) : ')*1000;
%--------------------------------------------------------------------------
%==========================================================================
%% system parameters
N=1;    % No. of IRSs
L=32;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient

%% variables
gamma_bar=logspace(-1,3,20); 

%% path-loss - correlation,large scale fading with shadowing
[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D,Path_loss_S_R]=PathLoss_model(M,K,N,Nr);

%% simulation
rate_sim_exct=[];
rate_sim_exct_sr=[];
for pp=1:length(gamma_bar)
    gamma_out_count_exct=0;
    gamma_out_count_exct_sr=0;
    for ns=1:NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
        alpha_sr=abs(sqrt(Path_loss_S_R).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------
        gamma_R=gamma_bar(pp).*sum(eta.*alpha_h.*alpha_g,2).^2;
        gamma_D=gamma_bar(pp).*alpha_f.^2;
        gamma_star=(gamma_R.*gamma_D)./(gamma_R+gamma_D+1);
        gamma_out_count_exct=gamma_out_count_exct+log2(1+gamma_star);
        
        gamma_R_sr=gamma_bar(pp).*alpha_sr.^2;
        gamma_D_sr=gamma_bar(pp).*alpha_f.^2;
        gamma_star_sr=(gamma_R_sr.*gamma_D_sr)./(gamma_R_sr+gamma_D_sr+1);
        gamma_out_count_exct_sr=gamma_out_count_exct_sr+log2(1+gamma_star_sr);
    end
%     average_snr_inst_ub=gamma_out_count_ub/NS;
% rate calculation---------------------------------------------------------
    rate_sim_exct=[rate_sim_exct gamma_out_count_exct/NS];
    rate_sim_exct_sr=[rate_sim_exct_sr gamma_out_count_exct_sr/NS];
    
%% analysis
% % %     Mu_X=0;
% % %     Var_X=0;
% % %     for ll=1:L
% % %         lambda_l_sqr=eta^2.*Path_loss_S_IRS.*Path_loss_IRS_R./4;
% % %         Mu_X=Mu_X+pi.*sqrt(lambda_l_sqr)./2;
% % %         Var_X=Var_X+lambda_l_sqr.*(16-pi^2)./4;
% % %     end
% % %     Mu_R=sqrt(gamma_bar(pp))*Mu_X;
% % %     Var_R=gamma_bar(pp)*Var_X;
% % %     sigma_D=gamma_bar(pp)*Path_loss_R_D;
% % %     psi=1/(qfunc(-Mu_R/sqrt(Var_R)));
% integral calculation-----------------------------------------------------
% % %     fun_int=@(u) (sqrt(Var_R).*u+Mu_R).*qfunc(u).*exp(-(sqrt(Var_R).*u+Mu_R).^2./sigma_D);
% % %     I_value_int=integral(fun_int,-Mu_R/sqrt(Var_R),Inf);
% % %     average_snr_ana_ub_int=2*sqrt(Var_R)*psi*I_value_int;
% % %     rate_Analytical_ub_int=[rate_Analytical_ub_int log2(1+average_snr_ana_ub_int)];
% % % % analytical calculation---------------------------------------------------
% % %     average_snr_ana_ub_apx=1/2*sigma_D*psi*(1-erf(-Mu_R/(sqrt(2*Var_R))));
% % %     rate_Analytical_ub_apx=[rate_Analytical_ub_apx log2(1+average_snr_ana_ub_apx)];  
% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(gamma_bar))]);
end

%% plotting the results
gamma_bar_db=10*log10(gamma_bar);
figure(3)
plot(gamma_bar_db,rate_sim_exct,'--','color',[1 0 1],'LineWidth',2,'MarkerSize',9);hold on;
plot(gamma_bar_db,rate_sim_exct_sr,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar_db,rate_sim_ub,':s','color',[0 0 1],'LineWidth',2,'MarkerSize',10);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_int,'-','color',[1 0.5 1],'LineWidth',2);hold on;
% plot(gamma_bar_db,rate_Analytical_ub_apx,'--','color',[0 0.5 0],'LineWidth',2);