%==========================================================================

%==========================================================================
clear all;
clc;
disp('_________________IRS and Relay Combination________________________');
% taking input for the number of iterations
NS=input('input number of samples for estimation (1000s) : ')*1000;
%--------------------------------------------------------------------------
%==========================================================================
%% system parameters
N=1;    % No. of IRSs
L=64;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
q_bits=4;
tau=pi/2^q_bits;

%% variables
gamma_bar=logspace(-3,3,20); 

%% path-loss - correlation,large scale fading with shadowing
[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);

%% simulation
rate_sim_exct=[];
rate_sim_ub=[];
rate_sim_lb=[]; 
rate_Analytical_ub_int=[];
rate_Analytical_ub_apx=[];
rate_sim_exct_q=[];
rate_sim_exct_q_apx=[];
rate_Analytical_ub_apx_q=[];
vec_out_all=[];
vec_out_ana=[];
for pp=1:length(gamma_bar)
    gamma_out_count_exct=0;
    gamma_out_count_ub=0;
    gamma_out_count_exct_q=0;
    gamma_out_count_exct_q_apx=0;
    vec_out=[];
    gamma_R_q_vec=[];
    for ns=1:NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------
        gamma_R=gamma_bar(pp).*sum(eta.*alpha_h.*alpha_g,2).^2;
        gamma_D=gamma_bar(pp).*alpha_f.^2;
        gamma_star=(gamma_R.*gamma_D)./(gamma_R+gamma_D+1);
        gamma_out_count_exct=gamma_out_count_exct+log2(1+gamma_star);
        gamma_out_count_ub=gamma_out_count_ub+gamma_star;
        
        q_error=2*tau*rand(1,L)-tau;
        gamma_R_q=gamma_bar(pp).*abs(sum(eta.*alpha_h.*alpha_g.*exp(1j*q_error),2)).^2;
        gamma_D_q=gamma_bar(pp).*alpha_f.^2;
        gamma_star_q=(gamma_R_q.*gamma_D_q)./(gamma_R_q+gamma_D_q+1);
        gamma_out_count_exct_q=gamma_out_count_exct_q+log2(1+gamma_star_q);
        
        gamma_R_q_vec=[gamma_R_q_vec gamma_R_q];
        
%         gamma_star_q_apx=gamma_D_q-gamma_D_q^2/gamma_R_q;%(gamma_R_q.*gamma_D_q)./(gamma_R_q+gamma_D_q);
% % %         if gamma_R_q>gamma_D_q
% % %             gamma_star_q_apx=gamma_D_q;
% % %             vec_out=[vec_out 1];
% % %         elseif gamma_R_q<gamma_D_q
% % %             gamma_star_q_apx=gamma_R_q;
% % %             vec_out=[vec_out 2];
% % %         else
% % %             gamma_star_q_apx=(gamma_R_q+gamma_D_q)/2;
% % %             vec_out=[vec_out 3];
% % %         end
%         gamma_out_count_exct_q_apx=gamma_out_count_exct_q_apx+log2(1+gamma_star_q_apx);        

    end
%     vec_out_all=[vec_out_all; vec_out];
%     average_snr_inst_ub=gamma_out_count_ub/NS;
% rate calculation---------------------------------------------------------
    rate_sim_exct=[rate_sim_exct gamma_out_count_exct/NS];
%     rate_sim_ub=[rate_sim_ub log2(1+average_snr_inst_ub)];
    
    rate_sim_exct_q=[rate_sim_exct_q gamma_out_count_exct_q/NS];
%     rate_sim_exct_q_apx=[rate_sim_exct_q_apx gamma_out_count_exct_q_apx/NS];
    
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
% analytical calculation---------------------------------------------------
% % %     average_snr_ana_ub_apx=1/2*sigma_D*psi*(1-erf(-Mu_R/(sqrt(2*Var_R))));
% % %     rate_Analytical_ub_apx=[rate_Analytical_ub_apx log2(1+average_snr_ana_ub_apx)];  
% quantized pahse calculation----------------------------------------------
% % %     E_gamma_D_q=gamma_bar(pp)*Path_loss_R_D;
% % %     E_YR=0;
% % %     Var_YR=0;
% % %     Var_YI=0;
% % %     for ll=1:N
% % %         E_YR=E_YR+pi*eta*sqrt(Path_loss_S_IRS/2)*sqrt(Path_loss_IRS_R/2)*sin(tau)/(2*tau);
% % %         Var_YR=Var_YR+eta^2*(Path_loss_S_IRS/2)*(Path_loss_IRS_R/2)*(sin(2*tau)/tau+2-pi^2*sin(tau)^2/(4*tau^2));
% % %         Var_YI=Var_YI+eta^2*(Path_loss_S_IRS/2)*(Path_loss_IRS_R/2)*(2-sin(2*tau)/tau);
% % %     end
% % %     E_YR_2=Var_YR+E_YR^2;
% % %     E_YI_2=Var_YI;
% % %     E_gamma_R_q=gamma_bar(pp)*(E_YR_2+E_YI_2);
% % %     Var_gamma_R_q=var(gamma_R_q_vec);
% % %     E_gamma_R_q_inv=1/E_gamma_R_q+Var_gamma_R_q/E_gamma_R_q^3;  
% % %     E_gamma_D_q2=gamma_bar(pp)^2*Path_loss_R_D^2*gamma(3);
% % %     E_gamma_star_apx_q=E_gamma_D_q-E_gamma_D_q2*E_gamma_R_q_inv;
    
% % %     if E_gamma_R_q>E_gamma_D_q
% % %         E_gamma_star_apx_q=E_gamma_D_q;
% % %         vec_out_ana=[vec_out_ana; 1];
% % %     elseif E_gamma_R_q<E_gamma_D_q
% % %         E_gamma_star_apx_q=E_gamma_D_q;
% % %         vec_out_ana=[vec_out_ana; 2];
% % %     else
% % %         E_gamma_star_apx_q=(E_gamma_R_q+E_gamma_D_q)/2;
% % %         vec_out_ana=[vec_out_ana; 3];
% % %     end
% % %     E_gamma_star_apx_q1=0;
% % % 
% % %     for ni=1:NS
% % %         if vec_out(1,ni)==1
% % %             E_gamma_star_apx_q1=E_gamma_star_apx_q1+E_gamma_D_q;
% % %         elseif vec_out(1,ni)==2
% % %             E_gamma_star_apx_q1=E_gamma_star_apx_q1+E_gamma_R_q;
% % %         else
% % %             E_gamma_star_apx_q1=E_gamma_star_apx_q1+(E_gamma_R_q+E_gamma_D_q)/2;
% % %         end
% % %     end
% % %     E_gamma_star_apx_q=E_gamma_star_apx_q1/NS;
%     rate_Analytical_ub_apx_q=[rate_Analytical_ub_apx_q log2(1+E_gamma_star_apx_q)];
    
    
% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(gamma_bar))]);
end

%% plotting the results
gamma_bar_db=10*log10(gamma_bar);

figure(2)
plot(gamma_bar_db,(rate_sim_exct_q./rate_sim_exct).*100,'--','color',[1, 0, 1]	,'LineWidth',2,'MarkerSize',9);hold on;
% plot(gamma_bar_db,(rate_Analytical_ub_apx_q./rate_Analytical_ub_apx).*100,'-','color',[0, 0, 1]	,'LineWidth',2,'MarkerSize',10);hold on;
