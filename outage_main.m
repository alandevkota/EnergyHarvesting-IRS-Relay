%==========================================================================
% Date: 02/17/2020

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
% L=512;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
% gamma_th=0.1;

%% variables
% gamma_bar=logspace(-1,3,20); 

%% path-loss - correlation,large scale fading with shadowing
% [Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);

%% simulation
outageProb_sim=[];  
outageProb_ana=[];
for pp=1:length(gamma_bar)
    gamma_out_count=0;
    for ns=1:NS
% actual channels----------------------------------------------------------
        alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
        alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
        alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------
        gamma_R=gamma_bar(pp).*sum(eta.*alpha_h.*alpha_g,2).^2;
        gamma_D=gamma_bar(pp).*alpha_f.^2;
        gamma_star=(gamma_R.*gamma_D)./(gamma_R+gamma_D+1);
        gamma_out_count=gamma_out_count+(gamma_star<gamma_th);
    end
    outage_inst=gamma_out_count/NS;
    outageProb_sim=[outageProb_sim outage_inst];

%% analysis
    Mu_X=0;
    Var_X=0;
    for ll=1:L
        lambda_l_sqr=eta^2.*Path_loss_S_IRS.*Path_loss_IRS_R./4;
        Mu_X=Mu_X+pi.*sqrt(lambda_l_sqr)./2;
        Var_X=Var_X+lambda_l_sqr.*(16-pi^2)./4;
    end
    Mu_R=sqrt(gamma_bar(pp))*Mu_X;
    Var_R=gamma_bar(pp)*Var_X;
    sigma_D=gamma_bar(pp)*Path_loss_R_D;
    psi=1/(qfunc(-Mu_R/sqrt(Var_R)));
    cdf_ana=1-psi.*qfunc((sqrt(gamma_th)-Mu_R)./sqrt(Var_R)).*exp(-gamma_th./sigma_D);
    outageProb_ana=[outageProb_ana cdf_ana];

% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(gamma_bar))]);
end

%% plotting the results
gamma_bar_db=10*log10(gamma_bar);
figure(2)
semilogy(gamma_bar_db,outageProb_sim,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
semilogy(gamma_bar_db,outageProb_ana,'-','color',[0 0 1],'LineWidth',2);hold on;
