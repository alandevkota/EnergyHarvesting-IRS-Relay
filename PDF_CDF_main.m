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
% L=512;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
gamma_bar=10;

%% path-loss - correlation,large scale fading with shadowing
% [Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);

%% simulation
gamma_star=[];
for pp=1:NS
% actual channels----------------------------------------------------------
    alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
    alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
    alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------
    gamma_R=gamma_bar.*sum(eta.*alpha_h.*alpha_g,2).^2;
    gamma_D=gamma_bar.*alpha_f.^2;
    gamma_star=[gamma_star  (gamma_R.*gamma_D)./(gamma_R+gamma_D+1)];
% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(NS)]);
end
%==========================================================================
%histrogram code
step=0.04;
nbins=0:step:5;
[ns_r1,x1]=hist(gamma_star,nbins);
pdf_sim=ns_r1/(NS*(step));
cdf_sim=cumsum(pdf_sim)*step;

%% analysis
Mu_X=0;
Var_X=0;
for ll=1:L
    lambda_l_sqr=eta^2.*Path_loss_S_IRS.*Path_loss_IRS_R./4;
    Mu_X=Mu_X+pi.*sqrt(lambda_l_sqr)./2;
    Var_X=Var_X+lambda_l_sqr.*(16-pi^2)./4;
end
Mu_R=sqrt(gamma_bar)*Mu_X;
Var_R=gamma_bar*Var_X;
sigma_D=gamma_bar*Path_loss_R_D;
psi=1/(qfunc(-Mu_R/sqrt(Var_R)));

pdf_ana=psi/(2*sigma_D).*exp(-x1./sigma_D) ...
    +psi./(2*sqrt(2*pi*Var_R.*x1)).*exp(-x1./sigma_D-(sqrt(x1)-Mu_R).^2./(2*Var_R)) ...
    -psi/(2*sigma_D).*erf((sqrt(x1)-Mu_R)./sqrt(Var_R)).*exp(-x1./sigma_D);

cdf_ana=1-psi.*qfunc((sqrt(x1)-Mu_R)./sqrt(Var_R)).*exp(-x1./sigma_D);

%% plotting the results
% figure
% plot(x1,pdf_sim,'--o','color',[1 0 0],'LineWidth',1); hold on;
% plot(x1,pdf_ana,'--','color',[0 0 1],'LineWidth',1); hold on;

figure(1)
plot(x1,cdf_sim,':','color',[0 0 1],'LineWidth',2);hold on;
plot(x1,cdf_ana,'--','color',[1 0 0],'LineWidth',2);
xlim([0 4])
ylim([0 1])
