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
L=32;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
gamma_bar=10;
q_bits=1;
tau=pi/2^q_bits;

%% path-loss - correlation,large scale fading with shadowing
[Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);

%% simulation
gamma_star=[];
for pp=1:NS
% actual channels----------------------------------------------------------
    alpha_h=abs(sqrt(Path_loss_S_IRS).*(randn(M,L)+1j*randn(M,L))/sqrt(2));
    alpha_g=abs(sqrt(Path_loss_IRS_R).*(randn(Nr,L)+1j*randn(Nr,L))/sqrt(2));
    alpha_f=abs(sqrt(Path_loss_R_D).*(randn(K,Nr)+1j*randn(K,Nr))/sqrt(2));
% SNRs calculation---------------------------------------------------------

    q_error=2*tau*rand(1,L)-tau;
    gamma_R_q=gamma_bar.*abs(sum(eta.*alpha_h.*alpha_g.*exp(1j*q_error),2)).^2;
    gamma_D_q=gamma_bar.*alpha_f.^2;
    gamma_star=[gamma_star  gamma_R_q];
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
Mu_YR=0;
Var_YR=0;
Mu_YI=0;
Var_YI=0;
for ll=1:L
    Mu_YR=Mu_YR+pi*eta*sqrt(Path_loss_S_IRS/2)*sqrt(Path_loss_IRS_R/2)*sin(tau)/(2*tau);
    Var_YR=Var_YR+eta^2*(Path_loss_S_IRS/2)*(Path_loss_IRS_R/2)*(sin(2*tau)/tau+2-pi^2*sin(tau)^2/(4*tau^2));
    Var_YI=Var_YI+eta^2*(Path_loss_S_IRS/2)*(Path_loss_IRS_R/2)*(2-sin(2*tau)/tau);
end
psi_R=1/(qfunc(-Mu_YR/sqrt(Var_YR)));
psi_I=1/(qfunc(-Mu_YI/sqrt(Var_YI)));
Mu_R=sqrt(gamma_bar)*Mu_YR;
Var_R=gamma_bar*Var_YR;
Mu_I=sqrt(gamma_bar)*Mu_YI;
Var_I=gamma_bar*Var_YI;
lambda=psi_R*psi_I/(8*pi*sqrt(Var_R*Var_I));
aa=-1/2.*(1/Var_R-1/Var_I);
bb=-1/2.*(2*Mu_R/Var_R);
cc=-1/2.*(Mu_R^2/Var_R+x1./Var_I);
dd=cc-bb.^2./(4*aa);
lambda_dash=lambda.*exp(dd);
alpha=bb./(2*aa);

N_large=10;
pdx_int=[];
pdx_int_asym=[];
for pp=1:length(x1)
%     fun_1=@(t) (1./sqrt(t.*(x1(pp)-t))).*exp(aa.*(sqrt(t)-alpha).^2);
    fun_1=@(t) (1./sqrt(x1(pp)-t.^2)).*exp(aa.*(t-alpha).^2);
    qq=integral(fun_1,0,Inf);
    pdx_int=[pdx_int 2*lambda_dash(pp).*qq];
    
    temp1=0;
    for kk=0:N_large
        fun_1_asym=@(t) t.^(2*kk).*exp(aa.*(t-alpha).^2);
        II=integral(fun_1_asym,0,Inf);
        temp1=temp1+(x1(pp).^(-kk-1/2)./4^kk).*nchoosek(2*kk,kk).*II;
    end
    pdx_int_asym=[pdx_int_asym 2*lambda_dash(pp).*temp1];
end

%% plotting the results
figure
plot(x1,pdf_sim,'--o','color',[1 0 0],'LineWidth',1); hold on;
plot(x1,pdx_int,'-','color',[0 0 1],'LineWidth',1); hold on;
% plot(x1,pdx_int_asym,'--','color',[0 0 0],'LineWidth',1); hold on;

% figure
% plot(x1,cdf_sim,':','color',[0 0 1],'LineWidth',2);hold on;
% plot(x1,cdf_ana,'--','color',[1 0 0],'LineWidth',2);
% xlim([0 4])
% ylim([0 1])
