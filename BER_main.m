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
% L=64;   % No. of elements
M=1;    % No. of BS antennas
K=1;    % No. of users 
Nr=1;   % No. of relays
sigma_N=1;  % noise power
eta=1;  % reflection coefficient
omega=1;
vartheta=2;

%% variables
gamma_bar=logspace(-3,3,20); 

%% path-loss - correlation,large scale fading with shadowing
% [Path_loss_S_IRS,Path_loss_IRS_R,Path_loss_R_D]=PathLoss_model(M,K,N,Nr);


%% simulation
ber_sim=[];
ber_ana=[];
I_3_int=[];
I_3_ana=[];
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
        gamma_out_count=gamma_out_count+omega*qfunc(sqrt(vartheta.*gamma_star));
    end
    ber_inst=gamma_out_count/NS;
    ber_sim=[ber_sim ber_inst];

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
    A=omega*psi*sqrt(vartheta)/(2*sqrt(2*pi));
    aa=vartheta/2+1/sigma_D;
    B=2*sqrt(Var_R);
    alpha=sqrt(aa*Var_R);
    beta=sqrt(aa)*Mu_R;
    l_min=-Mu_R/sqrt(Var_R);
    C=B/(2*sqrt(2)*alpha);
    l_min_dash=alpha*l_min+beta;
    D=exp(-beta^2/(2*alpha^2))/(alpha*sqrt(pi));
    rho=1+1/(2*alpha^2);
    gg=beta/alpha^2;
    E=D*exp(gg^2/(4*rho));
    qq=gg/(2*rho);
    b_min=l_min_dash-qq;
    
    I_1=-B*sqrt(pi)/(2*alpha)*qfunc(l_min)*erf(alpha*l_min+beta);
    I_dash1=sqrt(pi/2)*(1-erf(l_min/sqrt(2)));
    
%     fun_I_dash=@(tt) exp(-tt.^2./2)-exp(-tt.^2./2).*exp(-(alpha.*tt+beta).^2)./(sqrt(pi).*(alpha.*tt+beta));%erf(alpha.*tt+beta); 
%     I_dash=integral(fun_I_dash,l_min,inf);

%     fun_I_dash=@(tt) E.*(1./(tt+qq)).*exp(-rho.*tt.^2);
%     I_dash2=integral(fun_I_dash,b_min,inf);     
    
    N_large=4;
    I_dash2=0;
    for kk=0:N_large
        I_dash2=I_dash2+E*(-1)^kk*qq^(-1-kk)/(2*rho^(kk/2+1/2))*igamma(kk/2+1/2,rho*b_min^2);
    end
    
    
    I_dash=I_dash1-I_dash2;
    
    ber_ana=[ber_ana omega/2-A*(I_1+C.*I_dash)];
    
% output simulation progress-----------------------------------------------
    disp(['Simulation:solve for ' num2str(pp) ' out of ' num2str(length(gamma_bar))]);
end

%% plotting the results
gamma_bar_db=10*log10(gamma_bar);
figure(1)
semilogy(gamma_bar_db,ber_sim,'--o','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
semilogy(gamma_bar_db,ber_ana,'-','color',[0 0 1],'LineWidth',2);hold on;
xlim([-10 30])

