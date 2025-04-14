clc; 
clear; 
addpath('./algorithm'); 
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 50; 
N = 3000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 

% TMSDg = zeros(1,N); 
% TMSDb = zeros(1,N); 
% TMSDr = zeros(1,N); 
% TMSDk = zeros(1,N); 
% TMSDm = zeros(1,N); 

mu1 =0.05;
mu2 =0.1;
mu3 =0.15;
mu4 =0.2;
lamada = 0.4;

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for k = 1:K  
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
L = 4; 
unknown_w = unifrnd(-0.5,0.5,2*L,1) + 1j*unifrnd(-0.5,0.5,2*L,1);
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
       u2 = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N)); 
%       u2 = sqrt(0.4)*randn(1,N) + sqrt(0.6)*1j*randn(1,N); 
%       u2 = Noncircular_CG(N);

      Input = u2;   
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D = filter(conj(unknown_w(1:L)), 1, Input) + filter(conj(unknown_w(L+1:2*L)), 1, conj(Input)); %output
     noise = Complex_CG(N,0.008,0.002,0,8,2);
     Output = D + noise;
     Lw =  length(unknown_w); 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
   
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w,Lw,lamada,1,mu4,-10);
    MSDg = MSDg + Da.MSD;  
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w,Lw,lamada,0.1,mu4,-10);
    MSDb = MSDb + Da.MSD;
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w,Lw,lamada,0.01,mu4,-10);
    MSDr = MSDr + Da.MSD;
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w,Lw,lamada,0.001,mu4,-10);
    MSDk = MSDk + Da.MSD;

    
%    [Da]= T_WLCHS(Input,noise,unknown_w,L,mu4,lamada);
%    TMSDg = TMSDg + Da.TMSD;  
%    [Da]= T_WLCHS(Input,noise,unknown_w,L,mu4,lamada);
%    TMSDb = TMSDb + Da.TMSD;  
%    [Da]= T_WLCHS(Input,noise,unknown_w,L,mu4,lamada);
%    TMSDr = TMSDr + Da.TMSD;  
%    [Da]= T_WLCHS(Input,noise,unknown_w,L,mu4,lamada);
%    TMSDk = TMSDk + Da.TMSD;  

end

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
MSDg = MSDg/K;  
MSDb = MSDb/K;  
MSDr = MSDr/K;  
MSDk = MSDk/K;  
MSDm = MSDm/K;  

% TMSDg = TMSDg/K;  
% TMSDb = TMSDb/K;  
% TMSDr = TMSDr/K;  
% TMSDk = TMSDk/K;  
% TMSDm = TMSDm/K;  

n = 1:N; 
figure(1);
plot(n,10*log10(MSDg(n)),'g','linewidth',1);
hold on;
plot(n,10*log10(MSDb(n)),'b','linewidth',1);
hold on;
plot(n,10*log10(MSDr(n)),'r','linewidth',1);
hold on;
plot(n,10*log10(MSDk(n)),'k','linewidth',1);
hold on;


% plot(10*log10(TMSDg),'g','linewidth',1.5);
% hold on;
% plot(10*log10(TMSDb),'b','linewidth',1.5);
% hold on;
% plot(10*log10(TMSDr),'r','linewidth',1.5);
% hold on;
% plot(10*log10(TMSDk),'k','linewidth',1.5);
% hold on;

xlabel('iteration','FontSize',15,'FontName','Times New Roman');
ylabel('MSD(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('Simulation (\beta=1)','Simulation (\beta=0.1)','Simulation (\beta=0.01)','Simulation (\beta=0.001)','Theory (\beta=1)','Theory (\beta=0.1)','Theory (\beta=0.01)','Theory (\beta=0.001)');
set(h,'FontSize',10,'FontName','Times New Roman');