clc; 
clear; 
addpath('./algorithm'); 
addpath('./signal'); 
% ¡­¡­¡­¡­¡­¡­¡­¡­parameter¡­¡­¡­¡­¡­¡­¡­¡­¡­
K = 100; 
N = 3000; 
MSDg = zeros(1,N); 
MSDb = zeros(1,N); 
MSDr = zeros(1,N); 
MSDk = zeros(1,N); 
MSDm = zeros(1,N); 
MSDy = zeros(1,N); 

lamada = 0.5;
alpha = 5;
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
for k = 1:K  
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
L = 8; 
unknown_w = unifrnd(-0.5,0.5,L,1) + 1j*unifrnd(-0.5,0.5,L,1);
  
      u2 = sqrt(0.5)*(randn(1,N) + 1j*randn(1,N)); 
%      u2 = Noncircular_CG(N);
      Input = u2;   
% ¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­     
     D = filter(conj(unknown_w),1,Input); 
     noise = Complex_CG(N,0.008,0.002,0.01,8,2);
     Output = D + noise;

unknown_w2= [unknown_w;zeros(L,1)];
Lw = length(unknown_w2); 
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w2,Lw,0.5,0.8,0.2,5);
    MSDg = MSDg + Da.MSD;  
    [Da]= WLCHS(Input,Output,D,unknown_w2,Lw,0.2,lamada);
    MSDb = MSDb + Da.MSD;
    [Da]= NCCHS2(Input,Output,noise,D,unknown_w2,Lw,0.5,0.2,0.1,5);
    MSDr = MSDr + Da.MSD;
    [Da]= WLCHS(Input,Output,D,unknown_w2,Lw,0.1,lamada);
    MSDk = MSDk + Da.MSD;
    [Da] = NCCHS2(Input,Output,noise,D,unknown_w2,Lw,0.5,0.8,0.05,-10);
    MSDm = MSDm + Da.MSD;
    [Da] = WLCHS(Input,Output,D,unknown_w2,Lw,0.05,lamada);
    MSDy = MSDy + Da.MSD;

  
end

%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
%¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­¡­
MSDg = MSDg/K;  
MSDb = MSDb/K;  
MSDr = MSDr/K;  
MSDk = MSDk/K;  
MSDm = MSDm/K;  
MSDy = MSDy/K;  


n = 1:1:N; 
figure(1);
plot(n,10*log10(MSDg(n)),'--g','linewidth',1);
hold on;
plot(n,10*log10(MSDb(n)),'g','linewidth',1);
hold on;
plot(n,10*log10(MSDr(n)),'--b','linewidth',1);
hold on;
plot(n,10*log10(MSDk(n)),'b','linewidth',1);
hold on;
plot(n,10*log10(MSDm(n)),'--r','linewidth',1);
hold on;
plot(n,10*log10(MSDy(n)),'r','linewidth',1);
hold on;
xlabel('iteration','FontSize',15,'FontName','Times New Roman');
ylabel('MSD(dB)','FontSize',15,'FontName','Times New Roman');
h=legend('WL-NCAFHS \mu=0.02','WL-CHSAF \mu=0.02','WL-NCAFHS \mu=0.05','WL-CHSAF \mu=0.05 ','WL-NCAFHS \mu=0.1','WL-CHSAF \mu=0.1');
set(h,'FontSize',10,'FontName','Times New Roman');