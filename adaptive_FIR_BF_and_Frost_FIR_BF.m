clc;
clear all;
close all;
%% 阵列参数
M=12;
f0=1000;
fl=f0/2;
fu=f0;
fs=3.125*f0;
Ts=1/fs;
c=1500;
d=c/fu/2*[0:M-1];

% SNR_all=[-20:10];
% for iii=1:length(SNR_all)
    %% 频带参数
    Fpb=(0.16:0.005:0.32);                   %通带
    fpb=Fpb*fs;
    N_fk=length(fpb); 
    %% 滤波器参数
    L=8;
    D=L/2;
    h_morm_restraint=1; 
    %% 方位参数
    theta=(-90:2:90);
    thetas=10;
    thetai=-40;
    tau_sig = d/c*sin(thetas/180*pi);  
    tau_inter=d/c*sin(thetai/180*pi);
    INR=30;   % 干噪比30dB
%   SNR=SNR_all(iii);
    SNR=0;
    p_white_noise=1;
    p_inter=10^(INR/10)*p_white_noise;
    p_signal=10^(SNR/10)*p_white_noise;

     for ii=1:M    
         Tm(ii)=-round(tau_sig(ii)/Ts)*Ts;    
     end

    %% 构造理想数据协方差矩阵

    Bs=fu-fl;      
    fc=(fu+fl)/2;
    for ii=1:M
        for jj=1:L
            tau_ML_sig(ii,jj)=Tm(ii)+tau_sig(ii)+(jj-1)*Ts;
            tau_ML_inter(ii,jj)=Tm(ii)+tau_inter(ii)+(jj-1)*Ts;
        end
    end

    for ii=1:M
        for jj=1:L
            for kk=1:M
                for qq=1:L
                    Rd_sig(ii+(jj-1)*M,kk+(qq-1)*M)=p_signal/Bs*sinc(Bs*(tau_ML_sig(ii,jj)-tau_ML_sig(kk,qq)))*...
                        exp(-1i*2*pi*fc*(tau_ML_sig(ii,jj)-tau_ML_sig(kk,qq)));
                    Rd_inter(ii+(jj-1)*M,kk+(qq-1)*M)=p_inter/Bs*sinc(Bs*(tau_ML_inter(ii,jj)-tau_ML_inter(kk,qq)))*...
                        exp(-1i*2*pi*fc*(tau_ML_inter(ii,jj)-tau_ML_inter(kk,qq)));
                    if ii==kk
                        R_white_noise(ii+(jj-1)*M,kk+(qq-1)*M)=p_white_noise/Bs*sinc(Bs*(jj-qq)*Ts)*...
                        exp(-1i*2*pi*fc*(jj-qq)*Ts);
                    else
                        R_white_noise(ii+(jj-1)*M,kk+(qq-1)*M)=0;
                    end
                end
            end
        end
    end

    R=Rd_sig+Rd_inter+R_white_noise;
    V=chol(R);
     %% 自适应FIR波束形成器设计

      cvx_begin 
      variable hh(M*L) complex
      variable s_s(1)
      expression B_fk_thetas(N_fk)
      minimize (s_s)

      subject to

       for ii=1:N_fk

         tau_m=d/c*sin(thetas/180*pi);         
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetas(ii)=hh.'*vec_u;

       end

         p_out=norm(conj(V)*hh,2);   
         B_fk_thetas==1
         p_out<= s_s
         norm(hh,2)<=h_morm_restraint

       cvx_end
 
        SINR_out=10*log10(abs(hh.'*Rd_sig*conj(hh)/(hh.'*(Rd_inter+R_white_noise)*conj(hh))));
       
% end

% figure();
% plot(SNR_all,10*log10(abs(SINR_out)),'k-o');
% xlabel('SNR(in)/dB')
% ylabel('SINR(out)/dB')
% title('不同输入SNR时波束输出SINR')
% ylim([-15 25])

 %%  全局设计法波束图计算

for ii=1:N_fk
    
      for jj=1:length(theta)
          
       tau_m=d/c*sin(theta(jj)/180*pi);
       U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
       vec_u= U_fk_theta(:);
       Beam_pb(ii,jj)=hh.'*vec_u; 
       
      end   
end

figure();
for ii=1:N_fk
    
   plot(theta,20*log10(abs(Beam_pb(ii,:))));
   hold on
end
xlim([-90 90])  
% ylim([-60 5])
 title('重叠显示');
ylabel('波束/dB')
xlabel('方位/(°)')

% PLOT 三维图
 figure();
 [degree,normalized_freq]=meshgrid(theta,fpb/fs);
 surf(degree,normalized_freq,20*log10(abs(Beam_pb)));
 zlim([-100 0])
 xlim([-90 90])
 xlabel('方位/(^o)')
 ylabel('归一化频率')
 zlabel('波束/dB')
 title('3D显示');
