clc;
clear all;
close all;
%%
M=12;
f0=1000;
fl=f0/2;
fu=f0;
fs=3.125*f0;
Ts=1/fs;
c=1500;
d=c/fu/2*[0:M-1];
%%
Fpb=(0.16:0.005:0.32); %通带
Fb=0:0.01:0.5;
fb=Fb*fs;
fpb=Fpb*fs;

pb_index=round([Fpb(1)/0.005+1: Fpb(end)/0.005+1]);
theta=(-90:2:90);
thetaML=(-8:2:28);
thetaSL_left=(-90:2:-12);
thetaSL_right=(32:2:90);
thetaSL=[thetaSL_left thetaSL_right];

thetaML_index=[find(theta==thetaML(1)):find(theta==thetaML(end))];
thetaSL_index=[find(theta==thetaSL_left(1)):find(theta==thetaSL_left(end))...
    ,find(theta==thetaSL_right(1)):find(theta==thetaSL_right(end))];
thetas=10;
tau_sig = d/c*sin(thetas/180*pi);
thetai_1=-50;
thetai_2=-30;
tau_inter_1=d/c*sin(thetai_1/180*pi);
tau_inter_2=d/c*sin(thetai_2/180*pi);
INR=30;   % 干噪比30dB
SNR=0;
p_white_noise=1;
p_inter=10^(INR/10)*p_white_noise;
p_signal=10^(SNR/10)*p_white_noise;
%% 滤波器系数
 L=15;
 N_fk=length(Fpb);
 D=(L-1)/2;
 for ii=1:M
     Tm(ii)=-round(tau_sig(ii)/Ts+D)*Ts;    %为了对照书本写成floor ， 感觉round好点
 end
    %% 构造理想数据协方差矩阵
    Bs=fu-fl;
    fc=(fu+fl)/2;
    for ii=1:M
        for jj=1:L
            tau_ML_sig(ii,jj)=Tm(ii)+tau_sig(ii)+(jj-1)*Ts;
            tau_ML_inter_1(ii,jj)=Tm(ii)+tau_inter_1(ii)+(jj-1)*Ts;
            tau_ML_inter_2(ii,jj)=Tm(ii)+tau_inter_2(ii)+(jj-1)*Ts;
        end
    end

    for ii=1:M
        for jj=1:L
            for kk=1:M
                for qq=1:L
                    Rd_sig(ii+(jj-1)*M,kk+(qq-1)*M)=p_signal/Bs*sinc(Bs*(tau_ML_sig(ii,jj)-tau_ML_sig(kk,qq)))*...
                        exp(-1i*2*pi*fc*(tau_ML_sig(ii,jj)-tau_ML_sig(kk,qq)));
                    Rd_inter_1(ii+(jj-1)*M,kk+(qq-1)*M)=p_inter/Bs*sinc(Bs*(tau_ML_inter_1(ii,jj)-tau_ML_inter_1(kk,qq)))*...
                        exp(-1i*2*pi*fc*(tau_ML_inter_1(ii,jj)-tau_ML_inter_1(kk,qq)));
                    Rd_inter_2(ii+(jj-1)*M,kk+(qq-1)*M)=p_inter/Bs*sinc(Bs*(tau_ML_inter_2(ii,jj)-tau_ML_inter_2(kk,qq)))*...
                        exp(-1i*2*pi*fc*(tau_ML_inter_2(ii,jj)-tau_ML_inter_2(kk,qq)));
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
    
    Rd_inter=real(Rd_inter_1)+real(Rd_inter_2);
    Rd_sig=real(Rd_sig);
    R_white_noise=real(R_white_noise);
    R=Rd_sig+Rd_inter+R_white_noise;
 %% 全局法恒定主瓣优化设计 选择频率 f0 的常规波束
 reference_fre_index=find(fpb==(fu+fl)/2);
 theta_s_index=find(thetaML==thetas);

 omega=zeros(M*L,M*L);
 phita=zeros(M*L,M*L);
 for ii=1:N_fk
     
     for jj=1:length(thetaML)
         
         tau_m=d/c*sin(thetaML(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         U_f0_theta=exp(-1i*2*pi*fpb(reference_fre_index)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u_fk= U_fk_theta(:);
         vec_u_f0= U_f0_theta(:);
         omega=omega+((vec_u_fk-vec_u_f0)*(vec_u_fk-vec_u_f0)')./(N_fk*length(thetaML));
         
     end
     
      for jj=1:length(thetaSL)
         
         tau_m=d/c*sin(thetaSL(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u_fk=U_fk_theta(:);
         phita=phita+(vec_u_fk*vec_u_fk')./(N_fk*length(thetaSL));
         
     end
 end
 omega=real(omega);
 phita=real(phita);
 
%%
 RMS_SL=-33.2;
 RMS_ML_error=-51.2;
 cvx_begin
 variable hh(M*L)
 variable s_s(1)
 expression B_fk_thetaML(1)
 minimize (s_s)
 
 subject to
 
 tau_m=d/c*sin(thetas/180*pi);
 U_fk_theta=exp(-1i*2*pi*fpb(reference_fre_index)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
 vec_u_f0_thetas= U_fk_theta(:);
 B_fk_thetaML=hh.'*vec_u_f0_thetas;
         
 B_fk_thetaML==1
 hh.'*R*hh <= s_s
 hh.'*phita*hh <= 10^(RMS_SL/10)
 hh.'*omega*hh <= 10^(RMS_ML_error/10)

 cvx_end
 

 %%  全局设计法波束图计算
for ii=1:length(fpb)
    
      for jj=1:length(theta)
          
       tau_m=d/c*sin(theta(jj)/180*pi);
       U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
       vec_u= U_fk_theta(:);
       Beam_pb(ii,jj)=hh'*vec_u;   
      end

end

figure();
for ii=1:length(fpb)
    
    plot(theta,20*log10(abs(Beam_pb(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('FIR波束形成器各子带波束图');

% err_1=20*log10(max(max(abs(Beam_pb(:,thetaML_index)-Beam_pb(reference_fre_index,thetaML_index)))));
% err_2=20*log10(max(max(abs(Beam_pb(:,thetaSL_index)))));