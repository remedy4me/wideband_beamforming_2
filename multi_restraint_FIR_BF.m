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
%% 频带参数
Fpb=(0.16:0.005:0.32); %通带
fpb=Fpb*fs;
N_fk=length(fpb);
%% 方位参数
theta=(-90:2:90);
theta_ML=(-12:2:32);
thetaSL_left=(-90:2:-16);
thetaSL_right=(36:2:90);
theta_SL=[thetaSL_left thetaSL_right];

thetaML_index=[find(theta==theta_ML(1)):find(theta==theta_ML(end))];
thetaSL_index=[find(theta==thetaSL_left(1)):find(theta==thetaSL_left(end))...
    ,find(theta==thetaSL_right(1)):find(theta==thetaSL_right(end))];

SL=-30;
thetas=10;
thetai_1=-30;
thetai_2=-50;
tau_sig = d/c*sin(thetas/180*pi);
tau_inter_1=d/c*sin(thetai_1/180*pi);
tau_inter_2=d/c*sin(thetai_2/180*pi);
INR=30;   % 干噪比30dB
SNR=0;
p_white_noise=1;
p_inter=10^(INR/10)*p_white_noise;
p_signal=10^(SNR/10)*p_white_noise;
%% 滤波器参数
L=32;
D=L/2;
h_morm_restraint=1;
for ii=1:M
    Tm(ii)=-round(tau_sig(ii)/Ts+D)*Ts;
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

%% 自适应FIR波束形成器
reference_fre_index=find(fpb==(fu+fl)/2);

theta_s_index=find(theta_ML==thetas);

omega=zeros(M*L,M*L);
 for ii=1:N_fk
     
     for jj=1:length(theta_ML)
         
         tau_m=d/c*sin(theta_ML(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         U_f0_theta=exp(-1i*2*pi*fpb(reference_fre_index)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u_fk= U_fk_theta(:);
         vec_u_f0= U_f0_theta(:);
         omega=omega+((vec_u_fk-vec_u_f0)*(vec_u_fk-vec_u_f0)')./(N_fk*length(theta_ML));
         
     end
    
 end
 omega=real(omega);
 
%%    
RMS_ML_error=-30;
cvx_begin
variable hh(M*L)
variable s_s(1)
expression B_fk_thetaML(N_fk,length(theta_ML))
expression B_fk_thetaSL(N_fk,length(theta_SL))
minimize (s_s)

subject to

for ii=1:N_fk
    
    for jj=1:length(theta_ML)
        tau_m=d/c*sin(theta_ML(jj)/180*pi);
        U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts));
        vec_u=U_fk_theta(:);
        B_fk_thetaML(ii,jj)=hh.'*vec_u;
    end
    
    for jj=1:length(theta_SL)
        tau_m=d/c*sin(theta_SL(jj)/180*pi);
        U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts));
        vec_u=U_fk_theta(:);
        B_fk_thetaSL(ii,jj)=hh.'*vec_u;
    end

end

B_fk_thetaML(reference_fre_index,theta_s_index)==1 
hh.'*R*hh<= s_s                         
hh.'*omega*hh <= 10^(RMS_ML_error/10)   %% 恒定主瓣误差
max(abs(B_fk_thetaSL))<=10^(-25/20)     %% 旁瓣
norm(hh,2)<=10^(-10/20)                 %% 滤波器系数范数

cvx_end

SINR_out=10*log10(abs(hh.'*Rd_sig*hh/(hh.'*(Rd_inter+R_white_noise)*hh)));

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
for ii=1:length(fpb)
    plot(theta,20*log10(abs(Beam_pb(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('FIR波束形成器各子带波束图')