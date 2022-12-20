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
Fpb=(0.16:0.005:0.32);                   %通带
fpb=Fpb*fs;
N_fk=length(fpb);

%% 
% theta_ii=[-50:0];
% for iii=1:length(theta_ii)
    %% 方位参数
    theta=(-90:2:90);

    SL_all_band_left=round(interp1([Fpb(1),Fpb(end)],[-15,-5],Fpb));
    SL_all_band_right=round(interp1([Fpb(1),Fpb(end)],[35,25],Fpb));

    for ii=1:N_fk
        length_pb_theta_SL(ii)=length( [-90:6:SL_all_band_left(ii) SL_all_band_right(ii):6:90]);
    end

    SL=-30;
    thetas=10;
    thetai_1=-50;
    thetai_2=-30;
    tau_sig = d/c*sin(thetas/180*pi);
    tau_inter_1=d/c*sin(thetai_1/180*pi);
    tau_inter_2=d/c*sin(thetai_2/180*pi);
    INR=30;   % 干噪比30dB
    SNR=0;
    p_white_noise=1;
    p_inter=10^(INR/10)*p_white_noise;
    p_signal=10^(SNR/10)*p_white_noise;
    

    %% 滤波器参数
    L=16;
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
    V=chol(R);
    %% 自适应FIR波束形成器设计

    cvx_begin
    variable hh_1(M*L) 
    variable s_s(1)
    expression B_fk_thetas(N_fk)
    minimize (s_s)

    subject to

    for ii=1:N_fk

        tau_m=d/c*sin(thetas/180*pi);
        U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
        vec_u= U_fk_theta(:);
        B_fk_thetas(ii)=hh_1.'*vec_u;

    end


    B_fk_thetas==1
    norm(V*hh_1,2)<= s_s
    norm(hh_1,2)<=h_morm_restraint

    cvx_end

    SINR_out_noSL=10*log10(abs(hh_1.'*Rd_sig*hh_1/(hh_1.'*(Rd_inter+R_white_noise)*hh_1)));
    %% 自适应FIR波束形成器+旁瓣控制设计

    cvx_begin
    variable hh_2(M*L) 
    variable s_s(1)
    expression B_fk_thetas(N_fk)
    expression B_fk_thetaSL(sum(length_pb_theta_SL))
    minimize (s_s)

    subject to

    for ii=1:N_fk

        theta_SL=[-90:6:SL_all_band_left(ii) SL_all_band_right(ii):6:90];

        for jj=1:length(theta_SL)
            tau_m=d/c*sin(theta_SL(jj)/180*pi);
            U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts));
            vec_u=U_fk_theta(:);
            B_fk_thetaSL(sum(length_pb_theta_SL(1:ii))-length_pb_theta_SL(1)+jj)=abs(hh_2.'*vec_u);
        end

        U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_sig'+[0:L-1]*Ts));
        vec_u= U_fk_theta(:);
        B_fk_thetas(ii)=hh_2.'*vec_u;

    end

    B_fk_thetas==1
    norm(V*hh_2,2)<= s_s
    B_fk_thetaSL<=10^(SL/20)
    norm(hh_2,2)<=h_morm_restraint

    cvx_end

    SINR_out_2=10*log10(abs(hh_2.'*Rd_sig*hh_2/(hh_2.'*(Rd_inter+R_white_noise)*hh_2)));
     %%  全局设计法波束图计算

    for ii=1:N_fk

          for jj=1:length(theta)

           tau_m=d/c*sin(theta(jj)/180*pi);
           U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
           vec_u= U_fk_theta(:);
           Beam_pb_1(ii,jj)=hh_1.'*vec_u; 
           Beam_pb_2(ii,jj)=hh_2.'*vec_u; 
          end   
    end

% end

figure();
for ii=1:N_fk
    
    plot(theta,20*log10(abs(Beam_pb_1(ii,:))));
    hold on
end
xlim([-90 90])
 ylim([-120 5])
title('重叠显示');
ylabel('波束/dB')
xlabel('方位/(°)')

% PLOT 三维图
figure();
[degree,normalized_freq]=meshgrid(theta,fpb/fs);
surf(degree,normalized_freq,20*log10(abs(Beam_pb_1)));
zlim([-100 0])
xlim([-90 90])
xlabel('方位/(^o)')
ylabel('归一化频率')
zlabel('波束/dB')
title('3D显示');

figure();
for ii=1:N_fk
    
    plot(theta,20*log10(abs(Beam_pb_2(ii,:))));
    hold on
end
xlim([-90 90])
 ylim([-120 5])
title('重叠显示');
ylabel('波束/dB')
xlabel('方位/(°)')

% PLOT 三维图
figure();
[degree,normalized_freq]=meshgrid(theta,fpb/fs);
surf(degree,normalized_freq,20*log10(abs(Beam_pb_2)));
zlim([-100 0])
xlim([-90 90])
xlabel('方位/(^o)')
ylabel('归一化频率')
zlabel('波束/dB')
title('3D显示');

trap_dB_thetai_1=20*log10(abs(Beam_pb_1(:,find(theta==thetai_1))));
trap_dB_thetai_2=20*log10(abs(Beam_pb_1(:,find(theta==thetai_2))));
trap_dB_const_SL_thetai_1=20*log10(abs(Beam_pb_2(:,find(theta==thetai_1))));
trap_dB_const_SL_thetai_2=20*log10(abs(Beam_pb_2(:,find(theta==thetai_2))));

figure();
plot(Fpb,trap_dB_thetai_1,'b--o','MarkerFaceColor','k');
hold on
plot(Fpb,trap_dB_thetai_2,'b--o');
plot(Fpb,trap_dB_const_SL_thetai_1,'k-s');
plot(Fpb,trap_dB_const_SL_thetai_2,'k-^');
xlabel('归一化频率')
legend('旁瓣未控制,\theta1=-50°','旁瓣未控制,\theta2=-30°','旁瓣控制,\theta1=-50°','旁瓣控制,\theta2=-30°')
ylabel('凹槽深度/dB')
ylim([-120,-55])
