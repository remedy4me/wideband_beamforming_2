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
Fpb=(0.16:0.01:0.32);                   %通带
Fsb_left=(0:0.01:0.13);
Fsb_right=(0.35:0.01:0.50);
Fsb=[Fsb_left Fsb_right];               %阻带
Ftb_left=(0.14:0.01:0.15);
Ftb_right=(0.33:0.01:0.34);
Ftb=[Ftb_left Ftb_right];               %过渡带
Fb=0:0.01:0.5;

fb=Fb*fs;
fpb=Fpb*fs;
fsb=Fsb*fs;
ftb=Ftb*fs;

pb_index=round([Fpb(1)/0.01+1: Fpb(end)/0.01+1]);
sb_index=round([Fsb_left(1)/0.01+1:Fsb_left(end)/0.01+1 Fsb_right(1)/0.01+1:Fsb_right(end)/0.01+1]);
tb_index=round([Ftb_left(1)/0.01+1:Ftb_left(end)/0.01+1 Ftb_right(1)/0.01+1:Ftb_right(end)/0.01+1]);

theta=(-90:2:90);
thetaML=(-8:2:28);
thetaSL_left=(-90:2:-12);
thetaSL_right=(32:2:90);
thetaSL=[thetaSL_left thetaSL_right];

thetaML_index=[find(theta==thetaML(1)):find(theta==thetaML(end))];
thetaSL_left_index=[find(theta==thetaSL_left(1)):find(theta==thetaSL_left(end))];
thetaSL_right_index=[find(theta==thetaSL_right(1)):find(theta==thetaSL_right(end))];
thetaSL_index=[thetaSL_left_index thetaSL_right_index];

thetas=10;
tau = d/c*sin(thetas/180*pi);
%%  %%%求期望波束%%%
SL=-25;
wd=exp(-1i*2*pi*fl*d'*sin(thetas/180*pi)/c)/M;
a=exp(-1i*2*pi*fl*d'*sin(theta/180*pi)/c);
cbf_p=wd'*a;
energy_cbf_P=20*log10(abs(cbf_p));

energy_cbf_PML=energy_cbf_P(:,thetaML_index);
energy_cbf_PSL_right=SL*ones(1,length(thetaSL_right));
energy_cbf_PSL_left=SL*ones(1,length(thetaSL_left));

figure(1)
hold on
plot(theta,energy_cbf_P,'k-');
hold on
scatter(thetaML,energy_cbf_PML,'*');
hold on
plot(thetaSL_right,energy_cbf_PSL_right,'r');
hold on
plot(thetaSL_left,energy_cbf_PSL_left,'r');
legend('常规波束','期望主瓣','期望旁瓣')
title('期望波束')
xlabel('方位/(^o)')
ylabel('波束/dB')
ylim([-60 3])
xlim([-90 90])
grid on

 %% 全局法恒定主瓣优化设计
 p_pb_ideal=cbf_p(:,thetaML_index);
 
 L=64;
 N_fk=length(Fpb);
 pb_ideal_matrix=repmat(p_pb_ideal,N_fk,1);
 D=L/2;
 SL=-25;
 for ii=1:M
     
     Tm(ii)=-floor(tau(ii)/Ts+D)*Ts;    %为了对照书本写成floor ， 感觉round好点
 end
 
 cvx_begin
 variable hh(M*L)
 variable s_s(1)
 minimize (s_s)
 
 subject to
 
 for ii=1:length(fb)
     
     for jj=1:length(theta)
         
         tau_m=d/c*sin(theta(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_theta(ii,jj)=hh'*vec_u;
     end
     
 end
 
 B_trans_fk_thetaSL=B_fk_theta(tb_index,thetaSL_index);  % 过渡带旁瓣
 B_stop_fk_theta=B_fk_theta(sb_index,:);               %  阻带
 B_pass_fk_thetaSL=B_fk_theta(pb_index,thetaSL_index);  %通带旁瓣
 error_matrix=B_fk_theta(pb_index,thetaML_index)-pb_ideal_matrix;
 square_error=error_matrix(:)'*error_matrix(:);
 
 square_error <= s_s
 abs(B_trans_fk_thetaSL)<=10^(SL/20)
 abs(B_stop_fk_theta)<=10^(SL/20)
 abs(B_pass_fk_thetaSL)<=10^(SL/20)
 norm(hh,2)<=0.25
 
 cvx_end
        
for ii=1:L
   h_global(:,ii)=hh((ii-1)*M+1:ii*M) ;
end

%% 恒定主瓣响应 -频域各子带加权值优化设计
for ii=1:length(fpb)
    
    p_theta_ML=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaML/180*pi)/c);  
    p_theta_SL=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaSL/180*pi)/c);
 
    cvx_begin
    variable w_f(M,1) complex
    minimize (norm(w_f'*p_theta_ML-p_pb_ideal,2))
       
    subject to
    Beam_SL=w_f'*p_theta_SL;
    abs(Beam_SL)<=10^(SL/20)
    norm(w_f,2)<=0.4217
    cvx_end
    w_f_all(:,ii)=w_f;

end

 %%  FIR 滤波器优化设计
error_constraint=0.01;  
e_f_pass = exp(-1i*2*pi*(0:L-1).'*Fpb);
e_f_stop = exp(-1i*2*pi*(0:L-1).'*Fsb);
e_f_full = exp(-1i*2*pi*(0:L-1).'*Fb);  

delay_pb_matrix=exp(1i*2*pi*Tm'*Fpb*fs); %预延迟通带相位矩阵
delay_fb_matrix=exp(1i*2*pi*Tm'*Fb*fs);  %预延迟全频带相位矩阵
Hd_pass=conj(w_f_all).*delay_pb_matrix;
    
for ii=1:M

       Hd_pass_m=Hd_pass(ii,:);
       cvx_begin
       variable h(1,L)
       minimize(norm(h*e_f_pass - Hd_pass_m,2))
       
       subject to
       max(abs(h*e_f_stop)) <= error_constraint
       cvx_end
       h_m(ii,:)=h;
end

W_FK_global=h_global*e_f_pass;    
W_FK=h_m*e_f_pass ; %FIRl滤波器通带频响（未预延迟）
H_m=h_m*e_f_full;  %FIR滤波器全频带响应（未预延迟）

 %%  全局设计法波束图计算
for ii=1:length(fb)
    
      for jj=1:length(theta)
          
       tau_m=d/c*sin(theta(jj)/180*pi);
       U_fk_theta=exp(-1i*2*pi*fb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
       vec_u= U_fk_theta(:);
       Beam_all_band(ii,jj)=hh'*vec_u;                   
   end
end
%% 子带设计法波束图计算
for ii=1:length(fpb)    
    w_ff=w_f_all(:,ii);
    a=exp(-1i*2*pi*fpb(ii)*d'*sin(theta/180*pi)/c);   
    Beam_fre_domain(ii,:)=w_ff'*a;
end
%% 分布设计法波束图计算
for ii=1:length(fb)
    w_ff=H_m(:,ii).*conj(delay_fb_matrix(:,ii));
    a=exp(-1i*2*pi*fb(ii)*d'*sin(theta/180*pi)/c);   
    Beam_time_domain(ii,:)=w_ff.'*a;
end
%% 子带等效加权向量范数计算
for ii=1:length(fpb)
   frq_only_norm(ii)=norm(w_f_all(:,ii),2) ; %子带设计法
   w_ff=W_FK(:,ii).*conj(delay_pb_matrix(:,ii));
   frq_time_comb_norm(ii)=norm(w_ff,2);  % 分步设计法
   w_ff_global=W_FK_global(:,ii).*conj(delay_pb_matrix(:,ii));
   frq_time_global(ii)=norm(w_ff_global,2); %全局设计法
end

%% 给定子带加权向量范数约束 ，子带设计法
for ii=1:length(fpb)
    
    p_theta_ML=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaML/180*pi)/c);  
    p_theta_SL=exp(-1i*2*pi*fpb(ii)*d'*sin(thetaSL/180*pi)/c);
 
    cvx_begin
    variable w_f_2(M,1) complex
    minimize (norm(w_f_2'*p_theta_ML-p_pb_ideal,2))
       
    subject to
    Beam_SL=w_f_2'*p_theta_SL;
    abs(Beam_SL)<=10^(SL/20)
    norm(w_f_2,2)<=frq_time_global(ii)
    cvx_end
    w_f_all_2(:,ii)=w_f_2;

end
%% 子带设计法波束图计算(加权向量进行具体约束)及设计向量范数巨酸
for ii=1:length(fpb)    
    w_ff=w_f_all_2(:,ii);
    frq_only_norm_2(ii)=norm(w_ff,2) ; %子带设计法
    a=exp(-1i*2*pi*fpb(ii)*d'*sin(theta/180*pi)/c);   
    Beam_fre_domain_2(ii,:)=w_ff'*a;
end

%% 主瓣逼近均方根误差计算
    fpb_MLbeam_error=Beam_all_band(pb_index,thetaML_index)-pb_ideal_matrix;
    fpb_ML_Beam_error_frq_domian=Beam_fre_domain(:,thetaML_index)-repmat(p_pb_ideal,length(fpb),1);
    fpb_ML_Beam_error_time_domian=Beam_time_domain(pb_index,thetaML_index)-repmat(p_pb_ideal,length(fpb),1);
    fpb_ML_Beam_error_frq_domian_2=Beam_fre_domain_2(:,thetaML_index)-repmat(p_pb_ideal,length(fpb),1);
    for ii=1:length(fpb)
        error_square_global(ii)=sqrt(fpb_MLbeam_error(ii,:)*fpb_MLbeam_error(ii,:)'/length(thetaML)); 
        error_square_fre(ii)=sqrt(fpb_ML_Beam_error_frq_domian(ii,:)*fpb_ML_Beam_error_frq_domian(ii,:)'/length(thetaML));
        error_square_time(ii)=sqrt(fpb_ML_Beam_error_time_domian(ii,:)*fpb_ML_Beam_error_time_domian(ii,:)'/length(thetaML)); 
        error_square_fre2(ii)=sqrt(fpb_ML_Beam_error_frq_domian_2(ii,:)*fpb_ML_Beam_error_frq_domian_2(ii,:)'/length(thetaML));
    end
%% PLOT 三维图
 figure();
 [degree,normalized_freq]=meshgrid(theta,fb/fs);
 surf(degree,normalized_freq,20*log10(abs(Beam_time_domain)));
 zlim([-60 0])
 xlim([-90 90])
 ylim([0 0.5])
 xlabel('方位/(^o)')
 ylabel('归一化频率')
 zlabel('波束/dB')
 title('恒定主瓣响应宽带波束图');
 
 
 figure();
 surf(degree,normalized_freq,20*log10(abs(Beam_all_band)));
 zlim([-60 0])
 xlim([-90 90])
 ylim([0 0.5])
 xlabel('方位/(^o)')
 ylabel('归一化频率')
 zlabel('波束/dB')
 title('恒定主瓣响应宽带波束图');
 
 figure();
 plot(1:length(fpb),error_square_fre,'k-s');
 hold on
 plot(1:length(fpb),error_square_time,'k:+');
 plot(1:length(fpb),error_square_global,'b-o');
 plot(1:length(fpb),error_square_fre2,'r-x');
 title('主瓣均方根误差')
 xlabel('频率序号')
 ylabel('主瓣逼近均方根误差')
 legend('子带设计法1','分步设计法1','全局设计法','子带设计法2')    
 ylim([0.01 0.045])
 
 figure();
 plot(1:length(fpb),frq_only_norm,'k-s');
 
 hold on
 plot(1:length(fpb),frq_time_comb_norm,'k:+');
 plot(1:length(fpb),frq_time_global,'b-o');
 plot(1:length(fpb),frq_only_norm_2,'r-x');
 title('子带（等效）加权向量范数')
 xlabel('频率序号k')
 ylabel('L2范数大小')
 legend('子带设计法1','分步设计法1','全局设计法','子带设计法2')   
 ylim([0.3 0.6])