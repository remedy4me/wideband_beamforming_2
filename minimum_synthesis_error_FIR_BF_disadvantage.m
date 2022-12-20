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
Fpb=(0.16:0.005:0.32);                   %通带
Fb=0:0.01:0.5;

fb=Fb*fs;
fpb=Fpb*fs;

theta=(-90:2:90);
thetaML=(-8:2:28);
thetaML_2=(-2:2:22);
thetaSL_left=(-90:2:-12);
thetaSL_right=(32:2:90);
thetaSL=[thetaSL_left thetaSL_right];

thetaML_index=[find(theta==thetaML(1)):find(theta==thetaML(end))];
thetaML_2_index=[find(theta==thetaML_2(1)):find(theta==thetaML_2(end))];

thetas=10;
tau = d/c*sin(thetas/180*pi);
%%  %%%求期望波束%%%
SL=-30;
wd_fl=exp(-1i*2*pi*fl*d'*sin(thetas/180*pi)/c)/M;
a_fl=exp(-1i*2*pi*fl*d'*sin(theta/180*pi)/c);
cbf_p1=wd_fl'*a_fl;
energy_cbf_P=20*log10(abs(cbf_p1));
energy_cbf_P1ML=energy_cbf_P(:,thetaML_index);


wd_f2=exp(-1i*2*pi*(fl+fu)/2*d'*sin(thetas/180*pi)/c)/M;
a_f2=exp(-1i*2*pi*(fl+fu)/2*d'*sin(theta/180*pi)/c);
cbf_p2=wd_f2'*a_f2;
energy_cbf_P=20*log10(abs(cbf_p2));
energy_cbf_P2ML=energy_cbf_P(:,thetaML_2_index);


energy_cbf_PSL_right=SL*ones(1,length(thetaSL_right));
energy_cbf_PSL_left=SL*ones(1,length(thetaSL_left));

figure(1)
hold on
plot(thetaML,energy_cbf_P1ML,'b--o');
hold on
plot(thetaML_2,energy_cbf_P2ML,'r-s');
hold on
plot(thetaSL_right,energy_cbf_PSL_right,'k--');
hold on
plot(thetaSL_left,energy_cbf_PSL_left,'k--');
legend('f0=f1','f0=(fl+fu)/2','期望旁瓣级')
title('期望主瓣响应与期望旁瓣')
xlabel('方位/(^o)')
ylabel('期望波束/dB')
ylim([-70 3])
xlim([-90 90])
grid on

 %% 全局法恒定主瓣优化设计 选择频率 f0 的常规波束
 p_pb_ideal=cbf_p1(:,thetaML_index);
 
 L=15;
 N_fk=length(Fpb);
 pb_ideal_matrix=repmat(p_pb_ideal,N_fk,1);
 D=(L-1)/2;
 for ii=1:M
     Tm(ii)=-round(tau(ii)/Ts+D)*Ts;    
 end
 
 cvx_begin
 variable hh(M*L)
 variable s_s(1)
 expression B_fk_thetaML(N_fk,length(thetaML))
 expression B_fk_thetaSL(N_fk,length(thetaSL))
 expression error_matrix(N_fk,length(thetaML))
 minimize (s_s)
 
 subject to
 
 for ii=1:length(fpb)
     
     for jj=1:length(thetaML)
         
         tau_m=d/c*sin(thetaML(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetaML(ii,jj)=hh'*vec_u;
     end
     
     for jj=1:length(thetaSL)
         
         tau_m=d/c*sin(thetaSL(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetaSL(ii,jj)=abs(hh'*vec_u);
     end
     
 end
 
 error_matrix=abs(B_fk_thetaML-pb_ideal_matrix);
 
 max(error_matrix) <= s_s
 B_fk_thetaSL <= 10^(SL/20)

 cvx_end
        
 
  %% 全局法恒定主瓣优化设计 选择频率 (fl+fu)/2 的常规波束
 p_pb_ideal=cbf_p2(:,thetaML_2_index);
 pb_ideal_matrix=repmat(p_pb_ideal,N_fk,1);

 cvx_begin
 variable hh_2(M*L)
 variable s_s(1)
 expression B_fk_thetaML(N_fk,length(thetaML_2))
 expression B_fk_thetaSL(N_fk,length(thetaSL))
 expression error_matrix(N_fk,length(thetaML_2))
 minimize (s_s)
 
 subject to
 
 for ii=1:length(fpb)
     
     for jj=1:length(thetaML_2)
         
         tau_m=d/c*sin(thetaML_2(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetaML(ii,jj)=hh_2'*vec_u;
     end
     
     for jj=1:length(thetaSL)
         
         tau_m=d/c*sin(thetaSL(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetaSL(ii,jj)=hh_2'*vec_u;
     end
     
 end
 
 error_matrix=abs(B_fk_thetaML-pb_ideal_matrix);
 
 max(error_matrix) <= s_s
 abs(B_fk_thetaSL)<=10^(SL/20)

 cvx_end
 %%  全局设计法波束图计算
for ii=1:length(fpb)
    
      for jj=1:length(theta)
          
       tau_m=d/c*sin(theta(jj)/180*pi);
       U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
       vec_u= U_fk_theta(:);
       Beam_pb(ii,jj)=hh'*vec_u;    
       Beam_pb_2(ii,jj)=hh_2'*vec_u;  
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


figure();
[degree,normalized_freq]=meshgrid(theta,fpb/fs);
surf(degree,normalized_freq,20*log10(abs(Beam_pb)));
zlim([-80 0])
xlim([-90 90])
xlabel('方位/(^o)')
ylabel('归一化频率')
zlabel('波束/dB')
title('3D显示');

figure();
for ii=1:length(fpb)
    
    plot(theta,20*log10(abs(Beam_pb_2(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('FIR波束形成器各子带波束图');

error_beam=abs(Beam_pb(:,thetaML_index)-Beam_pb(1,thetaML_index));
delta_max=20*log10(max(max(error_beam)));