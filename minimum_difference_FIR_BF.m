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

thetas=10;
tau = d/c*sin(thetas/180*pi);
 %% 全局法恒定主瓣优化设计 选择频率 f0 的常规波束
 reference_fre_index=find(fpb==(fu+fl)/2);
 theta_s_index=find(thetaML==thetas);
 SL=-30;
 L=15;
 N_fk=length(Fpb);
 D=(L-1)/2;
 for ii=1:M
     Tm(ii)=-round(tau(ii)/Ts+D)*Ts;    %为了对照书本写成floor ， 感觉round好点
 end
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
         B_fk_thetaSL(ii,jj)=hh'*vec_u;
         
     end
     
 end
 
 pb_ideal_matrix=repmat(B_fk_thetaML(reference_fre_index,:),N_fk,1);
 error_matrix=abs(B_fk_thetaML-pb_ideal_matrix);
 
 B_fk_thetaML(reference_fre_index,theta_s_index)==1
 error_matrix <= s_s
 abs(B_fk_thetaSL)<=10^(SL/20)

 cvx_end
 
 %%
 cvx_begin
 variable hh_2(M*L)
 variable s_ss(1)
 expression B_fk_thetaML(1)
 expression B_fk_thetaSL(N_fk,length(thetaSL))
 minimize (s_ss)
 
 subject to
 
         tau_m=d/c*sin(thetas/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(reference_fre_index)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u_f0_thetas= U_fk_theta(:);
         B_fk_thetaML=hh_2.'*vec_u_f0_thetas;
         
 for ii=1:length(fpb)

     for jj=1:length(thetaSL)
         
         tau_m=d/c*sin(thetaSL(jj)/180*pi);
         U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
         vec_u= U_fk_theta(:);
         B_fk_thetaSL(ii,jj)=hh_2'*vec_u;
         
     end

 end

 B_fk_thetaML==1
 hh_2.'*omega*hh_2 <= s_ss
 abs(B_fk_thetaSL)<=10^(SL/20)

 cvx_end
 ML_avg_square_error=s_ss;
  %%
 cvx_begin
 variable hh_3(M*L)
 variable s_sss(1)
 expression B_fk_thetaML(1)
 minimize (s_sss)
 
 subject to
 
 tau_m=d/c*sin(thetas/180*pi);
 U_fk_theta=exp(-1i*2*pi*fpb(reference_fre_index)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
 vec_u_f0_thetas= U_fk_theta(:);
 B_fk_thetaML=hh_3.'*vec_u_f0_thetas;
         
 B_fk_thetaML==1
 hh_3.'*phita*hh_3 <=s_sss
 hh_3.'*omega*hh_3 <= ML_avg_square_error

 cvx_end
 %%  全局设计法波束图计算
for ii=1:length(fpb)
    
      for jj=1:length(theta)
          
       tau_m=d/c*sin(theta(jj)/180*pi);
       U_fk_theta=exp(-1i*2*pi*fpb(ii)*(Tm'+tau_m'+[0:L-1]*Ts)) ;
       vec_u= U_fk_theta(:);
       Beam_pb(ii,jj)=hh'*vec_u;  
       Beam_pb_2(ii,jj)=hh_2'*vec_u; 
       Beam_pb_3(ii,jj)=hh_3'*vec_u;
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
for ii=1:length(fpb)
    
    plot(theta,20*log10(abs(Beam_pb_2(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('FIR波束形成器各子带波束图');

figure();
for ii=1:length(fpb)
    
    plot(theta,20*log10(abs(Beam_pb_3(ii,:))));
    hold on 
end
xlim([-90 90])
ylabel('波束/dB')
xlabel('方位/(°)')
title('FIR波束形成器各子带波束图');

error_beam_1=abs(Beam_pb(:,thetaML_index)-Beam_pb(reference_fre_index,thetaML_index));
delta_max_1=20*log10(max(max(error_beam_1)));
