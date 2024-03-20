clear,clc
tic

% 
% lambda=1:1:200/3
% len=length(lambda);
% Q1=zeros(len,4);
% N2RC1=zeros(len,4);
% NRC1=zeros(len,4);
% separation_p1=zeros(len,4);
% separation_n1=zeros(len,4);
% separation_c1=zeros(len,4);
% x{1}=[F MG MP MG
%     F MP  G  MP
%     G VG  F MG
%     VG G MG G];
% x{2}=[MG MP P F
%     F MG P G
%     P MP VP F
%     G VG VG MG];
% x{3}=[G F MP MP
%     F P MP G
%     F VP P G
%     G VG AG MG];

disp('Step 2: numerical matrices')
Q(:,:,1)=[6 4 5 6
    7 6  7  4
    3 8  5 6
    8 7 7 6];
Q(:,:,2)=[5 6 4 6
    5 4  7  4
    7 8  5 6
    8 7 6 7];
Q(:,:,3)=[6 4 3 5
    5 6 3 7
    3 4 2 5
    7 8 8 6];
Q(:,:,4)=[7 5 4 4
    5 3 4 7
    5 2 3 7
    7 8 9 6];
[m,n]=size(Q(:,:,1));
t=4;

Q_mean=median(Q,3);
% q_mean = zeros(m,n);
% for k=1:3
%     q_mean = q_mean+q(:,:,k);
% end
% q_mean = q_mean/3;
% temp=[reshape(q_mean',1,[])];
%   % eval(['Y_mean''=  
% q_mean = sprintf('%4.2f  %4.2f  %4.2f   %4.2f \n', temp);
% disp('算术平均矩阵为')
% q_mean  

disp('步骤3，Q的arithmetic mean matrix：')

sprintf('%4.2f  %4.2f  %4.2f  %4.2f\n', Q_mean')

%%% 第3步，求加权矩阵的熵

for k=1:t 
    P(:,:,k)=Q(:,:,k)./sum(sum(Q(:,:,k)));
    E(k)=-sum(sum(P(:,:,k).*log(P(:,:,k))))/log(m*n);   
end

disp('Qk的熵为')
E

P=Q_mean/sum(sum(Q_mean));
E_mean=-sum(sum(P.*log(P)))/log(m*n) 

for k=1:t
    D(k) =abs( E(k)-E_mean);
    R(k) =E_mean./(E_mean+D(k));
end
disp('第4步，求4个偏离度')
D, R

disp('步骤8，weights of DMs 为：')

 lambda=R/sum(R) 
    
%  f{1}=[(0.65, 0.40) (0.40, 0.65) (0.50, 0.50) (0.65, 0.40)
%      (0.75, 0.30) (0.65, 0.40)  (0.75, 0.30)  (0.40, 0.65)
%      (0.30, 0.75) (0.85, 0.20)  (0.50, 0.50) (0.65, 0.40)
%      (0.85, 0.20) (0.75, 0.30)  (0.75, 0.30) (0.65, 0.40)];      

% alpha = 0:0.01:0.9156;
% len=length(alpha);
% Q1=zeros(len,4);
% NGU1=zeros(len,4);

% for row=1:len
mu{1}=[0.65 0.40 0.50 0.65   
     0.75 0.65 0.75 0.40     
     0.30 0.85 0.50 0.65
     0.85 0.75 0.75 0.65];
      
nu{1}=[0.40 0.65 0.50 0.40   
     0.30 0.40  0.30  0.65      
     0.75 0.20  0.50  0.40
     0.20 0.30  0.30  0.40];
 %  f{2}=[(0.50, 0.50) (0.65, 0.40) (0.40, 0.65) (0.65, 0.40)
%      (0.50, 0.50) (0.40, 0.65)  (0.75, 0.30)  (0.40, 0.65) 
%      (0.75, 0.30) (0.85, 0.20)  (0.50, 0.50) (0.65, 0.40)
%      (0.85, 0.20) (0.75, 0.30) (0.65, 0.40) (0.75, 0.30)];


mu{2}=[0.50 0.65 0.40 0.65   
      0.50 0.40 0.75 0.40
      0.75 0.85 0.50 0.65 
      0.85 0.75 0.65 0.75];
  
nu{2}=[0.50 0.40  0.65 0.40 
      0.50  0.65  0.30  0.65
      0.30  0.20  0.50  0.40
      0.20  0.30 0.40  0.30];
% f{3}=[(0.65, 0.40) (0.40, 0.65) (0.30, 0.75) (0.50, 0.50)
%     (0.50, 0.50) (0.65, 0.40) (0.30, 0.75) (0.75, 0.30)
%     (0.30, 0.75) (0.40, 0.65) (0.20, 0.85) (0.50, 0.50)
%     (0.75, 0.30) (0.85, 0.20) (0.85, 0.20) (0.65, 0.40)];  

mu{3}=[0.65 0.40 0.30 0.50    
    0.50 0.65 0.30 0.75 
    0.30 0.40 0.20 0.50  
    0.75 0.85 0.85 0.65 ];
  
nu{3}=[0.40 0.65 0.75 0.50   
     0.50 0.40 0.75  0.30   
     0.75 0.65 0.85  0.50  
     0.30 0.20 0.20  0.40];
mu{4}=[0.75 0.50 0.40 0.40
      0.50  0.30 0.40 0.75 
      0.50  0.20 0.30 0.75 
      0.75  0.85 0.95 0.65];
nu{4}=[0.30 0.50 0.65 0.65
       0.50 0.75 0.65 0.30
       0.50 0.85 0.75 0.30
       0.30 0.20 0.10 0.40];
 
 [m,n]=size(mu{1});

F=zeros(4,4,4);

disp('Step 8，Pythagorean fuzzy matrices');

for k=1:4
    temp=[reshape(mu{k}',1,[]); reshape(nu{k}',1,[])];
   eval(['F',num2str(k),'=  sprintf(''(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)\n'', temp)'])
end
clear F1 F2 F3 

%%%% 第2步，求加权决策矩阵

disp('Step 9，加权到决策者')

% weights=[0.3348, 0.3327, 0.3325];
% weight(k)=ones(4,4)*lambda(k)
% 
% Y=X;
% 
for k=1:4 
    weight(:,:,k)=ones(m,4)*lambda(k);
    xi(:,:,k)=sqrt(1-(1-mu{k}.^2).^weight(:,:,k)); 
   % xi(:,:,k)=mu{k}.^weight(:,:,k);
   % o(:,:,k)=sqrt(1-(1-nu{k}.^2).^weight(:,:,k));  
     o(:,:,k)=nu{k}.^weight(:,:,k);
   tmp=[reshape(xi(:,:,k)',1,[]);  reshape(o(:,:,k)',1,[])];
   eval(['H',num2str(k),'=  sprintf(''(%3.2f,  %3.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)\n'', tmp)']);
end


%%%% 转换成方案矩阵
disp('第10步，群决策矩阵')
xii=permute(xi,[3,2,1])
oi=permute(o,[3,2,1])
for k=1:4 
%    xii(i,:,:)=permute(xi,[3,2,1])
%    oi(i,:,:)=permute(o,[3,2,1])    
   temp=[reshape(xii(:,:,k)',1,[]);  reshape(oi(:,:,k)',1,[])];
   eval(['G',num2str(k),'=  sprintf(''(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)\n'', temp)']);
end

disp('第11步，加权到属性')

weights=[0.3,0.3,0.2,0.2];
weighted=ones(4,1)*weights

%Y=G;
for k=1:4   
     for j=1:4
         weighted(:,3,:);
   % tau(:,j,k)=xii(:,j,k).^weighted(:,j);
    
     tau(:,j,k)=sqrt(1-(1-xii(:,j,k).^2).^weighted(:,j));
   % upsilon(:,j,k)=sqrt(1-(1-oi(:,j,k).^2).^weighted(:,j));
    upsilon(:,j,k)=oi(:,j,k).^weighted(:,j);
    
     end
   % weighted
   tmp=[reshape(tau(:,:,k)',1,[]);  reshape(upsilon(:,:,k)',1,[])];
   eval(['Y',num2str(k),'=  sprintf(''(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)\n'', tmp)']);    
 end
% 方案理想解

%%% 第2步，求理想解
tau_pos=max(tau,[],3);
    tau_neg=min(tau,[],3);

   
    upsilon_pos=min(upsilon,[],3);%%%正理想解   
    upsilon_neg=max(upsilon,[],3); %%%负理想解

pos=[];
neg=[];
com=[];
pos=[pos, tau_pos(:,1),upsilon_pos(:,1), tau_pos(:,2), upsilon_pos(:,2), tau_pos(:,3), upsilon_pos(:,3), tau_pos(:,4), upsilon_pos(:,4)];
neg=[neg, tau_neg(:,1),upsilon_neg(:,1), tau_neg(:,2), upsilon_neg(:,2), tau_neg(:,3), upsilon_neg(:,3), tau_neg(:,4), upsilon_neg(:,4)];
com=[com, upsilon_pos(:,1),tau_pos(:,1), upsilon_pos(:,2), tau_pos(:,2), upsilon_pos(:,3), tau_pos(:,3), upsilon_pos(:,4), tau_pos(:,4)];

%%% 输出理想解 

disp('第12步，理想决策')

disp('maximum decision matrix Y_{+}')
sprintf('(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f, %4.2f)\n', pos')

disp('minimum decision matrix Y_{-}')
sprintf('(%4.2f, %4.2f)  (%4.2f,  %4.2f)  (%4.2f, %4.2f)  (%4.2f, %4.2f)\n', neg')

disp('余理想决策为：')
sprintf('(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f, %4.2f)\n', com')

for i=1:4  
    for j=1:4   
       chipos(i,j)=sqrt(1-tau_pos(i,j).^2-upsilon_pos(i,j).^2);
     %   chipos(i,j)=1-tau_pos(i,j)-upsilon_pos(i,j);
        % chi_{kj}^{+}
        chineg(i,j)=sqrt(1-tau_neg(i,j).^2-upsilon_neg(i,j).^2);
      % chineg(i,j)=1-tau_neg(i,j)-upsilon_neg(i,j);
        % chi_{kj}^{-}
    end
end
chipos
chineg

disp('第13步，遗憾矩阵为') 

%R=Y;
for k=1:4  
   rho(:,:,k)=tau_pos.*tau(:,:,k);
   varrho(:,:,k)=sqrt(upsilon_pos.^2+upsilon(:,:,k).^2-(upsilon_pos.^2).*(upsilon(:,:,k).^2));
   
   tmp=[reshape(rho(:,:,k)',1,[]);reshape(varrho(:,:,k)',1,[])];
   eval(['R',num2str(k),'=  sprintf(''(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)\n'', tmp)']);  
end
 %  https://www.ilovematlab.cn/thread-317135-1-1.html  

disp('第14步，最大遗憾矩阵Rmax为')  

rho_pos=max(rho,[],3); %%%正理想解
    

    varrho_pos=min(varrho,[],3); %%%正理想解 

Rmax=[];

Rmax=[Rmax, rho_pos(:,1),varrho_pos(:,1), rho_pos(:,2), varrho_pos(:,2), rho_pos(:,3), varrho_pos(:,3), rho_pos(:,4), varrho_pos(:,4)];

sprintf('(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f, %4.2f)\n', Rmax')

disp('第15步，群满意矩阵S_neg为')  

for k=1:4 
    psi(:,:,k)=tau(:,:,k).*tau_neg;  
    omega(:,:,k)=sqrt(upsilon_neg.^2+upsilon(:,:,k).^2-(upsilon_neg.^2).*(upsilon(:,:,k).^2));
  
    tmp=[reshape(psi(:,:,k)',1,[]);reshape(omega(:,:,k)',1,[])];
    eval(['S_neg',num2str(k),'= sprintf(''(%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f)\n'', tmp)']);  
end
disp('第16步，最大群满意矩阵S_neg_max为')  

psi_pos=max(psi,[],3); %%%正理想解
    
omega_pos=min(omega,[],3); %%%正理想解 

S_neg_max=[];

S_neg_max=[S_neg_max, psi_pos(:,1),omega_pos(:,1), psi_pos(:,2), omega_pos(:,2), psi_pos(:,3), omega_pos(:,3), psi_pos(:,4), omega_pos(:,4)];

sprintf('(%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f,  %4.2f)  (%4.2f, %4.2f)\n', S_neg_max')

% disp('第17步，群满意矩阵S_com为')  
% 
% for k=1:4    
%     iota{k}=tau{k}.*upsilon_pos;
% %     iota1(:,:,k)=iota{k};
% %     iota1
%     kappa{k}=sqrt(upsilon{k}.^2+tau_pos.^2-(tau_pos.^2).*(upsilon{k}.^2));
%     
%      tmp=[reshape(iota{k}',1,[]);reshape(kappa{k}',1,[])];
%      eval(['Scom',num2str(k),'= sprintf(''(%4.3f, %4.3f) (%4.3f, %4.3f) (%4.3f, %4.3f) (%4.3f, %4.3f)\n'', tmp)']);  
% end
% 
% disp('第18步，最大群满意矩阵S_com_max为')  
% 
% temp=[];
% for k=1:4  
%     tem(:,:,k)=iota{k};
%     iota_pos=max(tem,[],3); %%%正理想解
%     
%     tp(:,:,k)=kappa{k};
%     kappa_pos=min(tp,[],3); %%%正理想解       
% end
% 
% S_com_max=[];
% 
% S_com_max=[S_com_max, iota_pos(:,1),kappa_pos(:,1), iota_pos(:,2), kappa_pos(:,2), iota_pos(:,3), kappa_pos(:,3), iota_pos(:,4), kappa_pos(:,4)];
% 
% sprintf('(%4.3f, %4.3f) (%4.3f, %4.3f) (%4.3f, %4.3f) (%4.3f, %4.3f)\n', S_com_max')

disp('第19步，计算group utility')  

for k=1:4 
    % Yi*Y+
    chi(:,:,k)=sqrt(1-tau(:,:,k).^2-upsilon(:,:,k).^2);
    chi_pos=sqrt(1-tau_pos.^2-upsilon_pos.^2);
    chi_neg=sqrt(1-tau_neg.^2-upsilon_neg.^2);
    Y_Ypos(k)=sum(sum(tau(:,:,k).*tau_pos+upsilon(:,:,k).*upsilon_pos ...
     +chi(:,:,k).*chi_pos));
   % Yi*Y-
    Y_Yneg(k)=sum(sum(tau(:,:,k).*tau_neg+upsilon(:,:,k).*upsilon_neg ...
     +chi(:,:,k).*chi_neg));
end
disp('Yi与Ypos内积 ')  

Y_Ypos

disp('Yi与Yneg内积 ') 

Y_Yneg


for k=1:4
 GU(k)=min(4*n,Y_Ypos(k))/max(4*n,Y_Ypos(k)); 
end
GU

 
NGU=(GU-min(GU))./(max(GU)-min(GU))

disp('第20步，计算group regret ')  

for k=1:4  
    % Ri*Rmax
    sigma(:,:,k)=sqrt(1-rho(:,:,k).^2-varrho(:,:,k).^2);
    sigma_pos=sqrt(1-rho_pos.^2-varrho_pos.^2);
    R_Rmax(k)=sum(sum(rho(:,:,k).*rho_pos+varrho(:,:,k).*varrho_pos+sigma(:,:,k).*sigma_pos));
end

disp('Ri与Rmax内积 ') 

R_Rmax

for k=1:4
 GR(k)=min(4*n,R_Rmax(k))/max(4*n,R_Rmax(k)); 
end
GR 
 
NGR=(max(GR)-GR)./(max(GR)-min(GR))


disp('第21步，计算group satisfaction GS- ')  

for k=1:4  
    % S-*Smax
    varsigma(:,:,k)=sqrt(1-psi(:,:,k).^2-omega(:,:,k).^2);
    varsigma_pos=sqrt(1-psi_pos.^2-omega_pos.^2);
    Sneg_Smax(k)=sum(sum(psi(:,:,k).*psi_pos+omega(:,:,k).*omega_pos ...
     +varsigma(:,:,k).*varsigma_pos));
end

disp('Si^{-}与Smax^{-}内积 ')
Sneg_Smax

disp('GS_{i}^{-} ')

for k=1:4
 GSneg(k)=min(4*n,Sneg_Smax(k))/max(4*n,Sneg_Smax(k));
end
GSneg 
 
NGSneg=(GSneg-min(GSneg))./(max(GSneg)-min(GSneg))

% disp('第22步，计算group satisfaction GSc ')
% 
% for k=1:4
%     % Sc*Smax
%     varpi(:,:,k)=sqrt(1-iota{k}.^2-kappa{k}.^2);
%     varpi_pos=sqrt(1-iota_pos.^2-kappa_pos.^2);
%     %  s_{kj}^{ci}=(\iota_{kj}^{i},\kappa_{kj}^{i})
%     %  s_{kj}^{cmax}=(\iota_{kj}^{max},\kappa_{kj}^{min})
%     Scom_Smax(k)=sum(sum(iota{k}.*iota_pos+kappa{k}.*kappa_pos ...
%         +varpi(:,:,k).*varpi_pos));
% end
% disp('Sc*Smax ')
% Scom_Smax
% 
% for k=1:4
%     GScom(k)=4*n/(4*n+abs(4*n-Scom_Smax(k)));
% end
% GScom
% 
% NGScom=(GScom-min(GScom))./(max(GScom)-min(GScom))
% 
% disp('第23步，计算 comprehensive group satisfaction CGS ')
% 
% for k=1:4
%     CGS(k)=0.5*NGSneg(k)+0.5*NGScom(k);
% end

%CGS

disp('第24步，计算 comprehensive VIKOR measure Q ') 

% lambda=0:0.01:1/3;
%len=length(lambda);
%Q1=zeros(len,4);

%for row =1:len


% Lambda12 
% Q=lambda(row)*NGU/3+(2/3-lambda(row))*NGR+NGSneg/3; 
% Lambda13
% Q=lambda(row)*NGU+NGR/3+(2/3-lambda(row))*NGSneg; 
% Lambda1
 Q=NGU/3+NGR/3+NGSneg/3; 
% Q=NGU;
% Lambda2
% Q=NGU/3+lambda(row)*NGR+NGSneg/3; 
% Lambda3
% Q=NGU/3+NGR/3+lambda(row)*NGSneg; 
%  Lambda23 
%
% Q=NGU/3+lambda(row)*NGR+(2/3-lambda(row))*NGSneg; 

Q

% Q1(row,:)=Q; 
%  
% end
% 
%   
% plot(alpha,Q1(:,1),'b'),hold on,
% plot(alpha,Q1(:,2),'r.'),plot(alpha,Q1(:,3),'g--'),
% plot(alpha,Q1(:,4),'k-.');
% legend('A_{1}','A_{2}','A_{3}','A_{4}','Location','Southeast');
% % axis([0,1,-0.1,0.7]);
% %axis([0,0.35,1.2]);
% %set(gca,'xtick',[0:0.2:1]);
% % title('Qi');
% xlabel('\mu_{11}^{l}');
% ylabel('Four software products');
