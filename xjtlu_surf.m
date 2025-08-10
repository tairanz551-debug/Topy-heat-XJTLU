%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function xjtlu_surf(nelx,nely,volfrac,penal,rmin);
% INITIALIZE
% nelx=500;      %x�᷽��ĵ�Ԫ��
% nely=80;      %y�᷽��ĵ�Ԫ��
% volfrac=0.4;     %�����
% penal=5.0;         %���ϲ�ֵģ�͵ĳͷ�����
% rmin=1.2;       %���ȹ��˵İ뾶
x(1:nely,1:nelx) = volfrac; 
%��ʼ��
loop = 0;    %��ŵ��������ı���
change = 1.;     %����������ε�������Ʊ����Ĳ�ֵ�����ֵ�������жϵ�����ʱ����
% ��ʼ����
%while change > 0.01    %�����ε���������������Ʊ���x�Ĳ�ֵ�����ֵС��0.01ʱ����������
while (change > 0.01)&(loop<=100)
  loop = loop + 1;      %ÿѭ��һ�Σ�����������һ
  xold = x;           %��ǰһ�ε���Ʊ�����ֵ��xold�����Ժ��������µ�change
% FE-ANALYSIS
  [U]=FE(nelx,nely,x,penal);       %���нṹ������Ԫ���������ؽṹ��ȫ��λ�ƾ���  
% Ŀ�꺯���������ȷ���
  [KE] = lk;      %���ü��㵥Ԫ�նȾ�����ӳ��򣬷��ص�Ԫ�նȾ���
  c = 0.;        %�������Ŀ�꺯���ı���������Ŀ�꺯���ǽṹ�ܸն���󣬼������С
  for ely = 1:nely       %y�������
    for elx = 1:nelx       %x�����ÿ����Ԫ������xy����������Ƕ����е�Ԫ����ѭ��ͬһ�ֲ���
      n1 = (nely+1)*(elx-1)+ely;    %��Ԫ���Ͻڵ���
      n2 = (nely+1)* elx   +ely;     %��Ԫ���Ͻڵ���
     % Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);     %��ȫ��λ�ƾ����У������ɶȱ����ȡ��Ԫλ�ƾ���
     % c = c + x(ely,elx)^penal*Ue'*KE*Ue;               %����Ŀ�꺯����ֵ(�������)
     % dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;  %�����������ÿ����Ԫ��Ӧ��������ֵ����Ŀ�꺯������Ʊ����ĵ�����ָ����Ʊ����������µķ���
    Ue = U([n1; n2; n2+1; n1+1],1);     %��ȫ��λ�ƾ����У������ɶȱ����ȡ��Ԫλ�ƾ���
    c = c + (0.001+0.999*x(ely,elx)^penal)*Ue'*KE*Ue;  %����Ŀ�꺯���ĵ�ֵ������Сɢ������(���ɢ��ǿ��)
    dc(ely,elx) = -0.999*penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;    %�����������ÿ����Ԫ��Ӧ��������ֵ����Ŀ�꺯������Ʊ����ĵ�����ָ����Ʊ����������µķ���
    end
  end
% �����ȹ���
 [dc]   = check(nelx,nely,rmin,x,dc);    %���������ȹ����ӳ��򣬶������Ƚ��й��ˣ�������̸�ʽ��ͬʱʹ�߽�⻬
% �����Ż�׼����±���
  [x]    = OC(nelx,nely,x,volfrac,dc);    %�����Ż�׼�򣬸�����Ʊ���
% PRINT RESULTS ��ӡ���
  change = max(max(abs(x-xold)));    %������Ʊ����������ε����Ĳ�ֵ�ľ���ֵ�����ֵ
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
 ' ch.: ' sprintf('%6.3f',change )])      %����Ļ����ʾ����������Ŀ�꺯��ֵ������ʡ���Ʊ���������ֵ�ȵ�����Ϣ
% PLOT DENSITIES  
%  colormap(gray); 
%  colormap(parula); 
%  colormap(jet);

	clear data;
	data=importdata('./purple_color_map.txt');
	purple=data(:,1:3);
	
  colormap(jet);  
  
  h=imagesc(-x); 
%  axis equal; 
%	axis([0 14 0 4])
%  axis tight; 
%	axis normal
%	axis auto
%  set(gcf,'Position',[100 100 260 220])
	set(gcf, 'unit', 'centimeters', 'position', [10 5 28 8]);
  axis off;pause(1e-6);   %�Ż������ͼ����ʾ��������ǰ���ܶȳ��ԻҶ�ͼ����ʽ���������
  saveas(gca,strcat('./good_1','.png'))
%  savefig(h,'good_1.png')
end 
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE�Ż�׼����� %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%OC�㷨�ӳ���
function [xnew]=OC(nelx,nely,x,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;  %%11��12�������Լ�����������ճ��ӣ�moveΪ�����ƶ�����
while (l2-l1 > 1e-4)    
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));   %�˴�ΪOC�㷨�ĺ��ģ��ǹ�ʽ��2���ľ��巴ӳ
  if sum(sum(xnew)) - volfrac*nelx*nely > 0;    %���ö��ַ������������ճ���
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%�����ȹ����ӳ���
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
       fac = rmin-sqrt((i-k)^2+(j-l)^2);      %������˰뾶rmin�ڵ�ĳ����Ԫx(1,k)�ľ������
        sum = sum+max(0,fac);       %������˰뾶rmin�����е�Ԫ�ľ������֮��
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);   %����
      end
    end
   dcn(j,i) = dcn(j,i)/(x(j,i)*sum);   %���㵥Ԫx(j,i)���˺��������
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����Ԫ����ӳ���
function [U]=FE(nelx,nely,x,penal)
[KE] = lk;     %���ü��㵥Ԫ�նȾ����ӳ��򣬷��ص�Ԫ�նȾ���
%K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));   %����նȾ����ϡ�����
%F = sparse(2*(nely+1)*(nelx+1),1); U = zeros(2*(nely+1)*(nelx+1),1);   %�غ���F����λ�ƾ���U��ϡ�����
K = sparse((nelx+1)*(nely+1), (nelx+1)*(nely+1));
F = sparse((nely+1)*(nelx+1),1); U = zeros((nely+1)*(nelx+1),1);
for elx = 1:nelx
  for ely = 1:nely
    n1 = (nely+1)*(elx-1)+ely; 
    n2 = (nely+1)* elx   +ely;
    %edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]; %����Ԫ�İ˸����ɶȱ��
    %K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;  %������SIMP�ͷ��õ��ĵ�Ԫ�նȾ��󣨼���x(ely,elx)^penal*KE���������ɶ��ۼӣ��õ�����նȾ���
    edof = [n1; n2 ; n2+1; n1+1];
    K(edof,edof) = K(edof,edof) + (0.001+0.999*x(ely,elx)^penal)*KE;
  end
end
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F(:,1) = 0.01;   %��ʾ������Դ��"F(:,1) = 0.01;"��ʾ������ȼ��ȣ�q=0.01W/m^2*K

% F(:,1) = 0;
% F((nelx/2+1)*(nely+1)-(nely/2+1),1) = 0.01;  
% F(((nelx+1)*(nely+1)+1)/2,1) = 0.01; 

% ----- setup the graph -------

% plot a logo for heat sink

% nx = 4;
% ny = 4;
n = (nelx+1)*(nely+1);

a = 1:n;
a = a';

% nx=ny=4
% a = 1,2,3,...,25

b=reshape(a,nelx+1,nely+1);

% b =
%      1     6    11    16    21
%      2     7    12    17    22
%      3     8    13    18    23
%      4     9    14    19    24
%      5    10    15    20    25

%ax = linspace(-2,2,nelx+1);
ax = linspace(0,25,nelx+1);
% ax = linspace(0,2*pi,nelx+1);
% ax = linspace(0,2,nelx+1);
% ax = 0    0.2500    0.5000    0.7500    1.0000
% ay = linspace(0,2,nely+1);
% ay = linspace(2*pi,0,nely+1);
%ay = linspace(2,-2,nely+1);
ay = linspace(4,0,nely+1);

[X, Y] = meshgrid(ax,ay);

ox= reshape(X,n,1);
oy= reshape(Y,n,1);
o = [ox,oy];

% o =
% 
%          0    1.0000 (x,y)-> node number
%          0    0.7500 
%          0    0.5000
%          0    0.2500
%          0         0
%     0.5000    1.0000
%     0.5000    0.7500
%     0.5000    0.5000
%     0.5000    0.2500
%     0.5000         0
%     1.0000    1.0000
%     1.0000    0.7500
%     1.0000    0.5000
%     1.0000    0.2500
%     1.0000         0
%     1.5000    1.0000
%     1.5000    0.7500
%     1.5000    0.5000
%     1.5000    0.2500
%     1.5000         0
%     2.0000    1.0000
%     2.0000    0.7500
%     2.0000    0.5000
%     2.0000    0.2500
%     2.0000         0

% F=a';
% F=a;
c =[]; 

% x: o(i,1)
% y: o(i,2)
	

for i=1:n
%    if (o(i,1)<=0.75) & (o(i,1)>=0.25) & (o(i,2)>=0.25) & (o(i,2)<=0.75)
%    if (o(i,1)<=1.5) & (o(i,1)>=0.5) & (o(i,2)>=0.25) & (o(i,2)<=0.75)
%		if (o(i,1)<=0.56) & (o(i,1)>=0.44) & (o(i,2)>=0.44) & (o(i,2)<=0.56)
%		if (o(i,1)<=0.6) & (o(i,1)>=0.4) & (o(i,2)>=0.2) & (o(i,2)<=0.8)
%		if ((o(i,1)-0.5)^2+(o(i,2)-0.5)^2<=0.25^2) 
%		if (((o(i,1)-pi)^2+(o(i,2)-pi)^2-1)^3-(o(i,1)-pi)^2*(o(i,2)-pi)^3<=0) 
%		if ((o(i,1)^2+o(i,2)^2-1)^2<=o(i,1)^3)
%		if (o(i,1)^4+o(i,2)^4<=1)
%		if ((o(i,1)^2+o(i,2)^2)^3<=16*o(i,1)^2*o(i,2)^2)
%		if (o(i,2)^4-o(i,2)^2+(o(i,1)-0.5)^2<=0)
		
		
		% plot J
		dx =0.35;
		dc = 3;
		dh = 1;
		df=-0.2;
		dg=1.2;
		
		xx = o(i,1);
		yy = o(i,2);
%		if ((-1<=xx)&(xx<=1) &(1.5-dx<=yy) &(yy<=1.5+dx))| ((-dx<=xx)&(xx<=dx) &(-1.5+dx<=yy) &(yy<=1.5-dx)) | ((-0.7<=xx)&(xx<=dx) &(-1.5-dx<=yy) &(yy<=-1.5+dx)) % J
%		if ((-1<=xx)&(xx<=1) &(1.5-dx<=yy) &(yy<=1.5+dx))| ((-dx<=xx)&(xx<=dx) &(-1.5+dx<=yy) &(yy<=1.5-dx))  % T
%		if ((-dx<=xx)&(xx<=dx) &(-1.5<=yy) &(yy<=1.5)) | ((dx<=xx)&(xx<=1.2) &(-1.5<=yy) &(yy<=-1.5+2*dx))   % L
%		if ((-1<=xx)&(xx<=-1+2*dx) &(-1.5+2*dx<=yy) &(yy<=1.5)) | ((1-2*dx<=xx)&(xx<=1) &(-1.5+2*dx<=yy) &(yy<=1.5)) | ((-1<=xx)&(xx<=1) &(-1.5<=yy) &(yy<=-1.5+2*dx)) % U
%		if ((-xx-dx<=yy) & (yy<=-xx+dx) &(-1<=yy) &(yy<=1)) | ((xx-dx<=yy) & (yy<=xx+dx) &(-1<=yy) &(yy<=1)) % X
		
		if ((dc<=xx) & (xx<=dc+2*dx) &(dh<=yy) &(yy<=4-dh)) | (((xx-(2*dc+dx))^2+(1.2*yy-2)^2-1)^3-(xx-(2*dc+dx))^2*(1.2*yy-2)^3<=0) ...
		| ((-(xx-(3*dc+dx))-dx<=yy-2) & (yy-2<=-(xx-(3*dc+dx))+dx) &(-1<=yy-2) &(yy-2<=1)) | (((xx-(3*dc+dx))-dx<=yy-2) & (yy-2<=(xx-(3*dc+dx))+dx) &(-1<=yy-2) &(yy-2<=1))...
		| ((-1+df<=xx-(4*dc+dx))&(xx-(4*dc+dx)<=1-df) &(4-dh-2*dx<=yy) &(yy<=4-dh))| ((-dx<=xx-(4*dc+dx))&(xx-(4*dc+dx)<=dx) &(dh<=yy) &(yy<=4-dh-2*dx)) | ((-1<=xx-(4*dc+dx))&(xx-(4*dc+dx)<=dx) &(dh<=yy) &(yy<=dh+2*dx))...
		| ((-1+df<=xx-(5*dc+dx))&(xx-(5*dc+dx)<=1-df) &(4-dh-2*dx<=yy) &(yy<=4-dh))| ((-dx<=xx-(5*dc+dx))&(xx-(5*dc+dx)<=dx) &(dh<=yy) &(yy<=4-dh-2*dx))...
		| ((-dx<=xx-(6*dc+dx))&(xx-(6*dc+dx)<=dx) &(dh<=yy) &(yy<=4-dh)) | ((dx<=xx-(6*dc+dx))&(xx-(6*dc+dx)<=1.4) &(dh<=yy) &(yy<=dh+2*dx))...
		| ((7*dc+dx-dg<=xx)&(xx<=7*dc+dx-dg+2*dx) &(dh+2*dx<=yy) &(yy<=4-dh)) | ((7*dc+dx+dg-2*dx<=xx)&(xx<=7*dc+dx+dg) &(dh+2*dx<=yy) &(yy<=4-dh)) | ((7*dc+dx-dg<=xx)&(xx<=7*dc+dx+dg) &(dh<=yy) &(yy<=dh+2*dx)) 
        c = [c,a(i,1)]; % size (1,6)
    end
    
end

%for i=1:length(c)
%    F(c(i),1)=0.01;
%end




%------------------------------

% fixeddofs   = [nely/2+1-(nely/20):2:nely/2+1+(nely/20)];   %ʩ���¶ȱ߽��λ��Լ��,��ʾ��λ���¶�T=0��
	
% fixeddofs = [1:nely+1,nely+1:nely+1:nelx*nely+nelx+nely+1,nelx*nely+nelx+1:nelx*nely+nelx+nely+1,1:nely+1:nelx*nely+nelx+1];
	
%fixeddofs = [1,nely+1,nelx*(nely+1)+1,(nelx+1)*(nely+1)];

%size(fixeddofs)

% size(c) % (1,216)

% fixeddofs = [((nelx+1)*(nely+1)+1)/2];

fixeddofs = c;
	
alldofs     = [1:(nely+1)*(nelx+1)];    %����������нڵ����ɶ�
freedofs    = setdiff(alldofs,fixeddofs);    %��������������ɽڵ�����ɶ�
% SOLVING
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);       %����Ԫ��⣬�õ���������������ɽڵ��λ�� 
U(fixeddofs,:)= 0;         %�����й̶��Ľڵ�λ�Ƹ�ֵΪ0
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���㵥Ԫ�նȾ����ӳ���
function [KE]=lk
%E = 1.; 
%nu = 0.3;
%k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ... 
%   -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
%KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
 %                 k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
 %                 k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
 %                 k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
 %                 k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
 %                 k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
 %                 k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
  %                k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
  %�����άƽ���Ľڵ��ȴ����նȾ����ӳ���
   KE = [ 2/3 -1/6 -1/3 -1/6
          -1/6 2/3 -1/6 -1/3
          -1/3 -1/6 2/3 -1/6
          -1/6 -1/3 -1/6 2/3];



%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


