clc;clear;close all
load Fe-CN-XRD
X = rad(:,1:end); 
i = 1;
%------------------HHTȥ��----------------------------
for i = 1
 x               = X(i,1:end); 
t               = linspace(0,1,size(x,2));
[imf,ort,nbits] = emd(x);
[m_imf,n_imf]   = size(imf); 
%------------------ÿ��IMF˲ʱƵ�ʵļ���----------------------------
 for j           = 1:m_imf
 xn          = hilbert(imf(j,:));
 xr          = real(xn); 
  xi          = imag(xn);
A           = sqrt(xr.^2+xi.^2);
 P           = atan2(xi,xr);
 xh          = unwrap(P);
 fs          = 1000;
        xhd(j,:)    = fs*diff(xh)/(2*pi);
 [A, fa, tt] = hhspectrum(imf(j,:));
 [E,tt1]     = toimage(A,fa,tt,length(tt));
 E           = flipud(E);              
    end    
end
xhd_mean = mean(xhd,2);
Z=sum(imf(6:end,:));
std=std(xhd');
%-----------------------f����--------------------------------
alpha = 0.05;
p = 1-2*alpha;
xhd_mean_mean=mean(xhd_mean,1);
for i=1:m_imf;
    S(i,:)=((xhd_mean(i,:)-xhd_mean_mean)^2)/(9-1);
end
 F1 = S(9,:)/S(8,:);
 F2 = S(8,:)/S(7,:);
 F3 = S(7,:)/S(6,:);
 F4 = S(6,:)/S(5,:);
 F5 = S(5,:)/S(4,:);
 F6 = S(4,:)/S(3,:);
 F7 = S(3,:)/S(2,:);
 F8 = S(2,:)/S(1,:);
F = [F1 F2 F3 F4 F5 F6 F7 F8];
F_table = finv(p,2,4);
factor = find(F > F_table);
factor = factor(1);
Fn = factor+1;
% %------------------SGȥ��----------------------------
ynames = {'Fe-CN-XRD'};
X = rad(:,1:end); 
i = 1;
column  = 1;
columnname   = ynames{column};
disp(columnname);
switch column
    case 1
         maxrank = 2;
         window0 = 5; window1 = 11; window2 = 13; window3 = 21;
end
 disp('SG-smoothing');
     SG                     = sgdiff(X,maxrank,window3,0);
for i = 1
 w               = X(i,11:2951); 
end
for u = 1
 q               = SG(u,1:end);
end
%-----------------------ͼ--------------------------------
figure;plot(X(:,1:3001)');%ԭʼ����
figure;%IMF����
for j    = 1:m_imf
    subplot(m_imf,1,j);plot(imf(j,:));
end
figure;%˲ʱƵ�ʻ���
for j    = 1:m_imf
    subplot(m_imf,1,j);plot(t(1:length(t)-1),xhd(j,:));
end
figure;plot(xhd_mean,'-o');%ƽ��˲ʱƵ�ʻ���
figure;plot(x);hold on;plot(sum(imf(6:end,:)),'r-');%HHTȥ��������
figure;plot(X(i,16:2956));hold on;plot(SG(u,1:end),'g-');%SGȥ��������
figure,errorbar(xhd_mean,std);