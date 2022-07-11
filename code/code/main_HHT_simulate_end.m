clc;clear;close all
x=-4*pi:pi/100:pi;
y=2*exp(-((x-2)/2+4).^2).*abs(cos(x)-1.2);
plot(x,y);
f = y+ 0.05*randn(size(y));
figure;plot(f);
for i = 1
    x               = f(i,:);
    t               = linspace(0,1,size(x,2));%������ƽ��
    [imf,ort,nbits] = emd(x); %�ֽ�ֹͣ����   
    [m_imf,n_imf]   = size(imf);       
    for j           = 1:m_imf
       %------------------ÿ��IMF˲ʱƵ�ʵļ���----------------------------
        xn          = hilbert(imf(j,:)); %Discrete-time analytic signal via Hilbert transform
        xr          = real(xn); %������ʵ�� Complex real part
        xi          = imag(xn); %�������鲿 Complex imaginary part
        A           = sqrt(xr.^2+xi.^2); %A1ΪIMF1��˲ʱ���
        P           = atan2(xi,xr); % PΪIMF1��˲ʱλ��
        xh          = unwrap(P);%չ����λ��
        fs          = 1000;
        xhd(j,:)    = fs*diff(xh)/(2*pi);%xhdΪ˲ʱƵ��                   
    end    
end
 xhd_mean = mean(xhd,2);
Z=sum(imf(5:end,:));
std=std(xhd');
var=var(xhd');
V1=var(:,7)/var(:,8);
V2=var(:,6)/var(:,7);
V3=var(:,5)/var(:,6);
V4=var(:,4)/var(:,5);
V5=var(:,3)/var(:,4);
V6=var(:,2)/var(:,3);
V7=var(:,1)/var(:,2);
V=[V1 V2 V3 V4 V5 V6 V7];
%-----------------------f����--------------------------------
alpha = 0.05;
p = 1-2*alpha;
F_table =  finv(p,2,4);
factor = find(V > F_table);
factor = factor(1);
Fn = factor+1;
%------------------SGȥ��----------------------------
ynames = {'test'}
 X = f(:,1:end); 
i = 1;
column  = 1;
columnname   = ynames{column};
disp(columnname);
switch column
    case 1
         maxrank = 2;
         window0 = 5; window1 = 7; window2 = 13; window3 = 41;
end
 disp('SG-smoothing');
     SG                     = sgdiff(X,maxrank,window2,0);
for i = 1
 w               = X(i,9:497); 
end
for u = 1
 q               = SG(u,1:end);
end
%-----------------------HHT��SNR����--------------------------------
I_HHT  = Z;% I :original signal
In_HHT = f;% In:noisy signal(ie. original signal + noise signal)
Ps_HHT=sum((I_HHT).^2);%signal power
Pn_HHT=sum((I_HHT-In_HHT).^2); %noise power
snr_HHT=10*log10(Ps_HHT/Pn_HHT); % ����I�Ǵ��źţ�In�Ǵ����źţ�snr�������
%-----------------------SG��SNR����--------------------------------
I_SG = w;
In_SG= q;
snr=0;
Ps_SG=sum((I_SG).^2);%signal power
Pn_SG=sum((I_SG-In_SG).^2); %noise power
snr_SG=10*log10(Ps_SG/Pn_SG); % ����I�Ǵ��źţ�In�Ǵ����źţ�snr�������
%-----------------------ͼ--------------------------------
figure;%IMF����
for j    = 1:m_imf
    subplot(m_imf,1,j);plot(imf(j,:));
end
figure;%˲ʱƵ�ʻ���
for j    = 1:m_imf
    subplot(m_imf,1,j);plot(t(1:length(t)-1),xhd(j,:));
end
figure;plot(xhd_mean,'-o');%ƽ��˲ʱƵ�ʻ���
figure;plot(x);hold on;plot(sum(imf(5:end,:)),'r-');%HHTȥ��������
figure,plot(X(i,9:497));hold on;plot(SG(u,1:end),'g-');%SGȥ��������
figure,errorbar(xhd_mean,std);

