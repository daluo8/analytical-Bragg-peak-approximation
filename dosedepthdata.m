clc;
clear;
%user input
r0=5.1;
E=81.4;
sigma=0.26;
%
p=4.141*exp(-1.56*E)-2.25*10^-7*E^2+1.773;
alpha=2.2*10^-3;
epi=0.0;
rho=1;
beta=0.012;
gamma=0.6;

m=csvread('nonuclear/2%data/81.4MeV2%dose.csv',8,2);
x=m(:,1);
y=m(:,2);
x=x.*0.1;
y=flipud(y)./max(y);
[z,index]=max(y);
depth=index*0.1
plot(x,y,'--');
hold on;
%proximal part
x1=0:0.1:4;
y1=(17.93*(r0-x1).^(-0.435)+(0.444+31.7*epi/r0)*(r0-x1).^0.565)/(1+0.012*r0);
% y1=y1./max(y1);
% plot(x1,y1,'-');
% hold on;
%distal part
%equa 26
x2=0:0.1:depth-0.1;
y2=zeros(1,length(x2));
for i=1:1:length(x2)
    %PCF
    zeta=(r0-x2(i))/sigma;

    %Gamma(1/p)
    fun1=@(t)t.^(1/p-1).*exp(-t);
    g1=integral(fun1,0,Inf);
    %Gamma(1+1/p)
    fun2=@(t)t.^(1/p).*exp(-t);
    g2=integral(fun2,0,Inf);
    
    %D1=U(-1/2+1/p,-zeta)
    a=-1/2+1/p;
    fun3=@(t)t.^(a-1/2).*exp(-1/2.*t.^2).*cos(-zeta.*t+(1/2*a+1/4)*pi);
    d1=integral(fun3,0,Inf);
    D1=sqrt(2/pi)*exp(1/4*zeta^2)*d1;
    
    %D2=U(1/2+1/p,-zeta)
    a=1/2+1/p;
    fun4=@(t)t.^(a-1/2).*exp(-1/2*t.^2+zeta.*t);
    d2=integral(fun4,0,Inf);
    D2=exp(-1/4*zeta.^2)*d2/g2;
    
    y2(i)=exp(-zeta^2/4)*(sigma^(1/p)*g1*(1/sigma)*D1+(beta/p+gamma*beta+epi/r0)*D2)/(sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*r0));
    if y2(i)<0
        y2(i)=0;
    end
end
y2=y2./max(y2);
% plot(x2,y2,'-');
% hold on;
%equa 29

x3=depth-0.1:0.1:10;
% x3=0:0.1:10;
y3=zeros(1,length(x3));
for i=1:length(x3)
    z=-(r0-x3(i))./sigma;
    %d1=U(0.065)
    a=-1/2+0.565;
    fun1=@(t)t.^(-a-1/2).*exp(-1/2.*t.^2).*cos(z.*t+(1/2*a+1/4)*pi);
    d1=sqrt(2/pi)*exp(1/4*z^2)*integral(fun1,0,Inf);
    %d2=U(1.065)
    a=-1/2+1.565;
    fun2=@(t)t.^(a-1/2).*exp(-1/2*t.^2-z.*t);
    %gamma
    fun=@(t)t.^(1/2+a-1).*exp(-t)
    g=integral(fun,0,Inf);
    d2=exp(-1/4*z^2)/g*integral(fun2,0,Inf);
    y3(i)=(exp(-(r0-x3(i))^2/(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*d1+(0.157+11/26*epi/r0)*d2))/(1+0.012*r0);
% y2=(exp(-(r0-z2).^2./(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*
end
y3=y3./max(y3);
% x3=x3(2:length(x3));
% y3=y3(2:length(y3));
x2=x2(1:length(x2)-1);
y2=y2(1:length(y2)-1);
x=[x2 x3]
y=[y2 y3];
y=y./max(y);
save('nofit/81.4.mat','y');


plot(x,y,'-');
% y3=y3./max(y3);
% plot(x3,y3,'-');
hold on;
%MC data
% m=csvread('withoutnuclear.csv',8,2);
% x=m(:,1);
% y=m(:,2);
% x=flipud(x).*0.1;
% y=y./max(y);
% plot(x,y,'.-','Color',[0 1 0]);

xlabel('depth(cm)');
ylabel('D/Dmax');
legend('MC sim','analytical fit');