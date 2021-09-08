clc;
clear;
alpha=2.2*10^-3;
epi=0.;
rho=1;
beta=0.012;
gamma=0.6;
E=[40.0,48.2,59.6,70.8,76.2,78.6,81.4];
spread=[.10,.10,.10,.08,.04,.02,.02];
n=length(E);
dose=zeros(1,81);
Emax=E(n);
% data=zeros(n,100);
for i=1:n
    p=4.141*exp(-1.56*E(i))-2.25*10^-7*E(i)^2+1.773;
    r0=round(alpha*E(i)^p,1);
    sigma1=0.012*r0^0.935;
    sigma=sqrt(sigma1^2+(spread(i)*10)^2*alpha^2*p^2*E(i)^(2*p-2));
    depth1=round(r0-sigma,1);
    
    if depth1>0%proximal part
        x1=0:0.1:depth1;
        y1=(17.93*(r0-x1).^(-0.435)+(0.444+31.7*epi/r0)*(r0-x1).^0.565)/(1+0.012*r0);
%         y1=y1./max(y1);
    else%distal part
        depth1=-0.1;
    end
    %equa 26
    x2=depth1+0.1:0.1:r0-sigma;
    y2=zeros(1,length(x2));
    for j=1:1:length(x2)
        %PCF
        zeta=(r0-x2(j))/sigma;

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

        y2(j)=exp(-zeta^2/4)*(sigma^(1/p)*g1*(1/sigma)*D1+(beta/p+gamma*beta+epi/r0)*D2)/(sqrt(2*pi)*rho*p*alpha^(1/p)*(1+beta*r0));
        if y2(j)<0
            y2(j)=0;
        end
    end
%     y=y2./max(y2);
%     plot(x2,y,'-');
%     hold on;
    %equa 29
    x3=r0-sigma+0.1:0.1:8;
    y3=zeros(1,length(x3));
    for k=1:length(x3)
        z=-(r0-x3(k))./sigma;
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
        y2(k)=(exp(-(r0-x3(k))^2/(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*d1+(0.157+11/26*epi/r0)*d2))/(1+0.012*r0);
    % y2=(exp(-(r0-z2).^2./(4*sigma^2))*sigma^0.565*(11.26*sigma^-1*
    end
    if depth1<0
        y=[y2 y3];
    else
        y=[y1 y2 y3];
    end
    %weight
    w=rho*p*sin(pi/p)*alpha^(1/p)/(pi*alpha*(Emax^p-E(i)^p)^(1/p))*spread(i)*E(i);
    dose=dose+y./max(y)*w;
end
x=0:0.1:81;
dose=dose./max(dose);
plot(x,dose,'-');
hold on;
% %MC sim
% m=csvread('data1/10%data/40.0MeV10%dose.csv',8,2);
% x=m(:,1);
% y=m(:,2);
% x=flipud(x).*0.1;
% y=y./max(y);
% plot(x,y,'--');
% xlabel('depth(cm)');
% ylabel('D/Dmax');