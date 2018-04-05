%Calculates Cd found reported in "An empirical model to estimate the propagation
%of random breaking and nonbreaking waves over vegetation fields" by Mendez
%and Losada/Coastal Engineering 51(2004)

%Cd is for the nonbreaking random wave transformation model found in
%section 2.3

%Solves for Cd uding a Gauss–Newton algorithm for solving nonlinear least squares
%with a line search using the Armijo rule

clc
clear all
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))

miny=0.0;
maxy=1.2;
minx=0;
maxx=10.0;
ytick=(miny:0.2:maxy);
xtick=(minx:2.0:maxx);

markersize=15;
ticksize=17;
labelsize=17;
textsize=17;

base_dir=cd;
cd(base_dir)

files=dir('*.mat');

for i=1:numel(files)
    file1 = getfield(files,{i},'name');
    data=load(file1);
    
    N=data.N; 
    if N==0
        continue
    end
    
    ls=data.ls; 
    b=data.bv; 
    h=data.h;
    k=data.kp;
    
    if ls/h > 1.0
        s=1.0;
    else
        s=ls/h;
    end
    
    y=data.Hrms_veg/data.Hrms_veg(1); 
    
    %Hrms=data.Hmo_veg/1.416;
    %y=Hrms/Hrms(1);
    
    syms Ho B Cd x
    f=1./(1+B*Cd*x);
    grad=gradient(f,Cd);
    
    Ho=data.Hrms_veg(1);
    %Ho=Hrms(1);
    x=data.glocs(5:13);  
    term1=(1/(3*sqrt(pi)))*b*N*Ho*k;
    term2=(sinh(k*s*h)^3+3*sinh(k*s*h))/((sinh(2*k*h)+2*k*h)*sinh(k*h));
    B=term1*term2;
    a=1;
    
    [fitresult, gof] = check_Cd_cftool(x*B, y); %used to check implemented
    %algorithm (verified to work)
    
    count=1;
    
    Cd=0.0;
    step=100;
    while step>1*10^(-5)
        F=eval(f)-y;
        gradFT=eval(grad);
        
        df=gradFT'*F;
        gradF2=gradFT'*gradFT;
        step=linsolve(gradF2,-df);
        
        Cd_orig=Cd;
        fc1=eval(f);
        fc2=0.0001*step*eval(grad);
        Cd=Cd_orig+step;
        newf=eval(f);
        while newf(2:end)>fc1(2:end)+a*fc2(2:end);
            a=0.8*a;
            Cd=Cd_orig+a*step;
            newf=eval(f);
            count=count+1;
            if count>100
                break
            end
        end
    end
    
    eq=1./(1+B*Cd.*x);
    
    SSerr=sum((y-eq).^2);
    SStot=sum((y-mean(y)).^2);
    
    R2=1-(SSerr/SStot);
    
    Cd_str=num2str(Cd,'%6.2f');
    R2_str=num2str(R2,'%6.3f');
    
    hfig=figure(1);
    plot(x,y,'k+',x,eq,'-k','LineWidth',1.75,'MarkerSize',markersize);
    axis([minx maxx miny maxy])
    set(gca,'YTick',ytick,'FontWeight','Bold','FontSize',ticksize)
    set(gca,'XTick',xtick,'FontWeight','Bold','FontSize',ticksize)
    xlabel('x (m)','FontWeight','Bold','FontSize',labelsize)
    ylabel('H_{rms}/H_{rms,i}','FontWeight','Bold','FontSize',labelsize)
    str={strcat('C_D=',Cd_str);strcat('R^2=',R2_str)};
    text(8,1.1,str,'FontSize',textsize)
    grid on
    
    figout=strcat('bulk_Cd\',data.file,'_Cd.tif');
    
    print(hfig,'-dtiff','-r600',figout)
    clf
    
    %check=fitresult.a; %Cd provided by cftool
    
    save(strcat(data.file,'.mat'),'Cd','-append')
    
    clear Cd R2 x y Ho N L hveg k Cd_str R2_str Cd_ML
end

%%For Dubi Dataset...
%N=1200
%ls=0.2
%bv=0.025
%hveg=0.4
%y=[0.096;0.084;0.076;0.075;0.074;0.070;0.055;0.046] %Dubi dataset (top left in Mendez paper) 
%x=[0;1.15;2.15;3.15;4.15;5.15;6.15;7.15]%Dubi dataset (top left in Mendez paper) 

% %Example problem provided in pdf
% y=[3.2939;4.2699;7.1749;9.3008;20.259];
% 
% syms x1 x2 t
% f=x1*exp(x2*t);
% grad=gradient(f,[x1,x2]);
% 
% x1=0.1;
% x2=0.1;
% t=[1;2;4;5;8];
% 
% step=100;
% while abs(step)>1*10^(-5)
%     F=eval(f)-y;
%     ans1=eval(grad(1));
%     ans2=eval(grad(2));
%     gradFT=[ans1 ans2];
% 
%     df=gradFT'*F;
%     gradF2=gradFT'*gradFT;
%     step=linsolve(gradF2,-df);
%     x1=x1+step(1)*0.1;
%     x2=x2+step(2)*0.1;
% 
%     diff=(1/2)*(F')*F;
% end
% 
% f=x1*exp(x2*t);
% 
% plot(t,y,'bo',t,f,'r-')



