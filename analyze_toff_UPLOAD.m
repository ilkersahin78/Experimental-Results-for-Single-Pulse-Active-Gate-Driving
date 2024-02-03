
Vdspeak=max(Vds_lower);
Idpeak=min(Id);

% I and V will not be precise. Should be zeroed for fair comparison
% temp2=find(Vds_lower>600);
% bsl2=temp2(1);
%Voffsetilk=mean(Vds_lower(bsl2:end));
%Voffsetilk=mean(Vds_lower(100000:end));
%Vosc=Vds_lower(bsl2:end)-Voffsetilk;
% Oscillations and probe data die after 20000 samples
% Vosc=Vds_lower(bsl2:bsl2+20000)-mean(Vds_lower(bsl2:bsl2+20000));
% temp2=find(Vds_lower>600);
% bsl2=temp2(1);
% Vosc=Vds_lower-Vfit;
% Vosc=Vosc(bsl2:bsl2+20000);


%Ioffsetilk=mean(Id(80000:end));
%Iosc=Id(bsl2:bsl2+20000)-mean(Id(bsl2:bsl2+20000));
% temp2=find(Id<0);
% bsl2=temp2(1);
% Iosc=Id-Ifit;
% Iosc=Iosc(bsl2:bsl2+20000);


% Rise and Fall times are definition dependent (95% to target value)
temp2=find(Id<0.95*max(Id));
Ifall1=t(temp2(1));
temp2=find(Id<0);
Ifall2=t(temp2(1));
Idfalltime=Ifall2-Ifall1;
%Iring=Id(temp2(1):temp2(1)+50000)-mean(Id(temp2(1):temp2(1)+50000));
temp2=find(Id<0);
Iringtemp=zeros(1,length(Id));
Iringtemp(temp2(1):end)=Id(temp2(1):end);
Ideneme=smoothdata(Iringtemp,'gaussian',3000);
Iring=Iringtemp-Ideneme;


temp2=find(Vds_lower>0.05*max(Vds_lower));
Vrise1=t(temp2(1));
temp2=find(Vds_lower>600);
Vrise2=t(temp2(1));
Vdsrisetime=Vrise2-Vrise1;
%Vring=Vds_lower(temp2(1):temp2(1)+50000)-mean(Vds_lower(temp2(1):temp2(1)+50000));
temp2=find(Vds_lower>600);
Vringtemp=zeros(1,length(Vds_lower));
Vringtemp(temp2(1):end)=Vds_lower(temp2(1):end)-600;
Vdeneme=smoothdata(Vringtemp,'gaussian',3000);
Vring=Vringtemp-Vdeneme;

%Find integrals
bsl=find(t==0);
integral_EMI=0;
integral_EMI2=0;
integral_Vosc=0;
integral_Iosc=0;


% readings should be zeroed for fair comparison
EMIright=EMI(bsl:end);
EMIoffset=mean(EMIright);
EMIrightzeroed=EMIright-EMIoffset;

for temp1=1:20000 %length(EMIrightzeroed)
    integral_EMI=integral_EMI+abs(EMIrightzeroed(temp1));
end
integral_EMI=t(bsl+20000)*integral_EMI;

EMI2right=EMI2(bsl:end);
EMI2offset=mean(EMI2right);
EMI2rightzeroed=EMI2right-EMI2offset;

for temp1=1:20000 %length(EMI2rightzeroed)
    integral_EMI2=integral_EMI2+abs(EMI2rightzeroed(temp1));
end
integral_EMI2=t(bsl+20000)*integral_EMI2;

Vringright=Vring(bsl:end);
Vringoffset=mean(Vringright);
Vringrightzeroed=Vringright-Vringoffset;

for temp1=1:20000 %length(EMIrightzeroed)
    integral_Vosc=integral_Vosc+abs(Vringrightzeroed(temp1));
end
integral_Vosc=t(bsl+20000)*integral_Vosc;

Iringright=Iring(bsl:end);
Iringoffset=mean(Iringright);
Iringrightzeroed=Iringright-Iringoffset;

for temp1=1:20000 %length(EMIrightzeroed)
    integral_Iosc=integral_Iosc+abs(Iringrightzeroed(temp1));
end
integral_Iosc=t(bsl+20000)*integral_Iosc;


% %%%% PLOT BASE RESULTS %%%%
% 
% plot(t,Vds_lower/100)
% hold on
% plot (t,Vgs/3)
% hold on
% plot (t,Id/100)
% grid on
% title ('Vds Id Vgs pu 100 and 3')
% legend ('Vds/100','Vgs/3','Id/100')
% xlim([-2e-6 4e-6])
% savefig('results_Vds_Id_Vgs')


%%%% ANTENNA THING %%%%

% [f, P_EMI]=CalcFFT(t, EMI);
% [f, P_EMI2]=CalcFFT(t, EMI2);
% figure
% semilogx(f,P_EMI)
% hold on
% semilogx(f,P_EMI2)
% xlim([1e6 1e9])
% grid on
% legend ('H Sensor','E Sensor')
% title ('H and E Field Meas.')
% savefig('results_antenna_fft')

%%%% OSCILLATIONS THING %%%%

%create new time series
% temp1=find(t==0);
% temp2=length(Iring);
% tcurrent=t(temp1:temp1+temp2);
% 
% temp2=length(Vring);
% tvoltage=t(temp1:temp1+temp2);
% 
% [fcurosc, Pcurrent]=CalcFFT(tcurrent, Iring);
% [fvoltosc, Pvoltage]=CalcFFT(tvoltage, Vring);
% figure
% semilogx(fcurosc,Pcurrent)
% hold on
% semilogx(fvoltosc,Pvoltage)
% xlim([1e6 1e9])
% grid on
% legend ('Current','Voltage')
% title ('I&V Oscillations')
% savefig('results_VIoscillations_fft')

% temp1=find(t==0);
% temp2=length(Iosc);
% tcurrent=t(temp1:temp1+temp2);
% 
% temp2=length(Vosc);
% tvoltage=t(temp1:temp1+temp2);
% 
% [fcurosc, Pcurrent]=CalcFFT(tcurrent, Iosc);
% [fvoltosc, Pvoltage]=CalcFFT(tvoltage, Vosc);
% figure
% semilogx(fcurosc,Pcurrent)
% hold on
% semilogx(fvoltosc,Pvoltage)
% xlim([1e6 1e9])
% grid on
% legend ('Current','Voltage')
% title ('I&V Oscillations')
% savefig('results_VIoscillations_fft')


%save('off200run25','Vds_lower','Vgs','Igate','Id','EMI','EMI2','t', 'Vosc','Iosc','Vdspeak', 'Idpeak', 'integral_EMI','integral_EMI2','integral_Vosc','integral_Iosc', 'Vdsrisetime', 'Idfalltime', 'Eswitch' )
    
%save('off_21septpos','EMI','EMI2','t')

%plot(t, Igate, t, fitresultIgate(t))

%List Key results:

% Vdspeak
% Idpeak
% Bsensor=1000*integral_EMI
% Esensor=1000*integral_EMI2
% integral_Vosc
% integral_Iosc
% Vdsris_ns=1e9*Vdsrisetime
% Idfall_ns=1e9*Idfalltime
% LossmJ=Eswitch*1000


function [DataReturn] = CorrectOffset(Data, ExpectedValue, FracStart, Start)

if Start == 1
    offset = mean(Data(1:ceil(length(Data)*FracStart)));
else
    offset = mean(Data(end-ceil(length(Data)*FracStart):end));
end    

DataReturn=Data+(ExpectedValue-offset);

end

function [f,P1] = CalcFFT(tdata,ydata)
%CALCFFT Summary of this function goes here
%   Detailed explanation goes here


X=ydata;
L=length(tdata);

if mod(L,2) == 1
    X=X(1:(end-1));
    L=L-1;
end

T=tdata(2)-tdata(1);
Fs=1/T;


Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;


end

