%% PPM plotter
clc
clear all
files = dir('**/*.mat');
c=0;
h = waitbar(0,'Hello world');
janine = [];
Fnew = [];
for i = 1:length(files)
    waitbar(i/length(files),h);
    fid = open(files(i).name);
    if isfield(fid,'F_osc_acq')
        c = c+1;
        F_osc(c) = fid.F_osc_acq;
        dat(c).raw(:,1) = fid.f_dig;
        dat(c).raw(:,2) = fid.s_dig;
        dat(c).Fs = fid.fs;
        dat(c).OS = round(fid.fs/fid.F_osc_acq);
        dat(c).ncyc = round(length(fid.f_dig)/dat(c).OS);
        %phase(c,:) = fid.phase;
        [Phase(c,:),E(c)]=Qc(1/dat(c).Fs,dat(c).raw(:,1),dat(c).raw(:,2),F_osc(c),dat(c).ncyc,1);

    end
end
close(h)
yyaxis left
scatter(F_osc,Phase(:,1))
ylim([-0.1 0.1])
set(gca,'xscal','log')
yyaxis right
scatter(F_osc,E)
ylim([30 100])


function [Phase,Eout]=Qc(dt_us,f0,s0,f_Hz,ncyc,nstk)
Phase=[nan,nan,nan];
Sng=[nan nan nan nan];
status=0;
%% get info
ldat = length(f0);
dt_s = (ones(ldat,1).*dt_us);
data.R0(:,1)    =   [cumsum(dt_s)];
data.R0(:,2)    =   f0;
data.R0(:,3)    =   s0;
dt              =   mean(dt_s);                % calculate dt for acquired data
fss             =   1/dt;                       % calculate sampling rate 
scyc            =   round(fss/f_Hz);            % calculate oversampling
nsmpR           =   length(data.R0(:,1));       % number of samples (raw)
nsmpE           =   ncyc*nstk*scyc;             % numebr of samples (elaborated)
nsstk           =   scyc*ncyc;                  % samples in a stack
Status          =   [dt fss scyc nsmpR nsmpE nsstk];
%% interp on fixed dt
data.RT(:,1)    =   linspace(0,dt*nsmpR,nsmpR);
data.R1         =   interp1(data.R0(:,1),data.R0(:,[2 3]),data.RT,'linear','extrap');
%% detrend data
it              =   [round(scyc/2):scyc:nsmpR-round(scyc/2)];
data.RDTr(:,1)  =   data.RT(it);
sm              =   [fix(it-round(scyc/2)+1)',fix(it+round(scyc/2)-1)'];
for n=1:length(sm)
    data.RDTr(n,[2 3]) = mean(data.R1(sm(n,1):sm(n,2),:)); 
end
RDTrV = interp1(data.RDTr(:,1),data.RDTr(:,[2 3]),data.RT,'pchip','extrap');
data.R2=data.R1-RDTrV;
% %% delete pre and post
% synT    =   linspace(0,nsmpE*dt,nsmpE)'; synV=sin(2*pi()*f_Hz.*synT+1.5*pi());
% R3      =   data.R2-repmat(min(data.R2),nsmpR,1);
% R3      =   data.R2./repmat(max(data.R2),nsmpR,1);
% R3corr  =   xcorr(R3(:,1),synV); [m pf] =   max(R3corr);
% R3corr  =   xcorr(R3(:,2),synV); [m ps] =   max(R3corr);
Estr    =   1;%round((pf+ps)/2-nsmpR);
data.E0(:,1)        =   linspace(0,nsmpE*dt,nsmpE)';
data.E0(:,[2 3])    =   data.R1([Estr:Estr+nsmpE-1],:);
ForceM              =   mean([data.R0(:,2)]);
ShortM              =   mean([data.R0(:,3)]);
%% detrend data
it              = [round(scyc/2):scyc:nsmpE-round(scyc/2)];
data.EDTr(:,1)  = data.E0(it,1);
sm              = [fix(it-round(scyc/2)+1)',fix(it+round(scyc/2))'];
for n=1:length(sm)
    data.EDTr(n,[2 3]) = mean(data.E0(sm(n,1):sm(n,2),[2 3])); 
end
EDTrV = interp1(data.EDTr(:,1),data.EDTr(:,[2 3]),data.E0(:,1),'pchip','extrap');
data.E1=data.E0(:,[2 3])-EDTrV;
%% cut and overlap data
data.ET         =   data.E0(1:nsstk,1); 
tmp             =   reshape(data.E1(:,1),nsstk,nstk);
data.E2(:,1)    =   tmp;
tmp             =   reshape(data.E1(:,2),nsstk,nstk); 
data.E2(:,2)    =   tmp;
%% moving average
data.E3     =   data.E2;
% it          =   [round(nsav/2)+1:nsstk-round(nsav/2)];
% for n=1:length(it)
%     data.E3(n+round(nsav/2),:)=mean(data.E2(it(n)-round(nsav/2):it(n)+round(nsav/2)));
% end
%% calculate 
id      =   int32(nsstk*1/ncyc):int32(nsstk*(1-(1/ncyc)));
nsmpC   =   length(id);
Force   =   data.E3(id,1)-mean(data.E3(id,1));
Short   =   data.E3(id,2)-mean(data.E3(id,2));
%% fit the data
F1      =   sin(2*pi*f_Hz.*data.ET(id));
F2      =   cos(2*pi*f_Hz.*data.ET(id));
F       =   [F1 F2];
a       =   F\Force;
ForceA  =   norm(a);
ForceP  =   atan2(a(2),a(1));
a       =   F\Short;
ShortA  =   norm(a);
ShortP  =   atan2(a(2),a(1));
Phase(1)=   atan(tan(ForceP-ShortP));
%% AREA calculation
spent=0;
lost=0;
Force_a=(Force-min(Force))/(max(Force)-min(Force)).*2;
Short_a=(Short-min(Short))/(max(Short)-min(Short)).*2;
for c=1:nsmpC-1
    area=((Short_a(c+1)-Short_a(c))*(Force_a(c+1)+Force_a(c)))/2;
    lost=lost+area;
    if area>0
        spent=spent+abs(area);
    end
end
lost=lost/(ncyc-2);
spent=spent/(ncyc-2);
store=spent-lost;
Phase(2)=atan(lost/pi);
%% fft
NFFT = 2^nextpow2(nsmpC); % Next power of 2 from length of y
Yf = fft(Force_a-1,NFFT)/nsmpC;
Ys = fft(Short_a-1,NFFT)/nsmpC;
f = fss/2*linspace(0,1,NFFT/2+1);
pfftF=angle(Yf);
pfftS=angle(Ys);
pos=find(f>f_Hz);
Phase(3)=pfftF(pos(1))-pfftS(pos(1));
% dataOUT.ET=single(data.ET);
% dataOUT.E3=single(data.E3);
% dataOUT.R0=single(data.R0);
% dataOUT.RT=single(data.RT);
% Sng=[ForceM ForceA ShortM ShortA];
FA_kN = ForceA*4462.2;
SA = (0.0275^2-0.015^2)*pi();
L_Pa = FA_kN/SA;
U_m = ShortA*2.5e-6*1.494;
Ep = U_m/0.09;
Eout = 3*(L_Pa/Ep)/1e9;
end