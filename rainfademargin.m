%RAINFADEMARGIN Fade margin introduced by rain for given availability
%
%
%fadedepth_dB=rainfademargin(f_Hz,rain_mmh,pathelevation_deg,
%       polarizationtilt_deg,length_m,availability_perc,mode);
%f_Hz:                  frequency (Hz)
%rain_mmh:              for the location where the link will work,
%                       rain rate (mm/hour) occurring with probab. 0.01%
%                       obtainable via maps or computations in
%                       https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.837-7-201706-I!!PDF-E.pdf
%pathelevation_deg:     elevation of path from horizontal (degrees)
%polarizationtilt_deg:  polarization degrees (0 for horizontal, 90 for vertical)
%length_m:              horizontal length of link (m)
%availability_perc:     percentage of availability required (99.995 typ.)
%mode:                  recommendation: 'ITU-R P.530-16'(default) or 'ITU-R P.530-8'
%
%fadedepth_dB: dB of additional loss introduced by rain.
%
%SEE ALSO rainfadedistribution, rainfademargin
%
%CHANGELOG
%   29/10/2019: updated help and references
%   01/12/2016: created by Matteo Oldoni based on spreadsheet by Goran
%       Biscevic, using model ITU-R P.530 and P.838 (2016)
%       https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.838-3-200503-I!!PDF-E.pdf
%       https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.530-17-201712-I!!PDF-E.pdf
%       https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.837-7-201706-I!!PDF-E.pdf
%
function fadedepth_dB=rainfademargin(f_Hz,rain_mmh,pathelevation_deg,polarizationtilt_deg,length_m,availability_perc,mode)
if nargin<1
    f_Hz=80e9;
    rain_mmh=40;
    pathelevation_deg=0;
    polarizationtilt_deg=90;
    length_m=4.8e3;
    availability_perc=99;
end
fadedepth_dB=0;
if nargin<7 || isempty(mode)
    mode='ITU-R P.530-16';
end
G7=f_Hz/1e9;
G8=rain_mmh;
kh=10.^(log10kh(f_Hz)); G9=kh;
kv=10.^(log10kv(f_Hz)); G10=kv;
G11=alphah(f_Hz);
G12=alphav(f_Hz);
G15=pathelevation_deg*pi/180;
G16=polarizationtilt_deg*pi/180;
k=(G9+G10+(G9-G10).*POWER(COS(G15),2).*COS(2*G16))/2; G17=k;
alpha=(G9.*G11+G10.*G12+(G9.*G11-G10.*G12).*POWER(COS(G15),2).*COS(2*G16))./(2*G17); G18=alpha;
specificatt_dBkm=G17.*POWER(G8,G18);  G19=specificatt_dBkm;

G21=abs(length_m)/1000;
G22=availability_perc;
probability=(100-G22); G23=probability;
if strcmp(mode,'ITU-R P.530-8')
    %From ITU-R P.530-8 (and paper by Owolawi et al)
    dizero=35*EXP(-0.015*G8); G24=dizero;
    effectivelength_m=G21/(1+G21/G24); G25=effectivelength_m;
    totalattenexceeded001perc_dB=G19*G25; G26=totalattenexceeded001perc_dB;
    fadedepth_dB=G26*(0.12*POWER(G23,-(0.546+0.043*LOG10(G23)))); G27=fadedepth_dB;
else
    %From ITU-R P.530-16
    r_distancefactor=1./(0.477*POWER(G21,0.633).*POWER(G8,0.073*G18).*POWER(G7,0.123)-10.579*(1-(EXP(-0.024*G21)))); G29=r_distancefactor;
    effectivepathlength_km=G21.*G29; G30=effectivepathlength_km;
    totalattenexceeded001perc_dB=G30.*G19; G31=totalattenexceeded001perc_dB;
    cizero=0.12+0.4*LOG10(POWER(G7/10,0.8)); G32=cizero;
    ciuno=POWER(0.07,G32).*POWER(0.12,1-G32); G33=ciuno;
    cidue=0.855*G32+0.546*(1-G32); G34=cidue;
    citre=0.139*G32+0.043*(1-G32); G35=citre;
    fadedepth_dB=G31.*G33.*POWER(G23,-(G34+G35.*LOG10(G23))); G36=fadedepth_dB;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dalla ITU-R P.838-3
function lkh=log10kh(f_Hz)
G7=f_Hz/1e9;
a1=-5.33980; R10=a1;
a2=-0.35351; R11=a2;
a3=-0.23789; R12=a3;
a4=-0.94158; R13=a4;
b1=-0.10008; R14=b1;
b2=1.26970;  R15=b2;
b3=0.86036;  R16=b3;
b4=0.64552;  R17=b4;
c1=1.13098;  R18=c1;
c2=0.45400;  R19=c2;
c3=0.15354;  R20=c3;
c4=0.16817;  R21=c4;
mk=-0.18961; R22=mk;
ck=0.71147;  R23=ck;

lkh=R10*EXP(-(POWER((LOG10(G7)-R14)/R18,2)))+R11*EXP(-(POWER((LOG10(G7)-R15)/R19,2)))+R12*EXP(-(POWER((LOG10(G7)-R16)/R20,2)))+R13*EXP(-(POWER((LOG10(G7)-R17)/R21,2)))+R22*LOG10(G7)+R23;

function ah=alphah(f_Hz)
G7=f_Hz/1e9; 
a1=-0.14318;  U10=a1;
a2=0.29591;   U11=a2;
a3=0.32177;   U12=a3;
a4=-5.37610;  U13=a4;
a5=16.17210;  U14=a5;
b1=1.82442;   U15=b1;
b2=0.77564;   U16=b2;
b3=0.63773;   U17=b3;
b4=-0.96230;  U18=b4;
b5=-3.29980;  U19=b5;
c1=-0.55187;  U20=c1;
c2=0.19822;   U21=c2;
c3=0.13164;   U22=c3;
c4=1.47828;   U23=c4;
c5=3.43990;   U24=c5;
malpha=0.67849; U25=malpha;
calpha=-1.95537;U26=calpha;
ah=U10*EXP(-(POWER((LOG10(G7)-U15)/U20,2)))+U11*EXP(-(POWER((LOG10(G7)-U16)/U21,2)))+U12*EXP(-(POWER((LOG10(G7)-U17)/U22,2)))+U13*EXP(-(POWER((LOG10(G7)-U18)/U23,2)))+U14*EXP(-(POWER((LOG10(G7)-U19)/U24,2)))+U25*LOG10(G7)+U26;

function lkv=log10kv(f_Hz)
G7=f_Hz/1e9;
a1=-3.80595; R29=a1;
a2=-3.44965; R30=a2;
a3=-0.39902; R31=a3;
a4=0.50167;  R32=a4;
b1=0.56934;  R33=b1;
b2=-0.22911; R34=b2;
b3=0.73042;  R35=b3;
b4=1.07319;  R36=b4;
c1=0.81061;  R37=c1;
c2=0.51059;  R38=c2;
c3=0.11899;  R39=c3;
c4=0.27195;  R40=c4;
mk=-0.16398; R41=mk;
ck=0.63297;  R42=ck;
lkv=R29*EXP(-(POWER((LOG10(G7)-R33)/R37,2)))+R30*EXP(-(POWER((LOG10(G7)-R34)/R38,2)))+R31*EXP(-(POWER((LOG10(G7)-R35)/R39,2)))+R32*EXP(-(POWER((LOG10(G7)-R36)/R40,2)))+R41*LOG10(G7)+R42;

function av=alphav(f_Hz)
G7=f_Hz/1e9; 
a1=-0.07771; U29=a1;
a2=0.56727;  U30=a2;
a3=-0.20238; U31=a3;
a4=-48.29910; U32=a4;
a5=48.58330; U33=a5;
b1=2.33840;  U34=b1;
b2=0.95545;  U35=b2;
b3=1.14520;  U36=b3;
b4=0.791669; U37=b4;
b5=0.791459; U38=b5;
c1=-0.76284; U39=c1;
c2=0.54039;  U40=c2;
c3=0.26809;  U41=c3;
c4=0.116226; U42=c4;
c5=0.116479; U43=c5;
malpha=-0.053739; U44=malpha;
calpha=0.83433;   U45=calpha;
av=U29*EXP(-(POWER((LOG10(G7)-U34)/U39,2)))+U30*EXP(-(POWER((LOG10(G7)-U35)/U40,2)))+U31*EXP(-(POWER((LOG10(G7)-U36)/U41,2)))+U32*EXP(-(POWER((LOG10(G7)-U37)/U42,2)))+U33*EXP(-(POWER((LOG10(G7)-U38)/U43,2)))+U44*LOG10(G7)+U45;





function e=EXP(ex)
e=exp(ex);

function p=POWER(b,ex)
p=b.^ex;
function c=COS(a)
c=cos(a);

function l=LOG10(ex)
l=log10(ex);