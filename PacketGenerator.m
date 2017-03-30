function [y,NSYM,Mod]=PacketGenerator(t,itx,MCS,HTLENGTH,CarrierFreq)

NumHTLTF=1;
STBC=0;
HTLENGTH=HTLENGTH;

state=0;
LSTF=0;
HTLTF1=0;
HTSIG=0;
HTLTF2=0;
HTLTF3=0;
HTLTF4=0;
DATA=0;

%tHT_LTF1=8*10^-6;
%tHT_SIG=tHT_LTF1+8*10^-6;
%tHT_LTFs=tHT_SIG+8*10^-6;
%tHT_DATA=tHT_LTFs+(NumHTLTF-1)*4*10^-6;
tHT_LTF1=8*10^-6;
tHT_DATA=tHT_LTF1+8*10^-6;

[LSTF,~,Mod,state]=FieldGenerator('L-STF',t,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
[HTLTF1,~,Mod,state]=FieldGenerator('HT-LTF1',t-tHT_LTF1,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
[DATA,NSYM,Mod,state]=FieldGenerator('DATA',t-tHT_DATA,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
%[HTSIG,state]=FieldGenerator('HT-SIG',t-tHT_SIG,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
%[HTLTF2,state]=FieldGenerator('HT-LTF2',t-tHT_LTFs-(2-2)*4*10^-6,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
%[HTLTF3,state]=FieldGenerator('HT-LTF3',t-tHT_LTFs-(3-2)*4*10^-6,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
%[HTLTF4,state]=FieldGenerator('HT-LTF4',t-tHT_LTFs-(4-2)*4*10^-6,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);
%[DATA,state]=FieldGenerator('DATA',t-tHT_DATA,itx,MCS,NumHTLTF,STBC,HTLENGTH,state);

%y=LSTF+HTLTF1+HTSIG+HTLTF2+HTLTF3+HTLTF4+DATA;
%LSTF=0;
%HTLTF1=0;

FreqOffset=0;%Carrier frequency offset
y=(LSTF+HTLTF1+DATA).*exp(j*2*pi*CarrierFreq*(1+FreqOffset)*t);