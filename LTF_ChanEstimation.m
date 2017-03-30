function [H,Done,nstatus]=LTF_ChanEstimation(status,x,dt,i)
nstatus=1;
if status==0
   nstatus=0;
   Done=0;
   H=zeros(1,57);
   return;
end

pts=floor(8*10^-6/dt);
x=x(i-pts+1:i);

TGI=0.8*10^-6;
idxoffset=2*TGI/dt;
T=8*10^-6;
pts=((T-2*TGI)/2/dt);
X=fftshift(fft(x(idxoffset+1:idxoffset+pts)));

%SamplingFreq=1/dt;
%FreqReso=(SamplingFreq/2)/(pts/2);
%ptsShift=floor(abs(FreqOffset)/FreqReso)*sign(FreqOffset);
%X=circshift(X',ptsShift*-1)';

lx=length(X);
r=([X(lx/2-27:lx/2) X(lx/2+1:lx/2+29)]);
%[V,~]=max(abs(r));

%idx=(r>V/2)+(r<-V/2);
%[D,I]=max(idx);

HT_LTF1=[1,1,1,1,-1, -1,1,1,-1,1, -1,1,1,1,1, 1,1,-1,-1,1 ,1,-1,1,-1,1, 1,1,1,0,...
            1,-1,-1,1,1, -1,1,-1,1,-1, -1,-1,-1,-1,1, 1,-1,-1,1,-1, 1,-1,1,1,1, 1,-1,-1];

H=r.*HT_LTF1;
Done=1;