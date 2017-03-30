function [FreqOffset,CoarseFreqOffset_Done] = CoarseFreqOffset(CoarseDet_En,x,i,fsMHz) 

nstatus=1;
if status==0
   freqoffset=0;
   nstatus=0;
   finish=0;
   return;
end

pts=floor(0.8*10^-6/dt);
x=x(i-pts+1:i);

pts=floor(0.8*10^-6/2/dt);

z=sum(x(1:pts).*conj(x(pts+1:2*pts)));
freqoffset=angle(z)/(-2*pi*pts*dt);

finish=1;
fprintf('Maximum Freq Error Can be estimated : %d \n',1/2/pts/dt)
