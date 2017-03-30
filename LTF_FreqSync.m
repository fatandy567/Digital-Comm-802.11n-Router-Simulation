function[freqoffset,finish,nstatus]=LTF_FreqSync(status,x,dt,i)

nstatus=1;
if status==0
   freqoffset=0;
   nstatus=0;
   finish=0;
   return;
end
pts=floor(8*10^-6/dt);
x=x(i-pts+1:i);

TGI=0.8*10^-6;
pts=floor((8*10^-6-2*TGI)/2/dt);

idxoffset=floor(pts/2);%use the middle point of the time signal;
z=sum(x(idxoffset+1:pts+idxoffset).*conj(x(pts+idxoffset+1:2*pts+idxoffset)));
freqoffset=angle(z)/(-2*pi*pts*dt);

finish=1;
fprintf('Maximum Freq Error Can be estimated : %d \n',1/2/pts/dt)