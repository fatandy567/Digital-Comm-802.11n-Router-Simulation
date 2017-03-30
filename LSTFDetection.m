function [iscoming,nstatus]=LSTFDetection(En,x,dt,i)
nstatus=1;
iscoming=0;

if En ==0 % if the Block is turned off, then leave the block with returning 0;
    nstatus=0;
    iscoming=0;
    return 
end
points=floor(0.8*10^-6/dt);%points per STF symbol (0.8us each)
x_ext=[zeros(1,2*points) x zeros(1,points)];

x_delayed=x_ext(i:1:points+i-1);
x_now=x_ext(points+i:1:2*points+i-1);

c=sum(x_now.*conj(x_delayed));
p=sum(x_delayed.*conj(x_delayed));
m=abs(c)^2/p.^2;

if isnan(m)==true || m==inf %deal with 0/0;
    m=0;
end

if m>=0.8
    iscoming=1;
else
    iscoming=0;
end