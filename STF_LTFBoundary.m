function [finish,nstatus]=STF_LTFBoundary(En,x,dt,i)
nstatus=1;
m=0;
if En ==0
    nstatus=0;
    finish=0;
    return 
end
total_pts=floor(8*10^-6/dt);%points per LTF symbol (8us each)
cal_pts=floor(2*0.8*10^-6/dt);%points for calculations
x_ext=[zeros(1,total_pts) x];

x_seq=x_ext(i:1:total_pts+i-1);
x_1=x_seq(1:cal_pts);
x_2=x_seq(end-cal_pts+1:end);


c=sum(x_1.*conj(x_2));
p=sum(abs([x_1 x_2]).^2);
m=abs(c)/p;

if isnan(m)==true || m==inf %deal with 0/0;
    m=0;
end

if m>=0.48
    finish=1;
else
    finish=0;
end