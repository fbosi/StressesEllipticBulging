function [ dy2 ] = derivate2( y,x ) %differentiating y wrt x
% Function for calculating 2nd order derivatve
% taking input as arrays and output as derivate of same lenth array

n=length(y);

if length(x) == 1
    dx = x;
elseif length(x) > 1
    dx = x(2)-x(1);
end

for i=1:n
    if i==1
        dy2(i)=(y(i+1)-2*y(i)+y(i+1))/dx^2;     %((x(i+1)-x(i))^2+(x(i)-x(i+1))^2);      %forward difference   % symmetric assumption
    elseif i==n
        dy2(i)=dy2(i-1);%(y(i)-y(i-1))/(x(i)-x(i-1));      %backward difference
    else
        dy2(i)=(y(i+1)-2*y(i)+y(i-1))/dx^2;     %((x(i+1)-x(i))^2+(x(i)-x(i-1))^2);      %central difference
    end
end

