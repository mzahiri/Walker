function [V,W]= Nnet(Input, Output)


if max(abs(Input(:)))> 1 

Inputnor = Input / max(abs(Input(:)));

else

Inputnor = Input;

end

if max(abs(Output(:))) >1
Outputnor = Output / max(abs(Output(:)));
else
Outputnor = Output;

end
m = 6; %# hidden layers
[l,b] = size(Input);
[n,a] = size(Output);

V = rand(l,m); 
W = rand(m,n); 

count = 0;

[error delta_V delta_W] = Train(Inputnor,Outputnor,V,W);
errortot=[];

while count <=45000
count = count + 1;

W=W+delta_W;
V=V+delta_V;
count;
[error delta_V delta_W]=Train(Inputnor,Outputnor,V,W,delta_V,delta_W);
errortot=[errortot;error];
end
error
plot(errortot);
%Error_Mat(count)=error;
%Error_Rate=sum(Error_Mat)/count;
%figure;
%y=[1:count];
%plot(y, Error_Mat);

end