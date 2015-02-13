function [errorValue delta_V delta_W] = Train(Input, Output, V, W, delta_V, delta_W)



InputL1 = Input;
InputH = V' * InputL1;
[m n] = size(InputH);
OutputH = 1./(1+exp(-InputH));
Inputo = W'*OutputH;
clear m n;
[m n] = size(Inputo);
Outputo = 1./(1+exp(-Inputo));
errorValue = norm(Output - Outputo);
clear m n
[n p] = size(Output);
for i = 1 : n
for j = 1 : p
d(i,j) =(Output(i,j)-Outputo(i,j))*Outputo(i,j)*(1-Outputo(i,j));
end
end
Y = OutputH * d'; 
if nargin == 4
delta_W=zeros(size(W));
delta_V=zeros(size(V));
end
alpha=0.2;

delta_W= alpha.*Y;
error = W*d;
clear m n

[m n] = size(error);
for i = 1 : m

for j = 1 :n
d_star(i,j)= error(i,j)*OutputH(i,j)*(1-OutputH(i,j));
end

end
X = Input * d_star';
delta_V=alpha*X;

end