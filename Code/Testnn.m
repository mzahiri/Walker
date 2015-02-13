clear
clc
%Input= 10.* [0 0.4   0.5  0.7 0.8 0.9 1 ];

%Output= 1/100.*[2.0000    4.6250    8.0000   12.8750   20.0000   30.1250   34 ];
%Input=[0 1 0 1;0 0 1 1];
%Output=[0 1 1 0];
Input=linspace(0.2,1,30);
Output=(5.*Input.^3+Input);


[V,W]=Nnet(Input, Output);

Inputsam=[0.2 0.3 0.4 0.5 0.6 0.7 0.9 ];
%Inputsam=[0 1;1 1];

if max(abs(Inputsam(:)))> 1 
Norm_Input = Inputsam / max(abs(Input(:)));
else
Norm_Input = Inputsam;
end

Input=Norm_Input

Output_of_InputLayer = Input;
Input_of_HiddenLayer = V' * Output_of_InputLayer;
[m n] = size(Input_of_HiddenLayer);
Output_of_HiddenLayer = 1./(1+exp(-Input_of_HiddenLayer));
Input_of_OutputLayer = W'*Output_of_HiddenLayer;
clear m n;
[m n] = size(Input_of_OutputLayer);
Output_of_OutputLayer = 1./(1+exp(-Input_of_OutputLayer));
Output_of_OutputLayer=Output_of_OutputLayer*max(Output)
Output=(5.*Inputsam.^3+Inputsam)

%plot(Input,Output,'o')
%hold on
%plot(Input,Output_of_OutputLayer,'x')