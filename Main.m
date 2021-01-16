%Authors: Ahmed A.Allam,Marina Elgohary & Amr Mamdouh
%Date:15/1/2021

%%
%1) Generate the BPSK symbols (�1) of N = 5 users, each composed of randomly generated K = 10 symbols.
N=5;
K=10;
S=zeros(N,K);
for i=1:1:N
    S(i,:)=randi([0 1],K,1);
end
S(S==0)=-1;
%%
%2)Generate the maximal length spreading codes for the N users.
%The following genrator polynomials are obtained from the following link
%https://www.mathworks.com/help/comm/ref/commsrc.pn.html
NGenPolynomial=7;
C=zeros(N,NGenPolynomial);
switch NGenPolynomial
    case 7
      Genrator_Polynomial=[3 2 0];
    case 15
      Genrator_Polynomial=[4 3 0];
    case 31
      Genrator_Polynomial=[5 3 0];
end

 
 for i=1:1:N
      h=commsrc.pn('GenPoly',Genrator_Polynomial,'InitialStates',de2bi(i,floor(log2(NGenPolynomial+1))),'CurrentStates',de2bi(i,floor(log2(NGenPolynomial+1))),'Shift',0,'NumBitsOut',NGenPolynomial);
     C(i,:)=generate(h);
 end
C(C==0)=-1;
%%
%3)Spread the signal by multiplying each BPSK symbol with the spreading code
Spreaded_Code =zeros(N,K*NGenPolynomial);
for i=1:1:N
    Spreaded_Code(i,:)=kron(S(i,:),C(i,:));
end
%%
%4) Convolute the spreaded signal with the channel impulse response.
