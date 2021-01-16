%Authors: Ahmed A.Allam,Marina Elgohary & Amr Mamdouh
%Date:15/1/2021

%%
%1) Generate the BPSK symbols (±1) of N = 5 users, each composed of randomly generated K = 10 symbols.
N=5;
K=10;
S=zeros(N,K);
for i=1:1:N
    S(i,:)=randi([0 1],K,1);
end
S(S==0)=-1;
%figure(1);

% for i=1:N
%     subplot(N,1,i)
%     stairs(S(i,:));
%     ylim([-2 2]);
%     title(['Sent Message of user ',num2str(i)]);
% end

%%
%2)Generate the maximal length spreading codes for the N users.
%The following genrator polynomials are obtained from the following link
%https://www.mathworks.com/help/comm/ref/commsrc.pn.html
NGenPolynomial=7;
%for NGenPolynomial= [7 15 31] %case of plotting only
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
Spreaded_Signal =zeros(N,K*NGenPolynomial);
for i=1:1:N
    Spreaded_Signal(i,:)=kron(S(i,:),C(i,:));
end
%%
%4) Convolute the spreaded signal with the channel impulse response.
%we have 3 casses of channel needed in the report
ChannelCase=1;
%for ChannelCase=[1 2 3] in case of plotting only
 h=zeros(N,3);
h(1:N,:)=repmat([1 0 0],N,1);
switch ChannelCase
    case 2
      h(1,:)=[1 0.8 0.1];
    case 3
        h(1,:)=[1 0.8 0.1];
        h(2,:)=[1 0.2 0.35];
end
TransmittedSignal=zeros(N,K*NGenPolynomial+2);
for i=1:1:N
    TransmittedSignal(i,:)=conv(Spreaded_Signal(i,:),h(i,:));
end
RecievedSignal=sum(TransmittedSignal);
%%
%5) At the base station, apply correlator for each of the N users.
DespreadedMessage=zeros(N,K);
for i=1:1:N
   for j=1:NGenPolynomial:NGenPolynomial*K
       DespreadedSymbol=RecievedSignal(j:j+NGenPolynomial-1).*C(i,:);
       DespreadedMessage(i,ceil(j/NGenPolynomial))=(1/NGenPolynomial)*sum(DespreadedSymbol);
   end
end
%%
%6)Finally, apply hard decision decoding (threshold = 0) to estimate the transmitted BPSK symbols.
DecodedMessage=DespreadedMessage;
DecodedMessage(DecodedMessage>=0)=1;
DecodedMessage(DecodedMessage<0)=-1;
%%
%7)Plotting the needed figures
% figure();
% 
% stairs(RecievedSignal);
% ylim([-8 8]);
% title(['Recieved Message before Despreading ',' with PN sequence of Length: ',num2str(NGenPolynomial),'  case: ',num2str(ChannelCase)]);
% 
% figure();
% for i=1:N
%     subplot(N,1,i)
%     stairs(DespreadedMessage(i,:));
%     ylim([-2 2]);
%     title(['Despreaded Message of user: ',num2str(i),' with PN sequence of Length: ',num2str(NGenPolynomial),'  case: ',num2str(ChannelCase)]);
% end
% figure();
% for i=1:N
%     subplot(N,1,i)
%     stairs(DecodedMessage(i,:));
%     ylim([-2 2]);
%     title(['Estimated information of user: ',num2str(i),' with PN sequence of Length: ',num2str(NGenPolynomial),'  case: ',num2str(ChannelCase)]);
% end
%end
%end


