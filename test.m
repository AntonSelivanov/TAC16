% This MATLAB program checks the feasibility of LMIs from Theorems 1, 2 and Remarks 1, 3 of the paper 
% A. Selivanov and E. Fridman, “Event-Triggered H-infinity Control: a Switching Approach,” IEEE Transactions on Automatic Control, vol. 61, no. 10, pp. 3221–3226, 2016.

%% Parameters from Example 1
A=[0 1; 0 -3]; B=[0; 1]; C=[1 0]; K=3; 

% Theorem 1
h=.899; epsilon=.554; delta=.24; 
display(['Theorem 1: Omega=' num2str(LMI_TAC16_th1(A,B,C,K,h,epsilon,delta))]); 

% Remark 1
h=1.115; epsilon=4.6e-3; delta=.24; 
display(['Remark 1: Omega=' num2str(LMI_TAC16_rem1(A,B,C,K,h,epsilon,delta))]); 

%% Parameters from Example 2
A=[0 1 0 0; 0 0 -1 0; 0 0 0 1; 0 0 10/3 0]; 
B1=ones(4,1); 
B2=[0; .1; 0; -1/30]; 
C1=ones(1,4); 
D1=.1; 
C2=eye(4); 
D2=zeros(4,1); 
K=[2.9129 10.4357 287.9029 160.3271]; 

% Theorem 2
gamma=100; h=.065; etaM=.1; epsilon=.044; delta=0; 
display('Theorem 2: Omega='); 
disp(LMI_TAC16_th2(A,B1,B2,C1,D1,C2,D2,K,gamma,h,etaM,epsilon,delta)); 

% Remark 3
gamma=100; h=.036; etaM=.1; epsilon=.033; delta=0; 
display('Remark 3: Omega='); 
disp(LMI_TAC16_rem3(A,B1,B2,C1,D1,C2,D2,K,gamma,h,etaM,epsilon,delta)); 
