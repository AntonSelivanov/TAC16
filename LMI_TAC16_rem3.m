function OmegaVal=LMI_TAC16_rem3(A,B1,B2,C1,D1,C2,D2,K,gamma,h,etaM,epsilon,delta)
% This MATLAB program checks the feasibility of LMIs from Remark 3 of the paper 
% A. Selivanov and E. Fridman, “Event-Triggered H-infinity Control: a Switching Approach,” IEEE Transactions on Automatic Control, vol. 61, no. 10, pp. 3221–3226, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and LMILAB solver of Robust Control Toolbox (http://www.mathworks.com/products/robust/)

% Input: 
% A,B1,B2,C1,
% D1,C2,D2      - the parameters of the system (15); 
% K             - the static controller gain; 
% gamma         - the L_2-gain from (20); 
% h             - the sampling time in (5); 
% etaM          - the upper bound for the overall network-induced delay; 
% epsilon       - the event-triggering parameter from (5); 
% delta         - a desired decay rate; 

% Output: 
% OmegaVal - the value of Omega from (5). If OmegaVal is empty, the LMIs are not feasible. 

n=size(A,1); 
l=size(C2,1); 

%% Decision variables 
P=sdpvar(n); 
S0=sdpvar(n); 
S1=sdpvar(n); 
R0=sdpvar(n); 
R1=sdpvar(n); 
G0=sdpvar(n,n,'f'); 
G1=sdpvar(n,n,'f'); 
Omega=sdpvar(l); 

%% Notations 
tauM=h+etaM; 
H=etaM^2*R0+h^2*R1; 

%% The LMI for Psi (tau(t)\in(etaM,tauM])
Psi=blkvar; 
Psi(1,1)=A'*P+P*A+2*delta*P+S0-exp(-2*delta*etaM)*R0+C1'*C1; 
Psi(1,2)=exp(-2*delta*etaM)*R0; 
Psi(1,4)=P*B2*K*C2+C1'*D1*K*C2; 
Psi(1,5)=P*B2*K+C1'*D1*K; 
Psi(1,6)=P*B1; 
Psi(1,7)=P*B2*K*D2+C1'*D1*K*D2; 
Psi(1,8)=A'*H; 
Psi(2,2)=exp(-2*delta*etaM)*(S1-S0-R0)-exp(-2*delta*tauM)*R1; 
Psi(2,3)=exp(-2*delta*tauM)*G1; 
Psi(2,4)=exp(-2*delta*tauM)*(R1-G1); 
Psi(3,3)=-exp(-2*delta*tauM)*(R1+S1);
Psi(3,4)=exp(-2*delta*tauM)*(R1-G1'); 
Psi(4,4)=exp(-2*delta*tauM)*(G1+G1'-2*R1)+(D1*K*C2)'*D1*K*C2+epsilon*C2'*Omega*C2; 
Psi(4,5)=(D1*K*C2)'*D1*K; 
Psi(4,7)=(D1*K*C2)'*D1*K*D2+epsilon*C2'*Omega*D2; 
Psi(4,8)=(B2*K*C2)'*H; 
Psi(5,5)=(D1*K)'*D1*K-Omega; 
Psi(5,7)=(D1*K)'*D1*K*D2; 
Psi(5,8)=(B2*K)'*H; 
Psi(6,6)=-gamma^2*eye(size(B1,2)); 
Psi(6,8)=B1'*H; 
Psi(7,7)=(D1*K*D2)'*D1*K*D2+epsilon*D2'*Omega*D2-gamma^2*eye(size(D2,2)); 
Psi(7,8)=(B2*K*D2)'*H; 
Psi(8,8)=-H; 
Psi=sdpvar(Psi); 

%% The LMI for Phi (tau(t)\in[0,etaM])
Phi=blkvar; 
Phi(1,1)=A'*P+P*A+2*delta*P+S0-exp(-2*delta*etaM)*R0+C1'*C1; 
Phi(1,2)=exp(-2*delta*etaM)*G0; 
Phi(1,4)=P*B2*K*C2+exp(-2*delta*etaM)*(R0-G0)+C1'*D1*K*C2; 
Phi(1,5)=P*B2*K+C1'*D1*K; 
Phi(1,6)=P*B1; 
Phi(1,7)=P*B2*K*D2+C1'*D1*K*D2; 
Phi(1,8)=A'*H; 
Phi(2,2)=exp(-2*delta*etaM)*(S1-S0-R0)-exp(-2*delta*tauM)*R1; 
Phi(2,3)=exp(-2*delta*tauM)*R1; 
Phi(2,4)=exp(-2*delta*etaM)*(R0-G0'); 
Phi(3,3)=-exp(-2*delta*tauM)*(R1+S1);
Phi(4,4)=exp(-2*delta*etaM)*(G0+G0'-2*R0)+(D1*K*C2)'*D1*K*C2+epsilon*C2'*Omega*C2; 
Phi(4,5)=(D1*K*C2)'*D1*K; 
Phi(4,7)=(D1*K*C2)'*D1*K*D2+epsilon*C2'*Omega*D2; 
Phi(4,8)=(B2*K*C2)'*H; 
Phi(5,5)=(D1*K)'*D1*K-Omega; 
Phi(5,7)=(D1*K)'*D1*K*D2; 
Phi(5,8)=(B2*K)'*H; 
Phi(6,6)=-gamma^2*eye(size(B1,2)); 
Phi(6,8)=B1'*H; 
Phi(7,7)=(D1*K*D2)'*D1*K*D2+epsilon*D2'*Omega*D2-gamma^2*eye(size(D2,2)); 
Phi(7,8)=(B2*K*D2)'*H; 
Phi(8,8)=-H; 
Phi=sdpvar(Phi); 

%% Park's conditions
Park0=[R0 G0; G0' R0]; 
Park1=[R1 G1; G1' R1]; 

%% Solution of LMIs
LMIs=[P>=0, S0>=0, S1>=0, R0>=0, R1>=0, Omega>=0, Phi<=0, Psi<=0, Park0>=0, Park1>=0]; 
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    primal=check(LMIs); 
    if min(primal)>=0 && primal(1)>0
        OmegaVal=value(Omega); 
    end
else
    yalmiperror(sol.problem) 
end
