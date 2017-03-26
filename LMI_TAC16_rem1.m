function OmegaVal=LMI_TAC16_rem1(A,B,C,K,h,epsilon,delta)
% This MATLAB program checks the feasibility of LMIs from Remark 1 of the paper 
% A. Selivanov and E. Fridman, “Event-Triggered H-infinity Control: a Switching Approach,” IEEE Transactions on Automatic Control, vol. 61, no. 10, pp. 3221–3226, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B,C         - the parameters of the system (1); 
% K             - the static controller gain; 
% h             - the sampling time in (5); 
% epsilon       - the event-triggering parameter from (5); 
% delta         - a desired decay rate; 

% Output: 
% OmegaVal - the value of Omega from (5). If OmegaVal is empty, the LMIs are not feasible. 

n=size(A,1); 
l=size(C,1); 

%% Decision variables 
P=sdpvar(n); 
U=sdpvar(n); 
X=sdpvar(n,n,'f'); 
X1=sdpvar(n,n,'f'); 
P2=sdpvar(n,n,'f'); 
P3=sdpvar(n,n,'f'); 
Y1=sdpvar(n,n,'f'); 
Y2=sdpvar(n,n,'f'); 
Y3=sdpvar(n,n,'f'); 
Omega=sdpvar(l); 

%% The LMI for Xi
Xi=blkvar; 
Xi(1,1)=P+h*(X+X')/2; 
Xi(1,2)=h*X1-h*X; 
Xi(2,2)=-h*X1-h*X1'+h*(X+X')/2;
Xi=sdpvar(Xi);

%% The LMI with Psi0
Psi0=blkvar; 
Psi0(1,1)=A'*P2+P2'*A+2*delta*P-Y1-Y1'-(1/2-delta*h)*(X+X')+epsilon*C'*Omega*C; 
Psi0(1,2)=P-P2'+A'*P3-Y2+h*(X+X')/2; 
Psi0(1,3)=Y1'-P2'*B*K*C-Y3+(1-2*delta*h)*(X-X1); 
Psi0(1,4)=-P2'*B*K; 
Psi0(2,2)=-P3-P3'+h*U; 
Psi0(2,3)=Y2'-P3'*B*K*C-h*(X-X1); 
Psi0(2,4)=-P3'*B*K;
Psi0(3,3)=Y3+Y3'-(1/2-delta*h)*(X+X'-2*X1-2*X1'); 
Psi0(4,4)=-Omega; 
Psi0=sdpvar(Psi0);

%% The LMI with Psi1
Psi1=blkvar; 
Psi1(1,1)=A'*P2+P2'*A+2*delta*P-Y1-Y1'-(X+X')/2+epsilon*C'*Omega*C; 
Psi1(1,2)=P-P2'+A'*P3-Y2;  
Psi1(1,3)=Y1'-P2'*B*K*C-Y3+X-X1; 
Psi1(1,4)=h*Y1'; 
Psi1(1,5)=-P2'*B*K; 
Psi1(2,2)=-P3-P3'; 
Psi1(2,3)=Y2'-P3'*B*K*C; 
Psi1(2,4)=h*Y2'; 
Psi1(2,5)=-P3'*B*K;
Psi1(3,3)=Y3+Y3'-1/2*(X+X'-2*X1-2*X1'); 
Psi1(3,4)=h*Y3'; 
Psi1(4,4)=-h*U*exp(-2*delta*h); 
Psi1(5,5)=-Omega; 
Psi1=sdpvar(Psi1);

%% Solution of LMIs
LMIs=[P>=0, U>=0, Omega>=0, Xi>=0, Psi0<=0, Psi1<=0]; 
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 && all(primal([1,2,4])>0)
        OmegaVal=value(Omega); 
    end
else
    yalmiperror(sol.problem) 
end
