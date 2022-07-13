clear all
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MONTE CARLO
%   BASKET DOUBLE BARRIER CALL OPTION
%   2 asset sottostanti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 0.7;    % coefficiente di correlazione asset sottostanti    
N = 10000;    % numero traiettorie generate
r = 0.05;     % interest rate 
T = 1;        % maturity in anni

sigma1 = 0.25;% volatilita asset1
sigma2 = 0.25;% volatilita asset2

M = 100;      % numero di intervalli temporali 

k = T/M;

K = 1;       % strike price
BD = 1;      % lower barrier L
Bu = 2;      % upper barrier U

B_down = BD * exp(-0.5826*sigma1*sqrt(k)); %correzione barriera 
B_up   = BU * exp(+0.5826*sigma1*sqrt(k)); %correzione barriera 

sx1 = linspace(0,2,41); 
sx2 = linspace(0,2,41);

%matrice della soluzione
prezzoMC = zeros(length(sx2),length(sx1)) ;
prezzoMC_BC = zeros(length(sx2),length(sx1)) ;
prezzoBB = zeros(length(sx2),length(sx1)) ;
prezzoBB_BC = zeros(length(sx2),length(sx1)) ;

for xx = 1:length(sx1)
    for yy = 1:length(sx2)
   
    
S1_0 = sx1(xx) ;
S2_0 = sx2(yy) ;

 sum_payoffMC = 0;           % MC standard
 sum_payoffVA = 0;  % MC variabili antitetiche
 sum_payoffBB = 0;        % Brownian bridge
 sum_payoffVABB = 0; % brownian bridge + VA

 Z1=0;
 Z2=0;
 Z3=0;
 Z4=0;
 
if sx1(xx) + sx2(yy) < BU
    if sx1(xx) + sx2(yy) > BD
    
for i=1:N
   
    k=T/M;
    S1a=[];
    S1a=zeros(1,M+1); 
    S2a=[];
    S2a=zeros(1,M+1); 
    S1b=[];
    S1b=zeros(1,M+1); 
    S2b=[];
    S2b=zeros(1,M+1); 
    
    t=zeros(1,M+1);
    
    S1a(1)=S1_0;
    S2a(1)=S2_0;
    S1b(1)=S1_0;
    S2b(1)=S2_0;

 for i=1:M
    
    g=randn;
    f=randn;
    S1a(i+1)=S1a(i)*exp( (r - 0.5*(sigma1)^2)*k + sigma1*sqrt(k)*g);
    S1b(i+1)=S1b(i)*exp( (r - 0.5*(sigma1)^2)*k + sigma1*sqrt(k)*(-g));
    
    S2a(i+1)=S2a(i)*exp( (r - 0.5*(sigma2)^2)*k + sigma2*sqrt(k)*(...
        (sqrt(1 - (rho)^2))*f + (rho)*g) );
    S2b(i+1)=S2b(i)*exp( (r - 0.5*(sigma2)^2)*k + sigma2*sqrt(k)*(...
        (sqrt(1 - (rho)^2))*(-f) + (rho)*(-g)) );
    
 end
    Pa=S1a(M+1);
    Pb=S2a(M+1);
    Pc=S1b(M+1);
    Pd=S2b(M+1);

% CONDIZIONE:la somma del prezzo delle azioni deve rimanere nell'intervallo
% definito dalle 2 barriere

if min(S1a + S2a) > BD
    if max(S1a + S2a) < BU
       Z1 = max((Pa + Pb - K) , 0);
    else
        Z1 = 0;
    end    
else
        Z1 = 0;
end


if min(S1a + S2a) > B_down
    if max(S1a + S2a) < B_up
       Z2 = max((Pa + Pb - K) , 0);
    else
        Z2 = 0;
    end
else
    Z2 = 0;
end

%brownian bridge
    ProbUP = exp( - 2 .* log( BU./ (S1a(1)+S2a(1)) ) .*log(...
        BU./(Pa+Pb) ) ./ (T.*(sigma1^2)) ) ; 
    ProbDOWN = exp( -2.*log( BD./(S1a(1)+S2a(1)) ) .* log(...
        BD./(Pa+Pb) ) ./ (T.*(sigma1^2)) ) ;
    
    zeta1 = 1 - ProbUP ;  
    zeta2 = 1 - ProbDOWN ;
    
if min(S1a + S2a) > BD      %B_down  %BD
   if max(S1a + S2a) < BU   %B_up    %BU      
      Z3 = zeta1.*zeta2.*max(Pa + Pb - K, 0);
   else
      Z3=0;
   end  
else
    Z3 = 0; 
end
    
    
     %brownian bridge barriere corrette
ProbUP_meno = exp( -2.*log(B_up./(S1a(1)+S2a(1))).*log(...
    B_up./(Pa+Pb)) ./ (T.*(sigma1^2)) ) ; 
ProbDOWN_meno = exp( -2.*log(B_down./(S1a(1)+S2a(1))).*log(...
    B_down./(Pa+Pb)) ./ (T.*(sigma1^2)) ) ;
    
zeta1b = 1 - ProbUP_meno ;  
zeta2b = 1- ProbDOWN_meno ;
    
if min(S1a + S2a) >  B_down     %BD
   if max(S1a + S2a) < B_up     %BU   
      Z4 = max(Pa + Pb - K, 0) * zeta1b * zeta2b;
   else
      Z4 = 0;
   end 
else
    Z4 = 0;
end



    sum_payoffMC = sum_payoffMC + Z1;      % MC standard
    sum_payoffVA = sum_payoffVA + Z2;      % MC B corrette
    sum_payoffBB = sum_payoffBB +Z3;       % Brownian bridge
    sum_payoffVABB = sum_payoffVABB + Z4;  % brownian bridge B corrette
    
end
    
prezzoMC(xx,yy) = exp(-r*T)*(sum_payoffMC./N);
prezzoMC_BC(xx,yy) = exp(-r*T)*(sum_payoffVA./(N));
prezzoBB(xx,yy) = exp(-r*T)*(sum_payoffBB./N);
prezzoBB_BC(xx,yy) = exp(-r*T)*(sum_payoffVABB./(N));


    else
prezzoMC(xx,yy) = 0;
prezzoMC_BC(xx,yy) = 0;
prezzoBB(xx,yy) = 0;
prezzoBB_BC(xx,yy) = 0;
    end
    
else
prezzoMC(xx,yy) = 0;
prezzoMC_BC(xx,yy) = 0;
prezzoBB(xx,yy) = 0;
prezzoBB_BC(xx,yy) = 0;
end

    end
end
toc

% comandi di plot
% figure
% hold on
% surf(sx1,sx2,prezzoMC)
% title({'Two asset Duoble Barrier option ','K = 1, BU = 2, BD = 1,N = $5 \cdot 10^5$, M = 1000, T = 1, $\sigma_1=\sigma_2$ = 0.25, r = 0.05;'},'interpreter','latex')
% xlabel({'$S_1$'},'interpreter','latex')
% ylabel({'$S_2$'},'interpreter','latex')
% zlabel({'option price'},'interpreter','latex')
% hold off
 
% figure
% hold on
% surf(sx1,sx2,prezzoBB)
% title({'Two asset Duoble Barrier option Brownian Bridge ','K = 1, BU = 2, BD = 1,N = $5 \cdot 10^5$, M = 1000, T = 1, $\sigma_1=\sigma_2$ = 0.25, r = 0.05;'},'interpreter','latex')
% xlabel({'$S_1$'},'interpreter','latex')
% ylabel({'$S_2$'},'interpreter','latex')
% zlabel({'option price'},'interpreter','latex')
% hold off
 
% figure
% hold on
% surf(sx1,sx2,prezzoMC_BC)
% title({'Two asset Duoble Barrier option barriere corrette ','K = 1, BU = 2, BD = 1,N = $ 5 \cdot 10^5$, M = 1000, T = 1, $\sigma_1=\sigma_2$ = 0.25, r = 0.05;'},'interpreter','latex')
% xlabel({'$S_1$'},'interpreter','latex')
% ylabel({'$S_2$'},'interpreter','latex')
% zlabel({'option price'},'interpreter','latex')
% hold off