
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   DUOBLE BARRIER OPTION MONTE CARLO METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

tic;
%% parametri simulazione

S0 = linspace(0,3,301) ;
K = 1;
BU = 2;
BD = 1;
N = 10000;
M = 5000;
T = 1;
sigma = 0.25;
r = 0.05;
%% inizializzazione variabili
sum_payoff = 0;
sum_payoffVA = 0;
sum_payoffBB = 0;
sum_payoffVABB = 0;
sum_analitica1 = 0;
sum_analitica2 = 0;
Z1 = 0;
Z2 = 0;
Z3 = 0;
Z4 = 0;
k=T/M;
price1= zeros(1,length(S0));% prezzo opzione attualizzato
priceMC_BC = zeros(1,length(S0));
priceBB = zeros(1,length(S0)); % prezzo brownian bridge
priceBB_BC = zeros(1,length(S0));  % prezzo con metodi combinati
PriceA = zeros(1,length(S0));

%% correzione per la barriera da osservazione discreta

B_up = BU*exp(+0.5826*sigma*sqrt(k)); 
B_down = BD*exp(-0.5826*sigma*sqrt(k)); 

%% ciclo

for j=1:length(S0)
    
    sum_payoff = 0;
    sum_payoffVA = 0;
    sum_payoffBB = 0;
    sum_payoffVABB = 0;
    
    if S0(j) > BD
        if S0(j) < BU
for i=1:N
    
    ProbUP = 1 ; 
    ProbDOWN = 1 ;
    ProbUP1 = 1 ; 
    ProbDOWN1 = 1 ;
    
    [Sa]=traiettoria1_barrier_option_1d(S0(j),T,sigma,r,M);
    % uso r:intensita istantanea, sto facendo pricing
    %P=S(M+1) corrisponde al valore dell'optione alla scadenza
    Pa = Sa(M+1);
    %Pb = Sb(M+1);
    
    if min(Sa) > BD        %B_down %
        if max(Sa) < BU    %  B_up %
            Z1=max(Pa - K, 0);
        else
            Z1=0;
        end  
    else 
        Z1 = 0;
    end
    
    % MC barriere corrette
    
    if min(Sa) > B_down %BD
        if max(Sa) < B_up  %BU
           Z2=max(Pa - K, 0);
        else
           Z2=0;
        end 
    else 
        Z2 = 0;
    end
    
    %brownian bridge
    ProbUP = exp( -2.*log(BU./S0(j)).*log(BU./Pa) ./ (T.*(sigma^2)) ) ; 
    ProbDOWN = exp( -2.*log(BD./S0(j)).*log(BD./Pa) ./ (T.*(sigma^2)) ) ;
    
    zeta1 = 1 - ProbUP ;              
    zeta2 = 1 - ProbDOWN ;
    
    if min(Sa) > BD %B_down %BD
        if max(Sa) < BU   %B_up  %BU
            
      Z3 = zeta1.*zeta2.*max(Pa - K, 0);
    else
      Z3 = 0;
        end 
    else
        Z3 =0;
    end
    
     %brownian bridge + barriere corrette
    ProbUP1 = exp( -2.*log(B_up./S0(j)).*log(B_up./Pa) ./ (T.*(sigma^2)) ); 
    ProbDOWN1=exp(-2.*log(B_down./S0(j)).*log(B_down./Pa)./(T.*(sigma^2)));
    
    zeta1b = 1 - ProbUP1; 
                 
    zeta2b = 1- ProbDOWN1 ;
    
    if min(Sa) >  B_down %BD
        if max(Sa) < B_up  %BU
            
      Z4=max(Pa - K, 0) * zeta1b * zeta2b;
    else
      Z4 = 0;
        end 
    else
        Z4 = 0;
    end
    
    sum_payoff = sum_payoff + Z1;
    sum_payoffVA = sum_payoffVA  + Z2;
    sum_payoffBB = sum_payoffBB +Z3;
    sum_payoffVABB = sum_payoffVABB + Z4;
    
   
end
        

%% SOLUZIONE ANALITICA

sum_analtica1 = 0;
sum_analtica2 = 0;
kk1 = 0;
kk2 = 0;
kk3 = 0;
kk4 = 0;
kk5 = 0;
kk6 = 0;
d1=0;
d2=0;
d3=0;
d4=0;
sum_analitica1 = 0;
sum_analitica2 = 0;

% sol. analitica

for ww = -12:12
    nn = 0;
    nn = ww;
FF = BU;
d1 = (log((S0(j).*BU.^(2*nn))./(K*BD.^(2*nn))) + T*0.5*sigma.^2)./(...
    sigma*sqrt(T));
d2 = (log((S0(j).*BU.^(2*nn))./(FF*BD.^(2*nn))) + T*0.5*sigma.^2)./(...
    sigma*sqrt(T));
d3 = (log((BD.^(2*nn+2))./(K.*S0(j).*BU.^(2*nn))) + T*0.5*sigma.^2)./(...
    sigma*sqrt(T));
d4 = (log((BD.^(2*nn+2))./(FF.*S0(j).*BU.^(2*nn))) + T*0.5*sigma.^2)./(...
    sigma*sqrt(T)) ;

mu1 = 1 ;
mu2 = 0 ;
mu3 = 1 ;
   
kk1 = ( (BU.^(nn) ) ./ ( BD.^(nn) ) ).^(mu1);
kk2 = (BD./S0(j)).^(mu2);
kk3 = ( ( BD .^ (nn + 1) ) ./ ( (BU .^ (nn) ) .* S0(j) ) ).^(mu3);
kk4 = (( BU.^(nn) ) ./ ( BD.^(nn) )) .^ (mu1 - 2);
kk5 = (( BD./S0(j) ).^(mu2));
kk6 = (( BD.^(nn+1) )./( (BU.^(nn)) .* S0(j) ) ) .^ (mu3 - 2);
    
SolA =  kk1 * kk2 *( normcdf(d1) - normcdf(d2) ) - kk3 *( normcdf(d3) ...
        - normcdf(d4) ) ; 
SolB = kk4 * kk5 *( normcdf(d1 - sigma*sqrt(T) ) - normcdf(d2 ...
    - sigma*sqrt(T))) - kk6 *(normcdf(d3 - sigma*sqrt(T)) - normcdf(d4 ...
    - sigma*sqrt(T)));
    
    sum_analitica1 = sum_analitica1 + SolA;
    sum_analitica2 = sum_analitica2 + SolB;
end
        else 
        sum_analitica1 = 0;
        sum_analitica2 = 0;
        end
         
    else
    sum_analitica1 = 0;
    sum_analitica2 = 0;
        
    end

price1(j)=exp(-r*T)*(sum_payoff/N);% prezzo opzione attualizzato
priceMC_BC(j)=exp(-r*T)*(sum_payoffVA/(N));
priceBB(j)=exp(-r*T)*(sum_payoffBB/N); % prezzo brownian bridge
priceBB_BC(j)=exp(-r*T)*(sum_payoffVABB/(N));  % prezzo con metodi combinati
PriceA(j) = S0(j)*exp(-r*T)*sum_analitica1 - K*exp(-r*T)*sum_analitica2 ;

end
