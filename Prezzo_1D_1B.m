%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Barrier Option misure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
tic;

price1=[]; 
price2=[]; 
error1=[]; 
error2=[]; 

sum_payoff1 = 0;
sum_payoff2 = 0;

for W=1:100

 time_int(W)= 50*W;
 
 M = time_int(W);
 N = 100000;
 S0 = 50;
 K = 50;
 B = 30;
 
 T=1;
 k=T/M;
 sigma=0.2;
 r=0.1;  
 
 Ba=B*exp(-0.5826*sigma*sqrt(k));   
    
    
b=0; % da attivare se si considerano dividendi per sol analitica
    
sum_payoff1 = 0;
sum_payoff2 = 0;

Z1=0;
Z2=0;


parfor j=1:N


Sa=[];
Sa=zeros(1,M+1); 
Sb=[];
Sb=zeros(1,M+1);

t=zeros(1,M+1);
Sa(1)=S0;
Pa=0;
% richiamo funzione per generare traiettorie
[Sa,Pa]=traiettoria1_barrier_option_1d(S0,T,sigma,r,M); 

     if min(Sa) > B
      Z1=max(K-Pa,0);
    else
      Z1=0;
     end  
    
    if min(Sa) > Ba
      Z2=max(K-Pa,0);
    else
      Z2=0;
    end
    
    sum_payoff1 = sum_payoff1 + Z1 ; 
    sum_payoff2 = sum_payoff2 + Z2 ; 
    
end

price2(W)=exp(-r*T)*(sum_payoff2/(N));  % prezzo  con correzione a B
price1(W)=exp(-r*T)*(sum_payoff1/(N));  % prezzo Monte Carlo standard
% % 
% % %%%%%%%%%%%%%%%%%%%%%
% %         %HULL       %
% % %%%%%%%%%%%%%%%%%%%%%
% % 
lambda1=(r+0.5*sigma^2)/(sigma^2);
y9=(log((Ba^2)/(K*S0)))/(sigma*sqrt(T))+lambda1*sigma*sqrt(T);

d91=(T*(r+0.5*sigma^2) + log(S0/K))/(sigma*sqrt(T));
d92= d91 -sigma*sqrt(T) ;

x91=(log(S0/Ba))/(sigma*sqrt(T))+(sqrt(T)*sigma*lambda1);
y91=(log(Ba/S0))/(sigma*sqrt(T))+(sqrt(T)*sigma*lambda1);

PUTDI= -S0*(normcdf(-x91)) +(K)*(exp(-r*T))*(normcdf(sigma*sqrt(T)-x91))...
    - (K)*(exp(-r*T))*((Ba/S0)^(-2+(lambda1)*2))*(normcdf(...
    -sigma*sqrt(T)+y9)-normcdf(-sigma*sqrt(T)+y91))...
    + (S0)*((Ba/S0)^(2*(lambda1)))*(normcdf(y9)-normcdf(y91));

PUT= K*(exp(-r*T))*(normcdf(-d92))-(S0)*(normcdf(-d91));

% prezzo fornito dalla soluzione analitica libro HULL
price9(W) = PUT - PUTDI;  

% errori relativi rispetto sol. analitica
error1(W) = abs(1-((price1(W))./price9(W)));
error2(W) = abs(1-((price2(W))./price9(W)));

end

toc