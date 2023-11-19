%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a simple matlab implementation of the CMA-ES as discussed in 
% http://www.scholarpedia.org/article/Evolution_Strategies                  
% The code presented below should be regarded as a skeleton only            
% Note, the code presented is to be used under GNU General Public License   
% Author: Hans-Georg Beyer                                                  
% Email: Hans-Georg.Beyer_AT_fhv.at                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of the specific strategy and problem size:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu = 3;                 % number of parents
lambda = 12;            % number of offspring
yInit = ones(30,1);     % initial parent vector 
sigmaInit = 1;          % initial global mutation strength sigma 
sigmaMin = 1e-10;       % CMA-ES stops when sigma is smaller than sigmaMin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(yInit);      % determine search space dimensionality n   
tau = sqrt(n); tau_c = n^2; tau_sigma = sqrt(n); % time constants
Cov = eye(n);           % initial covariance matrix 
sigma = sigmaInit;      % initial sigma
s = zeros(n,1);         % set cumulation vector to zero
s_sigma = zeros(n,1);   % set cumulation vector to zero
% initializing individual population:
Individual.y = yInit; 
Individual.w = 0;
Individual.std = 0;
Individual.F = fitness(Individual.y);
for i=1:mu; ParentPop{i} = Individual; end;
yParent = yInit;        % initial centroid parent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evolution loop of the CMA-ES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)
 SqrtCov = chol(Cov)';                    % "square root" of covariance matrix
 for l = 1:lambda;                        % generate lambda offspring
  OffspringIndividual.std = randn(n,1);   % line (L1a)
  OffspringIndividual.w = sigma*(SqrtCov*OffspringIndividual.std); % line (L1a)
  OffspringIndividual.y = yParent + OffspringIndividual.w;         % line (L1b)
  OffspringIndividual.F = fitness(OffspringIndividual.y);  % determine fitness (L1c)
  OffspringPop{l} = OffspringIndividual;                   % offspring complete
 end;
 ParentPop = SortPop(OffspringPop, mu);   % sort population and take mu best
 disp(ParentPop{1}.F);                    % display best fitness in population
 Recombinant = CMArecomb(ParentPop);      % (L2) perform recombination 
 yParent = yParent + Recombinant.w;       % (L2) calculate new centroid parent
 s = (1-1/tau)*s + sqrt(mu/tau*(2-1/tau))*Recombinant.w/sigma;   % line (L3)
 Cov = (1-1/tau_c)*Cov + (s/tau_c)*s';                           % line (L4)
 Cov = (Cov + Cov')/2;                    % enforce symmetry of cov matrix
 s_sigma = (1-1/tau_sigma)*s_sigma + sqrt(mu/tau_sigma*(2-1/tau_sigma))*Recombinant.std; % line (L5)
 sigma = sigma*exp((s_sigma'*s_sigma - n)/(2*n*sqrt(n)));        % line (L6)
 if (sigma < sigmaMin ) break; end;       % termination condition
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yParent contains an approximation of the optimizer of the function "fitness(y)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

