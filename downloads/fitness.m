%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of matlab implementation of the CMA-ES as discussed in
% http://www.scholarpedia.org/article/Evolution_Strategies   
% The code presented below should be regarded as a skeleton only 
% Note, the code presented is to be used under GNU General Public License
% Author: Hans-Georg Beyer   
% Email: Hans-Georg.Beyer_AT_fhv.at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective function to be minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = fitness(x); 
% this is Schwefel's ellipsoid test function, a moderately conditioned 
% ellipsoid with a dominating isolated eigenvalue
 out = 0; 
 for i = 1:length(x); out = out + sum(x(1:i))^2; end;
end
