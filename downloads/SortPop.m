%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of matlab implementation of the CMA-ES as discussed in
% http://www.scholarpedia.org/article/Evolution_Strategies   
% The code presented below should be regarded as a skeleton only 
% Note, the code presented is to be used under GNU General Public License
% Author: Hans-Georg Beyer   
% Email: Hans-Georg.Beyer_AT_fhv.at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sorts population w.r.t. the individuals fitness in ascending order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sorted_pop = SortPop(pop, mu);
 for i=1:length(pop); fitnesses(i) = pop{i}.F; end;
 [sorted_fitnesses, index] = sort(fitnesses);
 for i=1:mu; sorted_pop{i} = pop{index(i)}; end
end
