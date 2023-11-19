%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of matlab implementation of the CMA-ES as discussed in
% http://www.scholarpedia.org/article/Evolution_Strategies   
% The code presented below should be regarded as a skeleton only 
% Note, the code presented is to be used under GNU General Public License
% Author: Hans-Georg Beyer   
% Email: Hans-Georg.Beyer_AT_fhv.at
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs recombination, i.e., calculating centroids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = CMArecomb(pop);
 r.w = 0; r.std = 0;
 for i=1:length(pop); r.w = r.w + pop{i}.w; r.std = r.std + pop{i}.std; end;
 r.w = r.w/length(pop); r.std = r.std/length(pop);
end
