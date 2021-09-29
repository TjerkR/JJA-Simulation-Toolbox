function  [pNodes,pDir] = JJAsim_2D_network_pathproperties(junctionIslands)
%extracts ordered nodes and junction directions of a path defined by juncion nodes.
%
%input: 
%junctionIslands  Nj by 2  Islands corresponding to the Nj junctions in a single path.
%
%output: 
%pNodes           Nj by 1  list of nodes in the path, in order.
%pDir             Nj by 1  direction of each junction with respect to the path. (value +1 or -1).
%
Nj = size(junctionIslands,1);
pNodes = zeros(Nj,1);
pDir = zeros(Nj,1);
if junctionIslands(1,1) == junctionIslands(2,1) || junctionIslands(1,1) == junctionIslands(2,2)
    pNodes(1) = junctionIslands(1,1);
    pDir(1) = -1;
else
    pNodes(1) = junctionIslands(1,2);
    pDir(1) = 1;
end
for j = 2:Nj
   pNodes(j) = setdiff(junctionIslands(j,:),pNodes(j-1)); 
   if pNodes(j) == junctionIslands(j,1)
       pDir(j) = -1;
   else
       pDir(j) = 1;
   end
end
end