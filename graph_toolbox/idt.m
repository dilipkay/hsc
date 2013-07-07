function [L Li] = idt(vertex, face, do_idt)

% Perform Intrinsic Delaunay Triangulation of a 
% mesh by edge flipping using the algorithm and formulas
% given in Fisher et. al's paper: An Algorithm for the Contruction
% of Intrinsic Delaunay Triangulations with Applications
% to Digital Geometry Processing

% find all edges in the mesh
% edge(1,i) and edge(2, i) are the 2 vertices
% sharing edge i; v2e(a, b) gives the edge index 
% of (a, b)
[edge v2e] = compute_edges(face);

% For edge (i, j) e2f(i, j) and e2f(j, i)
% are the indices of the faces adjacent to edge (i, j)
v2f = compute_edge_face_ring(face);

% ne = number of edges; nv = number of vertices
ne = size(edge, 2);
nv = size(v2f, 1);

% faces listed in a particular order
f1list = full(v2f(edge(1, :) + (edge(2, :) - 1)*nv));
% "reverse" ordering
f2list = full(v2f(edge(2, :) + (edge(1, :) - 1)*nv));

% prune to remove any boundary edges or "weird" edges
I1 = find(f1list > 0); 
I2 = find(f2list > 0); 
		        
% find vertices corresponding to one ordering of the faces
% f1v and f2v are of size 3 x ne
f1v = zeros(3, ne);
f2v = zeros(3, ne);
f1v(:, I1) = face(:, f1list(I1));
% and the othe ordering
f2v(:, I2) = face(:, f2list(I2));

%edge = edge(:, I);
%n = length(I);

% find the 3'rd vertex on this face 
edge(3, I2) = find_other(f2v(:, I2), edge(1, I2), edge(2, I2));
% and the 4th vertex
edge(4, I1) = find_other(f1v(:, I1), edge(1, I1), edge(2, I1));
% and angles 
[ignore edge(5, :)] = check_delaunay(vertex, edge(1:4, :));

idx = find(edge(1, :) > 0);
edge = edge(:, idx);
ne = length(idx);

% finally, check Delaunay criterion for each edge
% ic == 0 implis Delaunay is not satisfied so 
% edge must be flipped; for each edge (i, j)
% lengths{1:5} contains the lengths of (i, j)
% and the 4 other edges of the quadrilateral

max_edges = 2*ne;
stack = zeros(max_edges, 1);
mark = ones(max_edges, 1);
stack_top = 1;
stack_bottom = ne;
stack(1:stack_bottom) = [1:ne];

% to avoid repeated allocation, pre-allocate extra space for new edges
edge = [edge zeros(5, ne)];
tot_edge = 0;
num_nd = 0;

% save off the original data
face_orig = face;
v2f_orig = v2f;
v2e_orig = v2e;
ne_orig = ne;
edge_orig = edge;
f1v_orig = f1v;
f2v_orig = f2v;

if do_idt
	while stack_top < stack_bottom
		% pop an edge from the stack
		edge_idx = stack(stack_top);
		stack_top = stack_top + 1;
		% unmark it
		mark(edge_idx) = 0;
		% check for Delaunay
		a = edge(1, edge_idx);
		b = edge(2, edge_idx);
		c = edge(3, edge_idx);
		d = edge(4, edge_idx);
		%c = find_other(f2v(:, edge_idx), a, b);
		%d = find_other(f1v(:, edge_idx), a, b);

		if (a == 0)
			continue;
		end;
		if (a == b | a == c | a == d | b == c | b == d | c == d)
			fprintf('Error - possibly non-planar mesh');
			keyboard;
		end;

		ss = stack_bottom - stack_top;
		%fprintf('ss %d\n', ss);

		[is_delaunay cots] = check_delaunay(vertex, [a b c d]');
		edge(5, edge_idx) = cots;
		tot_edge = tot_edge + 1;
		if ~is_delaunay
			num_nd = num_nd + 1;

			% find the index of the 2 faces in the original triangulation
			%flocal1 = v2f(a, b); 
			%flocal2 = v2f(b, a);
			%if (flocal1 == 412 || flocal2 == 412)
				%keyboard;
			%end;

			% update with the new vertices (flipped edge)
			%face(:, flocal1) = [c; d; a];
			%face(:, flocal2) = [d; c; b];

			% update v2f - vertex to face mapping
			% and f1v/f2v - edge to vertex mapping
			%if (v2f(a, d) == flocal2)
				%v2f(a, d) = flocal1;
				%f1v(:, v2e(a, d)) = replace_other(f1v(:, v2e(a, d)), a, d, b, c);
			%else
				%v2f(d, a) = flocal1;
				%f1v(:, v2e(d, a)) = replace_other(f1v(:, v2e(d, a)), a, d, b, c);
			%end;

			%if (v2f(c, b) == flocal1)
				%v2f(c, b) = flocal2;
				%f2v(:, v2e(c, b)) = replace_other(f2v(:, v2e(c, b)), c, b, a, d);
			%else
				%v2f(b, c) = flocal2;
				%f2v(:, v2e(b, c)) = replace_other(f2v(:, v2e(b, c)), c, b, a, d);
			%end;

			% add the new edge to the list of edges
			ne = ne + 1;
			edge(1, ne) = c;
			edge(2, ne) = d;
			edge(3, ne) = a;
			edge(4, ne) = b;
			edge(5, ne) = 0;

			% the id's of the 4 edges
			%      a
			%    / | \
			%   c  |  d
			%   \  |  /
			%      b 
	  
			edge_ids = [v2e(a, c) v2e(c, b) v2e(b, d) v2e(d, a)];

			% adjust each of their entries i.e. replace neighbors a/b with c/d appropriately
			% a-c: replace b with d
			if (edge(3, edge_ids(1)) == b)
				edge(3, edge_ids(1)) = d;
			elseif (edge(4, edge_ids(1)) == b)
				edge(4, edge_ids(1)) = d;
			end;
			% c-b: replace a with d
			if (edge(3, edge_ids(2)) == a)
				edge(3, edge_ids(2)) = d;
			elseif (edge(4, edge_ids(2)) == a)
				edge(4, edge_ids(2)) = d;
			end;
			% b-d: replace a with c
			if (edge(3, edge_ids(3)) == a)
				edge(3, edge_ids(3)) = c;
			elseif (edge(4, edge_ids(3)) == a)
				edge(4, edge_ids(3)) = c;
			end;
			% d-a: replace b with c
			if (edge(3, edge_ids(4)) == b)
				edge(3, edge_ids(4)) = c;
			elseif (edge(4, edge_ids(4)) == b)
				edge(4, edge_ids(4)) = c;
			end;

			% invalidate the dropped edge
			edge(:, edge_idx) = 0;
			
			% add the 2 faces to the faces list
			%f1v(:, ne) = face(:, flocal1);
			%f2v(:, ne) = face(:, flocal2);
			
			% remove the dropped edge
			%f1v(:, edge_idx) = 0;
			%f2v(:, edge_idx) = 0;
			
			v2e(c, d) = ne;
			v2e(d, c) = ne;

			for j = 1:length(edge_ids)
				if ~mark(edge_ids(j))
					mark(edge_ids(j)) = 1;
					stack_bottom = stack_bottom + 1;
					stack(stack_bottom) = edge_ids(j);
				end
			end;
		end
	end
end;
L = sparse(nv, nv);

for i = 1:ne_orig
	a = edge_orig(1, i);
	b = edge_orig(2, i);
	if (a == 0 | b == 0)
		continue;
	end;
	L(a, b) = -edge_orig(5, i);
	L(b, a) = -edge_orig(5, i);
end;
L = L - spdiags(sum(L, 2), 0, size(L, 1), size(L, 2));

Li = sparse(nv, nv);
if (do_idt)
	for i = 1:ne
		a = edge(1, i);
		b = edge(2, i);
		c = edge(3, i);
		d = edge(4, i);
		% skip flipped edges
		if (a == 0)
			continue;
		end;
		Li(a, b) = -edge(5, i);
		Li(b, a) = -edge(5, i);
	end;
	Li = Li - spdiags(sum(Li, 2), 0, size(Li, 1), size(Li, 2));
end;

if (0)
	[Vg Dg] = eig(full(L), full(Li));
	[V D] = eig(full(L));
	D = diag(D);
	for i = 1:size(L, 1)
		Di(i,1) = V(:, i)'*Li*V(:, i);	
	end
	%figure; plot(D(2:end)./Di(2:end), 'r');
	ratio = D(2:end)./Di(2:end);
	Dg = sort(diag(Dg));
	ratiog = Dg(end)/Dg(1);
	fprintf('Energy CN %f; generalized CN %f\n', max(ratio)/min(ratio), ratiog);
end;


% given a 3 x N vector f, 1 x N vectors i and j, 
% find the complemetary 1 x N vector k i.e.
% k = f \ (i U j)
function k = find_other( f, i,j )
	n = size(f, 2);
	idx1 = zeros(3, n);
	for d = 1:3
		idx = find(f(d, :)  == i);
		idx1(d, idx) = 1;
		idx = find(f(d, :)  == j);
		idx1(d, idx) = 1;
	end;

	k = zeros(1, n);
	for d = 1:3
		idx = find(idx1(d, :) == 0);
		k(idx) = f(d, idx);
	end;

% replace entry c in vertices array with d
function vertices = replace_other(vertices, a, b, c, d)
	idx = find(vertices == c);
	vertices(idx) = d;
