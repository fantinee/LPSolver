function [flg, x, y, newB, t, s, r, obj] = Phase2(A,b,c,x0,B)
	[m,n] = size(A); 
	% Creating tableau
	T = A(:,B)\[b eye(m)];
	y = T(:,2:end)'*c(B);
	T = [T;[c'*x0,y']];
	% objective
	obj = c'*x0;

	% Start Simplex
	simplex = 1;
	ITER = 0;

	while (simplex == 1)
		% determine the next s and r values.
		y = T(end,2:end)';
		[zmin,s] = min(c-A'*y); 
		% check for convergence.
		if (abs(zmin) < 1e-14)
			% FOUND OPTIMAL
			simplex = 0;
			x       = zeros(n,1);
			x(B)    = T(1:end-1,1);
			obj     = c'*x;
			continue;
		end
		% here we define the t value
		t = T(1:end-1,2:end)*A(:,s);
		
		[flg,r] = Revisedgetr(n,s,B,T,t);

		if (flg == 1)
			info.run = 'Failure'; %% Failure case
			info.msg = 'Failure due degeneracy'; 
			simplex = 0;
			continue;
		end
		if (r < 1)
			info.run = 'Failure'; % Failure case
			info.msg = 'The problem has no lower bound'; 
			simplex = 0;
			continue;
		end
		x   = zeros(n,1);
		x(B)= T(1:end-1,1);
		ITER = ITER + 1;
		obj1 = c'*x;

        % update the revised simplex tableau.
        [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T);      
		if (flg == 1)
			info.run = 'Failure'; % Failure case
			info.msg = 'Failure due degeneracy'; 
			simplex = 0;
			continue;
		end
		newB   = B1;
		obj    = obj1;
		pause(1);
	end
end