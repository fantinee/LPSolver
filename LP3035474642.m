function [data, info] = LP3035474642(A, b,c)
%
% setup LP
%

flg = 0; %% initialize flg. 

[m,n] = size(A); %% find dimensions, m rows n columns

%
% STARTING SIMPLEX METHOD:
%% STARTING PHASE I SIMPLEX METHOD
% Computes basic feasible solution

A_phaseI = [A eye(m)]; %% A extended with I(m)
e_phaseI = ones(1,m); %% e = (1 .. 1)
c_phaseI = [zeros(1,n) e_phaseI]; %% new obj for phase I 
x_phaseI = [zeros(n,1); b]; %% x = [0 ... 0 b]'
B_phaseI = find(x_phaseI~=0); %% index from phase first Basis

%% STARTING PHASE II SIMPLEX METHOD

obj = c'*x;
disp(['Initial Objective = ', num2str(obj)]);
%disp('Displaying Initial solution x, c-A^T*y and their componentwise product');
%disp([x c-A'*y x.*(c-A'*y)]);
simplex = 1;
ITER = 0;
%pause(2);
while (simplex == 1)
%
% determine the next s and r values.
%
   y        = T(end,2:end)';
   [zmin,s] = min(c-A'*y); 
%
% check for convergence.
%
   if (abs(zmin) < 1e-14)
       disp('Simplex Method has converged');
       simplex = 0;
%       disp('Displaying Optimal Basis');
%       disp(B');
       x   = zeros(n,1);
       x(B) = T(1:end-1,1);
       obj  = c'*x;
       disp(['Optimal Objective = ', num2str(obj),' after ', num2str(ITER), ' iterations']);
       disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
       disp([x c-A'*y x.*(c-A'*y)]);
       continue;
   end

   t        = T(1:end-1,2:end)*A(:,s);
   [flg,r] = Revisedgetr(n,s,B,T,t);
   if (flg == 1)
       info.run = 'Failure'; % Failure case
       info.msg = 'Failure due degeneracy'; 
       %disp('LP is degenerate');
       simplex = 0;
       continue;
   end
   if (r < 1)
       disp('LP has no lower bound');
       simplex = 0;
       continue;
   end
   x   = zeros(n,1);
   x(B)= T(1:end-1,1);
   ITER = ITER + 1;
   f = ['Iteration ', num2str(ITER), ' Obj ', num2str(c'*x), '. Smallest component in c-A^T*y: ', ... 
         num2str(zmin), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
%   disp(f);
   obj1 = c'*x;
%
% update the revised simplex tableau.
%
   [T,B1,flg]=RevisedSimplexTableau(B,r,s,t,zmin,T);      
   if (flg == 1)
       %disp('LP is degenerate');
       info.run = 'Failure'; % Failure case
       info.msg = 'Failure due degeneracy'; 
       simplex = 0;
       continue;
   end
   B   = B1;
   obj = obj1;
%   disp('Current Basis is');
%   disp(B');
%   pause(1);
end
clear B1 f obj1 t zmin
end