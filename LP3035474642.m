function [data, info] = LP3035474642(A, b, c)
%
% setup LP
%

flg = 0; %% initialize flg. 

[m,n] = size(A); %% find dimensions, m rows n columns

%
% STARTING SIMPLEX METHOD:
% STARTING PHASE I SIMPLEX METHOD
% Computes basic feasible solution

A_phaseI = [A eye(m)]; %% A extended with I(m)
e_phaseI = ones(1,m); %% e = (1 .. 1)
c_phaseI = [zeros(1,n) e_phaseI]'; %% new obj for phase I 
x_phaseI = [zeros(n,1); b]; %% x = [0 ... 0 b]'
B_phaseI = find(x_phaseI~=0); %% index from phase first Basis

% getting feasible solution 

[flg, x0, y0, B0, t0, s0, r0, obj0] = Phase2(A_phaseI,b,c_phaseI,x_phaseI,B_phaseI); %% Solution phase I

if (flg==1)
  info.run = "Failure";
  info.msg = "Unable to compute a feasible solution in Phase I";
end

data.PhaseI.obj = obj0;
data.PhaseI.x = x0;

if (obj0 > 0) %% case where is not possible to form a feasible solution with x variables
  info.run = "Failure";
  info.msg = "Sum over Z!=0, not able to get a feasible solution."; 
end

%if (obj0 == 0) we move on to Phase II
% STARTING PHASE II SIMPLEX METHOD

x_phaseII = x0(1:n);

info.case = 1;
[flg, x, y, B, t, s, r, obj] = Phase2(A,b,c,x_phaseII,B0); %% Solution phase II
data.PhaseII.PrimalObj = obj;
data.PhaseII.x = x;
data.PhaseII.y = y;
data.PhaseII.z =c - (A')*y;

end