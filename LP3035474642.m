function [data, info] = LP3035474642(A, b, c)
%
% setup LP
%

flg = 0; %% initialize flg. 

[m,n] = size(A); %% find dimensions, m rows n columns

if(m>n)
  data = 0;
  info.run = "Failure";
  info.msg = "m > n inconsistent dimensions";
  return
end

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

[flg, x0, y0, B0, t0, s0, r0, obj0, info.PhaseI.loop] = Phase2(A_phaseI,b,c_phaseI,x_phaseI,B_phaseI); %% Solution phase I

if ((flg==1)||(flg==2))
  data = 0;
  info.run = "Failure";
  info.msg = "Unable to compute a feasible solution in Phase I";
  return
end

data.PhaseI.obj = obj0;
data.PhaseI.x = x0;

if (obj0 > 0) %% case where is not possible to form a feasible solution with x variables
  info.run = "Failure";
  info.msg = "Sum over Z!=0, not able to get a feasible solution."; 
  info.case = 3;
end

%if (obj0 == 0) we move on to Phase II
% STARTING PHASE II SIMPLEX METHOD
if (obj0 == 0)
  x_phaseII = x0(1:n); %% New feasible solution for phase II
  info.case = 1;
  [flg, data.PhaseII.x, data.PhaseII.y, B, data.PhaseII.t, s, r, data.PhaseII.Primalobj, info.PhaseII.loop] = Phase2(A,b,c,x_phaseII,B0); %% Solution phase II
  if (flg == 1)
    data = 0;
    info.run = "Failure";
    info.msg = "Degeneracy case during Phase II";
    return
  end

  data.PhaseII.z = c - (A')*data.PhaseII.y;

  if (flg == 2) %% case 2
    info.case = 3;
    n = len(data.PhaseII.x);
    data.PhaseII.lambda = data.PhaseII.x.*[zeros(1,s-1) 1 zeros(1,n-s)]; %% x_j zero for index!=s
  end

  data.PhaseII.Dualobj = b'*(data.PhaseII.y); %% Dual Objective
  if (abs(data.PhaseII.Dualobj - data.PhaseII.Primalobj)< 1e-7) %% Checks that Primal and Dual are equal
    info.run = "Success";
  else
    info.run = "Success";           
    info.case = 3;
  end
end

clear A b c
end
