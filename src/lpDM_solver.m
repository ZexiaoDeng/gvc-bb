function [ x, fval, output ] = lpDM_solver( Theta, ATilde, bTilde, U ) 
% LPDM_SOLVER 用于求解 LP( D, M ) 的线性规划问题
%
%       max      sum( ti/thetai )
%       s.t.     ATilde*U*t <= bTilde
%                ti >= 0, i = 1, ..., n
%
%   输入:
%       Theta  :     M 的锥方向 gamma-扩张倍数
%                    Theta = [ theta1 ; theta2 ; ... ; thetan ; ] ;
%       ATilde :     系统矩阵
%                    ATilde = [ A ; -I ; ] ;
%       U      :     M 的锥方向
%       bTidle :     右手向量
%                    bTidle = [ b - A*x0 ; x0 ; ] ;
%
n     = size( Theta, 1 ) ;

c     = -1./Theta ;
Aineq = ATilde*U ;
bineq = bTilde   ;
Aeq   = [] ;
beq   = [] ;
lb    = zeros( n, 1 ) ;
ub    = [] ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'interior-point-legacy', ...
                    'display'  , 'none' ) ;
                
[ x       , ...
  fval    , ...
  exitflag, ...
  output  , ...
  lambda  ] = linprog( c    ,        ...
                       Aineq, bineq, ...
                       Aeq  , beq  , ...
                       lb   , ub   , ops ) ;

fval = -fval ;
output.exitflag = exitflag ;
output.lambda   = lambda   ;



end