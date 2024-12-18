% ===========================================
% 最优化计算方法及其MATLAB程序实现
% 马昌凤, 柯艺芬, 谢亚君 著, 
% 国防工业出版社, 
% P225
% =========================================
%   min     f(x) = -3*x1^2 - x2^2 - 2*x3^2
%   s.t.    x1^2 + x2^2 + x3^2 = 3
%           x1 - x2 <= 0
%           x1 >= 0
%
function test

clc ;
clear ;

fprintf( 'test 01====================================\n' ) ;

Aineq = [ 1, -1, 0 ; ] ;
bineq = 0 ;
Aeq   = [] ;
beq   = [] ;
lb    = [ 0 ; -inf ; -inf ; ] ;
ub    = [ inf ; inf ; inf ; ] ;

% x0  = [ 0 ; 0 ; sqrt(3) ; ] ;
% x0  = [ sqrt( 1.5 ) ; sqrt( 1.5 ) ; 0 ] ;

x0    = [ 1.5 ; 1.5 ; 0 ; ] ;
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','sqp-legacy' ) ;
[ x       , ...
  fval    , ...
  exitflag, ...
  output  , ...
  lambda  , ...
  grad    , ...
  hessian ] = fmincon( @oracle ,        ...
                       x0     ,        ...
                       Aineq  , bineq, ...
                       Aeq    , beq  , ...
                       lb     , ub   , ...
                       @nonlcon,        ...
                       options       )

% 目标函数
function f = oracle( x )
    % 目标函数
    f = -3*x(1)^2 - x(2)^2 - 2*x(3)^2 ;
    
end

% 非线性约束函数
function [ c, ceq ] = nonlcon( x )
    % 非线性等式和不等式约束函数
    c = [] ;        % c(x) <= 0
    
    ceq = x(1)^2 + x(2)^2 + x(3)^2 - 3 ;    % ceq(x) = 0

end


% ============================================
% 全局优化引论
% R. Horst, P.M. Pardalos, N.V. Thoai 著
% 黄红选 译
% 梁治安 校
% P201
% =============================================
%   min f( x1, x2 ) = -( x1 - 20 )^2 - ( x2 - 10 )^2
%   s.t. x2 - x1/2 <= 10
%        ( x2 - 10 )^2 + ( x1 )^2 <= 500
%        x1 >= 0, x2 >= 0
%


clear ;
fprintf( 'test 02====================================\n' ) ;

Aineq = [ -1/2, 1 ; ] ;
bineq =  10 ;
Aeq   = [] ;
beq   = [] ;
lb    = zeros( 2, 1 ) ;
ub    = inf*ones( 2, 1 ) ;

x0    = [ 5 ; 6 ; ] ;

options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'iter', ...
                        'Algorithm','sqp-legacy' ) ;
[ x       , ...
  fval    , ...
  exitflag, ...
  output   ] = fmincon( @oracle2 ,        ...
                        x0       ,        ...
                        Aineq    , bineq, ...
                        Aeq      , beq  , ...
                        lb       , ub   , ...
                        @nonlcon2,        ...
                        options       )

% 目标函数
function f = oracle2( x )
    % 目标函数
    f = -( x(1) - 20 )^2 - ( x(2) - 10 )^2 ;
    
end

% 非线性约束函数
function [ c, ceq ] = nonlcon2( x )
    % 非线性等式和不等式约束函数
    c   = x(1)^2 + ( x(2) - 10 )^2  - 500 ;        % c(x) <= 0
    ceq = [] ;    % ceq(x) = 0

end


end





