function [ y, theta ] = gamma_extension( gamma, x, d, oracle )
%   GAMMA_EXTENSION  计算 gamma-扩张( gamma-extension )
%   GAMMA_EXTENSION : 以 x 为起点, d 为方向的射线 
%                   rho( x, d ) = { y: y = x + t*d, t >= 0 }
%                   上找到最长线段 [ x, x + theta*d ]
%                   使得对于所有的 y in [ x, x + theta*d ]
%                   满足 f( y ) >= gamma
%
%       y = x + theta*d
%
%   输入:
%       gamma  : 极小化过程中的某一阶段所得到的最好目标函数值
%       x      : 扩张起始点
%       d      : 扩张方向
%       oracle : 目标函数方程
%   输出:
%       y      : gamma-扩张
%       theta  : 射线方向延伸倍数
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P131
%

% 若 L( f; gamma ) 无界时, 置 theta1 = 1e9
theta1  = 1e9 ;
% options = optimoptions( 'fmincon'  , ...
%                         'Display'  ,'off', ...
%                         'Algorithm','interior-point' ) ;
% options = optimoptions( 'fmincon'  , ...
%                         'Display'  ,'off', ...
%                         'Algorithm','sqp' ) ;

% 寻找 theta 可以用托架-二分技术( bracketing-bisection technique )
t0      = 0 ;       % 初值
options = optimoptions( 'fmincon'  , ...
                        'Display'  ,'off', ...
                        'Algorithm','sqp-legacy' ) ;
theta2  = fmincon( @fun, t0, ...
                   [] , [], ...
                   [] , [], ...
                   [] , [], ...
                   @nonlcon, options ) ;

theta = min( theta1, theta2 ) ;

y = x + theta*d ;



function f = fun( t )
    f = -t ;
end

function [ c, ceq ] = nonlcon( t )
    
    c = - feval( oracle, x + t*d ) + gamma ;    % f( x + t*d ) >= gamma
    ceq = [] ;
    
end


end
