% ============================================
% 全局优化引论
% R. Horst, P.M. Pardalos, N.V. Thoai 著
% 黄红选 译
% 梁治安 校
% P178
% =============================================
%   min f( x )
%       x in D
%
%   D = { x: A*x <= b, x >= 0 } in R+^3
%
%   f(x) = -( abs( x1 + x2/2 + 2*x3/3) )^( 3/2 ) - x1^2
%   A = [  1,  1,  1   ; ...
%          1,  1, -1/4 ; ...
%         -2, -2,  1   ; ...
%          0,  0,  1   ; ] ;
%   b = [ 2 ; 1 ; 1 ; 3 ] ;
% 

clc ;
clear ;
close all ;

path = './bt-1.3' ;
addpath( path ) ;

BCP.oracle = @oracle ;
BCP.Aineq = [  1,  1,  1   ; ...
               1,  1, -1/4 ; ...
              -2, -2,  1   ; ...
               0,  0,  1   ; ...
              -1,  0,  0   ; ...
               0, -1,  0   ; ...
               0,  0, -1   ; ] ;
BCP.bineq = [ 2 ; ...
              1 ; ...
              1 ; ...
              3 ; ...
              0 ; ...
              0 ; ...
              0 ; ] ;

D.B = BCP.Aineq ;
D.b = BCP.bineq ;

P   = eval( polyh( D, 'h' ) ) ;   % 求出多胞体 P 的 H-rep, V-rep, P-rep   
CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
Adj = adj( P ) ;                  % 获取顶点对应的链表( 邻接表表示形式 )

plot( P ) ;
axis equal ;
axis( [ 0, 1.2, 0, 1.5, 0, 2 ] ) ;
grid on
hold on
view( [ 35, 27 ] ) ;
view( [ 77, 30 ] ) ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f\t%8.4f', ...
        CH.V( 1, i ), ...
        CH.V( 2, i ), ...
        CH.V( 3, i ) ) ; % 顶点
    fprintf( '%8d\t%8d\t%8d', ...
        Adj{i}(1), ...
        Adj{i}(2), ...
        Adj{i}(3) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end
          
xBar = CH.V( :, 4 ) ;

[ x, fval, output ] = cbb_solver( BCP )

BCP.oracle( xBar )

function f = oracle( x )
% 目标函数

    f = -( abs( x(1) + x(2)/2 + 2*x(3)/3) )^( 3/2 ) - x(1)^2 ;
end

