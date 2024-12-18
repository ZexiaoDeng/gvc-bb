function [ x, fval, output ] = cbb_solver( BCP )
%   CBB_SOLVER  基于锥形细分的分支定界求解器算法
%   CBB_SOLVER  用于求解基本凹规划问题 ( basic concave programming, BCP ), 形如
%
%       min     f( x )
%       s.t.    Aineq*x <= bineq
%               x >= 0
%
%   CBB_SOLVER  实现了如下算法:
%       (1) 锥形细分 ( conical subdivision )
%       (2) gamma-有效割平面 ( gamma-valid cuts )
%       (3) 分支定界框架 ( branch-and-bound )
%   输入:
%       BCP : 基本凹规划问题
%
%   输出:
%       x       :   最优解
%       fval    :   最优目标函数值
%       output  :   迭代 verbose 信息记录
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P148
%       Convex Analysis and Global Optimization, 2nd Edition, Hoang Tuy,
%       Spring, 2016, P180
%
% Copyright (C) 2021-2025 Zexiao Deng
%

path = './bt-1.3' ;
addpath( path ) ;

% ========================
% 阶段 1 ( phase 1 )
% ========================
% 计算当前最好解目标值 ( uppper_bound )
upper_bound = 1e1 ;
epsilon     = 1e-8 ;

D.B    = BCP.Aineq ;
D.b    = BCP.bineq ;
oracle = BCP.oracle ;

P   = polyh( D, 'h' ) ;           % 求出多面体 D 的 H-rep, V-rep, P-rep   
CH  = vrep( eval( P ) ) ;         % 获取多面体 D 的 V-rep
Adj = adj( eval( P ) ) ;          % 获取顶点对应的链表( 邻接表表示形式 )
V   = CH.V ;                      % 获取多面体 D 的顶点集

% fxBar = upper_bound
for idx = 1: size( V, 2 )
    if oracle( V( :, idx ) ) <= upper_bound
        % 局部最好解作为上界
        xBar  = V( :, idx ) ;
        fxBar = oracle( xBar ) ;
        break ;
    end
end
upper_bound = fxBar ;

% ========================
% 阶段 2 ( phase 2 )
% ========================
% 初始化根节点
for idx = 8: size( V, 2 )
    if oracle( V( :, idx ) ) >= upper_bound - epsilon
        % 任意给定初始点
        x0  = V( :, idx ) ;
        fx0 = oracle( x0 ) ;
        
        plot3( x0(1), x0(2), x0(3), '*', 'LineWidth', 2 )
        
        % 计算锥 M0 = { x | x = x0 + U0*t, t >= 0 }
        Z0  = V( :, Adj{idx} ) ;            % x0 的邻接点
        U0  = Z0  ;                         % n 条边为极方向
        for jdx = 1: size( Z0, 2 )
            U0( :, jdx ) = U0( :, jdx ) - x0 ;  
        end
        M0.V = x0 ;
        M0.D = U0 ;
        
        P = ( polyh( M0, 'v' ) ) ;   % 求出多胞体 P 的 H-rep, V-rep, P-rep
        opt.color = [ 1, 0, 0 ] ;
        plot( P, opt ) ;

        break ;
    end
end

% 计算局部上边界 gamma
upper_bound = fxBar ;

root.model = BCP ;     % 求解模型
root.x0    = x0 ;
root.M     = M0 ;

% 候选活跃节点 List
candidate_list(1)   = { root } ;

% 记录使用
LB = [] ;
UB = [] ;
k  = 0 ;

while ( ~isempty( candidate_list ) )
    
    % 选择节点, 深度优先搜索, 后进先出
    [ node, candidate_list ] = node_choice_DFS( candidate_list ) ;
    
    % =========================
    % 定界操作( bounding )
    % =========================
    M     = node.M  ;

    % 通过 gamma-epsilon 扩张计算 thetai
    gamma_eps = upper_bound - epsilon ;
    Theta     = zeros( size( M.D, 2 ), 1 ) ;
    Z_eps     = zeros( size( x0, 1 ), size( M.D, 2 ) ) ;
    for idx = 1: size( M.D, 2 )
        [ zi, thetai ] = gamma_extension( gamma_eps, M.V, M.D( :, idx ), oracle ) ;
        Z_eps( :, idx ) = zi ;
        Theta( idx, 1 ) = thetai ;
    end
    
    % 求解线性规划问题 LP( D, M )
    x0    = node.x0 ;
    Aineq = node.model.Aineq ;
    bineq = node.model.bineq ;
    
    n      = size( Theta, 1 ) ;
    ATilde = [ Aineq ; -eye( n ) ] ;            % LP( D, M ) 的系统矩阵
    bTilde = [ bineq - Aineq*x0 ; x0 ; ] ;      % LP( D, M ) 的右向量
    U      = M.D ;           % 锥 M 的极方向矩阵 U = [ u1, ..., un ]
    [ tBar_M , ...           % tBar_M = [ t1 ; ... ; tn ] ;
      mu_M   , ...           % mu( M )
      node.output ] = lpDM_solver( Theta , ...
                                   ATilde, ...
                                   bTilde, ...
                                   U      ) ;
    
    if ( node.output.exitflag <= 0 )
        % 剪枝操作( pruning )
        % 节点的问题是一个不可行问题
        lower_bound = inf ;
        fprintf( 'prune by infeasiblity!\n' ) ;
        continue ;
    end
    
    % 计算 f( D cap M ) 的下边界 beta_M = lower_bound
    if mu_M <= 1
        lower_bound = upper_bound - epsilon ;
    else
        fY = zeros( size( Theta, 1 ), 1 ) ;
        for idx = 1: size( Theta, 1 )
            thetai = Theta( idx ) ;
            ui     = U( :, idx ) ;
            yi     = x0 + mu_M*thetai*ui ;
            fyi    = oracle( yi ) ;
            fY( idx, 1 ) = fyi ;
        end
        lower_bound = min( fY ) ;
    end
    
    % =======================
    % 选择
    % =======================
    omega_M = x0 + U*tBar_M  ;
%     plot3( omega_M(1), omega_M(2), omega_M(3), 'r*', 'LineWidth', 2 )
    fomega_M = oracle( omega_M ) ;
    if fomega_M < upper_bound
        % 如果对于某个 M in mathcal_P
        % 并且好于当前最好解
        % 即f( omega( M ) ) < upper_bound
        % 执行下面任意一项
        
        % (2) 更新上界( 当前最好解 )操作
        for idx = 1: size( V, 2 )
            x_hat = V( :, idx ) ;
            fx_hat = oracle( x_hat ) ;
            if fx_hat <= fomega_M
                xBar        = x_hat  ;
                upper_bound = fx_hat ;
                
                break ;
            end
        end
    end
    
    if ( lower_bound >= upper_bound - epsilon )
        % 剪枝操作( pruning )
        % 节点的下边界大于当前最好解( 全局上边界 ), 剪枝
        fprintf( 'prune by bound!\n' ) ;
        continue ;
    end
    
    % ====================================
    % 分支操作( branching )
    % ====================================
    % (1) 单纯 omega-细分
    % 以 x0 为起点, 通过 omega( Mk )点, 进行 omega-细分
    sub_M = conical_subdivision( x0, omega_M - x0, tBar_M, U ) ;
    k = k + 1 ;
    
    figure
    for idx = 1: length( sub_M )
        len = length( candidate_list ) ;
        
        child_node   = node ;
        child_node.M = sub_M{ idx } ;
        
        P   = eval( polyh( sub_M{idx}, 'v' ) ) ;
        opt.color = [ 0, 0, idx+1 ]*0.1 ;
        plot( P, opt ) ;
        axis equal ;
        axis( [ 0, 1.2, 0, 1.5, 0, 2 ] ) ;
        grid on
        hold on
        view( [ 35, 27 ] ) ;
        view( [ 77, 30 ] ) ;
        
        candidate_list( len + 1 ) = { child_node } ;
    end
    
    LB = [ LB ; lower_bound ] ;
    UB = [ UB ; upper_bound ] ;
    
end

LB = [ LB ; lower_bound ] ;
UB = [ UB ; upper_bound ] ;

if abs( lower_bound - upper_bound ) <= 1e-7
    fprintf( 'Optimal solution found!\n' ) ;
    output.message = 'Optimal solution found!' ;
end
            
x        = xBar ;
fval     = upper_bound ;
output   = node.output ;
output.k = k ;

figure
plot( LB, '-^' ), hold on ;
plot( UB, '-s' ) ;



end


function [ node, candidate_list ] = node_choice_DFS( candidate_list )
    % 选择节点策略, 深度优先搜索, 后进先出
    len                   = length( candidate_list ) ;
    node                  = candidate_list{ len } ;
    candidate_list( len ) = []                    ;
    return ;
end






