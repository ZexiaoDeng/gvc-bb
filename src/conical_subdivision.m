function M = conical_subdivision( x0, w, lambda, U )
%   CONICAL_SUBDIVISION  计算锥形剖分
%
%       M0 = Union M(w) i in I
%       w in S
%       w = sum( lambda_i*di ), lambda_i >= 0, i = 1, ..., n
%       sum( lambda_i ) = 1, i = 1, ..., n
%
%   输入:
%       w      : 单纯形 S 的内方向
%       lambda : 一个有限指标集, int( Pi )非空
%       U      : 原多面锥的极方向
%                U = [ u1, u2, ..., un ]
%
%   输出:
%       M      : 子多面锥 cell
%                M = [ u1, u2, ..., ui-1, w, ui+1, ..., un ]
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P148
%

n       = size( U, 2 ) ;
flag    = lambda > 1e-8 ;
M       = cell( sum( flag ), 1 ) ;          % 锥集 M
counter = 0 ;                               % 累加器
Mi.V    = x0 ;                              % 锥 Mi 顶点
for idx = 1: n
    if flag( idx )
        Mi.D = U ;
        Mi.D( :, idx ) = w ;                % 锥 Mi 的 n 条边
        M( counter + 1, 1 ) = { Mi } ;
        counter = counter + 1 ;
    else
        continue ;
    end
end

return ;





end



