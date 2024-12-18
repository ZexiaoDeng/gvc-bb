function U = canonical_cone( x0, V )

% 计算
U  = V  ;
for idx = 1: size( V, 2 )
    U( :, idx ) = U( :, idx ) - x0 ;
end

M.V = x0 ;
M.D = U  ;

end


