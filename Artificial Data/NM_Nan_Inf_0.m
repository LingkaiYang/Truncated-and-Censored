function ind = NM_Nan_Inf_0(ID, L)
% ID: IDs of instances
% L: Likelihood of instances

% indices that are nan, inf or 0

ind_nan = ID( isnan(L) );
ind_inf = ID( isinf(L) );
ind_0 = ID( L==0 );

ind = [ind_nan, ind_inf, ind_0];

