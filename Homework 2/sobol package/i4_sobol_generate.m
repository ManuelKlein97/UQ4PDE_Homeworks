function r = i4_sobol_generate ( m, n, skip )

%*****************************************************************************80
%
%% i4_sobol_generate() generates a Sobol dataset.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 December 2009
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer M, the spatial dimension.
%
%    integer N, the number of points to generate.
%
%    integer SKIP, the number of initial points to skip.
%
%  Output:
%
%    real R(M,N), the points.
%
  for j = 1 : n
    key = skip + j - 1;
    [ r(1:m,j), key ]  = i4_sobol ( m, key );
  end

  return
end
