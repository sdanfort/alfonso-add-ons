function L = pdwtsFFT(n,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdwtsFFT".
%
% L = pdwtsFFT(n)
% L = pdwtsFFT(n,xyrange)
%
% Compute the cubature weights so that, if Pad is the
% matrix of Padua points computed through the call
% Pad = pdpts(n) or Pad = pdpts(n,xyrange), then the cubature 
% of the function f is given by L'*f(Pad(:,1),Pad(:,2)).
%-------------------------------------------------------------------------------
% INPUT.    
%
% n       : interpolation degree
% xyrange : an optional vector [a,b,c,d] defining the rectangle 
%           [a,b] x [c,d]. By default, xyrange = [-1,1,-1,1]
%
% OUTPUT.   
%
% L       : cubature weights
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% no external function is used.
%-------------------------------------------------------------------------------

% Copyright (C) 2008-2009 
% Marco Caliari, Stefano De Marchi, Alvise Sommariva, Marco Vianello.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%
% Author:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <stefano.demarchi@univr.i>   
%          Alvise Sommariva  <alvise@euler.math.unipd.it>
%          Marco Vianello    <marcov@euler.math.unipd.it>   
%
% Date: March 10, 2009.
%-------------------------------------------------------------------------------

if (nargin == 1)
  xyrange = [-1,1,-1,1];
else
  xyrange = varargin{1};
end
if (n == 0)
% degree 0
  L = (xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3));
else
% even,even moments matrix
  k = [2:2:n]';
  mom = zeros(n+1,1);
  mom(1) = 2*sqrt(2);
  mom(3:2:n+1) = 4./(1-k.^2);
  [M1,M2] = meshgrid(mom);
  M = M1.*M2;
  M(1,n+1) = M(1,n+1)/2;
  M(1,:) = M(1,:)/sqrt(2);
  M(:,1) = M(:,1)/sqrt(2);
  Mmom = fliplr(triu(fliplr(M)));
  Mmomhat = real(fft(Mmom,2*(n+1)));
  Mmomhat = Mmomhat(1:n+2,:);
  Mmomhathat = real(fft(Mmomhat,2*n,2));
  TMT = Mmomhathat(:,1:n+1);
% interpolation weights
  W = 2*ones(n+2,n+1)/(n*(n+1));
  W(1,:) = W(1,:)/2;
  W(:,1) = W(:,1)/2;
  W(n+2,:) = W(n+2,:)/2;
  W(:,n+1) = W(:,n+1)/2;
  L = W.*TMT;
  [M1,M2] = meshgrid([0:n],[0:n+1]);
  findM = find(mod(M1+M2,2));
% cubature weights
  L = L(findM);
  L = L*(xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))/4;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdwtsFFT
%
% at the Octave prompt
%
%!test
%! disp('Degree 0 weight')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(0,xyrange);
%! expected = 1;
%! assert(L/expected,1,10*eps);
%!test
%! disp('Degree 1 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(1,xyrange);
%! expected = [0.5;0.25;0.25];
%! assert(L,expected,10*eps);
%!test
%! disp('Degree 2 weights')
%! xyrange = [0,1,0,1];
%! L = pdwtsFFT(2,xyrange);
%! expected = [1/6;0;1/9;5/9;1/6;0];
%! assert(L,expected,10*eps);
