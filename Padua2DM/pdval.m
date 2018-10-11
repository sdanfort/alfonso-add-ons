function LnfX = pdval(C0f,xyrange,X,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdval".
%
% LnfX = pdval(C0f,xyrange,X)
% LnfX = pdval(C0f,xyrange,X1,X2)
%
% Evaluate the interpolation polynomial defined through its coefficient
% matrix C0f at the target points X(:,1),X(:,2) or at the meshgrid(X1,X2)
%-------------------------------------------------------------------------------
% INPUT.    
%
% C0f     : coefficient matrix
% xyrange : a vector [a,b,c,d] defining the rectangle [a,b] x [c,d]
% X       : a matrix with the abscissas of the target points in the
%           first column and the ordinates in the second one
% X1,X2   : two vectors defining the (mesh)grid X1 x X2 of the 
%           target points
%
% OUTPUT.   
%
% LnfX    : evaluation of the interpolation polynomial at the target points
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
% Authors:  
%          Marco Caliari     <marco.caliari@univr.it>
%          Stefano De Marchi <stefano.demarchi@univr.i>   
%          Alvise Sommariva  <alvise@euler.math.unipd.it>
%          Marco Vianello    <marcov@euler.math.unipd.it>   
%
% Date: March 10, 2009.
%-------------------------------------------------------------------------------

n = length(C0f)-1;
if (nargin == 3)
% scattered target points X(:,1),X(:,2)
% scale the target points exactly in the square [-1,1] x [-1,1]
  X1 = min(max(2*(X(:,1)'-xyrange(1))/(xyrange(2)-xyrange(1))-1,-1),1);
  X2 = min(max(2*(X(:,2)'-xyrange(3))/(xyrange(4)-xyrange(3))-1,-1),1);
  TX1 = cos([0:n]'*acos(X1));
  TX2 = cos([0:n]'*acos(X2));
  TX1(2:n+1,:) = TX1(2:n+1,:)*sqrt(2);
  TX2(2:n+1,:) = TX2(2:n+1,:)*sqrt(2);
  LnfX = sum((TX1'*C0f).*TX2',2);
else
% meshgrid points (X1,X2)
% scale the target points exactly in the square [-1,1] x [-1,1]
  Y = varargin{1};
  X1 = min(max(2*(X(:)'-xyrange(1))/(xyrange(2)-xyrange(1))-1,-1),1);
  X2 = min(max(2*(Y(:)'-xyrange(3))/(xyrange(4)-xyrange(3))-1,-1),1);
  TX1 = cos([0:n]'*acos(X1));
  TX2 = cos([0:n]'*acos(X2));
  TX1(2:n+1,:) = TX1(2:n+1,:)*sqrt(2);
  TX2(2:n+1,:) = TX2(2:n+1,:)*sqrt(2);
  LnfX = (TX1'*C0f*TX2)';
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdval
%
% at the Octave prompt
%
%!test
%! disp('Scattered target points')
%! xyrange = [0,1,0,1];
%! Pad = pdpts(20,xyrange);
%! C0f = pdcfsFFT(Pad,@funct);
%! X = [0,0;1/2,1/2;1,1];
%! Z = pdval(C0f,xyrange,X);
%! expected = [7.664205912849228e-01;3.2621734202884815e-01;...
%! 3.587865112678535e-02];
%! assert(Z,expected,10*eps)
%!test
%! disp('Meshgrid target points')
%! xyrange = [0,1,0,1];
%! Pad = pdpts(21,xyrange);
%! C0f = pdcfsFFT(Pad,@funct);
%! X1 = linspace(xyrange(1),xyrange(2),2);
%! X2 = linspace(xyrange(3),xyrange(4),2);
%! Z = pdval(C0f,xyrange,X1,X2);
%! expected = [7.664205912849229e-01,1.0757071952145181e-01;
%! 2.703371615911344e-01,3.5734971024838565e-02];
%! assert(Z,expected,10*eps)
