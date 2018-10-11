function [LnfX,errest] = pdint(n,xyrange,funct,X,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdint".
%
% [LnfX,errest] = pdint(n,xyrange,funct,X)
% [LnfX,errest] = pdint(n,xyrange,funct,X1,X2)
% [LnfX,errest] = pdint(n,xyrange,funct,X,[],opt1,opt2,...)
% [LnfX,errest] = pdint(n,xyrange,funct,X1,X2,opt1,opt2,...)
%
% Compute the interpolation polynomial of degree n on the 
% Padua points defined in the rectangle 
% [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)] of the function
% funct, evaluated at the target points X(:,1),X(:,2) or at the 
% meshgrid(X1,X2) and the interpolation error estimate
%-------------------------------------------------------------------------------
% INPUT.
%
% n       : interpolation degree
% xyrange : a vector [a,b,c,d] defining the rectangle [a,b] x [c,d]
% funct   : function to be interpolated in the form 
%           funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%           optional arguments for f
% X       : a matrix with the abscissas of the target points in the
%           first column and the ordinates in the second one
% X1,X2   : two vectors defining the (mesh)grid X1 x X2 of the 
%           target points
%
% OUTPUT.
%
% LnfX    : interpolation polynomial at X(:,1),X(:,2) or 
%           at meshgrid(X1,X2)
% errest  : interpolation error estimate
%-------------------------------------------------------------------------------
% EXAMPLES.
%
% 1)
% Compute the interpolation polynomial of degree 20 of the 
% Franke's function defined in funct.m and the interpolation error
% estimate, and plot the interpolant on a mesh
%
% xyrange = [0,1,0,1];
% X1 = linspace(xyrange(1),xyrange(2),100);
% X2 = linspace(xyrange(3),xyrange(4),100);
% [LnfX,errest] = pdint(20,xyrange,@funct,X1,X2);
% [XX1,XX2] = meshgrid(X1,X2);
% mesh(XX1,XX2,LnfX)
%
% 2)
% Compute the interpolation polynomial of degree 20 of the 
% Gaussian function defined in funct.m and the interpolation error
% estimate, and plot the interpolant on a set of scattered data
%
% xyrange = [-1,1,-1,1];
% X = 2*rand(1000,2)-1;
% [LnfX,errest] = pdint(20,xyrange,@funct,X,[],11);
% plot3(X(:,1),X(:,2),LnfX,'o')
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% 1. pdpts
% 2. pdcfsFFT
% 3. pdval
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

% Compute the Padua points in the rectangle defined by xyrange
Pad = pdpts(n,xyrange);
if (nargin == 4)
% Target points as X(:,1),X(:,2)    
  [C0f,errest] = pdcfsFFT(Pad,funct,xyrange);
  LnfX = pdval(C0f,xyrange,X);
elseif (nargin >= 5)
  Y = varargin{1};
  if (size(X) == size(Y))
% Target points as meshgrid(X,Y)
    [C0f,errest] = pdcfsFFT(Pad,funct,xyrange,varargin{2:end});
    LnfX = pdval(C0f,xyrange,X,Y);
  else
% Target points as X(:,1),X(:,2)    
    [C0f,errest] = pdcfsFFT(Pad,funct,xyrange,varargin{2:end});
    LnfX = pdval(C0f,xyrange,X);
  end
else
  error('Wrong number of input arguments')
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdint
%
% at the Octave prompt
%
%!test
%! disp('Degree 0')
%! xyrange = [-1,1,-1,1];
%! X1 = 0;
%! X2 = 0;
%! LnfX = pdint(0,xyrange,@funct,X1,X2,13);
%! expected = 1;
%! assert(LnfX/expected,1,10*eps);
%!test
%! disp('Degree 20')
%! xyrange = [0,1,0,1];
%! X1 = linspace(xyrange(1),xyrange(2),100);
%! X2 = linspace(xyrange(3),xyrange(4),100);
%! LnfX = pdint(20,xyrange,@funct,X1,X2);
%! expected = 1.219525109531196e+00;
%! assert(max(max(abs(LnfX)))/expected,1,10*eps); 
%!test
%! disp('Degree 21, meshgrid target points')
%! xyrange = [0,1,0,1];
%! X1 = linspace(xyrange(1),xyrange(2),100);
%! X2 = linspace(xyrange(3),xyrange(4),100);
%! LnfX = pdint(21,xyrange,@funct,X1,X2);
%! expected = 1.219534733107023e+00;
%! assert(max(max(abs(LnfX)))/expected,1,10*eps); 
%!test
%! disp('Degree 21, scattered target points')
%! xyrange = [-1,1,-1,1];
%! X = [-1,-1;0,0;1,1];
%! LnfX = pdint(21,xyrange,@funct,X,[],11);
%! expected = 9.99999999980831e-01;
%! assert(max(max(abs(LnfX)))/expected,1,10*eps); 
%!test
%! disp('Error estimate')
%! xyrange = [0,1,0,1];
%! X = [0,0;1,1];
%! [LnfX,errest] = pdint(21,xyrange,@funct,X);
%! expected = 2.2486904061664328e-03;
%! assert(errest/expected,1,200*eps); 
