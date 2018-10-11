function [varargout] = pdcub(n,xyrange,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdcub".
%
% PadL = pdcub(n,xyrange)
% [cubature,PadL] = pdcub(n,xyrange,funct,opt1,opt2,...)
%
% This function computes the Padua points defined in the rectangle 
% [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)] and the 
% corresponding cubature weights as a matrix of abscissas (first
% column), ordinates (second column) and weights (third column).
% Optionally, it computes the cubature of the function f (optional 
% input argument).
%-------------------------------------------------------------------------------
% INPUT.
%
% n        : interpolation degree
% xyrange  : a vector [a,b,c,d] defining the rectangle [a,b] x [c,d]
% funct    : function to be interpolated in the form 
%            funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%            optional arguments for funct 
%
% OUTPUT.
%
% PadL     : a matrix with the abscissas of Padua points in the 
%            first column, the ordinates in the second and the 
%            cubature weights in the third
% cubature : cubature of the integrand funct
%-------------------------------------------------------------------------------
% EXAMPLES.
%
% 1)
% Compute the Padua points and the cubature weights of degree 50
%
% xyrange = [0,1,0,1];
% PadL = pdcub(50,xyrange);
%
% 2)
% Compute the cubature of the Franke's function defined in funct.m
%
% xyrange = [0,1,0,1];
% [cubature,PadL] = pdcub(50,xyrange,@f);
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS CODE:
% 1. pdpts
% 2. pdwtsMM
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
% Date: March 14, 2009.
%-------------------------------------------------------------------------------

% Compute the Padua points in the rectangle defined by xyrange
PadL = pdpts(n,xyrange);
% Compute the cubature weights
PadL(:,3) = pdwtsMM(n,xyrange);


if (nargin < 2)
  error('Too few input arguments')
elseif (nargin >= 3)
  funct = varargin{1};
  varargout{1} = PadL(:,3)'*feval(funct,PadL(:,1),PadL(:,2),varargin{2:end});
  varargout{2} = PadL;
else
  varargout{1} = PadL;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdcub
%
% at the Octave prompt
%
%!test
%! disp('Degree 0 cubature')
%! xyrange = [0,1,0,1];
%! [cubature,PadL] = pdcub(0,xyrange,@funct,13);
%! expected = 1;
%! assert(cubature/expected,1,10*eps);
%!test
%! disp('Franke''s function cubature')
%! xyrange = [0,1,0,1];
%! [cubature,PadL] = pdcub(500,xyrange,@funct);
%! expected = 4.069695894915573e-01;
%! assert(cubature/expected,1,20*eps);  
