function [varargout] = pdpts(n,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdpts".
%
% Pad = pdpts(n)
% Pad = pdpts(n,xyrange)
% [X1,Y1,X2,Y2] = pdpts(n)
% [X1,Y1,X2,Y2] = pdpts(n,xyrange)
%
% Compute the (first family of) Padua points, either as a matrix Pad 
% with their abscissas in the first column and their ordinates in the
% second, or as two subgrids X1,Y1 and X2,Y2, respectively.
%-------------------------------------------------------------------------------
% INPUT.    
%
% n           : interpolation degree
% xyrange     : an optional vector [a,b,c,d] defining the rectangle 
%               [a,b] x [c,d]. Otherwise, xyrange = [-1,1,-1,1]
%
% OUTPUT.  
%
% Pad         : matrix of size ((n+1)*(n+2)/2) x 2 such
%               that (Pad(:,1),Pad(:,2)) defines the Padua points in the
%               rectangle [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)].  
% X1,Y1,X2,Y2 : the two subgrids X1,Y1 and X2,Y2 defining the Padua points
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
% Date: March 09, 2009.
%-------------------------------------------------------------------------------

if (nargin == 1)
% standard square [-1,1] x [-1,1]
  xyrange = [-1,1,-1,1];
else
% rectangle [xyrange(1),xyrange(2)] x [xyrange(3),xyrange(4)]
  xyrange = varargin{1};
end
if (n == 0)
% degree 0
  if (nargout ~= 4)
% points as a single matrix
    Pad = [xyrange(1),xyrange(3)];
    varargout = {Pad};
  else
% points as two subgrids
    X1 = xyrange(1);
    Y1 = xyrange(3);
    X2 = zeros(1,0);
    Y2 = zeros(1,0);
    varargout = {X1,Y1,X2,Y2};
  end  
else
% degree > 0
  zn = (xyrange(1)+xyrange(2)+(xyrange(2)-xyrange(1))*...
       cos(linspace(0,1,n+1)*pi))/2;
  zn1 = (xyrange(3)+xyrange(4)+(xyrange(4)-xyrange(3))*...
	cos(linspace(0,1,n+2)*pi))/2;
  if (nargout ~= 4)
% points as a single matrix
    [Pad1,Pad2] = meshgrid(zn,zn1);
    [M1,M2] = meshgrid([0:n],[0:n+1]);
    findM = find(mod(M1+M2,2));
    Pad = [Pad1(findM),Pad2(findM)];
    varargout = {Pad};
  else
% points as two (mesh)grids  
    En = zn(1:2:n+1);
    On = zn(2:2:n+1);
    En1 = zn1(1:2:n+2);
    On1 = zn1(2:2:n+2);
    [X1,Y1] = meshgrid(En,On1);
    [X2,Y2] = meshgrid(On,En1);
    varargout = {X1,Y1,X2,Y2};
  end  
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdpts
%
% at the Octave prompt
%
%!test
%! disp('Degree 0')
%! Pad = pdpts(0);
%! expected = [-1,-1];
%! assert(Pad,expected,10*eps); 
%!test
%! disp('Degree 1')
%! Pad = pdpts(1);
%! expected = [cos([0;1;1]*pi),cos([1;0;2]*pi/2)];
%! assert(Pad,expected,10*eps); 
%!test
%! disp('Degree 2 with xyrange')
%! Pad = pdpts(2,[0,1,0,2]);
%! expected = [(cos([0;0;1;1;2;2]*pi/2)+1)/2,cos([1;3;0;2;1;3]*pi/3)+1];
%! assert(Pad,expected,10*eps); 
