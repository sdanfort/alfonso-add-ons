function [C0f,errest,varargout] = pdcfsMM(X1,Y1,X2,Y2,funct,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdcfsMM".
%
% [C0f,errest] = padua2MM(X1,Y1,X2,Y2,funct,[],opt1,opt2,...)
% [C0f,errest,cubature] = padua2MM(X1,Y1,X2,Y2,funct,xyrange,opt1,opt2,...)
%-------------------------------------------------------------------------------
% INPUT.
%
% X1,Y1,X2,Y2 : Padua points, as computed through the call
%               [X1,Y1,X2,Y2] = pdpts(n) or 
%               [X1,Y1,X2,Y2] = pdpts(n,xyrange), being n the 
%               interpolation degree
% funct       : function to be interpolated in the form 
%               funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%               optional arguments for funct
% xyrange     : an optional vector [a,b,c,d] defining the rectangle 
%               [a,b] x [c,d] required only for cubature
%
% OUTPUT.
%
% C0f         : coefficient matrix
% errest      : interpolation error estimate
% cubature    : cubature through the coefficient matrix 
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

% recover the degree n
n = max(length(X2)+length(X1)-2,0);
if (n == 0)
% evaluate the function at the Padua point
  C0f = feval(funct,X1,Y1,varargin{2:end});
else
% evaluate the function at the Padua points
  GfT1 = feval(funct,X1,Y1,varargin{2:end});
  GfT2 = feval(funct,X2,Y2,varargin{2:end});
  GfT1 = GfT1*2;
  GfT2 = GfT2*2;
  GfT1 = GfT1/(n*(n+1));
  GfT2 = GfT2/(n*(n+1));
  GfT1(:,1) = GfT1(:,1)/2;
  GfT2(1,:) = GfT2(1,:)/2;
  if (mod(n,2) == 0)
    GfT1(:,n/2+1) = GfT1(:,n/2+1)/2;
    GfT1(n/2+1,:) = GfT1(n/2+1,:)/2;
  else
    GfT2((n+1)/2+1,:) = GfT2((n+1)/2+1,:)/2;
    GfT2(:,(n+1)/2) = GfT2(:,(n+1)/2)/2;
  end
% compute the coefficient matrix C0f by MM
  argn = linspace(0,pi,n+1);
  argn1 = linspace(0,pi,n+2);
  TE1 = cos([0:n]'*argn(1:2:n+1));
  TE1(2:n+1,:) = TE1(2:n+1,:)*sqrt(2);
  TO1 = cos([0:n]'*argn(2:2:n+1));
  TO1(2:n+1,:) = TO1(2:n+1,:)*sqrt(2);
  TE2 = cos([0:n]'*argn1(1:2:n+2));
  TE2(2:n+1,:) = TE2(2:n+1,:)*sqrt(2);
  TO2 = cos([0:n]'*argn1(2:2:n+2));
  TO2(2:n+1,:) = TO2(2:n+1,:)*sqrt(2);
  Gf1 = GfT1';
  Gf2 = GfT2';
  C0f1 = TE1*Gf1*TO2';
  C0f2 = TO1*Gf2*TE2';
  C0f = C0f1+C0f2;
  C0f = (fliplr(triu(fliplr(C0f))));
  C0f(n+1,1) = C0f(n+1,1)/2;
end
C0f2 = fliplr(C0f);
errest = sum(abs(diag(C0f2)));
if (n >= 1)
  errest = errest+sum(abs(diag(C0f2,-1)));
end
if (n >= 2)
  errest = errest+sum(abs(diag(C0f2,-2)));
end
errest = 2*errest;
if (nargout == 3)
% cubature required
  if (nargin < 3)
    error('xyrange required for cubature')
  else
    xyrange = varargin{1};
  end    
% compute the moments
  k = [0:2:n]';
  mom = 2*sqrt(2)./(1-k.^2);
  mom(1) = 2;
  [M1,M2] = meshgrid(mom);
  Mmom = M1.*M2;
  CM = C0f(1:2:n+1,1:2:n+1).*Mmom;
  cubature = (xyrange(2)-xyrange(1))*(xyrange(4)-xyrange(3))*sum(sum(CM))/4;
  varargout{1} = cubature;
end

%-------------------------------------------------------------------------------
% OCTAVE TESTS.
%-------------------------------------------------------------------------------
% Octave testing: type
%
% test pdcfsMM
%
% at the Octave prompt
%
%!test
%! disp('Constant function interpolation')
%! [X1,Y1,X2,Y2] = pdpts(10);
%! C0f = pdcfsMM(X1,Y1,X2,Y2,@funct,[],13);
%! expected = zeros(11);
%! expected(1,1) = 1;
%! assert(C0f,expected,10*eps)
%!test
%! disp('Degree 0 interpolation')
%! [X1,Y1,X2,Y2] = pdpts(0);
%! C0f = pdcfsMM(X1,Y1,X2,Y2,@funct,[],13);
%! expected = 1;
%! assert(C0f/expected,1,10*eps)
%!test
%! disp('Cubature through the coefficients')
%! xyrange = [0,1,0,1];
%! [X1,Y1,X2,Y2] = pdpts(500,xyrange);
%! [C0f,errest,cubature] = pdcfsMM(X1,Y1,X2,Y2,@funct,xyrange);
%! expected = 4.06969589491556e-01;
%! assert(cubature/expected,1,10*eps)
