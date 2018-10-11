function [C0f,errest,varargout] = pdcfsFFT(Pad,funct,varargin)

%-------------------------------------------------------------------------------
% USAGE of "pdcfsFFT".
%
% [C0f,errest] = pdcfsFFT(Pad,funct,[],opt1,opt2,...)
% [C0f,errest,cubature] = pdcfsFFT(Pad,funct,xyrange,opt1,opt2,...)
%
% Compute the coefficient matrix C0f by FFT and the error estimate.
% Optionally, compute the cubature through the coefficient matrix. 
%-------------------------------------------------------------------------------
% INPUT.
%
% Pad      : Padua points, as computed through the call
%            Pad = pdpts(n) or Pad = pdpts(n,xyrange), being n the 
%            interpolation degree
% funct    : function to be interpolated in the form 
%            funct(x,y,opt1,opt2,...), where opt1, opt2, ... are
%            optional arguments for funct
% xyrange  : an optional vector [a,b,c,d] defining the rectangle 
%            [a,b] x [c,d] required only for cubature
%
% OUTPUT.
%
% C0f      : coefficient matrix
% errest   : interpolation error estimate
% cubature : cubature through the coefficient matrix 
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
% Date: March 04, 2009.
%-------------------------------------------------------------------------------

N = size(Pad,1);
% recover the degree n from N = (n+1)(n+2)/2
n = round(-3+sqrt(1+8*N))/2;
if (n == 0)
% evaluate the function at the Padua point
  C0f = feval(funct,Pad(:,1),Pad(:,2),varargin{2:end});
else
% evaluate the function at the Padua points
  [M1,M2] = meshgrid([0:n],[0:n+1]);
  findM = find(mod(M1+M2,2));
  GfT = zeros((n+2)*(n+1),1);
  GfT(findM) = feval(funct,Pad(:,1),Pad(:,2),varargin{2:end});
  GfT = reshape(GfT,size(M1));
  GfT = GfT*2;
  GfT = GfT/(n*(n+1));
  GfT(1,:) = GfT(1,:)/2;
  GfT(n+2,:) = GfT(n+2,:)/2;
  GfT(:,1) = GfT(:,1)/2;
  GfT(:,n+1) = GfT(:,n+1)/2;
  Gf = GfT';
% compute the coefficient matrix C0f by FFT
  Gfhat = real(fft(Gf,2*n));
  Gfhat = Gfhat(1:n+1,:);
  Gfhathat = real(fft(Gfhat,2*(n+1),2));
  C0f = Gfhathat(:,1:n+1);
  C0f = 2*C0f;
  C0f(1,:) = C0f(1,:)/sqrt(2);
  C0f(:,1) = C0f(:,1)/sqrt(2);
  C0f = fliplr(triu(fliplr(C0f)));
  C0f(n+1,1) = C0f(n+1,1)/2;
end
C0f2 = fliplr(C0f);
% error estimate
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
% test pdcfsFFT
%
% at the Octave prompt
%
%!test
%! disp('Constant function interpolation')
%! Pad = pdpts(10);
%! C0f = pdcfsFFT(Pad,@funct,[],13);
%! expected = zeros(11);
%! expected(1,1) = 1;
%! assert(C0f,expected,10*eps)
%!test
%! disp('Degree 0 interpolation')
%! Pad = pdpts(0);
%! C0f = pdcfsFFT(Pad,@funct,[],13);
%! expected = 1;
%! assert(C0f/expected,1,10*eps)
%!test
%! disp('Cubature through the coefficients')
%! xyrange = [0,1,0,1];
%! Pad = pdpts(500,xyrange);
%! [C0f,errest,cubature] = pdcfsFFT(Pad,@funct,xyrange);
%! expected = 4.06969589491556e-01;
%! assert(cubature/expected,1,10*eps)
