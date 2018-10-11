%-------------------------------------------------------------------------------
% USAGE of "demo_pdint".
%
% A demo that computes the values of the interpolant of a function "f" 
% in a prescribed set of points, knowing its samples at Padua Points.
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdint
% 2. pdpts
% 3. funct
%-------------------------------------------------------------------------------
% TESTED ON:
%
% 1. GNU Octave 3.1.50.
% 2. Matlab 6.1.0.450. Release 12.1. 
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
% Date: March 16, 2009.
%-------------------------------------------------------------------------------

clear all; close all; more off;

n=100;                      % PADUA POINTS DEGREE.
xyrange=[0,pi,-pi,pi];      % DEFINITION OF THE RECTANGLE.

M=30;                       % EVALUATION POINTS IN THE SQUARE: DEGREE OF PADUA
                            % POINTS USED FOR TESTS. 

%-------------------------- SETTINGS END HERE ----------------------------------

X = pdpts(M,xyrange);       % [xyrange(1),xyrange(2)]x[xyrange(3),xyrange(4)]
                            % STORED IN A MATRIX HAVING SIZE  "M x 2".
                            % THE FIRST COLUMN CONTAINS THE ABSCISSAS,
                            % THE SECOND COLUMN CONTAINS THE ORDINATES.
                                        
LnfX = pdint(n,xyrange,@funct,X); 
                            % "LnfX" IS THE INTERPOLANT AT PADUA POINTS OF 
                            % THE FUNCTION, STORED IN A MATRIX HAVING SIZE 
                            % "M x 1", "M" BEING "(n+1)*(n+2)/2".
                                        
fX=feval(@funct,X(:,1),X(:,2)); 
                            % EVALUATION OF THE FUNCTION AT PADUA POINTS.

absmaxerr=norm(fX-LnfX);    % COMPUTING MAXIMUM ABSOLUTE ERROR.

aver_fvalue=mean(fX);

index=find(abs(fX)>0);      % COMPUTING MAXIMUM RELATIVE ERROR.
relmaxerr=norm(fX(index)-LnfX(index),inf)/aver_fvalue;        

fprintf('\n\t [n]: %5.0f [max.abs.err.]: %2.2e [rel.abs.err.]: %2.2e \n \n',...
n,absmaxerr,relmaxerr);   

