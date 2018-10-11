%-------------------------------------------------------------------------------
% USAGE of "demo_errors_pdcub".
%
% This demo computes the cubature of a function "f" on a rectangle by an 
% algebraic formula of degree "n" on Padua points. 
% A comparison of relative errors between FFT method and MM method is made.
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdpts
% 2. pdwtsFFT
% 3. pdwtsMM
% 4. funct
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
% Date: April 05, 2009.
%-------------------------------------------------------------------------------

clear all; close all; 
% more off;

n_vett=2:2:100;                            % PADUA POINTS DEGREES VECTOR.
xyrange=[0,1,0,1];                         % DEFINITION OF THE RECTANGLE.
n_exact=200;

%-------------------------- SETTINGS END HERE ----------------------------------

PadL_exact = pdpts(n_exact,xyrange);
fPadL_exact=feval(@funct,PadL_exact(:,1),PadL_exact(:,2));
w_exact = pdwtsFFT(n_exact,xyrange);
I_exact = w_exact'*fPadL_exact;                                           

for index=1:length(n_vett)
    
    n=n_vett(index);
    
    PadL = pdpts(n,xyrange);
    fPadL=feval(@funct,PadL(:,1),PadL(:,2));
    
    % FFT METHOD.
    wFFT = pdwtsFFT(n,xyrange);
    cubFFT(index) = wFFT'*fPadL;
    relerrFFT(index)=abs((cubFFT(index)-I_exact)/I_exact);
    
    % MM METHOD.
    wMM = pdwtsMM(n,xyrange);
    cubMM(index) = wMM'*fPadL;
    relerrMM(index)=abs((cubMM(index)-I_exact)/I_exact);
    
 fprintf('\n \n \t [DEGREE]: %4.0f [REL.ERR.FFT]: %2.5e [REL.ERR.MM]: %2.5e',...
        n,relerrFFT(index),relerrMM(index));
    
end

semilogy(n_vett,relerrFFT,'g-d',n_vett,relerrMM,'b-x');
legend('FFT','MM');
xlabel('Degree')
ylabel('Relative cubature error')

fprintf('\n \n');