%-------------------------------------------------------------------------------
% USAGE of "demo_errors_pdint".
%
% This demo compares the relative errors (w.r.t. maximum deviation from average)
% of the interpolant of "f" in Padua Points computed respectively by FFT and 
% MM method.
%-------------------------------------------------------------------------------
% FUNCTIONS CALLED BY THIS MAIN PROGRAM:
% 
% 1. pdcfsFFT
% 2. pdcfsMM
% 3. pdpts
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
% Date: March 16, 2009.
%-------------------------------------------------------------------------------

clear all; close all; 
% more off;

n_vett=2:2:100;                            % PADUA POINTS DEGREES VECTOR.
xyrange=[0,1,0,1];                         % DEFINITION OF THE RECTANGLE.

m=110;

%-------------------------- SETTINGS END HERE ----------------------------------

test_points=pdpts(m,xyrange);
f_test_points=feval(@funct,test_points(:,1),test_points(:,2));
aver_fvalue=mean(f_test_points);
max_dev_aver_fvalue=norm(f_test_points-aver_fvalue,inf);
    
for index=1:length(n_vett)
    
    n=n_vett(index);
    
    fprintf('\n \n \t [DEGREE]: %3.0f',n);
    
    % FFT METHOD.
    
    Pad = pdpts(n,xyrange);
    C0f_FFT = pdcfsFFT(Pad,@funct,xyrange);
    LnfX_FFT = pdval(C0f_FFT,xyrange,test_points);
  
    fft_err_inf=norm(f_test_points-LnfX_FFT,inf);
    fft_err_2rel(index)=fft_err_inf/max_dev_aver_fvalue;
    fprintf('\n \t [FFT] [REL.INTP.ERR.INF.]: %2.5e',fft_err_2rel(index));
    
    % MM METHOD.
    
    [X1,Y1,X2,Y2] = pdpts(n,xyrange);
    C0f_MM = pdcfsMM(X1,Y1,X2,Y2,@funct,xyrange);
    LnfX_MM = pdval(C0f_MM,xyrange,test_points);
    
    MM_err_inf=norm(f_test_points-LnfX_MM,inf);
    MM_err_2rel(index)=MM_err_inf/max_dev_aver_fvalue;
    
    fprintf('\n \t [MM ] [REL.INTP.ERR.INF.]: %2.5e',MM_err_2rel(index));
    
end

semilogy(n_vett,fft_err_2rel,'g-d',n_vett,MM_err_2rel,'b-x');
legend('FFT','MM');
xlabel('Degree');
ylabel('Interpolation relative error');

fprintf('\n');

