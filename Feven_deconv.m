function [zFeven, Feven] = Feven_deconv(zk, ktscap, A0, sgwin, sgdegree, tozero)
%% FEVEN_DECONV Implements the Sader-Jarvis algorithm for force calculation
%
%  Uses the Sader-Jarvis formula [1] to calculate the force from 
%  experimental AFM data. The pole at z=zmin is removed using the analytic 
%  solution of the integral. 
%  Savitzky-Golay filtering [2] can be used to perform a the 
%  numerical differentiation. 
%
%  Input data is given in the form of the cap-averaged force gradient
%  according to Söngen et al. [3]. This quantity is commonly calculated
%  from the frequency shift df, eigenfrequency f0 and cantilever 
%  stiffness k0 by
%  ktscap = k0.*(1-(((f0)+df(:))./(f0)).^2)
%
%  Arguments:
%  ------------
%   zk       - z axis
%   ktscap   - cap-averaged force gradient (measured; from df, in N/m)
%   A0       - oscillation amplitude (zero-peak, in m)
%   sgwin    - Window size for Savitzky-Golay filter
%   sgdegree - Degree of Savitzky-Golay filter polynom
%   tozero   - shift df data to zero at large z by averaging 
%              over tozero values 
%              (0: no shift along vertical axis)
%
%  Returns:
%  ------------
%   zFeven   - z axis
%   Feven    - even force
%
%  References:
%  [1] J.E. Sader and S.P. Jarvis, Appl. Phys. Lett. 84, 1801 (2004)
%      doi: 10.1063/1.1667267
%  [2] A. Savitzky, M.J.E. Golay, Anal. Chem. 36, 1627 (1964)
%      doi: 10.1021/ac60214a047
%  [3] H. Söngen, R. Bechstein, A. Kühnle, 
%      J. Phys. Cond. Matter 29, 274001 (2017)
%      doi: 10.1088/1361-648X/aa6f8b
%
%  Copyright 2020, Philipp Rahe
%
%  This is a script under the terms of the Creative Commons Attribution 
%  License (creativecommons.org/licenses/by/4.0), which permits 
%  unrestricted use, distribution, and reproduction in any medium, 
%  provided the original work is properly cited.
%
%  version 14.08.2020, Philipp Rahe (prahe@uos.de)

  % different checks on input parameters
  N = numel(zk);
  if N ~= numel(ktscap)
    error('Lengths of z (len %u) and df (len %u) do not match.', ...
          numel(zk), numel(ktscap));
  end
  if(~iscolumn(zk)), zk = zk';  end
  if(~iscolumn(ktscap)), ktscap = ktscap';  end
  
  % shift data to zero at large z (if requested)
  if(tozero > 0) 
    ktscap = ktscap - mean(ktscap((length(ktscap)-tozero):end));
  end
  
  % calculate the derivation in the Sader-Jarvis formula
  if( (sgwin>1) && (sgdegree>=1) && (sgwin>sgdegree) && mod(sgwin,2)>0 ) 
    % use Savitzky-Golay filter
    [~, ddfdz] = savgolfilter(zk, ktscap, sgwin, 1, sgdegree);
  else
    if( (sgwin~=0) )
      warning('Invalid filter parameter. Using central difference.');
    end
    % use central difference (left/right difference at edges)
    ddfdz = ktscap.*0;
    ddfdz(1)       = (ktscap(2)-ktscap(1))./(zk(2)-zk(1));
    ddfdz(2:end-1) = (ktscap(1:end-2) - ktscap(3:end) )./(zk(1:end-2) - zk(3:end));
    ddfdz(end)     = (ktscap(end)-ktscap(end-1))./(zk(end)-zk(end-1));
  end
  
  % calculate central difference for z (left/right difference at edges)
  ddz = zk.*0;
  ddz(1)       = zk(2)-zk(1);
  ddz(2:end-1) = (zk(1:end-2) - zk(3:end))/2.;
  ddz(end)     = zk(end)-zk(end-1);
  
  
  % prepare output array
  Feven = zeros(1,N-2);
  
  % iterate through data
  for j=1:(N-2)
    
    % integration
    int = -trapz(zk((j+1):N),...
                 (1 + sqrt(A0./(64*pi.*(zk((j+1):N)-zk(j))))).*...
                           ktscap((j+1):N) ...
	               - sqrt(A0.^3./(2.*(zk((j+1):N)-zk(j)))).* ...
                           ddfdz((j+1):N));

    % correction factors to remove pole of the integration
    corr0 = -ktscap(j)*abs(ddz(j));
    corr1 = -sqrt(A0/(16*pi))*ktscap(j)*sqrt(abs(ddz(j)));
    corr2 = sqrt(4*A0^3/2).*ddfdz(j).*sqrt(abs(ddz(j)));
    
    % resulting even force
    Feven(j) = int + corr0 + corr1 + corr2;
    
  end
  
  % return according z axis. 
  zFeven = zk(1:N-2);


end
