function [mag, phase] = mag_curve(struture, q, delta, omega, ystate, N)
%%       Transfer Function H(s) Curve for the Fractional-order Chaos 
%%                  Circuit Parameter Selection Problem

% struture : The structure of fractance, Desired/Chain/Tree/Ladder/ZPK
% q : the order of fractional-order chaos system
% delta : dB, the maximum discrepancy between desired H(s) and \hat{H}(s)
% omega :  vector of frequency sample points
% N :  the number of fractance system order

% Copyright (c) 2020, Kunpeng Wang.
% 

assert(ischar(struture)&&isrow(struture),'struture must be a 1xN char.')

tau0 = 100;  % the relaxation time consistant
% Select transfer function type :
switch struture
	case 'Desired'
		Hs = (1+tau0*(1j*omega)).^(-q);
        mag = abs(Hs);
        phase = unwrap(angle(Hs))*180/pi;
        
	case 'Chain'
        kappa = ystate(end);
        y = ystate(1:end-1)';
        yr = y(1:N);
        yc = y(N+1:end)*1e-12; % unit is pF
        
        Hs = 0;
        for i=1:N
            Hs = Hs + 1./(yc(i)*1j*omega + 1./yr(i));
        end
        
        Hs = kappa*Hs;
        mag = abs(Hs);
        phase = unwrap(angle(Hs))*180/pi;
        
	case 'Tree'
        kappa = ystate(end);
        y = ystate(1:end-1)';
        yr = y(1:N);
        yc = y(N+1:end)*1e-12;
 
        switch N
            case 1
                Hs = 1./(yc(1)*1j*omega + 1./yr(1));
            case 2
                Hs = 1./(1./(1./(yc(1)*1j*omega)) + ...
                     1./(yr(1)+1./(yc(2)*1j*omega+1./yr(2))));
            case 3
                Hs = 1./(1./(1./(yc(3)*1j*omega+1./yr(3))+ ...
                     1./(yc(1)*1j*omega)) + ...
                     1./(1./(yc(2)*1j*omega+1./yr(2))+yr(1)));
            case 4
                Hs = 1./(1./(1./(1./(1./(yc(4)*1j*omega+1./yr(4)) + ...
                     yr(2))+yc(2)*1j*omega)+yr(1)) + ...
                     1./(1./(yc(3)*1j*omega+1./yr(3)) + ...
                     1./(yc(1)*1j*omega)));
            case 5
                Hs = 1./(1./(1./(1./(1./(yc(5)*1j*omega+1./yr(5)) + ...
                     1./(yc(2)*1j*omega)) + ...
                     1./(1./(yc(4)*1j*omega+1./yr(4))+yr(2))) + ...
                     yr(1))+1./(1./(yc(3)*1j*omega+1./yr(3)) + ...
                     1./(yc(1)*1j*omega)));
            case 6
                Hs = 1./(1./(1./(1./(1./(yc(6)*1j*omega + ...
                     1./yr(6))+yr(3)) + yc(3)*1j*omega) + ...
                     1./(yc(1)*1j*omega)) + ...
                     1./(1./(1./(1./(yc(5)*1j*omega + ...
                     1./yr(5)) + 1./(yc(2)*1j*omega)) + ...
                     1./(1./(yc(4)*1j*omega + 1./yr(4))+yr(2)))+yr(1)));
            case 7
                Hs = 1./(1./(1./(1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega)) + ...
                     1./(1./(yc(6)*1j*omega + 1./yr(6))+yr(3))) + ...
                     1./(yc(1)*1j*omega)) + ...
                     1./(1./(1./(1./(yc(5)*1j*omega+1./yr(5)) + ...
                     1./(yc(2)*1j*omega)) + 1./(1./(yc(4)*1j*omega + ...
                     1./yr(4))+yr(2)))+yr(1)));
            case 8
                Hs = 1./(1./(1./(1./(1./(1./(1./(yc(8)*1j*omega + ...
                     1./yr(8))+yr(4))+yc(4)*1j*omega)+yr(2)) + ...
                     1./(1./(yc(5)*1j*omega+1./yr(5)) + ...
                     1./(yc(2)*1j*omega)))+yr(1)) + ...
                     1./(1./(1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega))+1./(1./(yc(6)*1j*omega + ...
                     1./yr(6))+yr(3))) + 1./(yc(1)*1j*omega)));
            case 9
                Hs = 1./(1./(1./(1./(1./(1./(1./(yc(9)*1j*omega + ...
                     1./yr(9))+1./(yc(4)*1j*omega)) + ...
                     1./(1./(yc(8)*1j*omega + ...
                     1./yr(8))+yr(4)))+yr(2)) + ...
                     1./(1./(yc(5)*1j*omega + ...
                     1./yr(5))+1./(yc(2)*1j*omega)))+yr(1)) + ...
                     1./(1./(1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega))+1./(1./(yc(6)*1j*omega + ...
                     1./yr(6))+yr(3))) + ...
                     1./(yc(1)*1j*omega)));
            case 10
                Hs = 1./(1./(1./(1./(1./(1./(1./(yc(10)*1j*omega + ...
                     1./yr(10))+yr(5))+yc(5)*1j*omega) + ...
                     1./(yc(2)*1j*omega)) + ...
                     1./(1./(1./(1./(yc(9)*1j*omega+1./yr(9)) + ...
                     1./(yc(4)*1j*omega))+1./(1./(yc(8)*1j*omega + ...
                     1./yr(8))+yr(4)))+yr(2)))+yr(1)) + ...
                     1./(1./(1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega)) + ...
                     1./(1./(yc(6)*1j*omega+1./yr(6))+yr(3))) + ...
                     1./(yc(1)*1j*omega)));
            case 11
                Hs = 1./(1./(1./(1./(1./(1./(1./(yc(11)*1j*omega + ...
                     1./yr(11))+1./(yc(5)*1j*omega)) + ...
                     1./(1./(yc(10)*1j*omega+1./yr(10))+yr(5))) + ...
                     1./(yc(2)*1j*omega)) + ...
                     1./(1./(1./(1./(yc(9)*1j*omega+1./yr(9)) + ...
                     1./(yc(4)*1j*omega)) + ...
                     1./(1./(yc(8)*1j*omega + ...
                     1./yr(8))+yr(4)))+yr(2)))+yr(1)) + ...
                     1./(1./(1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega)) + ...
                     1./(1./(yc(6)*1j*omega+1./yr(6))+yr(3))) + ...
                     1./(yc(1)*1j*omega)));
            case 12
                Hs = 1./(1./(1./(1./(1./(1./(1./(yc(12)*1j*omega + ...
                     1./yr(12))+yr(6))+yc(6)*1j*omega)+yr(3)) + ...
                     1./(1./(yc(7)*1j*omega+1./yr(7)) + ...
                     1./(yc(3)*1j*omega))) + ...
                     1./(yc(1)*1j*omega)) + ...
                     1./(1./(1./(1./(1./(1./(yc(11)*1j*omega + ...
                     1./yr(11))+1./(yc(5)*1j*omega)) + ...
                     1./(1./(yc(10)*1j*omega + 1./yr(10))+yr(5))) + ...
                     1./(yc(2)*1j*omega)) + ...
                     1./(1./(1./(1./(yc(9)*1j*omega+1./yr(9)) + ...
                     1./(yc(4)*1j*omega))+1./(1./(yc(8)*1j*omega + ...
                     1./yr(8))+yr(4)))+yr(2)))+yr(1)));
            case 13
                Hs = 1./(1./(1./(1./(1./(1./(1./(...
                     1./yr(9)+1j*omega*yc(9)) + ...
                     1./(1j*omega*yc(4))) + ...
                     1./(1./(1./yr(8)+1j*omega*yc(8))+yr(4)))+yr(2)) + ...
                     1./(1./(1./(yr(5)+1./(1./yr(10)+1j*omega*yc(10))) + ...
                     1./(1./(1./yr(11)+1j*omega*yc(11)) + ...
                     1./(1j*omega*yc(5)))) + ...
                     1./(1j*omega*yc(2))))+yr(1)) + ...
                     1./(1./(1./(1./(1./yr(7)+1j*omega*yc(7)) + ...
                     1./(1j*omega*yc(3)))+1./(1./(1./(yr(6) + ...
                     1./(1./yr(12)+1j*omega*yc(12))) + ...
                     1./(1./(1./yr(13)+1j*omega*yc(13)) + ...
                     1./(1j*omega*yc(6))))+yr(3)))+1./(1j*omega*yc(1))));
            case 14
                Hs = 1./(1./(1./(1./(1./(1./(1./(...
                     1./yr(9)+1j*omega*yc(9)) + ...
                     1./(1j*omega*yc(4))) + ...
                     1./(1./(1./yr(8)+1j*omega*yc(8))+yr(4)))+yr(2)) + ...
                     1./(1./(1./(yr(5) + ...
                     1./(1./yr(10)+1j*omega*yc(10))) + ...
                     1./(1./(1./yr(11)+1j*omega*yc(11)) + ...
                     1./(1j*omega*yc(5)))) + ...
                     1./(1j*omega*yc(2))))+yr(1)) + ...
                     1./(1./(1./(1./(1./(yr(7) + ...
                     1./(1./yr(14)+1j*omega*yc(14)))+1j*omega*yc(7)) + ...
                     1./(1j*omega*yc(3)))+1./(1./(1./(yr(6) + ...
                     1./(1./yr(12)+1j*omega*yc(12))) + ...
                     1./(1./(1./yr(13)+1j*omega*yc(13)) + ...
                     1./(1j*omega*yc(6))))+yr(3)))+1./(1j*omega*yc(1))));
            case 15
                Hs = 1./(1./(1./(1./(1./(1./(1./(...
                     1./yr(9)+1j*omega*yc(9)) + ...
                     1./(1j*omega*yc(4))) + ...
                     1./(1./(1./yr(8)+1j*omega*yc(8))+yr(4)))+yr(2)) + ...
                     1./(1./(1./(yr(5) + ...
                     1./(1./yr(10)+1j*omega*yc(10))) + ...
                     1./(1./(1./yr(11)+1j*omega*yc(11)) + ...
                     1./(1j*omega*yc(5)))) + ...
                     1./(1j*omega*yc(2))))+yr(1)) + ...
                     1./(1./(1./(1./(1./(yr(7) + ...
                     1./(1./yr(14)+1j*omega*yc(14))) + ...
                     1./(1./(1./yr(15)+1j*omega*yc(15)) + ...
                     1./(1j*omega*yc(7)))) + ...
                     1./(1j*omega*yc(3)))+1./(1./(1./(yr(6) + ...
                     1./(1./yr(12)+1j*omega*yc(12))) + ...
                     1./(1./(1./yr(13)+1j*omega*yc(13)) + ...
                     1./(1j*omega*yc(6))))+yr(3)))+1./(1j*omega*yc(1))));
            otherwise
                error('For tree system, order N="%s" is not supported.',N)
        end
 
        Hs = kappa*Hs;
        mag = abs(Hs);
        phase = unwrap(angle(Hs))*180/pi;
        
	case 'Ladder'
        kappa = ystate(end);
        y = ystate(1:end-1)';
        yr = y(1:N);
        yc = y(N+1:end)*1e-12;
        Hs = ladder_hs(N, omega, yr, yc);
        
        Hs = kappa*Hs;
        mag = abs(Hs);
        phase = unwrap(angle(Hs))*180/pi;
        
    case 'ZPK'
        p_0 = 1/tau0 * 10^(delta/(20*q));
        a = 10^(delta/(10*(1-q)));
        b = 10^(delta/(10*q));
        w_max = max(omega)/(2*pi);
        % the total number of the poles in H(s)
        tf_N = 1 + floor(log10(w_max/p_0)/log10(a*b)); 

        hsz = zeros(1, tf_N);
        hsp = zeros(1, tf_N+1);
        for i=0:tf_N
           if i < tf_N
                hsz(i+1) = ((a*b)^i*a*p_0); % zeros
           end
           hsp(i+1) = ((a*b)^i*p_0); % poles
        end
        hsk = b^tf_N * p_0; % gain

        Hs = zpk(hsz, hsp, hsk); % transfer function H(s)
        [mag, phase] = bode(Hs, omega);
        mag = squeeze(mag);
        phase = -(squeeze(phase)+180);
        %tf_N = tf_N + 1; % the order of transfer function
	otherwise
		error('Struture "%s" is not supported.',struture)
end

end

% subfunction of ladder structure
function Hs = ladder_hs(n, omega, yr, yc)
    if n == 1
        Hs = 1./(yc(n)*1j*omega + 1./yr(n));
    else
        Hs = 1./(yc(n)*1j*omega + 1./(yr(n) + ...
             ladder_hs(n-1, omega, yr, yc)));
    end
end

