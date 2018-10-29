function [ eD ] = disturbance_effect( D,tau,taue,taui,eoi )
eD = sum( D .* (eoi * exp(-tau/taue) - exp(-tau/taui)) );
end

