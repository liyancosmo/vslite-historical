function [ eD ] = disturbance_effects( D,taui,taue,eoi )
eD = zeros(size(D));
for i = 1:length(D)
    Di = D(1:i);
    tau = i - (1:i);
    eD(i) = disturbance_effect(Di,tau,taui,taue,eoi);
end
end

function [ eD ] = disturbance_effect( D,tau,taui,taue,eoi )
eD = sum( D .* (eoi * exp(-tau/taue) - exp(-tau/taui)) );
end

