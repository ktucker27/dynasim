function [a, b, w, eloc_evs] = nqs_learn(a, b, w, eloc, wave, eta, num_iterations, num_samps, num_steps)
% nqs_learn Learn the parameters (a, b, w) of a neural quantum state (NQS)
%           for a given local energy operator eloc

eloc_evs = zeros(num_iterations,1);
for i=1:num_iterations
    if mod(i,10) == 0
        disp(['i = ', num2str(i)]);
    end
    
    % Compute the local energy expectation
    eloc_evs(i,1) = nqs_ev(eloc, wave, a, b, w, num_samps, num_steps);
    
    % Compute the statistical average of the G operator for all parameters
    geva = 0;
    gevb = 0;
    gevw = 0;
    for j=1:num_samps
        sz = nqs_sample(wave, a, b, w, num_steps);
        theta = b + w*sz;
        elocval = feval(eloc, a, b, w, sz);
        geva = geva + 2*real((elocval - eloc_evs(i,1))*sz);
        gevb = gevb + 2*real((elocval - eloc_evs(i,1))*conj(tanh(theta)));
        gevw = gevw + 2*real((elocval - eloc_evs(i,1))*conj(tanh(theta))*sz');
    end
    geva = geva/num_samps;
    gevb = gevb/num_samps;
    gevw = gevw/num_samps;
    
    % Update parameters according to gradient descent
    a = a - eta/sqrt(i+1)*geva;
    b = b - eta/sqrt(i+1)*gevb;
    w = w - eta/sqrt(i+1)*gevw;
    
    if isnan(norm(a))
        disp(['Received NaN at i = ', num2str(i)]);
        break;
    end
end