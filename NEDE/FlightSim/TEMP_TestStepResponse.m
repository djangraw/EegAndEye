% TEMP_TestStepResponse.m
%
% Tests whether the impulse response is in fact the derivative of a step
% response and, if so, what the impulse response for our given step
% response looks like.
%
% Created 9/16/14 by DJ.

stepresp = [0 0 0 0 1 2 3 2 1 1 1 1 1];
impresp = [0 diff(stepresp)];
step_recon = conv(ones(1,length(stepresp)),impresp,'full');
step_recon = step_recon(1:length(stepresp));

figure(22);
subplot(3,1,1);
plot(stepresp);
ylabel('step response')

subplot(3,1,2);
plot(impresp);
ylabel('impulse response')

subplot(3,1,3);
plot(step_recon);
ylabel('reconstructed step response')
xlabel('samples')

isequal(stepresp,step_recon)