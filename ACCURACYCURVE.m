x = (1:30);
e1e2 = [.909 .818 .886 .841 .841 .818 .909 .864 .795 .818 .909 .818 .750 .795 .773 .795 .864 .818 .886 .773 .773 .864 .773 .841 .909 .886 .909 .864 .750 .909];
e1e3 = [.568 .705 .818 .750 .750 .727 .705 .773 .750 .705 .727 .818 .750 .795 .773 .727 .818 .818 .750 .795 .773 .773 .750 .727 .705 .750 .818 .818 .773 .705];
e2e3 = [.524 .422 .435 .537 .544 .510 .531 .544 .571 .517 .464 .437 .538 .571 .591 .523 .537 .527 .481 .473 .574 .568 .591 .544 .544 .537 .574 .473 .444 .534];

figure(1);
plot(x, e1e2*100, x, e1e3*100, x, e2e3*100)
title('Accuracy for Electrode Configurations')
xlabel('Classifier Trial #'), xlim([1 30])
ylabel('Accuracy (%)'), ylim([0 100])
legend('e1 & e2', 'e1 & e3', 'e2 & e3')

fprintf('The mean accuracy was %s%% and the highest accuracy was %s%% for electrode configuration e1 e2.\n', num2str(mean(e1e2)*100), num2str(max(e1e2)*100));

a = [.694 .714 .701 .653 .680 .633 .728 .701 .731 .714 .673 .653 .701 .727 .727 .701 .686 .673 .633 .731 .731 .714 .668 .641 .694 .644 .684 .714 .727 .633];
b = [.680 .639 .573 .573 .682 .674 .627 .636 .573 .659 .705 .611 .644 .626 .626 .654 .636 .611 .705 .682 .591 .591 .611 .611 .636 .653 .622 .682 .621 .705];
c = [.455 .523 .500 .500 .500 .591 .523 .477 .568 .523 .591 .568 .636 .455 .409 .591 .523 .423 .455 .568 .500 .477 .455 .568 .523 .523 .455 .636 .500 .500];

figure(2);
plot(x, a*100, x, b*100, x, c*100)
title('Accuracy')
xlabel('Classifier Trial #'), xlim([1 30])
ylabel('Accuracy (%)'), ylim([0 100])
legend('Spectral Power, Spectral Power Difference, DASDV, Max, Mean, Variance', '4 Features & Correlated (Overlapping) FFT Features', '4 Features: DASDV, RMS, Max, Variance')

fprintf('The mean accuracy was %s%% and the highest accuracy was %s%% for electrode configuration e1 e2.\n', num2str(mean(a)*100), num2str(max(a)*100));