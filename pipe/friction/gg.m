% Sample data
x = [1, 2, 3, 4, 5];
y = [2, 4, 5, 4, 6];

% Perform linear regression with zero intercept
degree = 1;
coefficients = polyfit(x, y, degree, 'Options', [0 0]);
a = coefficients(1);

% Create a scatter plot and overlay the fitted curve
scatter(x, y, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
fitted_y = a * x;

% Overlay the fitted curve without a constant term (zero intercept)
hold on;
plot(x, fitted_y, 'r', 'LineWidth', 2);
hold off;

% Add labels and a legend
xlabel('x');
ylabel('y');
legend('Data', 'Fitted Curve', 'Location', 'Best');

% Display the slope 'a'
disp(['The slope (a) is: ', num2str(a)]);