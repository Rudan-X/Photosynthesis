function ex_matlab()
x = linspace(0, 4*pi, 50);
y = cos(x);
plot(x,y)
xlabel('x')
ylabel('cos(x)')
title('cos plot')
grid on
save('res.mat','x')
quit
end
