x = 0.01:0.01:100;
y = zeros(1,length(x));
for i=1:length(x)
    y(i) = (1+x(i)/2)^(1/x(i));
end

plot(x,y)
xlabel('x')
ylabel('y')