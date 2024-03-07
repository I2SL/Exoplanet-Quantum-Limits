function [x,y] = UniformSpiralSampling(num_samples,num_turns,r0,th0)
% uniform spiral sampling for a spiral that starts at r0,th0, goes through
% num_turns, and ends at radial distance 1

% expansion rate
a = (1-r0)/(2*pi*num_turns);

% spacing between samples
c = 1/(num_samples-1)*(2*pi*num_turns)*(a/2 * (2*pi*num_turns) + r0);

% Generate angular samples
th = 0;
for i = 1:(num_samples-1)
    aa = a/2;
    bb = r0;
    cc = -a/2*th(i)^2 - r0*th(i) - c;
    th(i+1,1) = (-bb + sqrt(bb^2 - 4 * aa * cc ))/(2*aa);
end

r = a*th + r0;
th = th + th0;
[x,y] = pol2cart(th,r);

%{
scatter(x,y,'filled','b')

thp = linspace(min(th),max(th),10000);
rp = a*thp +b;
rp = rp/max(rp);
[xp,yp] = pol2cart(thp,rp);
hold on
scatter(0,0,'k','filled')
plot(xp,yp,'k')
hold off
xlim([-1,1])
ylim([-1,1])
axis square
%}
end