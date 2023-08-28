function y = gof(xtrue,xest)

if sum(size(xtrue) ~= size(xest))
    xest = xest';
end
y = (1 - sqrt(sum((xtrue - xest).^2)/sum((xtrue - mean(xtrue)).^2)))*100;