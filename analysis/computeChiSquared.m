function chiSquare = computeChiSquared(est, true, plotOn, figNum)

%counting_variance = cov(est);
counting_variance = ones(256,256);

chiSquareMap = 1 ./ counting_variance .* (est - true).^2;

chiSquare = sum(sum( chiSquareMap ));

if plotOn == 1
    figure(figNum); subplot(1,2,2);
    contourf(chiSquareMap);
    colorbar(); title(sprintf("\\chi^2 Map, \\chi^2=%.3e", chiSquare));
end

end