intJump = 1:10:120;
distStopBins = constVariables.distStop;
for ii = 1:(length(intJump)-1)
    distStopBins(constVariables.distStop>=intJump(ii) & constVariables.distStop<intJump(ii+1))=ii;
end

intJumpY = 1:10:240;
cyclicYBins = constVariables.cyclicY;
for ii = 1:(length(intJumpY)-1)
    cyclicYBins(constVariables.cyclicY>=intJumpY(ii) & constVariables.cyclicY<intJumpY(ii+1))=ii;
end

neuronRate = spkMat(28,:);
avgRatedistStop = [];
for ii = 1:(length(intJump)-1)
    avgRatedistStop(ii) = sum(neuronRate(distStopBins==ii));
end

avgRateY = [];
for ii = 1:(length(intJumpY)-1)
    avgRateY(ii) = sum(neuronRate(cyclicYBins==ii));
end

figure
subplot(1,2,1)
plot(avgRatedistStop)

subplot(1,2,2)
plot(avgRateY)