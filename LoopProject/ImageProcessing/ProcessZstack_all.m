function ProcessZstack_all(channel2, registration,data, rescale,rescale_factor, regstack);
global path
path = zStack_avgPerPlane(channel2, registration);
RegisterStackToFunctionalImaging_v2(data, rescale, rescale_factor, regstack);
PlotROIsMasks_2(data, rescale, rescale_factor);
end