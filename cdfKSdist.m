% --- Calculate CDF and KS distance from given PDF

function [CDFtheor, CDFemp, dKS] = cdfKSdist(PDFtheor, PDFemp)

CDFtheor = zeros(length(PDFtheor), 1);
CDFemp = zeros(length(PDFemp), 1);
for i = 1:length(PDFtheor)
    CDFtheor(i, 1) = sum(PDFtheor(i:end));
    CDFemp(i, 1) = sum(PDFemp(i:end));
end

dKS = max(abs(CDFtheor - CDFemp));