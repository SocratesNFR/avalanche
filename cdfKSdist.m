% --- Calculate CDF and KS distance from given PDF

function [CDF, cdfEMP, dKS] = cdfKSdist(PDF, pdfEMP)

CDF = zeros(length(PDF), 1);
cdfEMP = zeros(length(pdfEMP), 1);
for i = 1:length(PDF)
    CDF(i, 1) = sum(PDF(i:end));
    cdfEMP(i, 1) = sum(pdfEMP(i:end));
end

dKS = max(abs(CDF - cdfEMP));