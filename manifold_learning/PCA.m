function [Y] = PCA (X,d)
    opts.disp = 0;
    [Y,~] = eigs(X*X',d,'lm',opts);
end