function softPts = computeSoftCurve(q, curves)
    n       = numel(q)/2;
    samples = 100;
    T       = eye(4);
    softPts = zeros(3, samples*n);
    skelpts = zeros(3, n + 1);
    thetas = q(2:2:end);
    phis   = q(1:2:end);

    L = curves.L;
   

    for i = 1:n
        Xi  = position_vector(thetas(i), samples, phis(i), L, T);
        idx = (i-1)*samples + (1:samples);
        softPts(:,idx) = Xi;
        T = T * Transform(thetas(i), phis(i), L);
        skelpts(:,i+1) = softPts(:,samples*i);
    end


    
end
