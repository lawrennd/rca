function gridSearch(minpar, maxpar, samples, grad_limit)
% minpar = 1e-4;
% maxpar = totalvar;
% samples = 10;
% grad_limit = 1e-1;

grad_converged = false;
while ~grad_converged
    par = logspace(log10(minpar), log10(maxpar), samples);
    history_L = Ly_f(par, WWt_hat, S_f, n);
    history_GradL = Grad_Ly_f_sigma2_n(par, WWt_hat, S_f, n);
    
    semilogx(par, history_L, 'b.'), hold on
    semilogx(par, history_GradL, 'r.'), hold on
    
    [~,imax] = max(history_L);   maxLpar = par(imax);
    maxLgrad = Grad_Ly_f_sigma2_n(maxLpar, WWt_hat, S_f, n);
    
    grad_converged = abs(maxLgrad) < grad_limit;
    if ~grad_converged
        if maxLgrad < 0  % Shrink search bounds.
            maxpar = maxLpar;
            posgrad_pars = par( history_GradL(1:imax) > 0 );
            minpar = posgrad_pars(end);
        elseif maxLgrad > 0
            minpar = maxLpar;
            neggrad_pars = par([boolean(zeros(imax-1,1)); history_GradL(imax:end) < 0]);
            maxpar = neggrad_pars(1);
        else
            grad_converged = true;
        end
    end
end