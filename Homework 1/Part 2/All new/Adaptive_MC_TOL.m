function [solution] = Adaptive_MC_TOL(TOL, randomfield, param, start_M, start_I)

% Check if I is even number

I = 2*start_I;  % Double initial I because we work with half I's
M = start_M;
runn_tol = 2*TOL;
N = 10;
x_hat = linspace(0, 1, N + 1);

while runn_tol > TOL
    rng(1);
    Y = normrnd(0, 1, [N M]);
    x = linspace(0, 1, I+1);
    h = 1/I;
    x_half = zeros([I+1 1]);
    a = zeros([M I+1]);

    for k=1:I
        x_half(k) = (x(k+1) + x(k))/2;
    end

    if randomfield == 1
        for sample=1:M
            ind = zeros([length(x) N]);
            
            % Calculate for every node in the I mesh
            for l=1:length(x)
                % Check for every node in I mesh if its in [x^(k-1), x^(k)]
                for r=2:N+1
                    if (x(l) >= x_hat(r-1) & x(l) <= x_hat(r))
                        ind(l, r-1) = Y(r-1, sample);
                    else
                        ind(l, r-1) = 0;
                    end
                end
            end
                
            stairs = crop(ind);
            a(sample, :) = 1 + param*stairs;
            
            F = zeros([I+1 1]);
            for j=1:I+1
                F(j) = 4*pi^2*cos(2*pi*x_half(j));
            end
            u_h = zeros([M I+1]);
        
            for sample=1:M
                A = zeros([I+1 I+1]);
                for j=2:I
                    A(j, j-1) = -a(sample, j-1)/(h^2);
                    A(j, j) = (a(sample, j-1) + a(sample, j+1))/(h^2);
                    A(j, j+1) = -a(sample, j+1)/(h^2);
                end
                u_h(sample, :) = A(2:I, :)\F(2:I);
            end
            QuantOfInt = zeros([M 1]);
            QuantOfInt2 = zeros([M 1]);
            for sample=1:M
                QuantOfInt(sample) = h*sum(u_h(sample, :));
                QuantOfInt2(sample) = 2*h*sum(u_h(sample, 1:2:end));
            end
        end
    end

    if randomfield == 2
        a = zeros([M I+1]);
        a_coarse = zeros([M floor(I/2)+1]);
        u_h = zeros([M I+1]);
        u_h2 = zeros([M floor(I/2)+1]);
        A = zeros([I+1 I+1]);
        A2 = zeros([floor(I/2)+1 floor(I/2)+1]);
        F_coarse = zeros([floor(I/2)+1 1]);

        C = Matern_specialcases_only(x_half, param);
        [EVec, EVal] = eig(C);
        SV = zeros([I 1]);
        EVec = flip(EVec, 2);
        for k=1:I
            SV(k) = sqrt(EVal(I+1-k, I+1-k));
        end
        for sample=1:M
            inner = zeros([length(x) N]);
            for k=1:N
                inner(:, k) = SV(k)*Y(k).*EVec(:, k);
            end
            kappa = sum(inner, 2);
            a(sample, :) = exp(kappa);
    
            for i=1:I/2+1
                if i == 1
                    a_coarse(sample, i) = mean(a(sample, 1:2));
                elseif i == I/2+1
                    a_coarse(sample, i) = mean(a(sample, end-1:end));
                else
                    a_coarse(sample, i) =  mean(a(sample, 2*(i)-2:2*(i)));
                end
            end
    
            for j=1:I+1
                F(j) = 4*pi^2*cos(2*pi*x(j));
            end
            for j=1:floor(I/2)+1
                F_coarse(j) = 4*pi^2*cos(2*pi*x_half(j));
            end
            for j=2:I
                A(j, j-1) = -a(sample, j-1)/(h^2);
                A(j, j) = (a(sample, j-1) + a(sample, j+1))/(h^2);
                A(j, j+1) = -a(sample, j+1)/(h^2);
            end
            
            h2 = 1/(I/2);
            for j=2:floor(I/2)
                A2(j, j-1) = -a_coarse(sample, j-1)/(h2^2);
                A2(j, j) = (a_coarse(sample, j-1) + a_coarse(sample, j+1))/(h2^2);
                A2(j, j+1) = -a_coarse(sample, j+1)/(h2^2);
            end
            u_h(sample, :) = A(2:I, :)\F(2:I);
            u_h2(sample, :) = A2(2:floor(I/2), :)\F_coarse(2:floor(I/2));
            
            QuantOfInt(sample) = h*sum(u_h(sample, :));
            QuantOfInt2(sample) = h2*sum(u_h2(sample, :));
        end
    end

    % MC estimate after M samples using I+1 gridpoints
    MC_est = 1/M*sum(QuantOfInt2)
    MC_var = 1/(M-1)*sum((QuantOfInt2 - MC_est).^2)

    runn_tol = TOL/2;
end


end