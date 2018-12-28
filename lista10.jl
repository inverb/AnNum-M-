#Mateusz Basiak
#Lista 10

using Printf;

#Obilczanie rekurencyjne jednego współczynnika wielomianu
function coeff(a, b, i, n, k, w)
    if(a == n+1)
        if(b == 0)
            return w;
        else
            return Float64(0.0);
        end
    end
    if(a == k)
        return coeff(a+Float64(1.0), b, i, n, k, w);
    end
    wyn = coeff(a+1, b, i, n, k, w * Float64(-a));
    if(b > 0)
        wyn = wyn + coeff(a+1, b-1, i, n, k, w);
    end
    return wyn;
end

#Obliczanie całki
function integral(n, k)
    wyn = Float64(0.0);
    for i in 0:n
        N = Float64(1.0);
        for j in 1:(i+1)
            N = N * n;
        end
        N = N / (i+1);
        M = coeff(Float64(0.0), i, i, n, k, Float64(1.0));
        N = N * M;
        wyn = wyn + N;
    end
    return wyn;
end

function A_k(h, n, k)
    sign = Float64(1.0);
    for i in 1:(n-k)
        sign = sign * Float64(-1.0);
    end
    k1 = k2 = Float64(1.0);
    for i in 1:k
        k1 = k1 * Float64(i);
    end
    for i in 1:(n-k)
        k2 = k2 * Float64(i);
    end
    k1 = k1 * k2;
    h = (sign * h)/k1;
    I = integral(n,k);
    h = h * I;
    #@printf("%f : %f\n", k ,h);
    return h;
end

#Liczy sumy B_k
function suma(a, b, n)
    s = Float64(0.0);
    h = (b-a)/n;
    for i in 0:n
        B = A_k(h,n,i);
        B = B/(b-a);
        if(B < 0.0)
            @printf("Less than zero: %d %d %f\n", n, i, B);
        end
        B = abs(B);
        s = s+B;
    end
    @printf("Value for n = %d is: %f\n", n, s);
end

#Liczy kwadratury
function quadrature(f, n, a, b)
    h = (b-a)/n;
    q = 0.0;
    for k in 0:n
        A = A_k(h, n, k);
        x = f(a + k * h);
        q = q + A * x;
    end
    #println(q);
    return q;
end

#Zadanie 4
#for i in 1:20
#    suma(0, 1, i);
#end

#zadanie 5
function foo(x)
    wyn = x*x + 1.0;
    wyn = 1 / wyn;
    return wyn;
end

#@printf("Exact value is: %f\n", 2 * atan(4));
#for i in 1:5
#    q = quadrature(foo, 2*i, -4, 4);
#    @printf("Value of quadrature for n = %d is: %f\n", 2*i, q);
#end
