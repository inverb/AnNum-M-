#Pracownia z Analizy Numerycznej (M)
#Zadanie P.3.5
#Mateusz Basiak
#Wrocław, 13.01.2019r.

#Funkcja implementująca wzór Adamsa-Bashfortha
function Adams_Bashforth(h, F0, F1, F2, F3)
    wyn = Float64(55) * F0 - Float64(59) * F1 + Float64(37) * F2 - Float64(9) * F3;
    wyn = (wyn * h) / Float64(24);
    return wyn;
end

#Funkcja implementująca wzór Adamsa-Moulthona
function  Adams_Moulton(h, F_1, F0, F1, F2)
    wyn = Float64(9) * F_1 + Float64(19) * F0 - Float64(5) * F1 + F2;
    wyn = (wyn * h) / Float64(24);
    return wyn;
end

#Funkcja obliczająca wartość g w kolejnym punkcie
function next_point(F, h, s_1, F0, F1, F2, F3, g_0)
    wyn = Adams_Bashforth(h, F0, F1, F2, F3);
    F_1 = F(s_1, g_0 + wyn);
    wyn = Adams_Moulton(h, F_1, F0, F1, F2);
    return g_0 + wyn;
end

#Funkcja pomocnicza obliczająca wartość funkcji g w punkcie metodą Rungego-Kutty
function Runge_Kutta(F, s_0, g_0, kh)
    G1 = kh * F(s_0, g_0);
    G2 = kh * F(s_0 + kh/2, g_0 + G1/2);
    G3 = kh * F(s_0 + kh/2, g_0 + G2/2);
    G4 = kh * F(s_0 + kh, g_0 + G3);
    wyn = g_0 + (G1 + 2*G2 + 2*G3 + G4)/6;
    return wyn;
end

#Główny algorytm
function algorithm(F, s_0, g_0, l, n)
    h = l/n;

    #Tablica punktów równoodległych
    S = Array{Float64}(undef, n+1);
    S[1] = s_0;
    for i in 2:(n+1)
        S[i] = s_0 + l*(i-1)/n;
    end

    #Tablica wartości funkcji g w punktach S
    G = Array{Float64}(undef, n+1);
    G[1] = g_0;

    #Tablica wartości funkcji F w kolejnych punktach
    Ftab = Array{Float64}(undef, n+1);
    Ftab[1] = F(s_0, g_0);

    #Obliczanie trzech pierwszych wartości g metodą Rungego-Kutty
    for i in 2:4
        G[i] = Runge_Kutta(F, s_0, g_0, (i-1)*h);
        Ftab[i] = F(S[i], G[i]);
    end
    #Główna pętla algorytmu
    for i in 5:(n+1)
        G[i] = next_point(F, h, S[i], Ftab[i-1], Ftab[i-2], Ftab[i-3], Ftab[i-4], G[i-1]);
        Ftab[i] = F(S[i], G[i]);
    end

    return G;
end
