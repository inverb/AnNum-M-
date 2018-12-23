#Pracownia z Analizy Numerycznej (M)
#Zadanie P.2.20
#Mateusz Basiak
#Wrocław, 13.12.2018r.

#Funkcja wyliczająca współczynniki dla k=1
function podstawa_rekursji(N,a,b,c,a_,b_,c_,A)
    A_ = Array{Float64}(undef,N+1);
    A_[2] = (b[N] * A[N+1])/b_[1];
    A_[1] = A[N] - (a[N] * A[N+1]) + (a_[1] * A_[2]);
    return A_;
end

#Funkcja wykonująca jedną iterację rekursji
function rekursja(N,a,b,c,a_,b_,c_,A,k,A_2,A_1)
    Awyn = Array{Float64}(undef, N+1);

    Awyn[1] = A[N-k+1] - (a[N-k+1] * A_1[1]) + (a_[1] * b[N-k+1] * A_1[1])/b_[1] + (c_[2] * b[N-k+1] * A_1[2])/b_[2] - (c[N-k+2] * A_2[1]);

    Awyn[k] = (b[N-k+1] * A_1[k-1])/b_[k-1] + (a_[k] * b[N-k+1] * A_1[k])/b_[k] - (a[N-k+1] * A_1[k]);

    Awyn[k+1] = (b[N-k+1] * A_1[k])/b_[k];

    for n in 2:(k-1)
        Awyn[n] = -(a[N-k+1] * A_1[n]) + (b[N-k+1] * A_1[n-1])/b_[n-1] + (a_[n] * b[N-k+1] * A_1[n])/b_[n];
        Awyn[n] = Awyn[n] + (c_[n+1] * b[N-k+1] * A_1[n+1])/b_[n+1] - (c[N-k+2] * A_2[n]);
    end
    return Awyn;
end

#Główna funkcja wykonująca algorytm
#Wszystkie indeksy są przesunięte o 1 w stosunku do algorytmu w sprawozdaniu,
#ponieważ Julia indeksuje tablice od 1.
function algorytm(N,a,b,c,a_,b_,c_,A)
    #Tablice pomocnicze
    A_1 = Array{Float64}(undef, N+1);
    A_2 = Array{Float64}(undef, N+1);
    Awyn = Array{Float64}(undef, N+1);

    A_2[1] = A[N+1];
    A_1=podstawa_rekursji(N,a,b,c,a_,b_,c_,A);

    for k in 2:N
        Awyn=rekursja(N,a,b,c,a_,b_,c_,A,k,A_2,A_1);
        A_2=A_1;
        A_1=Awyn;
    end
    return Awyn;
end
