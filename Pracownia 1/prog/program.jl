#Pracownia z Analizy Numerycznej (M)
#Zadanie P.1.11
#Mateusz Basiak
#Wrocław, 28.10.2018r.

#Funkcja znajdująca n-ty element szeregu Taylora
#funkcji sin(x) w punkcie x
function element_szeregu_Taylora(x,n)
    licznik = Float64(1.0);
    mianownik = Float64(1.0);
    for i in 1:n
        licznik = licznik*x;
        mianownik = mianownik*Float64(i);
    end
    return licznik/mianownik
end

#Funkcja obliczająca sin(x) za pomocą
#szeregu Taylora, obliczając każdy wyraz osobno.
function sin_Taylor1(x)
    wynik = x;
    wynik1 = wynik;
    j=Float64(3.0);
    while(true)
        t=element_szeregu_Taylora(x,j); #znajdowanie i-tego wyrazu szeregu Taylora
        if(j%4==3) wynik1-=t;
        else wynik1+=t; end
        j=j+2;
        if(wynik1 == wynik)
            return wynik;
        end
        wynik=wynik1;
    end
end

#Funkcja obliczająca cos(x) za pomocą
#szeregu Taylora, obliczając każdy wyraz osobno.
function cos_Taylor1(x)
    wynik = Float64(1.0);
    j=2;
    wynik1 = wynik;
    while(true)
        t=element_szeregu_Taylora(x,j); #znajdowanie i-tego wyrazu szeregu Taylora
        if(j%4==2) wynik1-=t;
        else wynik1+=t; end
        j=j+2;
        if(wynik1 == wynik)
            return wynik;
        end
        wynik=wynik1;
    end
    return wynik;
end

#Funkcja obliczająca sin(x) za pomocą 1000 pierwszych
#wyrazów szeregu Taylora, korzystając z przenawiasowania działań.
function sin_Taylor2(x)
    wynik = Float64(1.0);
    i = 1000;
    while i>0
        wynik = wynik * x * x;
        mianownik = (Float64(2.0) * Float64(i)) * (Float64(2.0) * Float64(i) + Float64(1.0)); #(2i)*(2i+1)
        wynik = wynik / mianownik;
        wynik = Float64(1.0) - wynik;
        i = i - 1;
    end
    wynik = x * wynik;
    return wynik;
end

#Funkcja obliczająca cos(x) za pomocą 1000 pierwszych
#wyrazów szeregu Taylora, korzystając z przenawiasowania działań.
function cos_Taylor2(x)
    wynik = Float64(1.0);
    i = 1000;
    while i>0
        wynik = wynik * x * x;
        mianownik = (Float64(2.0) * Float64(i) - Float64(1.0)) * (Float64(2.0) * Float64(i)); #(2i-1)*(2i)
        wynik = wynik / mianownik;
        wynik = Float64(1.0) - wynik;
        i = i - 1;
    end
    return wynik;
end

#Funkcja obliczająca sin(x) wyliczając rekurencyjnie
#licznik i mianownik p_n i q_n korzystając z własności ułamka łańcuchowego.
function sin_fraction(x)
    p0=Float64(0.0);
    p1=x;
    q0=Float64(1.0);
    q1=Float64(1.0);
    i=2;
    while(p0/q0 != p1/q1)
        a = (Float64(2.0)*i-Float64(4.0))*(Float64(2.0)*i-Float64(3.0))*x*x;
        if(i == 2) a = x*x; end
        b = (Float64(2.0)*i-Float64(2.0))*(Float64(2.0)*i-Float64(1.0))-x*x;
        p = p1*b + p0*a;
        q = q1*b + q0*a;
        p0=p1;
        p1=p;
        q0=q1;
        q1=q;
        i = i + 1;
    end
    return p1/q1;
end

#Funkcja obliczająca sin(x) za pomocą iloczynu
function sin_mult(x)
    wynik = Float64(1.0);
    wynik1 = wynik;
    i=1;
    while(true)
        czynnik = x*x/(Float64(pi)*Float64(pi)*Float64(i)*Float64(i));
        czynnik = Float64(1.0) - czynnik;
        wynik1 = wynik1 * czynnik;
        if(wynik1 == wynik)
            return wynik * x;
        end
        wynik=wynik1;
        i=i+1;
    end
end

#Funkcja obliczająca cos(x) za pomocą iloczynu
function cos_mult(x)
    wynik = Float64(1.0);
    wynik1 = wynik;
    i=1;
    while(true)
        iprim = i - Float64(0.5);
        czynnik = x*x/(Float64(pi)*Float64(pi)*Float64(iprim)*Float64(iprim));
        czynnik = Float64(1.0) - czynnik;
        wynik1 = wynik1 * czynnik;
        if(wynik1 == wynik)
            return wynik;
        end
        wynik=wynik1;
        i = i+1;
    end
end
