'''Varianta 2 de implementare a metodei indexului'''

import math
import random
import numpy as np


def L(x):
    l = math.log(x)
    return np.exp(math.sqrt(l*math.log(l)))

def Ln(x, n = 1/3):
    l = math.log(x)
    return np.exp(l**n*math.log(l)**(1-n))

def cmmdc(a, b):
    # Functie ce se foloseste de Algoritmul lui Euclid pentru a calcula cel mai mare divizor comun dintre a si b
    while b != 0:
        a, b = b, a % b
    return a

def verificare_prim(p):
    if p <= 1:
        return False
    for i in range(2, int(math.sqrt(p)) + 1):
        if p % i == 0:
            return False
    return True

def Factorizare(a):
    # Functia calculeaza descompunerea in factori primi a numarului a si o returneaza sub forma unui dictionar ce are ca si chei intregii primi si puterile la care se gasesc acestia in descompunere, la valorile dictionarului
    d = {} 
    if verificare_prim(a): # Cazul in care a este prim
        d[a] = 1
        return d
    
    # Cazul in care a nu este prim:
    # Izolam cazul in care s-ar gasi 2 in descompunerea lui a, pentru ca mai apoi sa nu mai luam in considerare numerele pare
    if a % 2 == 0:
        # a contine 2 in descompunerea sa in factori primi
        # calculam la ce putere se regaseste 2
        putere = 1
        nr = a // 2
        r = 0
        while r == 0:
            r = nr % 2
            nr = nr // 2
            if r == 0:
                putere += 1
            else:
                r = 1
        d[2] = putere
    
    # Cautam alti factori primi din descompunerea lui a
    for i in range(3, a + 1):
        # Eliminam factorii pari, deoarece deja i-am luat in considerare
        if i % 2 == 0:
            continue  # Sare peste restul codului din bucla si trece la urmatoarea iteratie
            
        # Ne uitam numai la i prim
        if verificare_prim(i):
            #print(i)
            if a % i == 0:
                # a contine i in descompunerea sa in factori primi
                # calculam la ce putere se regaseste i
                putere = 1
                nr = a // i
                r = 0
                while r == 0:
                    r = nr % i
                    nr = nr // i
                    if r == 0:
                        putere += 1
                    else:
                        r = 1
                d[i] = putere
    
    return d


def exponentiere(x, e, p):
    # Functie care imi calculeaza x^e mod p (inlocuieste functia pow)
    # Daca x congruent cu 0 mod p:
    if x % p == 0: 
        if e < 0:
            print(f"Nu exista invers modular pentru {x} mod {p}.")
            return None
        else:
            return 0
    # Daca x diferit de 0 mod p si e este negativ:
    if e < 0:
        # Calculam mai intai inversul modular
        # Verificam daca exista inversul lui x mod p
        if cmmdc(x, p) == 1:
            # x^(-1) = x^(phi(p)-1)
            # Calculam phi(p), unde p este diferit de 1:
            phi = p
            factorizare_p = Factorizare(p)
            for i in factorizare_p.keys():
                phi = phi * (i - 1) // i # scriem asa ca sa ne asiguram ca obtinem un numar intreg
            e_nou = phi - 1
            nr = exponentiere(x, e_nou, p)
            nr = exponentiere(nr, np.abs(e), p)
            return nr
        
        else:
            print(f"Nu exista invers modular pentru {x} mod {p}.")
            return None
    # Daca x diferit de 0 mod p si e este pozitiv:
    r = e % 2
    e = e // 2
    exp = x # x^(2^i)
    nr = 1
    if r == 1:
        nr = exp
        nr = nr % p
        
    while e > 0:
        r = e % 2
        e = e // 2
        exp = exp**2
        exp = exp % p
        if r == 1:
            nr = nr * exp
            nr = nr % p
    return nr


def verificare_elem_primitiv(g, p):
    if g % p == 0 and g % p == 1:
        return False
    else:
        verif = g
        count = 1
        while verif % p != 1:
            verif = verif * g
            count += 1
        if count == p - 1:
            #print(f'ordinul lui g = {g} este {count}')
            #print(f'Elementul g = {g} este element primitiv in grupul F{p}* cu inmultirea')
            return True
        return False
    
def B_smooth(B, gi):
    # Functia verifica daca numarul introdus gi este B-smooth, adica in descompunerea sa se afla numai numerele prime din baza B la diferite puteri
    # Functia returneaza True si descompunerea in factori primi a numarului introdus daca e B-smooth, respectiv False si None, daca numarul nu este B-smooth
    d = {} 
    # Verificam pentru fiecare factor prim p din baza daca se regaseste in descompunerea lui gi si la ce putere
    for p in B:
        if gi % p == 0:
            # gi contine p in descompunerea sa in factori primi
            # calculam la ce putere se regaseste p
            putere = 1
            # De fiecare daca cand gi se imparte exact la p, il impartim pe gi la p, pentru a scoate pe rand numarul p din descompunerea lui gi, pana nu se mai gaseste deloc
            gi = gi // p 
            r = 0
            while r == 0:
                r = gi % p
                if r == 0:
                    putere += 1
                    gi = gi // p
                else:
                    r = 1
            d[p] = putere
        else:
            d[p] = 0
    # Daca la final, dupa toate impartirile, am ramas cu gi = 1, atunci numarul gi este B-smooth, altfel nu 
    if gi == 1:
        return True, d
    else:
        return False, None



def verificare_sistem(A, p):
    # Calculam determinantul matricii sistemului si ne asiguram ca este inversabil mod p
    det_A = int(round(np.linalg.det(A))) # folosim int si round pentru a nu obtine erori
    det_A = int(det_A % p)
    
    if cmmdc(det_A, p) != 1:
        #print('Sistemul nu este compatibil determinat. Cauta alt sistem de congruente.')
        return False

    return True



def invers_modulo(a, p):
    # Calculează inversul lui a mod p
    return exponentiere(int(a), -1, p)


def invers_matrice_modulo(A,p):
    # Imi calculeaza inversului matricei A mod p
    if verificare_sistem(A, p):
        n = len(A)
        A = A % p  # Asiguram ca lucram cu resturi mod p
        I = np.eye(n, dtype=int)  # Matricea identitate
        AI = np.hstack((A, I))  # Formam matricea [A | I]
        for coloana in range(n):
            # Gasim pivotul si il aducem pe diagonala
            linie_pivot = -1
            for linie in range(coloana, n):
                if AI[linie, coloana] % p != 0:
                    linie_pivot = linie
                    pivot = AI[linie, coloana]
                    break

            if linie_pivot < 0:
                print(f"Sistemul nu are soluții unice modulo {p}.")
                return None

            if linie_pivot != coloana:
                AI[[coloana, linie_pivot]] = AI[[linie_pivot, coloana]] # interschimbam liniile pentru a aduce pivotul pe pozitia corecta

            invers_pivot = invers_modulo(pivot, p)
            AI[coloana] = (AI[coloana] * invers_pivot) % p # AI[coloana] este linia pivotului
    
            # Eliminam toate celelalte elemente din coloana curenta
            for linie in range(n):
                if linie != coloana:
                    factor = AI[linie, coloana]
                    AI[linie] = (AI[linie] - factor * AI[coloana]) % p 
    
        return AI[:, n:]  # Partea dreapta este inversa lui A
    else:
        print(f'Matricea A = {A} nu este inversabila.')
        return None

def adaugare_congruenta(Ab, Ab_initial, g, p, B, mod):
    # Functia imi adauga o noua linie la matricea sistemului de congruente atat in matricea sistemului modulo factorul prim curent, cat si in matricea sistemului de congruente de la care am plecat
    nr_gasite = 0
    contor = 0 # numara cate numere au fost testate pentru a gasi unul B-smooth
    while nr_gasite < 1:
        contor += 1
        ti = random.randint(0, p-2) # alegem un exponent random intre 0 si p-2 inclusiv
        gi = exponentiere(g, ti, p)
        valoare_adevar, factorizare = B_smooth(B, gi)
        if valoare_adevar:
            nr_gasite += 1
            exponenti = np.zeros(len(B)).astype(int)
            for i in factorizare.keys():
                for prim in B:
                    if i == prim:
                        exponenti[B.index(prim)] = factorizare[i]
    linie_noua = np.append(exponenti, ti)
    Ab_initial = np.vstack([Ab_initial, linie_noua])
    linie_noua - linie_noua % mod
    Ab = np.vstack([Ab, linie_noua])
    
    return Ab, Ab_initial, contor

def ridicare_solutie(x, Ab, dimensiune, mod, putere):
    x = x.T
    p = mod
    # Ridicam solutia mod p^n (p = mod, n = putere):
    # Extragem liniile si coloanele ce contin pivoti din matricea Ab
    A = Ab[:dimensiune, :dimensiune]
    b = Ab[:dimensiune, -1]
    # Calculam inversa matricii sistemului modulo p
    A_inv = invers_matrice_modulo(A, p)

    for i in range(2, putere+1):
        F = A @ x
        F = F - b
        delta = A_inv @ F
        x = x - delta
        x = x % (p**i)       
    return x.flatten() # pentra a returna un vector clasic in ambele cazuri


def cautare_C_congruente(B, g, p):
    # Functia imi gaseste initial |B| congruente, dupa care adauga treptat cate una, pana obtinem un sistem de congruente cu solutie unica modulo fiecare factor prim in parte din descompunerea lui p-1
    nr_gasite = 0
    contor = 0 # imi numara cate numere prime este necasar sa testeze pentru a gasi indeajuns de multe numere B-smooth cat sa rezolve sistemul de congruente
    nr_b_smooth = []
    lista_exponenti = []
    lista_ti = []
    # Generam mai intai |B| congruente
    while nr_gasite < len(B): 
        contor += 1
        ti = random.randint(0, p-2) # alegem un exponent random intre 0 si p-2 inclusiv
        gi = exponentiere(g, ti, p)
        valoare_adevar, factorizare = B_smooth(B, gi)
        if valoare_adevar:
            nr_b_smooth.append((gi, factorizare))
            nr_gasite += 1
            exponenti = np.zeros(len(B)).astype(int)
            for i in factorizare.keys():
                for prim in B:
                    if i == prim:
                        exponenti[B.index(prim)] = factorizare[i]
            
            lista_exponenti.append(exponenti)
            lista_ti.append(ti)

    #Cream matricea sistemului de congruente
    A = np.zeros((nr_gasite, len(B))).astype(int)
    b = np.zeros((nr_gasite,1)).astype(int)
    for i in range(nr_gasite):
        b[i,0] = lista_ti[i]
        for j in range(len(B)):
            A[i,j] = lista_exponenti[i][j]
    Ab = np.hstack((A, b))
    
    # Aplicam metoda Gauss Jordan pentru fiecare factor prim din descompunearea lui p-1
    fact = Factorizare(p-1) # dictionar ce contine fiecare factor prim din descompunere si puterea la care se gaseste
    lista_solutii = []
    for prim in fact.keys():
        mod = prim
        putere = fact[prim]
        copie_Ab = np.copy(Ab) # O sa lucram local cu o copie a matricei Ab
        copie_Ab = copie_Ab % mod
        # Vrem sa facem pivoti pe fiecare coloana din A
        coloana = 0
        n = len(B)
        while coloana < n:
            linie_pivot = -1
            for linie in range(coloana, n):
                if copie_Ab[linie, coloana] % mod != 0:
                    linie_pivot = linie
                    pivot = copie_Ab[linie, coloana]
                    break

            if linie_pivot < 0:
                # Sistemul nu are solutie unica modulo factorul prim curent, asadar mai trebuie sa adaug o congruenta
                am_gasit_pivot = False
                while am_gasit_pivot == False:
                    copie_Ab, Ab, contor_plus = adaugare_congruenta(copie_Ab, Ab, g, p, B, mod)
                    contor = contor + contor_plus
                    
                    # Daca nu exista pivot pe prima coloana (adica coloana == 0)
                    if coloana == 0:
                        if copie_Ab[-1, 0] % mod != 0:
                            pivot = copie_Ab[-1, 0] % mod
                            am_gasit_pivot = True
                            linie_pivot = len(copie_Ab) - 1 # ultima linie
                    else:
                        # Coloana curenta nu este prima
                        # Transformam noua linie prin metoda Gauss Jordan
                        for i in range(coloana):
                            factor = copie_Ab[-1, i]
                            copie_Ab[-1] = (copie_Ab[-1] - factor * copie_Ab[i]) % mod
                        # Verificam daca avem pivot pe coloana curenta pe ultima linie din matricea copie_Ab
                        if copie_Ab[-1,coloana] % mod != 0:
                            pivot = copie_Ab[-1, coloana] % mod
                            am_gasit_pivot = True
                            linie_pivot = len(copie_Ab) - 1 # ultima linie
                            
            if linie_pivot != coloana:
                copie_Ab[[coloana, linie_pivot]] = copie_Ab[[linie_pivot, coloana]] # interschimbam liniile pentru a aduce pivotul pe pozitia corecta
                Ab[[coloana, linie_pivot]] = Ab[[linie_pivot, coloana]]

            invers_pivot = invers_modulo(pivot, mod)
            copie_Ab[coloana] = (copie_Ab[coloana] * invers_pivot) % mod # copie_Ab[coloana] este linia pivotului

            # Eliminam valorile de sub si de deasupra pivotului din coloana curenta
            for linie in range(len(copie_Ab)):
                if linie != coloana:
                    factor = copie_Ab[linie, coloana]
                    copie_Ab[linie] = (copie_Ab[linie] - factor * copie_Ab[coloana]) % mod
            coloana += 1
        # Am obtinut numarul necesar de pivoti (|B|), asadar putem extrage solutiile modulo factorul prim curent
        # Solutiile se afla pe ultima coloana a matricei copie_Ab
        sol = copie_Ab[:n, -1]
        # Daca factorul prim curent se afla la o putere > 1 in descompunerea lui p-1 trebuie sa ridicam solutia modulo factorul prim la puterea respectiva
        if putere > 1:
            sol = ridicare_solutie(sol, Ab, len(B), mod, putere)
        lista_solutii.append((sol, mod**putere))

    return lista_solutii, contor 


def LCR(lista_solutii, nr_necunoscute, p):
    x_vector = [] # vectorul de solutii
    # Extragen datele din lista
    for j in range(nr_necunoscute):
        # Calculam pentru fiecare necunoscuta in parte
        x = 0 # necunoscuta pe care vrem sa o aflam
        lista_nr_prime = []
        lista_intregi = []
        for solutie in lista_solutii:
            lista_nr_prime.append(solutie[1])
            lista_intregi.append(solutie[0][j])
        for i in range(len(lista_intregi)):
            pi = lista_nr_prime[i]
            ai = lista_intregi[i]
            bi = p // pi
            inv_bi = invers_modulo(bi, pi) # inversul lui bi mod pi
            x = x + ai * bi * inv_bi
        x = x % p
        x_vector.append(int(x))
    return x_vector
            
def etapa1(p, g, h, B):
    # Cautam C>= |B| numere B-smooth de forma g^t mod p si alcatuim sistemul de congruente
    # Functia imi gaseste initial |B| congruente, dupa care adauga treptat cate una, pana obtinem un sistem de congruente cu solutie unica modulo fiecare factor prim in parte din descompunerea lui p-1
    lista_solutii, contor = cautare_C_congruente(B, g, p)
    # lista_solutii este o lista ce contine perechi de solutii si modulul la care se regasesc
    
    # Ne folosim de LCR pentru a ajunge la solutia unica a sistemului
    nr_necunoscute = len(B)
    solutie_sistem = LCR(lista_solutii, nr_necunoscute, p-1)
    #print((A@solutie_sistem)%(p-1)) # verificare, trb sa dea b
    # Asadar, am aflat logaritmii elementelor din baza, marcand astfel finalul etapei 1.
    logaritmi_baza = solutie_sistem
    #print('Logaritmii bazei sunt: ', logaritmi_baza)   
    return logaritmi_baza, contor         



def Metoda_Indexului(p, g, h):
    # Verificam daca numerele p si g respecta conditiile metodei indexului 
    
    if not verificare_prim(p):
        print(f'p = {p} nu este prim.')
        return None
    
    if not verificare_elem_primitiv(g, p):
        print(f'Elementul g = {g} nu este primitiv.')
        return None
    
    # Transformam numerele g si h modulo p
    g = g % p
    h = h % p
    
    # Cream baza de numere prime
    
    b = round(L(p)**(1/math.sqrt(2)))
    B = []

    for i in range(b + 1):
            if verificare_prim(i):
                B.append(i)
                
    print("Baza: ", B)

    
    # ETAPA 1: Aflam log_g(pi), unde pi sunt elementele din baza B
    logaritmi_baza, contor = etapa1(p, g, h, B)

    # ETAPA 2: Aflarea logaritmului cautat log_g(h)
    # Cautam un numar B-smooth de forma h*g^s mod p, pentru s random, 1 <= s <= p-2
    count = 0
    while count < 1:
        s = random.randint(1, p-2) # alegem un exponent random intre 0 si p-2 inclusiv
        gs = exponentiere(g, s, p) # functie folosita pentru exponentiere modulara
        nr = (h * gs) % p
        valoare_adevar, factorizare = B_smooth(B, nr)
        if valoare_adevar:
            # Am gasit un numar B-smooth si avem si descompunerea sa
            count = 1
            exponenti = np.zeros(len(B)).astype(int)
            for i in factorizare.keys():
                for prim in B:
                    if i == prim:
                        exponenti[B.index(prim)] = factorizare[i]

    # Logaritmam expresia si extragem log_g(h)
    log = -s
    for i in range(len(B)):
        log = log + exponenti[i] * logaritmi_baza[i]
    log = log % (p-1)
    print("Solutia din Zp-1 a problemei logaritmului discret g^x = h mod p este: ", log)
    return log, contor, B



    
   
#g = int(input("Introduceti un numar intreg g : "))
#h = int(input("Introduceti un numar intreg h : "))
#p = int(input("Introduceti un numar intreg p : "))

p = 611953
g = 7
h = 78902


lista_contor = []
for i in range(1):
    log, contor, B = Metoda_Indexului(p, g, h)
    lista_contor.append(contor)
    print(f'Pentru rezolvarea sistemului de congruente cu {len(B)} necunoscute au fost testate {contor} numere.')

print("Numarul de teste efectuate la fiecare iteratie: ",lista_contor)
print("Numarul mediu de teste efectuate: ", sum(lista_contor)/len(lista_contor))
verificare = exponentiere(g, log, p)

print(verificare) 

#%%

def calcul_ordine(p):
    for i in range(p):
        if i % p == 0:
            print(f'ordinul lui {i} este infinit')
        else:
            verif = i
            count = 1
            while verif % p != 1:
                verif = verif * i
                count += 1
            print(f'ordinul lui {i} este {count}')
            

def gasire_element_primitiv(p):
    for i in range(p):
        if i % p != 0:
            if verificare_elem_primitiv(i, p):
                return i
   
#print(gasire_element_primitiv(611953))

# Nr. de generatoare ale unui grup ciclic este phi(phi(n)), pt Zn* grup ciclic












