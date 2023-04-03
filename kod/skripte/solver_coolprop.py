from funkcije_coolprop import brzina_ekspanzije, stanje_gustoca, zokg_brzina, zom_zokg_tlak, zoe_entalpija, \
    stanje_temperatura, temp_stijenke, alpha_unutarnje, eksplicitna, manhattan_rjesenja, citaj_stanje
from velicine import *
from uvjeti import h_visokotlacno, p_visokotlacno, p_niskotlacno, stanje_0
import os
import matplotlib.pyplot as plt

def vanjska_iteracija(trenutak, pretpostavka, vrijednost_proslog_trenutka, podrelaksacija):
    # plt.plot(pretpostavka[3], label=str(broj_iteracija))
    tlak_na_ulazu = pretpostavka[1, 0]
    brzina_na_ulazu = brzina_ekspanzije(p_visokotlacno[trenutak] - tlak_na_ulazu, pretpostavka[3,0])
    entalpija_na_ulazu = h_visokotlacno[trenutak] - brzina_na_ulazu ** 2 / 2
    gustoca_na_ulazu = stanje_gustoca(tlak_na_ulazu, entalpija_na_ulazu)

    brzina = zokg_brzina(pretpostavka[0], brzina_na_ulazu, pretpostavka[1], p_niskotlacno[trenutak], pretpostavka[2],
                         pretpostavka[3], gustoca_na_ulazu, dx, podrelaksacija[0])

    tlak = zom_zokg_tlak(brzina, pretpostavka[1], p_visokotlacno[trenutak], p_niskotlacno[trenutak], pretpostavka[2],
                         pretpostavka[3], vrijednost_proslog_trenutka[3], dx)
    tlak = eksplicitna(tlak, pretpostavka[1], podrelaksacija[1])

    entalpija = zoe_entalpija(brzina, vrijednost_proslog_trenutka[0], brzina_na_ulazu, pretpostavka[3],
                              vrijednost_proslog_trenutka[3], gustoca_na_ulazu, tlak, vrijednost_proslog_trenutka[1],
                              vrijednost_proslog_trenutka[2], entalpija_na_ulazu, pretpostavka[6], dx, pretpostavka[2],
                              podrelaksacija[2])

    gustoca = stanje_gustoca(tlak, entalpija)
    gustoca = eksplicitna(gustoca, pretpostavka[3], podrelaksacija[3])

    temperatura = stanje_temperatura(tlak, entalpija)
    temperatura = eksplicitna(temperatura, pretpostavka[4], podrelaksacija[4])

    alpha_unutarnje_izracunato = alpha_unutarnje(tlak, entalpija, brzina * gustoca, pretpostavka[6]) # ne zapisuje se

    temperatura_stijenke = temp_stijenke(vrijednost_proslog_trenutka[5], temperatura, alpha_unutarnje_izracunato, dx)
    temperatura_stijenke = eksplicitna(temperatura_stijenke, pretpostavka[5], podrelaksacija[5])

    specificni_toplinski_tok = alpha_unutarnje_izracunato * (temperatura_stijenke - temperatura)
    specificni_toplinski_tok = eksplicitna(specificni_toplinski_tok, pretpostavka[6], podrelaksacija[6])

    rjesenje = np.stack((brzina, tlak, entalpija, gustoca, temperatura, temperatura_stijenke, specificni_toplinski_tok))
    return rjesenje

def rjesavac_trenutka(trenutak, vrijednost_proslog_trenutka, podrelaksacija, rezidual_limit):
    global broj_iteracija
    broj_iteracija = 0
    rjesenje_1 = vrijednost_proslog_trenutka
    uvjet = True
    while uvjet:
        rjesenje_2 = vanjska_iteracija(trenutak, rjesenje_1, vrijednost_proslog_trenutka, podrelaksacija)
        rezidual = manhattan_rjesenja(rjesenje_1, rjesenje_2, 1)
        uvjet = any( np.greater(rezidual, rezidual_limit) )
        rjesenje_1 = rjesenje_2
        broj_iteracija = broj_iteracija + 1
        print(broj_iteracija)
        print(rezidual)
    return rjesenje_1

def rjesavac(pocetni_trenutak, podrelaksacija, rezidual_limit):
    for t in range(pocetni_trenutak, n_t):
        prosli_trenutak = citaj_stanje(t - 1)
        rjesenje_trenutka = rjesavac_trenutka(t, prosli_trenutak, podrelaksacija, rezidual_limit)
        snimi_trenutak(rjesenje_trenutka, t)
    return print('rjeseno')

def rjesavac_uz_moje_pretpostavke(podrelaksacija, rezidual_limit):
    rjesenje_trenutka = rjesavac_trenutka(0, stanje_0, podrelaksacija, rezidual_limit)
    snimi_trenutak(rjesenje_trenutka, 0)
    return print('rjeseno')

def snimi_trenutak(stanje, vremenski_trenutak):
    info = np.array([vremena[vremenski_trenutak], broj_iteracija, np.log10(zeljeni_rezidual), imp])
    putanja_direktorija = '/home/josip/Desktop/diplomski/kod/rezultati_coolprop/' + str(vremenski_trenutak) + '/'
    os.makedirs( os.path.dirname( putanja_direktorija ), exist_ok=True )
    putanja_info = putanja_direktorija + 'informacije.csv'
    np.savetxt(putanja_info, info, delimiter=",")
    for i in range( broj_velicina ):
        putanja_velicine = putanja_direktorija + raspored_velicina[i] + '.csv'
        np.savetxt(putanja_velicine, stanje[i], delimiter=",")
    return print('broj iteracija: ', broj_iteracija)

zeljeni_rezidual = 1e-5

imp = 0.1
eksp = 0.2
rjesavac(13, [imp, eksp, imp, eksp, eksp, eksp, eksp], np.repeat(zeljeni_rezidual, broj_velicina))


