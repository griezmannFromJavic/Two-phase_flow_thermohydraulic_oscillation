import numpy as np

### diskretizacije i svojstva materijala ###
rt = 'water'
n_l = 50
l_c = 20
tau = 10 * 60
dt = 0.2

def geometrijska_mreza(r):
    a = l_c * (1 - r) / (1 - r ** n_l)
    lokacije = np.zeros(n_l + 1)
    for i in range(1, n_l + 1):
        lokacije[i] = a * (1 - r ** i) / (1 - r)
    return lokacije


# lokacije_x = np.linspace(0, l_c, n_l + 1)
lokacije_x = geometrijska_mreza(0.9)

lokacije_x_p = lokacije_x[1:] + 0.5 * (lokacije_x[1] - lokacije_x[0])
lokacije_x_p = np.insert(lokacije_x_p, 0, lokacije_x[0], axis=0)
lokacije_x_p[-1] = lokacije_x[-1]
dx = lokacije_x[1:] - lokacije_x[:-1]
dx_p = lokacije_x_p[1:] - lokacije_x_p[:-1]
lokacije_sredista = (lokacije_x[:-1] + lokacije_x[1:]) / 2
# vremena = np.arange(0, tau, dt)
vremena = np.zeros( int(tau / dt) )
n_t = len(vremena)

c_c = 880 
lambda_c = 45 
rho_c = 7850
k = 0.01 / 1000
alpha_vanjsko = 2500
temperatura_ogrijevnog_medija = 273.15 + 750 + 100 # 285
g = 9.81
beta = np.pi / 2 # vertikalno
# beta = 0
sin_beta = np.sin(beta)

# dn = 25
dn = 40
dn_unutarnji = [10, 15, 20, 25, 32, 40, 50, 65, 80, 100, 125, 150, 200, 250, 300]
dn_vanjski = [17.2, 21.3, 26.9, 33.7, 42.4, 48.3, 60.3, 76.1, 88.9, 114.3, 139.7, 168.3, 219.1, 273.0, 323.9]
d_u = dn / 1000
d_v = dn_vanjski[ dn_unutarnji.index(dn) ] / 1000
epsilon = k / d_u

a_u = d_u ** 2 * np.pi / 4
a_v = d_v ** 2 * np.pi / 4

raspored_velicina = ['brzina', 'tlak', 'entalpija', 'gustoca', 'temperatura', 'temp_stijenke', 'toplinski_tok']
broj_velicina = len(raspored_velicina)

broj_iteracija = 0
