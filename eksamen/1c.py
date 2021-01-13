import numpy as np


def piecewise_func(c, Ne, psi, D, N):
    """
    Input parametere:
        c    : Koeffisient vektoren, oppnådd ved å løse matriselikningen Ac = b
        Ne   : Antall elementer
        psi  : En liste med basisfunksjonene som funksjoner, altså at de kan ta et argument x
        D    : Grensebetingelsen for u(1,t) = D
        N    : Antall punkter jeg vil evaluere i hver element
    """

    c = np.append(c, D)    # Legger til grensebetingelsen slik at vi unngår en if-test når vi finner løsningen

    N_p1 = int(Ne/2)       # Antall P1 elementer
    N_p2 = int(Ne/2)       # Antall P2 elementer


    x = np.linspace(0, 1, Ne+1)

    # Degrees of freedom map som holder de lokale nodene for P2 elementene
    dof_map_P2 = []
    node = 0
    for e in range(N_p2):
        dof_map_P2.append([])
        for i in range(node, 3+node):
            dof_map_P2[e].append(i)
        node += 2

    # Degrees of freedom map som holder de lokale nodene for P1 elementene
    dof_map_P1 = []
    for e in range(N_p1):
        dof_map_P1.append([])
        for i in range(node, 2+node):
            dof_map_P1[e].append(i)
        node += 1


    x_total = []    # Holder den totale x-arrayen
    u = []          # Holder den totale løsningen

    counter = 0
    for index, e in enumerate(dof_map_P2):
        start = x[index]
        slutt = x[index+1]

        x_lokal = np.linspace(start, slutt, N)

        u.append(c[dof_map_P2[e][0]]*psi[dof_map_P2[e][0]](x_lokal) + c[dof_map_P2[e][1]]*psi[dof_map_P2[e][1]](x_lokal) + c[dof_map_P2[e][2]]*psi[dof_map_P2[e][2]](x_lokal))

        x_total = np.concatenate((x_total, x_lokal))
        counter = index

    for index, e in enumerate(dof_map_P1):
        start = x[index+counter]
        slutt = x[index+1+counter]

        x_lokal = np.linspace(start, slutt, 10)

        u.append(c[dof_map_P1[e][0]]*psi[dof_map_P1[e][0]](x_lokal) + c[dof_map_P1[e][1]]*psi[dof_map_P1[e][1]](x_lokal) + c[dof_map_P1[e][2]]*psi[dof_map_P1[e][2]](x_lokal))

        x_total = np.concatenate((x_total, x_lokal))


    u = np.concatenate(u, axis=0)      # Siden u er en liste med arrayer, gjør dette u til én array

    # For å unngå at samme x-koordinat blir brukt to ganger, finner indeksen til de elementene i x_total som er like
    idx = np.unique(x_total, return_index=True)[0]

    X = x_total[idx]
    U = u[idx]

    return X, U
