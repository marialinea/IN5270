import numpy as np


def L2_norm(u_e, u, Ns):
    """
    Input parametere:
        u_e : den eksakte løsningen som er en funksjon av x
        u   : den tilnærmede løsningen i diskrete punkter
        Ns  : antall integrasjonspunkter
    """

    h = 1/(Ns-1)
    x = np.linspace(0, 1, Ns)


    # Regner ut arealet av trapesene som tilhører de indre punktene
    inder_punkter = 0
    for i in range(1,Ns):
        indre_punkter += (u_e(x[i]) - u[i])**2

    endepunkter = ((u_e(x[0]) - u[0])**2 + (u_e(x[-1]) - u[-1])**2)/2

    return np.sqrt(h*(endepunkter + indre_punkter))
