import numpy as np
import matplotlib.pyplot as plt

# геометрические высоты, м
H_tab = np.array([-2000, 0, 11019, 20063, 32162, 47350, 51412, 71802, 86152, 95411, 104128, 120000, 140000, 160000, 200000,250000])
# геопотенциальные высоты, м
Hg_tab = np.array([-2000, 0, 11000, 20000, 32000, 47000, 51000, 71000, 85000, 94000, 102450, 117777, 140000, 160000, 200000, 250000])

# температурные градиенты, К/м
beta_tab = 1e-3 * np.array([-6.5, -6.5, 0, 1, 2.8, 0, -2.8, -2, 0, 3, 11, 11.259, (6.8 / 10), (3.97 / 10), (1.75 / 10), (0.57 / 10)])

# температуры, К
T_tab = np.array([301.15, 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.65, 186.525, 211.752, 334.417, 559.6, 695.6, 854.4, 941.9])

# давления, Па
p_tab = np.array([127774, 101325, 22632, 5474.87, 868.014, 110.906, 66.9384, 3.95639, 0.35, 0.07, 0.0143661, 0.00266618, 0.000826375, 0.000303620, 0.0000853025, 0.0000247564])



def pressure(h):
    a, b, c, d = (-4.66327058e-11, 1.64892202e-05, -1.72987157e+00, 5.04367594e+04)     # All
    a, b, c, d = (-4.48341384e-09,  3.77384224e-04, -1.05079889e+01,  9.91625699e+04)   # :20
    # a, b, c, d = (-4.74859848e-10, 8.64404044e-05, -4.77222136e+00, 7.68979966e+04)     # :50
    return a * h**3 + b * h**2 + c * h + d

def air_density(h):
    a, b, c, d = (-3.85008277e-14,  3.60207951e-09, -1.13411053e-04,  1.22108259e+00)  # :20
    # a, b, c, d = (-4.74859848e-10, 8.64404044e-05, -4.77222136e+00, 7.68979966e+04)     # :50
    return a * h ** 3 + b * h ** 2 + c * h + d


def air_temperature(h):
    return T0 - 0.0065 * h



def atm440181(H):
    """
    :param H: высота над уровнем Земли
    :return: (давление, плотность, температура, скорость звука)
    """
    assert H is not None, "Given H is None"

    R = 287.05287 # газовая постоянная
    R1 = 8314.32 # и еще одна
    g = 9.80665 # среднее ускорение свободного падения
    mu_0 = 28.964420 # молярная масса воздуха на уровне моря

    H_g = H * 6356767 / ( H + 6356767 ) # пересчет высоты в геопотенциальную по среднему радиусу Земли
    p, T = 0., 0.

    if H < 0:
        p = p_tab[0]
        rho = p_tab[0] * mu_0 / (R1 * T_tab[0])
        T = T_tab[0]
        a = 20.046796 * np.sqrt(T)
        return p, rho, T, a


    if H >= H_tab[-1]:
        T = 186.65
        p = 0
        ro = 0
        a = 273

        return p, ro, T, a
    else:
        for i in range(len(H_tab)):
            if H >= H_tab[i] and H < H_tab[i+1]:
                T = T_tab[i] + beta_tab[i] * (H_g - Hg_tab[i])
                if H_g > 100000:
                    T = 186.65

                if beta_tab[i] == 0:
                    p = p_tab[i] * np.exp( -g * (H_g - Hg_tab[i]) / (R * T))
                else:
                    p = p_tab[i] * ( ( 1 + beta_tab[i]*( H_g - Hg_tab[i] )/T_tab[i] )**( -g /(R*beta_tab[i]) ) )
                break

        ro = p * mu_0 / (R1 * T)
        a = 20.046796*np.sqrt(T)
        return p, ro, T, a


p0, ro0, T0, a0 = atm440181(0)




if __name__ == "__main__":
    print(atm440181(-10))
    Nt = 100
    hmax = 200.e3
    heights = np.linspace(0, hmax, Nt)
    p = np.zeros(Nt)
    ro = np.zeros(Nt)
    T = np.zeros(Nt)
    a = np.zeros(Nt)
    for i, h in enumerate(heights):
        p[i], ro[i], T[i], a[i] = atm440181(h)

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,8))
    axes[0][0].plot(heights, p, 'r.-', label='p')
    axes[0][1].plot(heights, ro, 'b.-', label='$\\rho$')
    axes[1][0].plot(heights, T, 'g*-', label='T')
    axes[1][1].plot(heights, a, 'k*-', label='a')
    for row in axes:
        for ax in row:
            ax.legend()
            ax.grid()

    plt.show()
