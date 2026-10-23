import math
from CoolProp.CoolProp import PropsSI


class Pipe:
    def __init__(self, name, L, D_i=0.020, D_o=0.022, t_ins=0.010, k_ins=0.036, T_inf=20.0, P=101325):
        self.name = name
        self.L = L
        self.D_i = D_i
        self.D_o = D_o
        self.t_ins = t_ins
        self.k_ins = k_ins
        self.T_inf = T_inf
        self.P = P

        self.r_i = D_i / 2
        self.r_o = D_o / 2
        self.r_ins = self.r_o + t_ins
        self.D_outer = 2 * self.r_ins
        self.A_flow = math.pi * self.r_i**2
        self.A_cu = math.pi * (self.r_o**2 - self.r_i**2)

        self.k_cu = 401.0
        self.rho_cu = 8933.0
        self.Cp_cu = 385.0
        self.epsilon = 0.80
        self.sigma = 5.670374419e-8
        self.g = 9.81

        self.h_i_stag = 37.5

    def water_properties(self, T):
        T_K = T + 273.15
        rho = PropsSI("D", "T", T_K, "P", self.P, "Water")
        mu = PropsSI("VISCOSITY", "T", T_K, "P", self.P, "Water")
        Cp = PropsSI("CPMASS", "T", T_K, "P", self.P, "Water")
        Pr = PropsSI("PRANDTL", "T", T_K, "P", self.P, "Water")
        k = PropsSI("CONDUCTIVITY", "T", T_K, "P", self.P, "Water")
        return rho, mu, Cp, Pr, k

    def calculate_stagnant_UA(self, T_bulk):
        h_i = self.h_i_stag

        UA, T_surface = self.calculate_UA(
            T_bulk,
            h_i
        )

        A_i = 2 * math.pi * self.r_i * self.L
        R_i = 1 / (h_i * A_i)

        Q = UA * (T_bulk - self.T_inf)
        T_inner = T_bulk - Q * R_i

        return UA, T_surface, T_inner, None, None, h_i

    def flow_properties(self, T, v):
        rho, mu, Cp, Pr, k = self.water_properties(T)

        m_dot = rho * self.A_flow * v
        Re = rho * v * self.D_i / mu

        Nu = (
            3.66
            if Re < 2300
            else 0.023 * Re**0.8 * Pr**0.4
        )

        h_i = Nu * k / self.D_i

        return rho, Cp, m_dot, Re, Nu, h_i

    def calculate_UA(self, T_bulk, h_i):
        self.r_ins = self.r_o + self.t_ins
        self.D_outer = 2 * self.r_ins

        A_i = 2 * math.pi * self.r_i * self.L
        A_outer = 2 * math.pi * self.r_ins * self.L

        R_i = 1 / (h_i * A_i)

        R_cu = (
            math.log(self.r_o / self.r_i) /
            (2 * math.pi * self.k_cu * self.L)
        )

        R_ins = (
            math.log(self.r_ins / self.r_o) /
            (2 * math.pi * self.k_ins * self.L)
            if self.t_ins > 0
            else 0
        )

        T_surface = (T_bulk + self.T_inf) / 2
        T_sur_K = self.T_inf + 273.15

        for _ in range(20):
            T_film_K = (
                (T_surface + self.T_inf) / 2
            ) + 273.15

            rho_air = PropsSI(
                "D", "T", T_film_K,
                "P", self.P, "Air"
            )

            mu_air = PropsSI(
                "VISCOSITY", "T", T_film_K,
                "P", self.P, "Air"
            )

            k_air = PropsSI(
                "CONDUCTIVITY", "T", T_film_K,
                "P", self.P, "Air"
            )

            Cp_air = PropsSI(
                "CPMASS", "T", T_film_K,
                "P", self.P, "Air"
            )

            Pr_air = PropsSI(
                "PRANDTL", "T", T_film_K,
                "P", self.P, "Air"
            )

            nu_air = mu_air / rho_air
            alpha_air = k_air / (rho_air * Cp_air)
            beta = 1 / T_film_K

            Ra_D = (
                self.g *
                beta *
                abs(T_surface - self.T_inf) *
                self.D_outer**3 /
                (nu_air * alpha_air)
            )

            Nu_air = (
                0.60 +
                (0.387 * Ra_D**(1/6)) /
                (
                    (
                        1 +
                        (0.559 / Pr_air)**(9/16)
                    )**(8/27)
                )
            )**2

            h_o = Nu_air * k_air / self.D_outer

            T_surface_K = T_surface + 273.15

            h_r = (
                self.epsilon *
                self.sigma *
                (T_surface_K + T_sur_K) *
                (T_surface_K**2 + T_sur_K**2)
            )

            R_o = (
                1 /
                ((h_o + h_r) * A_outer)
            )

            UA = (
                1 /
                (
                    R_i +
                    R_cu +
                    R_ins +
                    R_o
                )
            )

            T_new = (
                self.T_inf +
                UA *
                (T_bulk - self.T_inf) *
                R_o
            )

            if abs(T_new - T_surface) < 0.001:
                T_surface = T_new
                break

            T_surface = T_new

        return UA, T_surface

    def steady_flow_temperature(self, T_in, v):
        rho, Cp, m_dot, Re, Nu, h_i = (
            self.flow_properties(T_in, v)
        )

        UA, T_surface = self.calculate_UA(
            T_in,
            h_i
        )

        T_out = (
            self.T_inf +
            (T_in - self.T_inf) *
            math.exp(
                -UA / (m_dot * Cp)
            )
        )

        return (
            T_out,
            UA,
            m_dot,
            Re,
            Nu,
            T_surface
        )

    def stagnant_temperature(self, T_init, t):
        rho, _, Cp, _, _ = (
            self.water_properties(T_init)
        )

        UA, _, _, _, _, _ = (
            self.calculate_stagnant_UA(T_init)
        )

        V_water = self.A_flow * self.L
        V_cu = self.A_cu * self.L

        m_water = rho * V_water
        m_cu = self.rho_cu * V_cu

        C_total = (
            m_water * Cp +
            m_cu * self.Cp_cu
        )

        return (
            self.T_inf +
            (T_init - self.T_inf) *
            math.exp(
                -(UA / C_total) * t
            )
        )

    def travel_time(self, T, v):
        rho, _, m_dot, _, _, _ = (
            self.flow_properties(T, v)
        )

        V = self.A_flow * self.L

        return rho * V / m_dot