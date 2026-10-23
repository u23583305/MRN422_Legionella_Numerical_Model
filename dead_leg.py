import math
from CoolProp.CoolProp import PropsSI


class DeadLeg:
    def __init__(self, name, L, TC_positions, D_i=0.020, D_o=0.022, t_ins=0.010, k_ins=0.036, T_inf=20.0, P=101325):
        self.name, self.L, self.TC_positions = name, L, TC_positions
        self.D_i, self.D_o = D_i, D_o
        self.t_ins, self.k_ins = t_ins, k_ins
        self.T_inf, self.P = T_inf, P

        self.r_i, self.r_o = D_i / 2, D_o / 2
        self.r_ins = self.r_o + t_ins
        self.D_outer = 2 * self.r_ins

        self.A_water = math.pi * self.r_i**2
        self.A_cu = math.pi * (
            self.r_o**2 - self.r_i**2
        )

        self.k_cu = 401.0
        self.rho_cu = 8933.0
        self.Cp_cu = 385.0

        self.epsilon = 0.80
        self.sigma = 5.670374419e-8
        self.g = 9.81

        self.h_i_stag = 37.5

        self.boundaries = (
            [0.0] +
            [
                (
                    TC_positions[i] +
                    TC_positions[i + 1]
                ) / 2
                for i in range(
                    len(TC_positions) - 1
                )
            ] +
            [L]
        )

        self.section_lengths = [
            self.boundaries[i + 1] -
            self.boundaries[i]
            for i in range(
                len(TC_positions)
            )
        ]

    def water_properties(self, T):
        T_K = T + 273.15

        rho = PropsSI(
            "D", "T", T_K,
            "P", self.P, "Water"
        )

        Cp = PropsSI(
            "CPMASS", "T", T_K,
            "P", self.P, "Water"
        )

        k = PropsSI(
            "CONDUCTIVITY", "T", T_K,
            "P", self.P, "Water"
        )

        return rho, Cp, k

    def calculate_UA(
        self,
        T_bulk,
        L_section
    ):
        self.r_ins = self.r_o + self.t_ins
        self.D_outer = 2 * self.r_ins

        A_i = (
            2 *
            math.pi *
            self.r_i *
            L_section
        )

        A_outer = (
            2 *
            math.pi *
            self.r_ins *
            L_section
        )

        h_i = self.h_i_stag

        R_i = (
            1 /
            (h_i * A_i)
        )

        R_cu = (
            math.log(
                self.r_o /
                self.r_i
            ) /
            (
                2 *
                math.pi *
                self.k_cu *
                L_section
            )
        )

        R_ins = (
            math.log(
                self.r_ins /
                self.r_o
            ) /
            (
                2 *
                math.pi *
                self.k_ins *
                L_section
            )
            if self.t_ins > 0
            else 0
        )

        T_surface = (
            T_bulk +
            self.T_inf
        ) / 2

        T_sur_K = (
            self.T_inf +
            273.15
        )

        for _ in range(50):
            T_film_K = (
                (
                    T_surface +
                    self.T_inf
                ) / 2
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
                "CONDUCTIVITY",
                "T", T_film_K,
                "P", self.P,
                "Air"
            )

            Cp_air = PropsSI(
                "CPMASS",
                "T", T_film_K,
                "P", self.P,
                "Air"
            )

            Pr_air = PropsSI(
                "PRANDTL",
                "T", T_film_K,
                "P", self.P,
                "Air"
            )

            nu_air = (
                mu_air /
                rho_air
            )

            alpha_air = (
                k_air /
                (
                    rho_air *
                    Cp_air
                )
            )

            beta = (
                1 /
                T_film_K
            )

            Ra_D = (
                self.g *
                beta *
                abs(
                    T_surface -
                    self.T_inf
                ) *
                self.D_outer**3 /
                (
                    nu_air *
                    alpha_air
                )
            )

            Nu_air = (
                0.60 +
                (
                    0.387 *
                    Ra_D**(1/6)
                ) /
                (
                    1 +
                    (
                        0.559 /
                        Pr_air
                    )**(9/16)
                )**(8/27)
            )**2

            h_o = (
                Nu_air *
                k_air /
                self.D_outer
            )

            T_surface_K = (
                T_surface +
                273.15
            )

            h_r = (
                self.epsilon *
                self.sigma *
                (
                    T_surface_K +
                    T_sur_K
                ) *
                (
                    T_surface_K**2 +
                    T_sur_K**2
                )
            )

            R_o = (
                1 /
                (
                    (h_o + h_r) *
                    A_outer
                )
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
                (
                    T_bulk -
                    self.T_inf
                ) *
                R_o
            )

            if (
                abs(
                    T_new -
                    T_surface
                ) < 0.001
            ):
                T_surface = T_new
                break

            T_surface = T_new

        return UA, T_surface

    def section_capacities(self, T):
        rho, Cp, _ = (
            self.water_properties(T)
        )

        capacities = []

        for L_section in self.section_lengths:
            m_water = (
                rho *
                self.A_water *
                L_section
            )

            m_cu = (
                self.rho_cu *
                self.A_cu *
                L_section
            )

            capacities.append(
                m_water * Cp +
                m_cu * self.Cp_cu
            )

        return capacities

    def simulate(
        self,
        T_init,
        duration,
        dt=1.0,
        T_main=None
    ):
        n = len(
            self.TC_positions
        )

        if isinstance(
            T_init,
            (list, tuple)
        ):
            T = list(T_init)
        else:
            T = [T_init] * n

        if len(T) != n:
            raise ValueError(
                "T_init must contain one temperature for each dead-leg thermocouple."
            )

        times = [0.0]

        curves = [
            [T[i]]
            for i in range(n)
        ]

        T_average = (
            sum(T) / n
        )

        _, _, k_water = (
            self.water_properties(
                T_average
            )
        )

        capacities = (
            self.section_capacities(
                T_average
            )
        )

        t = 0.0

        while t < duration:
            T_old = T.copy()
            T_new = T.copy()

            T_boundary = (
                T_main(t)
                if callable(T_main)
                else T_main
            )

            for i in range(n):
                Q_axial = 0.0

                if (
                    i == 0
                    and
                    T_boundary is not None
                ):
                    dx = (
                        self.TC_positions[0]
                    )

                    G = (
                        k_water *
                        self.A_water +
                        self.k_cu *
                        self.A_cu
                    ) / dx

                    Q_axial += (
                        G *
                        (
                            T_boundary -
                            T_old[i]
                        )
                    )

                if i > 0:
                    dx = (
                        self.TC_positions[i] -
                        self.TC_positions[i - 1]
                    )

                    G = (
                        k_water *
                        self.A_water +
                        self.k_cu *
                        self.A_cu
                    ) / dx

                    Q_axial += (
                        G *
                        (
                            T_old[i - 1] -
                            T_old[i]
                        )
                    )

                if i < n - 1:
                    dx = (
                        self.TC_positions[i + 1] -
                        self.TC_positions[i]
                    )

                    G = (
                        k_water *
                        self.A_water +
                        self.k_cu *
                        self.A_cu
                    ) / dx

                    Q_axial += (
                        G *
                        (
                            T_old[i + 1] -
                            T_old[i]
                        )
                    )

                UA, _ = (
                    self.calculate_UA(
                        T_old[i],
                        self.section_lengths[i]
                    )
                )

                Q_loss = (
                    UA *
                    (
                        T_old[i] -
                        self.T_inf
                    )
                )

                T_new[i] = (
                    T_old[i] +
                    (
                        (
                            Q_axial -
                            Q_loss
                        ) /
                        capacities[i]
                    ) *
                    dt
                )

            T = T_new
            t += dt

            times.append(t)

            for i in range(n):
                curves[i].append(
                    T[i]
                )

        return times, curves

    def risk_time(
        self,
        times,
        curve,
        T_risk_low,
        T_risk_high
    ):
        risk = [
            t
            for t, T in zip(
                times,
                curve
            )
            if (
                T_risk_low
                <= T
                <= T_risk_high
            )
        ]

        return (
            (
                max(risk) -
                min(risk)
            ) / 60
            if risk
            else 0
        )