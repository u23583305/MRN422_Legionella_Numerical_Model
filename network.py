from pipe import Pipe
from dead_leg import DeadLeg


class Network:

    def __init__(self, T_source=60.0, T_inf=20.0, v=0.5):

        self.T_source = T_source
        self.T_inf = T_inf
        self.v = v

        self.pipes = []
        self.dead_leg_connections = {}

        # Actual thermocouple positions measured from 15B
        self.main_TC_positions = {
            "15B": 0.000,
            "7B": 0.320,
            "5B": 0.450,
            "6B": 0.630,
            "4B": 0.830,
            "8A": 1.030,
            "11A": 1.485,
            "12A": 1.845
        }

        # Actual fitting positions measured from 15B
        self.fitting_positions = {
            "T-piece 1": 0.100,
            "Elbow": 0.390,
            "T-piece 2": 1.250,
            "T-piece 3": 1.725
        }


    # ======================================================
    # NETWORK CONSTRUCTION
    # ======================================================

    def add_pipe(self, pipe):

        self.pipes.append(pipe)


    def attach_dead_leg(self, after_pipe, dead_leg):

        if after_pipe not in self.dead_leg_connections:
            self.dead_leg_connections[after_pipe] = []

        self.dead_leg_connections[after_pipe].append(
            dead_leg
        )


    # ======================================================
    # STEADY FLOW
    # ======================================================

    def steady_flow(self):

        results = []

        T_current = self.T_source

        distance = 0.0
        total_travel_time = 0.0


        for i, pipe in enumerate(self.pipes):

            T_in = T_current

            (
                T_out,
                UA,
                m_dot,
                Re,
                Nu,
                T_surface
            ) = pipe.steady_flow_temperature(
                T_in,
                self.v
            )

            travel_time = pipe.travel_time(
                T_in,
                self.v
            )

            distance += pipe.L
            total_travel_time += travel_time

            result = {
                "index": i,
                "pipe": pipe.name,
                "length": pipe.L,
                "distance": distance,
                "T_in": T_in,
                "T_out": T_out,
                "UA": UA,
                "m_dot": m_dot,
                "Re": Re,
                "Nu": Nu,
                "surface_temperature": T_surface,
                "travel_time": travel_time,
                "cumulative_travel_time": total_travel_time,
                "dead_legs": []
            }

            if i in self.dead_leg_connections:

                for dead_leg in self.dead_leg_connections[i]:

                    result["dead_legs"].append({
                        "name": dead_leg.name,
                        "connection_temperature": T_out
                    })

            results.append(result)

            T_current = T_out


        return results


    # ======================================================
    # CURRENT SIMPLE RECIRCULATION MODEL
    # ======================================================

    def recirculation(
        self,
        initial_temperatures,
        duration,
        dt=1.0
    ):

        times = [
            i * dt
            for i in range(
                int(duration / dt) + 1
            )
        ]

        curves = [
            [T]
            for T in initial_temperatures
        ]

        arrival_times = []

        cumulative_time = 0.0


        for pipe in self.pipes:

            cumulative_time += pipe.travel_time(
                self.T_source,
                self.v
            )

            arrival_times.append(
                cumulative_time
            )


        steady_results = self.steady_flow()

        target_temperatures = [
            result["T_out"]
            for result in steady_results
        ]


        for t in times[1:]:

            for i in range(
                len(self.pipes)
            ):

                if t >= arrival_times[i]:

                    curves[i].append(
                        target_temperatures[i]
                    )

                else:

                    curves[i].append(
                        initial_temperatures[i]
                    )


        return (
            times,
            curves,
            arrival_times
        )


    # ======================================================
    # RESULT HELPERS
    # ======================================================

    def pipe_temperatures(self, results):

        distances = [0.0]

        temperatures = [
            self.T_source
        ]


        for result in results:

            distances.append(
                result["distance"]
            )

            temperatures.append(
                result["T_out"]
            )


        return (
            distances,
            temperatures
        )


    def connection_temperatures(
        self,
        results
    ):

        connections = {}


        for result in results:

            for dead_leg in result["dead_legs"]:

                connections[
                    dead_leg["name"]
                ] = dead_leg[
                    "connection_temperature"
                ]


        return connections


    def total_length(self):

        return sum(
            pipe.L
            for pipe in self.pipes
        )


    def total_travel_time(
        self,
        results
    ):

        if results:

            return results[-1][
                "cumulative_travel_time"
            ]

        return 0.0


    # ======================================================
    # MAIN-LINE STAGNATION
    # ======================================================

    def main_stagnation(
        self,
        T_init,
        duration,
        dt=1.0
    ):

        names = list(
            self.main_TC_positions.keys()
        )

        positions = list(
            self.main_TC_positions.values()
        )

        n = len(names)


        if isinstance(T_init, dict):

            T = [
                T_init[name]
                for name in names
            ]

        elif isinstance(
            T_init,
            (list, tuple)
        ):

            T = list(T_init)

        else:

            T = [T_init] * n


        if len(T) != n:

            raise ValueError(
                "T_init must contain one "
                "temperature for each "
                "main-line thermocouple."
            )


        boundaries = [0.0]

        boundaries += [
            (
                positions[i]
                + positions[i + 1]
            ) / 2
            for i in range(n - 1)
        ]

        boundaries += [
            positions[-1]
        ]


        section_lengths = [
            boundaries[i + 1]
            - boundaries[i]
            for i in range(n)
        ]


        times = [0.0]

        curves = {
            name: [T[i]]
            for i, name in enumerate(names)
        }


        t = 0.0


        while t < duration:

            T_old = T.copy()
            T_new = T.copy()


            for i in range(n):

                section = Pipe(
                    names[i],
                    L=section_lengths[i],
                    D_i=0.020,
                    D_o=0.022,
                    t_ins=0.0,
                    T_inf=self.T_inf
                )


                (
                    rho,
                    _,
                    Cp,
                    _,
                    k_water
                ) = section.water_properties(
                    T_old[i]
                )


                V_water = (
                    section.A_flow
                    * section_lengths[i]
                )

                V_cu = (
                    section.A_cu
                    * section_lengths[i]
                )


                m_water = (
                    rho
                    * V_water
                )

                m_cu = (
                    section.rho_cu
                    * V_cu
                )


                C_total = (
                    m_water * Cp
                    + m_cu * section.Cp_cu
                )


                Q_axial = 0.0


                if i > 0:

                    dx = (
                        positions[i]
                        - positions[i - 1]
                    )

                    G_left = (
                        k_water
                        * section.A_flow
                        + section.k_cu
                        * section.A_cu
                    ) / dx

                    Q_axial += (
                        G_left
                        * (
                            T_old[i - 1]
                            - T_old[i]
                        )
                    )


                if i < n - 1:

                    dx = (
                        positions[i + 1]
                        - positions[i]
                    )

                    G_right = (
                        k_water
                        * section.A_flow
                        + section.k_cu
                        * section.A_cu
                    ) / dx

                    Q_axial += (
                        G_right
                        * (
                            T_old[i + 1]
                            - T_old[i]
                        )
                    )


                # Calibrated effective stagnant internal convection
                UA, _, _, _, _, _ = (
                    section.calculate_stagnant_UA(
                        T_old[i]
                    )
                )


                Q_loss = (
                    UA
                    * (
                        T_old[i]
                        - self.T_inf
                    )
                )


                T_new[i] = (
                    T_old[i]
                    + (
                        (
                            Q_axial
                            - Q_loss
                        )
                        / C_total
                    ) * dt
                )


            T = T_new

            t += dt

            times.append(t)


            for i, name in enumerate(names):

                curves[name].append(
                    T[i]
                )


        return times, curves


    # ======================================================
    # COUPLED MAIN + DEAD-LEG STAGNATION
    # ======================================================

    def coupled_stagnation(
        self,
        T_init_main,
        T_init_dead_leg,
        duration,
        dt=1.0,
        t_ins=0.0,
        k_ins=0.036
    ):

        # --------------------------------------------------
        # MAIN-LINE NODES
        # --------------------------------------------------

        main_positions = {
            "15B": 0.000,
            "7B": 0.320,
            "5B": 0.450,
            "6B": 0.630,
            "4B": 0.830,
            "8A": 1.030,
            "T2": 1.250,
            "11A": 1.485,
            "12A": 1.845
        }

        main_names = list(
            main_positions.keys()
        )

        x_main = list(
            main_positions.values()
        )

        n_main = len(
            main_names
        )


        # --------------------------------------------------
        # INITIAL T-PIECE 2 TEMPERATURE
        # --------------------------------------------------

        x_8A = main_positions["8A"]
        x_T2 = main_positions["T2"]
        x_11A = main_positions["11A"]

        fraction = (
            (x_T2 - x_8A)
            / (x_11A - x_8A)
        )

        T_T2_initial = (
            T_init_main["8A"]
            + fraction
            * (
                T_init_main["11A"]
                - T_init_main["8A"]
            )
        )


        T_main = []

        for name in main_names:

            if name == "T2":

                T_main.append(
                    T_T2_initial
                )

            else:

                T_main.append(
                    T_init_main[name]
                )


        # --------------------------------------------------
        # DEAD-LEG NODES
        # --------------------------------------------------

        dead_leg = (
            self.dead_leg_connections[7][0]
        )

        # Apply the insulation condition passed to
        # coupled_stagnation.
        dead_leg.t_ins = t_ins
        dead_leg.k_ins = k_ins

        dead_names = [
            "7A",
            "9A"
        ]

        x_dead = [
            0.07,
            0.39
        ]

        T_dead = [
            T_init_dead_leg["7A"],
            T_init_dead_leg["9A"]
        ]


        # --------------------------------------------------
        # MAIN-LINE CONTROL VOLUMES
        # --------------------------------------------------

        main_boundaries = [
            0.0
        ]

        main_boundaries += [
            (
                x_main[i]
                + x_main[i + 1]
            ) / 2
            for i in range(
                n_main - 1
            )
        ]

        main_boundaries += [
            x_main[-1]
        ]


        main_lengths = [
            main_boundaries[i + 1]
            - main_boundaries[i]
            for i in range(
                n_main
            )
        ]


        # --------------------------------------------------
        # DEAD-LEG CONTROL VOLUMES
        # --------------------------------------------------

        dead_boundaries = [
            0.0,
            (
                x_dead[0]
                + x_dead[1]
            ) / 2,
            dead_leg.L
        ]

        dead_lengths = [
            dead_boundaries[i + 1]
            - dead_boundaries[i]
            for i in range(
                len(dead_names)
            )
        ]


        # --------------------------------------------------
        # STORAGE
        # --------------------------------------------------

        times = [
            0.0
        ]

        main_curves = {
            name: [T_main[i]]
            for i, name
            in enumerate(main_names)
        }

        dead_curves = {
            name: [T_dead[i]]
            for i, name
            in enumerate(dead_names)
        }


        t = 0.0


        # --------------------------------------------------
        # TIME LOOP
        # --------------------------------------------------

        while t < duration:

            main_old = (
                T_main.copy()
            )

            dead_old = (
                T_dead.copy()
            )

            main_new = (
                T_main.copy()
            )

            dead_new = (
                T_dead.copy()
            )


            # ==============================================
            # MAIN LINE
            # ==============================================

            for i in range(
                n_main
            ):

                section = Pipe(
                    main_names[i],
                    L=main_lengths[i],
                    D_i=0.020,
                    D_o=0.022,
                    t_ins=t_ins,
                    k_ins=k_ins,
                    T_inf=self.T_inf
                )


                (
                    rho,
                    _,
                    Cp,
                    _,
                    k_water
                ) = section.water_properties(
                    main_old[i]
                )


                V_water = (
                    section.A_flow
                    * main_lengths[i]
                )

                V_cu = (
                    section.A_cu
                    * main_lengths[i]
                )


                m_water = (
                    rho
                    * V_water
                )

                m_cu = (
                    section.rho_cu
                    * V_cu
                )


                C_total = (
                    m_water
                    * Cp
                    + m_cu
                    * section.Cp_cu
                )


                Q_axial = 0.0


                # ------------------------------------------
                # UPSTREAM AXIAL CONDUCTION
                # ------------------------------------------

                if i > 0:

                    dx = (
                        x_main[i]
                        - x_main[i - 1]
                    )

                    G_left = (
                        k_water
                        * section.A_flow
                        + section.k_cu
                        * section.A_cu
                    ) / dx

                    Q_axial += (
                        G_left
                        * (
                            main_old[i - 1]
                            - main_old[i]
                        )
                    )


                # ------------------------------------------
                # DOWNSTREAM AXIAL CONDUCTION
                # ------------------------------------------

                if i < (
                    n_main - 1
                ):

                    dx = (
                        x_main[i + 1]
                        - x_main[i]
                    )

                    G_right = (
                        k_water
                        * section.A_flow
                        + section.k_cu
                        * section.A_cu
                    ) / dx

                    Q_axial += (
                        G_right
                        * (
                            main_old[i + 1]
                            - main_old[i]
                        )
                    )


                # ------------------------------------------
                # DEAD-LEG COUPLING AT T-PIECE 2
                # ------------------------------------------

                if (
                    main_names[i]
                    == "T2"
                ):

                    dx_dead = (
                        x_dead[0]
                    )

                    G_dead = (
                        k_water
                        * section.A_flow
                        + section.k_cu
                        * section.A_cu
                    ) / dx_dead

                    Q_dead_connection = (
                        G_dead
                        * (
                            dead_old[0]
                            - main_old[i]
                        )
                    )

                    Q_axial += (
                        Q_dead_connection
                    )


                # ------------------------------------------
                # AMBIENT HEAT LOSS
                # ------------------------------------------

                (
                    UA,
                    _,
                    _,
                    _,
                    _,
                    _
                ) = section.calculate_stagnant_UA(
                    main_old[i]
                )


                Q_loss = (
                    UA
                    * (
                        main_old[i]
                        - self.T_inf
                    )
                )


                main_new[i] = (
                    main_old[i]
                    + (
                        (
                            Q_axial
                            - Q_loss
                        )
                        / C_total
                    )
                    * dt
                )


            # ==============================================
            # DEAD LEG
            # ==============================================

            for i in range(
                len(dead_names)
            ):

                T_current = (
                    dead_old[i]
                )

                (
                    rho,
                    Cp,
                    k_water
                ) = dead_leg.water_properties(
                    T_current
                )


                L_section = (
                    dead_lengths[i]
                )


                m_water = (
                    rho
                    * dead_leg.A_water
                    * L_section
                )

                m_cu = (
                    dead_leg.rho_cu
                    * dead_leg.A_cu
                    * L_section
                )


                C_total = (
                    m_water
                    * Cp
                    + m_cu
                    * dead_leg.Cp_cu
                )


                Q_axial = 0.0


                # ------------------------------------------
                # T2 TO 7A COUPLING
                # ------------------------------------------

                if i == 0:

                    dx = (
                        x_dead[0]
                    )

                    G_T2 = (
                        k_water
                        * dead_leg.A_water
                        + dead_leg.k_cu
                        * dead_leg.A_cu
                    ) / dx

                    T_T2 = (
                        main_old[
                            main_names.index(
                                "T2"
                            )
                        ]
                    )

                    Q_axial += (
                        G_T2
                        * (
                            T_T2
                            - dead_old[i]
                        )
                    )


                # ------------------------------------------
                # PREVIOUS DEAD-LEG NODE
                # ------------------------------------------

                if i > 0:

                    dx = (
                        x_dead[i]
                        - x_dead[i - 1]
                    )

                    G_left = (
                        k_water
                        * dead_leg.A_water
                        + dead_leg.k_cu
                        * dead_leg.A_cu
                    ) / dx

                    Q_axial += (
                        G_left
                        * (
                            dead_old[i - 1]
                            - dead_old[i]
                        )
                    )


                # ------------------------------------------
                # NEXT DEAD-LEG NODE
                # ------------------------------------------

                if i < (
                    len(dead_names) - 1
                ):

                    dx = (
                        x_dead[i + 1]
                        - x_dead[i]
                    )

                    G_right = (
                        k_water
                        * dead_leg.A_water
                        + dead_leg.k_cu
                        * dead_leg.A_cu
                    ) / dx

                    Q_axial += (
                        G_right
                        * (
                            dead_old[i + 1]
                            - dead_old[i]
                        )
                    )


                # ------------------------------------------
                # AMBIENT HEAT LOSS
                # ------------------------------------------

                UA, _ = (
                    dead_leg.calculate_UA(
                        dead_old[i],
                        L_section
                    )
                )


                Q_loss = (
                    UA
                    * (
                        dead_old[i]
                        - self.T_inf
                    )
                )


                dead_new[i] = (
                    dead_old[i]
                    + (
                        (
                            Q_axial
                            - Q_loss
                        )
                        / C_total
                    )
                    * dt
                )


            # ==============================================
            # UPDATE TEMPERATURES
            # ==============================================

            T_main = (
                main_new
            )

            T_dead = (
                dead_new
            )

            t += dt

            times.append(
                t
            )


            for i, name in enumerate(
                main_names
            ):

                main_curves[
                    name
                ].append(
                    T_main[i]
                )


            for i, name in enumerate(
                dead_names
            ):

                dead_curves[
                    name
                ].append(
                    T_dead[i]
                )


        return (
            times,
            main_curves,
            dead_curves
        )


        # ======================================================
    # COUPLED MAIN + DEAD-LEG RECIRCULATION
    # ======================================================

    def coupled_recirculation(
        self,
        T_init_main,
        T_init_dead_leg,
        inlet_time,
        inlet_temperature,
        flow_rate_L_min,
        duration,
        dt=0.05,
        t_ins=0.0,
        k_ins=0.036,
        dead_leg_cavity_velocity=0.0
    ):

        import numpy as np

        # --------------------------------------------------
        # MAIN-LINE NODES - unchanged from validated model
        # --------------------------------------------------
        main_positions = {
            "15B": 0.000,
            "7B": 0.320,
            "5B": 0.450,
            "6B": 0.630,
            "4B": 0.830,
            "8A": 1.030,
            "T2": 1.250,
            "11A": 1.485,
            "12A": 1.845
        }
        main_names = list(main_positions.keys())
        x_main = list(main_positions.values())
        n_main = len(main_names)

        x_8A = main_positions["8A"]
        x_T2 = main_positions["T2"]
        x_11A = main_positions["11A"]
        fraction = (x_T2 - x_8A) / (x_11A - x_8A)
        T_T2_initial = (
            T_init_main["8A"]
            + fraction * (T_init_main["11A"] - T_init_main["8A"])
        )

        T_main = []
        for name in main_names:
            if name == "T2":
                T_main.append(T_T2_initial)
            else:
                T_main.append(T_init_main[name])

        # --------------------------------------------------
        # DEAD LEG - spatially resolved control volumes
        # --------------------------------------------------
        dead_leg = self.dead_leg_connections[7][0]
        dead_leg.t_ins = t_ins
        dead_leg.k_ins = k_ins

        # Internal numerical nodes. 7A and 9A remain at their
        # measured locations and are used for validation output.
        # The first numerical node is placed at the measured 7A location.
        # This lets the near-junction control volume exchange directly with T2,
        # instead of forcing the T2 disturbance through an artificial 0.025 m
        # node before it can reach 7A.
        x_dead = [
            0.070,
            0.120,
            0.170,
            0.220,
            0.270,
            0.320,
            0.370,
            0.390
        ]
        n_dead = len(x_dead)

        # Only two temperatures are measured initially in the dead leg.
        # Use a linear profile between the measured 7A and 9A values.
        # The first control volume spans the T2-to-7A entrance region and
        # is represented by the measured 7A initial temperature.
        T_7A_initial = T_init_dead_leg["7A"]
        T_9A_initial = T_init_dead_leg["9A"]
        T_dead = []
        for x in x_dead:
            if x == 0.070:
                T_dead.append(T_7A_initial)
            else:
                f = (x - 0.070) / (0.390 - 0.070)
                T_dead.append(
                    T_7A_initial + f * (T_9A_initial - T_7A_initial)
                )

        # --------------------------------------------------
        # CONTROL-VOLUME LENGTHS
        # --------------------------------------------------
        main_boundaries = [0.0]
        main_boundaries += [
            (x_main[i] + x_main[i + 1]) / 2
            for i in range(n_main - 1)
        ]
        main_boundaries += [x_main[-1]]
        main_lengths = [
            main_boundaries[i + 1] - main_boundaries[i]
            for i in range(n_main)
        ]

        dead_boundaries = [0.0]
        dead_boundaries += [
            (x_dead[i] + x_dead[i + 1]) / 2
            for i in range(n_dead - 1)
        ]
        dead_boundaries += [dead_leg.L]
        dead_lengths = [
            dead_boundaries[i + 1] - dead_boundaries[i]
            for i in range(n_dead)
        ]

        # --------------------------------------------------
        # FLOW RATE
        # --------------------------------------------------
        flow_rate_m3_s = flow_rate_L_min / 1000.0 / 60.0
        D_i = 0.020
        A_flow = np.pi * D_i ** 2 / 4
        velocity = flow_rate_m3_s / A_flow

        # --------------------------------------------------
        # CREATE MAIN SECTIONS ONCE
        # --------------------------------------------------
        main_sections = []
        for i in range(n_main):
            section = Pipe(
                main_names[i],
                L=main_lengths[i],
                D_i=0.020,
                D_o=0.022,
                t_ins=t_ins,
                k_ins=k_ins,
                T_inf=self.T_inf
            )
            main_sections.append(section)

        # --------------------------------------------------
        # PROPERTY CACHES
        # --------------------------------------------------
        main_property_cache = {}
        main_UA_cache = {}
        dead_property_cache = {}
        dead_UA_cache = {}

        def main_properties(section, T):
            key = round(T, 1)
            if key not in main_property_cache:
                main_property_cache[key] = section.water_properties(key)
            return main_property_cache[key]

        def main_UA(section, section_index, T, h_i_flow):
            key = (section_index, round(T, 1), round(h_i_flow, 1))
            if key not in main_UA_cache:
                UA, _ = section.calculate_UA(round(T, 1), h_i_flow)
                main_UA_cache[key] = UA
            return main_UA_cache[key]

        def dead_properties(T):
            key = round(T, 1)
            if key not in dead_property_cache:
                dead_property_cache[key] = dead_leg.water_properties(key)
            return dead_property_cache[key]

        def dead_UA(T, L_section, section_index):
            key = (section_index, round(T, 1), round(L_section, 5))
            if key not in dead_UA_cache:
                UA, _ = dead_leg.calculate_UA(round(T, 1), L_section)
                dead_UA_cache[key] = UA
            return dead_UA_cache[key]

        # --------------------------------------------------
        # STORAGE
        # --------------------------------------------------
        times = [0.0]
        main_curves = {
            name: [T_main[i]]
            for i, name in enumerate(main_names)
        }
        dead_curves = {
            "7A": [T_7A_initial],
            "9A": [T_9A_initial]
        }
        t = 0.0

        # --------------------------------------------------
        # TIME LOOP
        # --------------------------------------------------
        while t < duration:
            main_old = T_main.copy()
            dead_old = T_dead.copy()
            main_new = T_main.copy()
            dead_new = T_dead.copy()

            T_inlet = np.interp(t, inlet_time, inlet_temperature)
            T_T2_old = main_old[main_names.index("T2")]

            # --------------------------------------------------
            # SPATIALLY DISTRIBUTED CAVITY-FLOW EXCHANGE
            # --------------------------------------------------
            # The effective cavity velocity remains the single calibrated
            # parameter. Exchange is applied conservatively between adjacent
            # control volumes rather than directly from 7A to 9A.
            cavity_interface_Q = [0.0] * n_dead

            if dead_leg_cavity_velocity > 0.0:
                # T2 -> near-junction control volume represented by 7A.
                # This is the direct cavity-exchange entrance term.
                rho_cav, Cp_cav, _ = dead_properties(
                    (T_T2_old + dead_old[0]) / 2
                )
                cavity_interface_Q[0] = (
                    rho_cav
                    * Cp_cav
                    * dead_leg.A_water
                    * dead_leg_cavity_velocity
                    * (T_T2_old - dead_old[0])
                )

                # Exchange between successive dead-leg control volumes.
                # The Shiina heat-balance attenuation is evaluated locally
                # over each axial interval.
                for j in range(1, n_dead):
                    dx_cav = x_dead[j] - x_dead[j - 1]
                    T_mean = (dead_old[j - 1] + dead_old[j]) / 2
                    rho_cav, Cp_cav, _ = dead_properties(T_mean)
                    UA_pen = dead_UA(T_mean, dx_cav, j)
                    K_pen = UA_pen / dx_cav
                    a_pen = (
                        4.0 * K_pen * dx_cav
                        / (
                            rho_cav
                            * Cp_cav
                            * dead_leg_cavity_velocity
                            * dead_leg.D_i ** 2
                        )
                    )
                    penetration_factor = 1.0 / (1.0 + a_pen)
                    cavity_interface_Q[j] = (
                        rho_cav
                        * Cp_cav
                        * dead_leg.A_water
                        * dead_leg_cavity_velocity
                        * penetration_factor
                        * (dead_old[j - 1] - dead_old[j])
                    )

            Q_cavity_T2 = cavity_interface_Q[0]

            # ==================================================
            # MAIN LINE - physics unchanged
            # ==================================================
            for i in range(n_main):
                section = main_sections[i]
                rho, mu, Cp, Pr, k_water = main_properties(
                    section, main_old[i]
                )

                V_water = section.A_flow * main_lengths[i]
                V_cu = section.A_cu * main_lengths[i]
                m_water = rho * V_water
                m_cu = section.rho_cu * V_cu
                C_total = (
                    m_water * Cp
                    + m_cu * section.Cp_cu
                )

                Re = rho * velocity * section.D_i / mu
                if Re < 2300:
                    Nu = 3.66
                else:
                    Nu = 0.023 * Re ** 0.8 * Pr ** 0.4
                h_i_flow = Nu * k_water / section.D_i

                UA = main_UA(section, i, main_old[i], h_i_flow)
                Q_loss = UA * (main_old[i] - self.T_inf)

                if i == 0:
                    T_upstream = T_inlet
                else:
                    T_upstream = main_old[i - 1]

                m_dot = rho * flow_rate_m3_s
                Q_flow = m_dot * Cp * (T_upstream - main_old[i])

                Q_axial = 0.0
                if i > 0:
                    dx = x_main[i] - x_main[i - 1]
                    G_left = (
                        k_water * section.A_flow
                        + section.k_cu * section.A_cu
                    ) / dx
                    Q_axial += G_left * (main_old[i - 1] - main_old[i])

                if i < n_main - 1:
                    dx = x_main[i + 1] - x_main[i]
                    G_right = (
                        k_water * section.A_flow
                        + section.k_cu * section.A_cu
                    ) / dx
                    Q_axial += G_right * (main_old[i + 1] - main_old[i])

                if main_names[i] == "T2":
                    # Axial conduction from T2 to the 7A entrance control volume.
                    dx_dead = x_dead[0]
                    G_dead = (
                        k_water * section.A_flow
                        + section.k_cu * section.A_cu
                    ) / dx_dead
                    Q_axial += G_dead * (dead_old[0] - main_old[i])

                    # Equal-and-opposite cavity-flow exchange.
                    Q_axial -= Q_cavity_T2

                dTdt = (Q_flow + Q_axial - Q_loss) / C_total
                main_new[i] = main_old[i] + dTdt * dt

            # ==================================================
            # DEAD LEG - spatial control volumes
            # ==================================================
            for i in range(n_dead):
                rho, Cp, k_water = dead_properties(dead_old[i])
                L_section = dead_lengths[i]

                m_water = rho * dead_leg.A_water * L_section
                m_cu = dead_leg.rho_cu * dead_leg.A_cu * L_section
                C_total = (
                    m_water * Cp
                    + m_cu * dead_leg.Cp_cu
                )

                Q_axial = 0.0

                if i == 0:
                    dx = x_dead[0]
                    G_T2 = (
                        k_water * dead_leg.A_water
                        + dead_leg.k_cu * dead_leg.A_cu
                    ) / dx
                    Q_axial += G_T2 * (T_T2_old - dead_old[i])

                if i > 0:
                    dx = x_dead[i] - x_dead[i - 1]
                    G_left = (
                        k_water * dead_leg.A_water
                        + dead_leg.k_cu * dead_leg.A_cu
                    ) / dx
                    Q_axial += G_left * (dead_old[i - 1] - dead_old[i])

                if i < n_dead - 1:
                    dx = x_dead[i + 1] - x_dead[i]
                    G_right = (
                        k_water * dead_leg.A_water
                        + dead_leg.k_cu * dead_leg.A_cu
                    ) / dx
                    Q_axial += G_right * (dead_old[i + 1] - dead_old[i])

                UA = dead_UA(dead_old[i], L_section, i)
                Q_loss = UA * (dead_old[i] - self.T_inf)

                # Conservative pairwise cavity exchange:
                # incoming from upstream interface, outgoing downstream.
                Q_cavity = cavity_interface_Q[i]
                if i < n_dead - 1:
                    Q_cavity -= cavity_interface_Q[i + 1]

                dead_new[i] = (
                    dead_old[i]
                    + ((Q_axial + Q_cavity - Q_loss) / C_total) * dt
                )

            # --------------------------------------------------
            # IMPOSE MEASURED 15B BOUNDARY AND UPDATE
            # --------------------------------------------------
            T_next = np.interp(
                t + dt, inlet_time, inlet_temperature
            )
            main_new[0] = T_next

            T_main = main_new
            T_dead = dead_new
            t += dt
            times.append(t)

            for i, name in enumerate(main_names):
                main_curves[name].append(T_main[i])

            # 7A is an exact numerical node; 9A is the closed-end node.
            T_7A = np.interp(0.070, x_dead, T_dead)
            T_9A = np.interp(0.390, x_dead, T_dead)
            dead_curves["7A"].append(T_7A)
            dead_curves["9A"].append(T_9A)

        return (
            times,
            main_curves,
            dead_curves,
            velocity
        )

# ==========================================================
# ACTUAL EXPERIMENTAL NETWORK
# ==========================================================

network = Network(
    T_source=60.0,
    T_inf=20.0,
    v=0.5
)


# Main flow path:
#
# 15B -> T1 -> 7B -> elbow -> 5B -> 6B ->
# 4B -> 8A -> T2 -> 11A -> T3 -> 12A


pipe_1 = Pipe(
    "15B to T-piece 1",
    L=0.10,
    t_ins=0.0
)

pipe_2 = Pipe(
    "T-piece 1 to 7B",
    L=0.22,
    t_ins=0.0
)

pipe_3 = Pipe(
    "7B to Elbow",
    L=0.07,
    t_ins=0.0
)

pipe_4 = Pipe(
    "Elbow to 5B",
    L=0.06,
    t_ins=0.0
)

pipe_5 = Pipe(
    "5B to 6B",
    L=0.18,
    t_ins=0.0
)

pipe_6 = Pipe(
    "6B to 4B",
    L=0.20,
    t_ins=0.0
)

pipe_7 = Pipe(
    "4B to 8A",
    L=0.20,
    t_ins=0.0
)

pipe_8 = Pipe(
    "8A to T-piece 2",
    L=0.22,
    t_ins=0.0
)

pipe_9 = Pipe(
    "T-piece 2 to 11A",
    L=0.235,
    t_ins=0.0
)

pipe_10 = Pipe(
    "11A to T-piece 3",
    L=0.24,
    t_ins=0.0
)

pipe_11 = Pipe(
    "T-piece 3 to 12A",
    L=0.12,
    t_ins=0.0
)


network.add_pipe(pipe_1)
network.add_pipe(pipe_2)
network.add_pipe(pipe_3)
network.add_pipe(pipe_4)
network.add_pipe(pipe_5)
network.add_pipe(pipe_6)
network.add_pipe(pipe_7)
network.add_pipe(pipe_8)
network.add_pipe(pipe_9)
network.add_pipe(pipe_10)
network.add_pipe(pipe_11)


# ==========================================================
# DEAD LEG
# ==========================================================

# T-piece 2 -> 7A = 0.07 m
# 7A -> 9A = 0.32 m
# 9A is the closed end.

dead_leg_1 = DeadLeg(
    "Dead Leg 1",
    L=0.39,
    TC_positions=[
        0.07,
        0.39
    ],
    t_ins=0.0
)


# pipe_8 is index 7.
# Its outlet is T-piece 2.

network.attach_dead_leg(
    after_pipe=7,
    dead_leg=dead_leg_1
)


# ==========================================================
# TEST NETWORK
# ==========================================================

if __name__ == "__main__":

    results = network.steady_flow()

    connections = (
        network.connection_temperatures(
            results
        )
    )


    print("\nMAIN PIPE NETWORK")


    for result in results:

        print(
            f'{result["pipe"]} | '
            f'T_in = '
            f'{result["T_in"]:.3f} °C | '
            f'T_out = '
            f'{result["T_out"]:.3f} °C | '
            f'travel = '
            f'{result["travel_time"]:.2f} s'
        )


    print("\nDEAD-LEG CONNECTION")


    for name, temperature in (
        connections.items()
    ):

        print(
            f"{name} | "
            f"T-piece temperature = "
            f"{temperature:.3f} °C"
        )


    print("\nNETWORK SUMMARY")

    print(
        f"Total length = "
        f"{network.total_length():.3f} m"
    )

    print(
        f"Total travel time = "
        f"{network.total_travel_time(results):.2f} s"
    )


    print("\nMAIN THERMOCOUPLE POSITIONS")


    for name, position in (
        network.main_TC_positions.items()
    ):

        print(
            f"{name} | "
            f"x = {position:.3f} m"
        )


    print("\nFITTING POSITIONS")


    for name, position in (
        network.fitting_positions.items()
    ):

        print(
            f"{name} | "
            f"x = {position:.3f} m"
        )


    # ======================================================
    # SIMPLE COUPLING TEST
    # ======================================================

    test_main = {
        "15B": 60.0,
        "7B": 60.0,
        "5B": 60.0,
        "6B": 60.0,
        "4B": 60.0,
        "8A": 60.0,
        "11A": 60.0,
        "12A": 60.0
    }

    test_dead = {
        "7A": 50.0,
        "9A": 50.0
    }


    (
        test_times,
        test_main_curves,
        test_dead_curves
    ) = network.coupled_stagnation(
        T_init_main=test_main,
        T_init_dead_leg=test_dead,
        duration=60.0,
        dt=1.0
    )


    print("\nCOUPLED STAGNATION TEST")

    print(
        f"T2 | "
        f"T(0) = "
        f"{test_main_curves['T2'][0]:.3f} °C | "
        f"T(60s) = "
        f"{test_main_curves['T2'][-1]:.3f} °C"
    )

    print(
        f"7A | "
        f"T(0) = "
        f"{test_dead_curves['7A'][0]:.3f} °C | "
        f"T(60s) = "
        f"{test_dead_curves['7A'][-1]:.3f} °C"
    )

    print(
        f"9A | "
        f"T(0) = "
        f"{test_dead_curves['9A'][0]:.3f} °C | "
        f"T(60s) = "
        f"{test_dead_curves['9A'][-1]:.3f} °C"
    )