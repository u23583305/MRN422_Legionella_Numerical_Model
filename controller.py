class Controller:
    def __init__(self, T_on=45.0, T_off=55.0, max_runtime=600, dTdt_min=0.01):
        self.T_on = T_on
        self.T_off = T_off
        self.max_runtime = max_runtime
        self.dTdt_min = dTdt_min

        self.pump_on = False
        self.runtime = 0.0

    def reset(self):
        self.pump_on = False
        self.runtime = 0.0

    def start_condition(self, T_control):
        return T_control <= self.T_on

    def stop_condition(self, T_control, dTdt=None):
        if T_control >= self.T_off:
            return True

        if self.runtime >= self.max_runtime:
            return True

        if dTdt is not None and dTdt <= self.dTdt_min:
            return True

        return False

    def update(self, T_control, dt, dTdt=None):
        if not self.pump_on:
            if self.start_condition(T_control):
                self.pump_on = True
                self.runtime = 0.0

        else:
            self.runtime += dt

            if self.stop_condition(T_control, dTdt):
                self.pump_on = False

        return self.pump_on