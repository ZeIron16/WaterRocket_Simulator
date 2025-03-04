import math
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from sympy.interactive import printing


_PI = math.pi
_RHO_EAU = 1000.0
_R = 8.314
_MOL_WEIGHT_AIR = 0.029 

class Rocket:
    def __init__(self, empty_volume, water_volume, radius, nozzle_radius, empty_mass, drag_coeff, pressure, surface):
        self.empty_volume = empty_volume
        self.water_volume = water_volume
        self.radius = radius
        self.nozzle_radius = nozzle_radius
        self.empty_mass = empty_mass
        self.drag_coeff = drag_coeff
        self.pressure = pressure
        self.surface = surface

class Environment:
    def __init__(self, g, rho, T_ext, P_ext, gamma):
        self.g = g
        self.rho = rho
        self.T_ext = T_ext
        self.P_ext = P_ext
        self.gamma = gamma

def pourcent_eau(r):
    return r.water_volume / r.empty_volume

def temp2(r, e):
    res = 1.0 - pourcent_eau(r)
    res = math.pow(res, e.gamma - 1)
    return res * e.T_ext

def beta(e):
    return 1.03 + (0.021 * e.gamma)

def press2(r, e):
    res = (r.empty_volume - r.water_volume) / r.empty_volume
    res = math.pow(res, e.gamma)
    return r.pressure * res

def p_trans(e):
    res = (e.gamma + 1) / 2.0
    res = math.pow(res, e.gamma / (e.gamma - 1))
    return e.P_ext * res

r_m = _R / _MOL_WEIGHT_AIR

def _c2(r, e):
    t2 = temp2(r, e)
    res = e.gamma * t2 * r_m
    return math.sqrt(res)

def imp(r, e):
    p2 = press2(r, e)
    b = beta(e)
    pt = p_trans(e)
    c2 = _c2(r, e)

    terme1 = p2 / (b * e.P_ext)
    terme1 = math.pow(terme1, (e.gamma - 1) / (2 * e.gamma))
    terme1 = 1 - terme1
    terme1 = pt * terme1 / (p2 * (e.gamma - 1))

    terme2 = b * e.P_ext / p2
    terme2 = math.pow(terme2, (e.gamma + 1) / (2 * e.gamma))
    terme2 = 1 - terme2

    terme3 = 8.0 / (e.gamma + 1)
    terme3 = math.sqrt(terme3)

    terme4 = p2 * r.empty_volume / c2

    res = terme3 * terme4 * (terme1 + terme2)
    return res / r.empty_mass



def _tau(r, e):
    c2 = _c2(r, e)
    a_star = r.nozzle_radius * _PI * r.nozzle_radius
    res = r.empty_volume / (a_star * c2)
    res = res * (2 / (e.gamma - 1))

    terme1 = (e.gamma + 1) / 2.0
    terme1 = math.pow(terme1, (e.gamma + 1) / (2 * (e.gamma - 1)))

    res = res * terme1
    return res


def t_gas(r, e):
    tau = _tau(r, e)
    p2 = press2(r, e)
    b = beta(e)

    term1 = p2 / (b * e.P_ext)
    term1 = math.pow(term1, (e.gamma - 1) / (2 * e.gamma))

    res = term1 - 1
    return tau * res


def ft_gas(r, e, t, v):
    tau = _tau(r,e)
    a_star = r.nozzle_radius * _PI * r.nozzle_radius
    terme1 = 2 / (e.gamma + 1)
    terme1 = terme1**(1/(e.gamma - 1))
    terme1 = terme1*2*a_star*press2(r,e)
    f_t = 1 + (t/tau)
    f_t = f_t**(2*e.gamma/(-e.gamma + 1))
    f_t = f_t*terme1 - (e.P_ext*a_star)
    f_t = f_t / r.empty_mass
    f = f_t - e.g - 0.5*r.surface*e.rho*r.drag_coeff
    return f




def runge_kutta_4(r, e, f, u0, t0, tf, dt):
    t = np.arange(t0, tf, dt)
    u = np.zeros(len(t))
    u[0] = u0
    
    for i in range(1, len(t)):
        k1 = dt * f(r,e,t[i-1], u[i-1])
        k2 = dt * f(r,e,t[i-1] + 0.5*dt, u[i-1] + 0.5*k1)
        k3 = dt * f(r,e,t[i-1] + 0.5*dt, u[i-1] + 0.5*k2)
        k4 = dt * f(r,e,t[i-1] + dt, u[i-1] + k3)
        u[i] = u[i-1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    
    return t, u



def v_e(r, e):
    p2 = press2(r, e)
    res = 2 * (p2 - e.P_ext)
    res = res / _RHO_EAU
    return math.sqrt(res)

def expulsion_eau_v1(T, r, e, time_intervals, velocities, dt):
    ve = v_e(r, e)
    a_star = r.nozzle_radius * _PI * r.nozzle_radius
    a = r.radius * _PI * r.radius

    debit = ve * a_star
    t_f = r.water_volume / debit
    T[0] = t_f

    m0 = r.empty_mass + (_RHO_EAU * r.water_volume)

    terme1 = _RHO_EAU * a_star * ve * t_f / m0
    terme1 = 1 - terme1

    vf = math.log(terme1)
    vf = -ve * vf - e.g * t_f
    T[1] = 0

    terme2 = math.log(terme1)
    terme3 = m0 / (_RHO_EAU * a)

    zf = ve * t_f * (1 - terme2)
    zf += terme3 * terme2
    zf -= 0.5 * e.g * math.pow(t_f, 2)
    T[2] = zf

    total_time = 0.0
    current_water_volume = r.water_volume

    while current_water_volume > 0:
        current_water_volume -= dt * ve * a_star
        if current_water_volume < 0:
            current_water_volume = 0
        r.water_volume = current_water_volume
        current_mass = r.empty_mass + _RHO_EAU * current_water_volume

        # Update velocity based on mass and thrust
        current_thrust = _RHO_EAU * a_star * ve * ve
        acceleration = current_thrust / current_mass - e.g
        T[1] += acceleration * dt

        time_intervals.append(total_time)
        velocities.append(T[1])
        total_time += dt
    T[1] = vf

def methode1(T, r, e):
    if r.pressure != 0 and r.water_volume != 0:
        time_intervals = []
        velocities = []
        dt = 0.001  

        expulsion_eau_v1(T, r, e, time_intervals, velocities, dt)

        t0 = 0.0
        tf = t_gas(r,e)
        dt = 0.001
        v = T[1]
        
        t_i, v_int = runge_kutta_4(r, e, ft_gas, v, t0, tf, dt)

        for i in range (len(t_i)) :
            t_i[i] += T[0]

        time_int = np.concatenate((time_intervals,t_i))
        velo = np.concatenate((velocities, v_int))

        T[1] += imp(r, e)
        T[0] += t_gas(r, e)
        v_moy = (v + T[1]) / 2.0
        T[2] += v_moy * t_gas(r, e)

        return time_int, velo
    return [], []


def simulate_rocket(r, e):
    T = [0.0, 0.0, 0.0]
    time_intervals, velocities = methode1(T, r, e)
    return time_intervals, velocities

r = Rocket(0.002, 0.0005, 0.05, 0.011, 0.1, 0.5, 700000, 0.0079)
e = Environment(9.81, 1.2, 293.15, 100000, 1.4)

time_intervals, velocities = simulate_rocket(r, e)


plt.plot(time_intervals, velocities)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Vitesse lors de la poussee')
plt.grid(True)
plt.show()
