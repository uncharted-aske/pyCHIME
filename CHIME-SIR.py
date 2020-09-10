import matplotlib.pyplot as plt

"""
This is a Python port of the CHIME-SIR.for program.
"""


def sir(s, i, r, beta, gamma, n):
    s_n = (-beta * s * i) + s
    i_n = (beta * s * i - gamma * i) + i
    r_n = gamma * i + r

    scale = n / (s_n + i_n + r_n)

    return s_n * scale, i_n * scale, r_n * scale


def sim_sir(s_n, i_n, r_n, gamma, i_day, N_p, policy_betas, policy_days, T, S, I, R):
    n = s_n + i_n + r_n
    d = i_day
    for p_idx in range(1, N_p+1):
        beta = policy_betas[p_idx-1]
        N_d = policy_days[p_idx-1]
        for d_idx in range(1, N_d+1):
            T.append(d)
            S.append(s_n)
            I.append(i_n)
            R.append(r_n)
            s_n, i_n, r_n = sir(s_n, i_n, r_n, beta, gamma, n)
            d += 1


def get_beta(intrinsic_growth_rate, gamma, susceptible, relative_contact_rate):
    inv_contact_rate = 1.0 - relative_contact_rate
    updated_growth_rate = intrinsic_growth_rate + gamma
    return updated_growth_rate / susceptible * inv_contact_rate


def get_growth_rate(doubling_time):
    if doubling_time == 0.0:
        return 0.0
    return 2.0 ** (1.0 / doubling_time) - 1.0


if __name__ == "__main__":
    s_n = 999.0
    i_n = 1.0
    r_n = 0.0

    beta = 0.0
    doubling_time = 0.0
    growth_rate = 0.0

    i_day = 17
    n_days = 20
    N_p = 3
    # N_t = 121
    infectious_days = 14
    relative_contact_rate = 0.05

    gamma = 1.0 / infectious_days

    policy_betas = []
    policy_days = []
    T = []
    S = []
    # E = []
    I = []
    R = []

    for p_idx in range(1, N_p+1):
        doubling_time = (p_idx - 1.0) * 5.0
        growth_rate = get_growth_rate(doubling_time)
        beta = get_beta(growth_rate, gamma, s_n, relative_contact_rate)
        policy_betas.append(beta)
        policy_days.append(n_days * p_idx)

    sim_sir(s_n, i_n, r_n, gamma, i_day, N_p, policy_betas, policy_days, T, S, I, R)

    plt.plot(T, S, label="Susceptible")
    plt.plot(T, I, label="Infected")
    plt.plot(T, R, label="Recovered")
    plt.xlabel('day')
    plt.ylabel('# of people')
    plt.legend()
    plt.show()