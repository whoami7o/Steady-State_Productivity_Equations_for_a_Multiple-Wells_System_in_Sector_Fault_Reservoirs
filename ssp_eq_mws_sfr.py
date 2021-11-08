# # function library for for Steady-State Productivity Equations for a Multiple-Wells System in Sector Fault Reservoirs

# import block
import numpy as np
from math import *


# set decimal places for np.arrays
np.set_printoptions(precision=3)


# # basic functions

def drawdown(wells_pressure: list[float], reservior_pressure: float) -> list[float]:
    # # calculating drawdawn vector
    return [reservior_pressure - p for p in wells_pressure]


def number(radial_permeability, payzone_thickness, oil_visc, form_vf, unit_conv_factor=86.4) -> float:
    # # calculating a number in Darcy's equation

    # radial permeability : K_r, m^2*e-6
    # payzone_thickness : H, m
    # oil_vics : mu, mPa*s
    # form_vf : B, m^3/m^3
    # unit_conv_factor : F_D = 86.4 for SI, 0.001127 for oil field units

    return (unit_conv_factor * 2*pi * radial_permeability * payzone_thickness) / (oil_visc * form_vf)


def get_avg_perm(radial_permeability: float, vertical_permeability: float) -> float:
    # # average permeability : K_a, m^2*e-6

    # radial permeability : K_r, m^2*e-6
    # vertical_permeability : K_r, m^2*e-6

    return radial_permeability**(2/3) * vertical_permeability**(1/3)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # Dimensionless Variables for fault drainage region

# x,y,z_coordinate : x,y,z , m
# drainage_radius : R_e, m
# average_perm : K_a, m^2*e-6
# radial permeability : K_r, m^2*e-6
# vertical_permeability : K_r, m^2*e-6
# r_coordinate : R_i, m
# wellbore_radius : R_w, m


# x_D
def get_x_d(x_coordinate, drainage_radius, average_perm, radial_permeability) -> float:
    return (x_coordinate/drainage_radius) * (average_perm/radial_permeability)**(1/2)


# y_D
def get_y_d(y_coordinate, drainage_radius, average_perm, radial_permeability) -> float:
    return (y_coordinate/drainage_radius) * (average_perm/radial_permeability)**(1/2)


# z_D
def get_z_d(z_coordinate, drainage_radius, average_perm, vertical_permeability) -> float:
    return (z_coordinate/drainage_radius) * (average_perm/vertical_permeability)**(1/2)


# H_D
def get_h_d(h, drainage_radius, average_perm, vertical_permeability) -> float:
    return (h/drainage_radius) * (average_perm/vertical_permeability)**(1/2)


# R_eD
def get_r_ed(average_perm: float, radial_permeability: float) -> float:
    return (average_perm/radial_permeability)**(1/2)


# R_jD
def get_r_jd(r_coordinate: float, drainage_radius: float, average_perm: float, radial_permeability: float) -> float:
    return (r_coordinate/drainage_radius) * (average_perm/radial_permeability)**(1/2)


# R_wD
def get_r_wd(wellbore_radius, drainage_radius, radial_permeability, vertical_permeability) -> float:
    return (wellbore_radius/drainage_radius) * ((radial_permeability/vertical_permeability)**(1/4) +
                                                (vertical_permeability/radial_permeability)**(1/4)) * (1/2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # calculating matrix's elements

# non-diagonal element
def get_a_ij(r_di, r_dj, phi_i, phi_j, tau) -> float:
    lambda_1 = ((r_di ** (2 * tau) + r_dj ** (2 * tau) - 2 * (r_di * r_dj) ** (tau) * cos(
        tau * (phi_i - phi_j) * pi / 180)) /
                (1 - 2 * ((r_di * r_dj) ** (tau)) * cos(tau * (phi_i - phi_j) * pi / 180) + (r_di * r_dj) ** (2 * tau)))

    lambda_2 = ((r_di ** (2 * tau) + r_dj ** (2 * tau) - 2 * (r_di * r_dj) ** (tau) * cos(
        tau * (phi_i + phi_j) * pi / 180)) /
                (1 - 2 * ((r_di * r_dj) ** (tau)) * cos(tau * (phi_i + phi_j) * pi / 180) + (r_di * r_dj) ** (2 * tau)))

    return (1 / 2) * log(lambda_1 * lambda_2)


# diagonal element
def get_a_ii(r_di, r_wdi, phi_i, tau) -> float:
    lambda_3 = 2 * tau * r_wdi * r_di ** (2 * tau - 1) * sin(tau * phi_i * pi / 180)

    lambda_4 = (1 - r_di ** (2 * tau)) * (1 - 2 * (r_di ** (2 * tau)) * cos(2 * tau * pi / 180) + r_di ** (4 * tau))

    return log(lambda_3 / lambda_4)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def get_a_matrix(wells_location: list[tuple], drainage_radius, average_perm, radial_perm, vertical_perm, tau):

    # # calculating [A] matrix aka influence matrix
    # well_location: [(r_i, phi_i, r_wi)]
    a_matrix = []

    for i in range(len(wells_location)):
        ai_row = []
        for j in range(len(wells_location)):
            if i == j:
                r_di = get_r_jd(wells_location[i][0], drainage_radius, average_perm, radial_perm)
                r_wdi = get_r_wd(wells_location[i][2], drainage_radius, radial_perm, vertical_perm)
                a = get_a_ii(r_di, r_wdi, wells_location[i][1], tau)
            else:
                r_di = get_r_jd(wells_location[i][0], drainage_radius, average_perm, radial_perm)
                r_dj = get_r_jd(wells_location[j][0], drainage_radius, average_perm, radial_perm)
                a = get_a_ij(r_di, r_dj, wells_location[i][1], wells_location[j][1], tau)

            ai_row.append(-a)

        a_matrix.append(ai_row)

    return np.array(a_matrix)


def get_ds_matrix(wells_skinf: list[float]):
    # # calculating [D_s] matrix aka diagonal matrix with skinfactors of each well as diagonal elements

    return np.diag(wells_skinf)


def get_dqinv_matrix(wells_production: list[float]):
    # # # calculating [D_q]^(-1) matrix

    return np.diag([1/q for q in wells_production])


def find_q(wells_info: list[tuple], reservoir_info: tuple) -> float and list:
    # # find production vector [Q] in case [S] is known
    # well properties should be : tuple(p_wf, r_i, phi_i, r_w, skin)
    # reservoir properties should be : tuple(p_e, r_e, PHI, h, k_r, k_z, mu, B)

    # create main vectors [wf pressure], [well coordinates], [skin], [wellbore radius]
    p_well_vec = []
    well_location_vec = []  # (r_i, phi_i, r_wi)
    skin_fac_vec = []

    for well in wells_info:
        p_well_vec.append(well[0])
        well_location_vec.append((well[1], well[2], well[3]))
        skin_fac_vec.append(well[4])

    # in-between calculations
    drawdown_vec = drawdown(p_well_vec, reservoir_info[0])
    darcys_number = number(reservoir_info[4], reservoir_info[3], reservoir_info[-2], reservoir_info[-1])
    avg_permeability = get_avg_perm(reservoir_info[-4], reservoir_info[-3])
    tau_number = 180/reservoir_info[2]

    influence_matrix = get_a_matrix(well_location_vec, reservoir_info[1], avg_permeability, reservoir_info[-4],
                                    reservoir_info[-3], tau_number)
    skin_matrix = get_ds_matrix(skin_fac_vec)

    # final answer
    q_vec = darcys_number * np.matmul(np.linalg.inv(np.add(influence_matrix, skin_matrix)), drawdown_vec)
    q_sum = sum(q_vec)

    return q_sum, q_vec.tolist()


def find_s(wells_info: list[tuple], reservoir_info: tuple) -> list:
    # # find production vector [S] in case [Q] is known
    # well properties should be : tuple(p_wf, r_i, phi_i, r_w, production)
    # reservoir properties should be : [p_e, r_e, PHI, h, k_r, k_z, mu, B]

    # create main vectors [wf pressure], [well coordinates], [Q], [wellbore radius]
    p_well_vec = []
    well_location_vec = []
    prod_vec = []

    for well in wells_info:
        p_well_vec.append(well[0])
        well_location_vec.append((well[1], well[2], well[3]))
        prod_vec.append(well[4])

    # in-between calculations
    drawdown_vec = drawdown(p_well_vec, reservoir_info[0])
    darcys_number = number(reservoir_info[4], reservoir_info[3], reservoir_info[-2], reservoir_info[-1])
    avg_permeability = get_avg_perm(reservoir_info[-4], reservoir_info[-3])
    tau_number = 180 / reservoir_info[2]

    influence_matrix = get_a_matrix(well_location_vec, reservoir_info[1], avg_permeability, reservoir_info[-4],
                                    reservoir_info[-3], tau_number)
    dqinv_matrix = get_dqinv_matrix(prod_vec)

    # final answer
    s = np.matmul(dqinv_matrix,
                  np.subtract(darcys_number * np.array(drawdown_vec), np.matmul(influence_matrix, prod_vec))
                  )

    return s.tolist()


if __name__ == '__main__':

    # # test full-solution functions
    # (p_wf, r_i, phi_i, r_w, skin)
    well_input = [
        (5.1, 300, 5, 0.1, 5),  # well_1
        (5.2, 400, 10, 0.1, 5),  # well_2
        (6.4, 800, 15, 0.1, 5),  # well_3
        (5.5, 1000, 20, 0.1, 5),  # well_4
        (4.8, 200, 25, 0.1, 5),  # well_5
        (8.7, 1500, 30, 0.1, 5),  # well_6
        (5.8, 600, 35, 0.1, 5),  # well_7
        (9.3, 2000, 40, 0.1, 5)  # well_8
    ]
    # [p_e, r_e, PHI, h, k_r, k_z, mu, B]
    reservoir_input = (18.0, 2500, 50, 20, 0.1, 0.025, 5.0, 1.25)

    sum_rate, rate_vec = find_q(well_input, reservoir_input)
    print(f'Q_sum = {sum_rate: .2f}\nQ = {[round(_, 2) for _ in rate_vec]}')

    # (p_wf, r_i, phi_i, r_w, skin)
    well_input_1 = [
        (5.1, 300, 5, 0.1, 17.91),  # well_1
        (5.2, 400, 10, 0.1, 10.80),  # well_2
        (6.4, 800, 15, 0.1, 39.74),  # well_3
        (5.5, 1000, 20, 0.1, 39.36),  # well_4
        (4.8, 200, 25, 0.1, 12.03),  # well_5
        (8.7, 1500, 30, 0.1, 32.59),  # well_6
        (5.8, 600, 35, 0.1, 34.03),  # well_7
        (9.3, 2000, 40, 0.1, 52.45)  # well_8
    ]
    skin_vec = find_s(well_input_1, reservoir_input)
    print(f'S = {[round(_, 2) for _ in skin_vec]}')

    # # # test for CASE 1# # #
    # print('[A] + [D_s]\n',
    #       np.add(
    #           get_a_matrix([(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)],
    #                        2500, get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50),
    #           get_ds_matrix([5 for _ in range(8)])), '\n'
    #       )
    #
    # print('([A] + [D_s])^(-1)\n',
    #       np.linalg.inv(np.add(
    #           get_a_matrix([(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)],
    #                        2500, get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50),
    #           get_ds_matrix([5 for _ in range(8)]))), '\n'
    #       )
    #
    # print('\n')
    # print('\t\t\t# # # SKIN # # #')
    # for skin in (0, 5, 10):
    #     print(f'S = {skin}', end='\t')
    #     print(
    #         f'Q_sum = {sum(number(0.1, 20, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix([(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500, get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50), get_ds_matrix([skin for _ in range(8)]))), np.transpose(np.array(drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0))))): .2f}')
    #
    #     print('\tQ = ',
    #           number(0.1, 20, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix(
    #               [(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500,
    #               get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50), get_ds_matrix([skin for _ in range(8)]))),
    #                                                np.transpose(np.array(
    #                                                    drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0)))), '\n'
    #           )
    #
    # print('\n')
    # print('\t\t\t# # # PHI # # #')
    # for phi in (60, 90, 120):
    #     print(f'PHI = {phi}', end=' ')
    #     print(
    #         f'Q_sum = {sum(number(0.1, 20, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix([(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500, get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / phi), get_ds_matrix([0 for _ in range(8)]))), np.transpose(np.array(drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0))))): .2f}')
    #
    #     print('\t Q = ',
    #           number(0.1, 20, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix(
    #               [(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500,
    #               get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / phi), get_ds_matrix([0 for _ in range(8)]))),
    #                                                np.transpose(np.array(
    #                                                    drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0)))), '\n'
    #           )
    #
    # print('\n')
    # print('\t\t\t# # # H # # #')
    # for h in (18, 20, 22):
    #     print(f'H = {h}', end='\t')
    #     print(
    #         f'Q_sum = {sum(number(0.1, h, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix([(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500, get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50), get_ds_matrix([0 for _ in range(8)]))), np.transpose(np.array(drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0))))): .2f}')
    #
    #     print('\tQ = ',
    #           number(0.1, h, 5, 1.25) * np.matmul(np.linalg.inv(np.add(get_a_matrix(
    #               [(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500,
    #               get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50), get_ds_matrix([0 for _ in range(8)]))),
    #                                               np.transpose(np.array(
    #                                                   drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0)))), '\n'
    #           )
    # # # # end # # #
    #
    # # # # test for CASE 2 # # #
    # print('S = ',
    #       np.matmul(get_dqinv_matrix([17.91, 10.80, 39.74, 39.36, 12.03, 32.59, 34.03, 52.45]), np.subtract(
    #           number(0.1, 20, 5, 1.25) * np.array(drawdown([5.1, 5.2, 6.4, 5.5, 4.8, 8.7, 5.8, 9.3], 18.0)), np.matmul(
    #               get_a_matrix(
    #                   [(300, 5), (400, 10), (800, 15), (1000, 20), (200, 25), (1500, 30), (600, 35), (2000, 40)], 2500,
    #                   get_avg_perm(0.1, 0.025), 0.1, 0.025, 0.1, 180 / 50),
    #               np.transpose(np.array([17.91, 10.80, 39.74, 39.36, 12.03, 32.59, 34.03, 52.45])))))
    #       )
    # # # # end # # #
