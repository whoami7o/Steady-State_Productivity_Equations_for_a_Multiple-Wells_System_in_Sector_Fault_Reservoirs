# import block
import ssp_eq_mws_sfr as sspeq
import pandas as pd
from tqdm import tqdm


def get_input_path():
    msg = 'Enter .xlsx absolute file path or .xlsx name if script and file in the same dir:\n'
    path = input(msg).strip().replace('"', '').replace("'", '')
    try:
        excel_file = pd.read_excel(path, header=None)
        return excel_file
    except OSError:
        print('\nNo such directory or Wrong path name\n')
        return get_input_path()


def get_case_bound(dataframe):
    start = dataframe.index[dataframe[0] == 'START'].to_list()
    end = dataframe.index[dataframe[0] == 'END'].to_list()
    return start, end


def split_properties(dataframe):
    # boundaries
    resv_bound1 = int(dataframe.index[dataframe[0] == 'RESV'].to_list()[0])
    resv_bound2 = int(dataframe.index[dataframe[0] == 'WELL'].to_list()[0])
    # slicing
    resv_prop = dataframe.iloc[[_ for _ in range(resv_bound1 + 1, resv_bound2)]].reset_index(drop=True)
    well_prop = dataframe.iloc[resv_bound2 + 1:].reset_index(drop=True)
    # finishing
    resv_prop.columns = resv_prop.iloc[0]
    resv_prop = resv_prop[1:]

    mode_marker = well_prop.iloc[0].to_list()[5].upper()

    well_prop.columns = well_prop.iloc[0]
    well_prop = well_prop[1:].dropna(axis=1).drop(columns='index')

    return resv_prop.values.tolist(), well_prop.values.tolist(), mode_marker


def do_calculation(file):

    start_i, end_i = get_case_bound(file)

    output_file = pd.DataFrame()

    for i in tqdm(range(len(start_i))):
        case_i = file.iloc[[_ for _ in range(start_i[i] + 1, end_i[i])]].reset_index(drop=True)

        resv_i, wells_i, mode = split_properties(case_i)
        if mode == 'SKIN':
            q_sum, q = sspeq.find_q(wells_i, *resv_i)
            # print(f'CASE {i}\t{mode} is known\nQ_sum : {q_sum: .3f}\nQ_vec : {[round(_, 2) for _ in q]}')
            well_numeration = [f'well {i + 1}' for i in range(len(q))] + ['SUM']
            q_vec = q + [q_sum]
            solution = pd.concat([pd.DataFrame([f'CASE {i}', 'Q (rate)']), pd.DataFrame([well_numeration, q_vec])])
        elif mode == 'Q':
            s_vec = sspeq.find_s(wells_i, *resv_i)
            # print(f'CASE {i}\t{mode} is known\nS_vec : {[round(_, 2) for _ in s_vec]}')
            well_numeration = [f'well {i + 1}' for i in range(len(s_vec))]
            solution = pd.concat([pd.DataFrame([f'CASE {i}', 'SKIN']), pd.DataFrame([well_numeration, s_vec])])
        else:
            # will append empty df in case smth will go wrong
            solution = pd.DataFrame()

        output_file = pd.concat([output_file, solution])

    return output_file


if __name__ == '__main__':
    input_file = get_input_path()
    solution = do_calculation(input_file)
    solution.to_excel('cmd_solution.xlsx', index=False, header=False)
    print('DONE')
