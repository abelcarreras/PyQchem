import numpy as np


# Simple order
# set higher dipole moment states first if energy gap is lower than eps_energy
def get_order_states_list(states, eps_moment=0.1, eps_energy=0.05):

    import itertools
    order = np.arange(len(states))

    for subset in itertools.combinations(range(len(states)), 2):
        i, j = subset

        if np.abs(states[i]['total energy'] - states[j]['total energy']) < eps_energy:
            tmi = np.linalg.norm(states[i]['transition moment'])
            tmj = np.linalg.norm(states[j]['transition moment'])
            if tmi - tmj < eps_moment:
                order[i], order[j] = order[j], order[i]

    return order


def correct_order_list(list, order):

    if np.array(order).shape != np.array(list).T.shape:
        print(np.array(order).shape, np.array(list).T.shape)
        raise Exception('Error in correcting order (not same shape)')

    alist = np.array(list)

    try:
        ordered_list = []
        for l, o in zip(alist.T, order):
            # print('l', l, 'o', o)
            ordered_list.append(l[o])
    except:
        ordered_list = []
        for l, o in zip(alist, order):
            # print('l', l, 'o', o)
            ordered_list.append(list[o])

    return np.array(ordered_list).T.tolist()
